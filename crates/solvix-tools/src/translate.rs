//! Translation rules: Abaqus script AST -> `solvix_core::Model`
//! (issue #35).
//!
//! v1 coverage (documented boundary):
//! * one TWO_D_PLANAR DEFORMABLE_BODY part per model; the sketch is
//!   reduced to its line/rectangle bounding box + `CircleByCenterRadius`
//!   holes; the mesh is the structured grid from [`crate::meshgen`]
//!   (CPS4), seeded by `seedPart(size=...)`;
//! * `Material`/`Elastic`/`Density`, `HomogeneousSolidSection`
//!   (material + thickness), `SectionAssignment` (region ignored —
//!   single body);
//! * `StaticStep`/`DynamicStep` (reduced to constant-dt Newmark),
//!   `DisplacementBC` (u1/u2/u3 -> dofs), `EncastreBC`,
//!   `ConcentratedForce` (cf1..cf3), `GravityLoad` (model-wide),
//!   `Pressure` (edge region -> element P-faces, best effort);
//! * regions: `vertices/edges.findAt(...)` geometric selections and
//!   named `a.Set`/`a.Surface` regions. Geometric edges resolve against
//!   the rectangle boundary sides; vertices to the nearest mesh node.
//!   Keyword arguments only — positional arguments beyond `findAt`
//!   data are warned and ignored.

use std::collections::{HashMap, HashSet};

use solvix_core::{
    BoundaryCondition, ConcentratedLoad, DistributedLoad, DynamicsSpec, ElementSet, Material,
    Model, Step,
};

use crate::meshgen::{self, Hole};
use crate::python_parse::{eval_expr, Expr, Script, Stmt, StmtKind, Value};

/// Fatal translation error (nothing solvable can be produced).
#[derive(Debug, Clone, PartialEq)]
pub struct TranslateError {
    pub line: usize,
    pub message: String,
}

impl std::fmt::Display for TranslateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

fn terr<T>(line: usize, message: impl Into<String>) -> Result<T, TranslateError> {
    Err(TranslateError {
        line,
        message: message.into(),
    })
}

/// Non-fatal notes: `warnings` = degradation we recovered from,
/// `unsupported` = constructs recognized but NOT translated.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Report {
    pub warnings: Vec<String>,
    pub unsupported: Vec<String>,
}

pub struct Translation {
    pub model: Model,
    pub report: Report,
}

pub fn translate(src: &str) -> Result<Translation, TranslateError> {
    let script: Script = crate::python_parse::parse_script(src).map_err(|e| TranslateError {
        line: e.line,
        message: e.message,
    })?;
    let mut t = Translator::default();
    for st in &script.stmts {
        t.stmt(st)?;
    }
    t.finish()
}

// ---------------------------------------------------------------------
// Intermediate objects
// ---------------------------------------------------------------------

/// Geometric region reference, resolved to node ids once the mesh
/// exists.
#[derive(Debug, Clone)]
enum Region {
    /// Vertex selections (`vertices.findAt`): nearest mesh node.
    Coords(Vec<[f64; 2]>),
    /// Edge selections (`edges.findAt`): rectangle-boundary sides by
    /// midpoint.
    EdgeMids(Vec<[f64; 2]>),
    /// Named `a.Set`/`a.Surface` region.
    Named(String),
    /// A whole-part geometry selection (`p.faces`, `p.Set(faces=...)`).
    AllNodes,
}

enum Obj {
    Sketch(SketchData),
    Part(usize),
    Assembly,
    Region(Region),
    Opaque,
}

#[derive(Default, Clone)]
struct SketchData {
    segments: Vec<([f64; 2], [f64; 2])>,
    circles: Vec<(f64, f64, f64)>,
}

struct PartData {
    name: String,
    line: usize,
    sketch: Option<SketchData>,
    seed: Option<f64>,
    meshed: bool,
    section: Option<String>,
}

struct MatAcc {
    name: String,
    elastic: Option<(f64, f64)>,
    density: Option<f64>,
}

#[derive(Clone)]
struct PendingBc {
    region: Region,
    dofs: Vec<usize>,
}

#[derive(Clone)]
struct StepAcc {
    name: String,
    dynamics: Option<DynamicsSpec>,
    bcs: Vec<PendingBc>,
    cloads: Vec<(Region, usize, f64)>,
    pressures: Vec<(Region, f64)>,
    gravity: Option<(f64, [f64; 3])>,
}

impl StepAcc {
    fn new(name: &str) -> Self {
        StepAcc {
            name: name.to_string(),
            dynamics: None,
            bcs: Vec::new(),
            cloads: Vec::new(),
            pressures: Vec::new(),
            gravity: None,
        }
    }
}

/// A decomposed call statement: resolved attribute/index path plus raw
/// arguments. `mdb.models['M'].StaticStep(name=...)` becomes
/// `path = [mdb, models, Model-1, StaticStep]`.
struct Call<'a> {
    path: Vec<String>,
    pos: &'a [Expr],
    kwargs: &'a [(String, Expr)],
    line: usize,
}

// ---------------------------------------------------------------------
// Translator
// ---------------------------------------------------------------------

#[derive(Default)]
struct Translator {
    data: HashMap<String, Value>,
    objs: HashMap<String, Obj>,
    mats: Vec<MatAcc>,
    sections: Vec<(String, String, Option<f64>)>,
    parts: Vec<PartData>,
    sets: Vec<(String, Region)>,
    steps: Vec<StepAcc>,
    initial_bcs: Vec<PendingBc>,
    model_name: Option<String>,
    report: Report,
    warned: HashSet<&'static str>,
}

impl Translator {
    fn warn(&mut self, line: usize, msg: impl AsRef<str>) {
        self.report
            .warnings
            .push(format!("line {line}: {}", msg.as_ref()));
    }

    fn warn_once(&mut self, key: &'static str, line: usize, msg: impl AsRef<str>) {
        if self.warned.insert(key) {
            self.warn(line, msg);
        }
    }

    fn unsupp(&mut self, line: usize, msg: impl AsRef<str>) {
        self.report
            .unsupported
            .push(format!("line {line}: {}", msg.as_ref()));
    }

    // -- path/expression helpers --------------------------------------

    fn flatten(&self, e: &Expr) -> Option<Vec<String>> {
        match e {
            Expr::Name(n) => Some(vec![n.clone()]),
            Expr::Attr(b, n) => {
                let mut v = self.flatten(b)?;
                v.push(n.clone());
                Some(v)
            }
            Expr::Index(b, i) => {
                let key = match eval_expr(i, &self.data) {
                    Ok(Value::Str(s)) => s,
                    Ok(Value::Num(x)) => format!("{x}"),
                    _ => return None,
                };
                let mut v = self.flatten(b)?;
                v.push(key);
                Some(v)
            }
            _ => None,
        }
    }

    fn as_call<'s>(&self, e: &'s Expr, line: usize) -> Option<Call<'s>> {
        let Expr::Call { func, pos, kwargs } = e else {
            return None;
        };
        Some(Call {
            path: self.flatten(func)?,
            pos,
            kwargs,
            line,
        })
    }

    fn kwarg<'k>(&self, kwargs: &'k [(String, Expr)], name: &str) -> Option<&'k Expr> {
        kwargs.iter().find(|(k, _)| k == name).map(|(_, e)| e)
    }

    fn kwarg_val(&self, kwargs: &[(String, Expr)], name: &str) -> Result<Option<Value>, String> {
        match self.kwarg(kwargs, name) {
            None => Ok(None),
            Some(e) => eval_expr(e, &self.data).map(Some),
        }
    }

    fn kwarg_str(&self, kwargs: &[(String, Expr)], name: &str) -> Result<Option<String>, String> {
        match self.kwarg_val(kwargs, name)? {
            Some(Value::Str(s)) => Ok(Some(s)),
            Some(_) => Err(format!("'{name}' must be a string")),
            None => Ok(None),
        }
    }

    fn kwarg_num(&self, kwargs: &[(String, Expr)], name: &str) -> Result<Option<f64>, String> {
        match self.kwarg_val(kwargs, name)? {
            Some(v) => v
                .as_f64()
                .ok_or_else(|| format!("'{name}' must be a number"))
                .map(Some),
            None => Ok(None),
        }
    }

    /// A `(x, y)` pair kwarg like `point1=(-0.5, -0.05)`.
    fn kwarg_point(&self, c: &Call, name: &str) -> Option<[f64; 2]> {
        match self.kwarg_val(c.kwargs, name) {
            Ok(Some(Value::Tuple(v)) | Some(Value::List(v))) => {
                let cc: Vec<f64> = v.iter().filter_map(Value::as_f64).collect();
                if cc.len() >= 2 {
                    Some([cc[0], cc[1]])
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    fn warn_pos(&mut self, c: &Call, what: &str) {
        if !c.pos.is_empty() {
            self.warn(
                c.line,
                format!("positional argument(s) in {what} ignored (v1: keyword arguments only)"),
            );
        }
    }

    /// True when `path` addresses the (single) model; warns on
    /// multi-model scripts.
    fn model_path_ok(&mut self, m: &str, line: usize) -> bool {
        match &self.model_name {
            None => {
                self.model_name = Some(m.to_string());
                true
            }
            Some(first) if first == m => true,
            Some(first) => {
                self.warn_once(
                    "multi-model",
                    line,
                    format!(
                        "model '{m}' ignored (v1: single model — everything lands in '{first}')"
                    ),
                );
                false
            }
        }
    }

    // -- statement dispatch --------------------------------------------

    fn stmt(&mut self, st: &Stmt) -> Result<(), TranslateError> {
        match &st.kind {
            StmtKind::Import | StmtKind::Noop => Ok(()),
            StmtKind::Block { keyword } => {
                self.unsupp(
                    st.line,
                    format!("'{keyword}' block skipped (v1: no control flow in the translatable subset)"),
                );
                Ok(())
            }
            StmtKind::Assign { target, value } => self.assign(target, value, st.line),
            StmtKind::Expr(e) => self.expr_stmt(e, st.line),
        }
    }

    fn expr_stmt(&mut self, e: &Expr, line: usize) -> Result<(), TranslateError> {
        if matches!(e, Expr::Name(_)) {
            return Ok(());
        }
        let Some(c) = self.as_call(e, line) else {
            self.unsupp(line, "statement ignored (not a supported call)");
            return Ok(());
        };
        self.action(&c)
    }

    fn assign(&mut self, target: &str, value: &Expr, line: usize) -> Result<(), TranslateError> {
        if let Some(c) = self.as_call(value, line) {
            self.bind(target, &c);
            return Ok(());
        }
        if let Some(path) = self.flatten(value) {
            if path.len() == 4
                && path[0] == "mdb"
                && path[1] == "models"
                && path[3].eq_ignore_ascii_case("rootAssembly")
                && self.model_path_ok(&path[2], line)
            {
                self.objs.insert(target.to_string(), Obj::Assembly);
                return Ok(());
            }
            if path.len() == 1 {
                if let Some(o) = self.objs.get(&path[0]) {
                    let o = clone_obj(o);
                    self.objs.insert(target.to_string(), o);
                    return Ok(());
                }
            }
            if let Some(r) = self.region_from_path(&path) {
                self.objs.insert(target.to_string(), Obj::Region(r));
                return Ok(());
            }
            self.unsupp(line, "assignment ignored (unsupported expression)");
            return Ok(());
        }
        match eval_expr(value, &self.data) {
            Ok(v) => {
                self.data.insert(target.to_string(), v);
            }
            Err(_) => self.unsupp(line, "cannot evaluate assignment — ignored"),
        }
        Ok(())
    }

    /// Bind the target of `x = <call>()` to the object the call
    /// creates.
    fn bind(&mut self, target: &str, c: &Call) {
        let p = &c.path;
        let last = p.last().map(|s| s.as_str()).unwrap_or_default();
        // Region-producing assignments (`region = ...findAt(...)`).
        if let Some(r) = self.region_from_call(c) {
            self.objs.insert(target.to_string(), Obj::Region(r));
            return;
        }
        if let Some(r) = self.region_from_path(p) {
            self.objs.insert(target.to_string(), Obj::Region(r));
            return;
        }
        // mdb.models['M'].ConstrainedSketch(...)
        if p.len() == 4
            && p[0] == "mdb"
            && p[1] == "models"
            && last.eq_ignore_ascii_case("ConstrainedSketch")
        {
            if self.model_path_ok(&p[2], c.line) {
                self.objs
                    .insert(target.to_string(), Obj::Sketch(SketchData::default()));
            }
            return;
        }
        // mdb.models['M'].Part(...)
        if p.len() == 4 && p[0] == "mdb" && p[1] == "models" && last.eq_ignore_ascii_case("Part") {
            let dim = self
                .kwarg(c.kwargs, "dimensionality")
                .and_then(|e| match e {
                    Expr::Name(n) => Some(n.clone()),
                    _ => None,
                })
                .unwrap_or_default();
            let name = self
                .kwarg_str(c.kwargs, "name")
                .unwrap_or_default()
                .unwrap_or_else(|| "Part-1".into());
            if dim.eq_ignore_ascii_case("TWO_D_PLANAR") {
                if !self.parts.is_empty() {
                    self.warn(c.line, "multiple parts unsupported (v1: first part wins)");
                    self.objs.insert(target.to_string(), Obj::Opaque);
                    return;
                }
                self.parts.push(PartData {
                    name,
                    line: c.line,
                    sketch: None,
                    seed: None,
                    meshed: false,
                    section: None,
                });
                self.objs
                    .insert(target.to_string(), Obj::Part(self.parts.len() - 1));
            } else {
                self.unsupp(
                    c.line,
                    format!("Part dimensionality '{dim}' unsupported (v1: TWO_D_PLANAR only)"),
                );
                self.objs.insert(target.to_string(), Obj::Opaque);
            }
            return;
        }
        // mesh.ElemType(...)
        if p.len() == 2 && p[0] == "mesh" && last.eq_ignore_ascii_case("ElemType") {
            self.unsupp(c.line, "mesh.ElemType ignored (v1: CPS4 structured grid)");
            self.objs.insert(target.to_string(), Obj::Opaque);
            return;
        }
        if p.len() >= 2 {
            // Object method call used as an assignment (`a.Set(...)`).
            let head = p[0].clone();
            if let Some(obj) = self.objs.get(&head).map(clone_obj) {
                match obj {
                    Obj::Assembly => {
                        let m = last.to_ascii_lowercase();
                        if m == "set" || m == "surface" {
                            self.assembly_set_or_surface(c);
                            return;
                        }
                        if m == "instance" {
                            self.assembly_instance(c);
                            self.objs.insert(target.to_string(), Obj::Opaque);
                            return;
                        }
                        if m == "datumcsysbydefault" {
                            self.unsupp(c.line, "assembly datums unsupported");
                            self.objs.insert(target.to_string(), Obj::Opaque);
                            return;
                        }
                    }
                    Obj::Part(_) if last.eq_ignore_ascii_case("Set") => {
                        self.warn_once(
                            "geom-set",
                            c.line,
                            "geometry set covers the whole part (v1: no geometry kernel)",
                        );
                        self.objs
                            .insert(target.to_string(), Obj::Region(Region::AllNodes));
                        return;
                    }
                    _ => {}
                }
            }
        }
        self.unsupp(
            c.line,
            format!(
                "assignment '{} = {}(...)' unsupported — ignored",
                target, last
            ),
        );
        self.objs.insert(target.to_string(), Obj::Opaque);
    }

    /// Top-level action statement.
    fn action(&mut self, c: &Call) -> Result<(), TranslateError> {
        let p = &c.path;
        match p.first().map(String::as_str) {
            Some("session") | Some("print") => Ok(()), // display bookkeeping
            Some("sys") => Ok(()),
            Some("mdb") => {
                if p.len() == 2 && matches!(p[1].as_str(), "save" | "saveAs" | "close") {
                    return Ok(());
                }
                self.mdb_action(c)
            }
            Some(v) => self.obj_action(v, c),
            None => Ok(()),
        }
    }

    fn mdb_action(&mut self, c: &Call) -> Result<(), TranslateError> {
        let p = &c.path;
        if p.len() < 4 || p[1] != "models" {
            self.unsupp(c.line, format!("mdb call '{}' unsupported", p.join(".")));
            return Ok(());
        }
        if !self.model_path_ok(&p[2], c.line) {
            return Ok(());
        }
        let rest = &p[3..];
        let method = rest.first().map(|s| s.as_str()).unwrap_or_default();
        let ml = method.to_ascii_lowercase();

        // materials['Steel'].Elastic / .Density
        if rest.len() == 3 && rest[0].eq_ignore_ascii_case("materials") {
            let mat = rest[1].clone();
            let prop = rest[2].to_ascii_lowercase();
            let table = match self.kwarg_val(c.kwargs, "table") {
                Ok(t) => t,
                Err(e) => {
                    self.warn(c.line, e);
                    return Ok(());
                }
            };
            let Some(Value::Tuple(rows) | Value::List(rows)) = table else {
                self.warn(c.line, "'table' must be a tuple of rows — ignored");
                return Ok(());
            };
            let flat: Vec<Vec<f64>> = rows
                .iter()
                .map(|row| match row {
                    Value::Tuple(r) | Value::List(r) => {
                        r.iter().filter_map(Value::as_f64).collect()
                    }
                    v => v.as_f64().map(|x| vec![x]).unwrap_or_default(),
                })
                .collect();
            let Some(m) = self.mats.iter_mut().find(|m| m.name == mat) else {
                self.unsupp(
                    c.line,
                    format!("'{prop}' for unknown material '{mat}' — ignored"),
                );
                return Ok(());
            };
            match prop.as_str() {
                "elastic" => {
                    if let Some(row) = flat.first() {
                        if row.len() >= 2 {
                            m.elastic = Some((row[0], row[1]));
                        } else {
                            self.warn(c.line, "*ELASTIC row needs (E, nu) — ignored");
                        }
                    }
                }
                "density" => {
                    if let Some(row) = flat.first() {
                        if let Some(&d) = row.first() {
                            m.density = Some(d);
                        }
                    }
                }
                _ => self.unsupp(
                    c.line,
                    format!("material property '{prop}' unsupported (v1: Elastic, Density)"),
                ),
            }
            return Ok(());
        }

        match ml.as_str() {
            "material" => {
                self.warn_pos(c, "Material");
                match self.kwarg_str(c.kwargs, "name") {
                    Ok(Some(n)) => {
                        if self.mats.iter().any(|m| m.name == n) {
                            self.warn(c.line, format!("material '{n}' redefined"));
                        }
                        self.mats.push(MatAcc {
                            name: n,
                            elastic: None,
                            density: None,
                        });
                    }
                    Ok(None) => self.warn(c.line, "Material without a name — ignored"),
                    Err(e) => self.warn(c.line, e),
                }
            }
            "homogeneoussolidsection" => {
                self.warn_pos(c, "HomogeneousSolidSection");
                match (
                    self.kwarg_str(c.kwargs, "name"),
                    self.kwarg_str(c.kwargs, "material"),
                    self.kwarg_num(c.kwargs, "thickness"),
                ) {
                    (Ok(Some(n)), Ok(Some(m)), Ok(t)) => {
                        self.sections.push((n, m, t));
                    }
                    _ => self.warn(
                        c.line,
                        "HomogeneousSolidSection needs name/material — ignored",
                    ),
                }
            }
            "staticstep" => {
                self.warn_pos(c, "StaticStep");
                let name = self
                    .kwarg_str(c.kwargs, "name")
                    .unwrap_or_default()
                    .unwrap_or_else(|| format!("Step-{}", self.steps.len() + 1));
                if let Some(prev) = self.kwarg_str(c.kwargs, "previous").unwrap_or_default() {
                    if let Some(last) = self.steps.last() {
                        if !prev.is_empty() && prev != last.name && prev != "Initial" {
                            self.warn(
                                c.line,
                                format!("previous='{prev}' is not the last step '{}' — steps kept in creation order", last.name),
                            );
                        }
                    }
                }
                if self.steps.iter().any(|s| s.name == name) {
                    self.warn(c.line, format!("duplicate step name '{name}'"));
                }
                self.steps.push(StepAcc::new(&name));
            }
            "dynamicstep" => {
                self.warn_pos(c, "DynamicStep");
                let name = self
                    .kwarg_str(c.kwargs, "name")
                    .unwrap_or_default()
                    .unwrap_or_else(|| format!("Step-{}", self.steps.len() + 1));
                let dur = self
                    .kwarg_num(c.kwargs, "timePeriod")
                    .unwrap_or_default()
                    .unwrap_or(1.0);
                let inc = self
                    .kwarg_num(c.kwargs, "initialInc")
                    .unwrap_or_default()
                    .unwrap_or(dur / 10.0);
                let mut acc = StepAcc::new(&name);
                if dur > 0.0 && inc > 0.0 {
                    let n = (dur / inc).round().max(1.0);
                    acc.dynamics = Some(DynamicsSpec {
                        dt: dur / n,
                        duration: dur,
                        rayleigh_alpha: 0.0,
                        rayleigh_beta: 0.0,
                    });
                    self.warn(
                        c.line,
                        format!("DynamicStep reduced to constant dt = {} Newmark", dur / n),
                    );
                } else {
                    self.warn(
                        c.line,
                        "DynamicStep with non-positive time — treated as static",
                    );
                }
                self.steps.push(acc);
            }
            "displacementbc" | "encastrebc" => {
                self.warn_pos(c, method);
                self.displacement_bc(c, ml == "encastrebc")?;
            }
            "concentratedforce" => {
                self.warn_pos(c, "ConcentratedForce");
                let create = self
                    .kwarg_str(c.kwargs, "createStepName")
                    .unwrap_or_default();
                if create.as_deref() == Some("Initial") {
                    self.warn(c.line, "loads cannot live in the Initial step — ignored");
                    return Ok(());
                }
                let forces: Vec<(usize, f64)> = ["cf1", "cf2", "cf3"]
                    .iter()
                    .enumerate()
                    .filter_map(|(i, k)| {
                        self.kwarg_num(c.kwargs, k)
                            .ok()
                            .flatten()
                            .filter(|v| *v != 0.0)
                            .map(|v| (i + 1, v))
                    })
                    .collect();
                if forces.is_empty() {
                    self.warn(c.line, "ConcentratedForce without cf1/cf2/cf3 — ignored");
                    return Ok(());
                }
                let Some(region) = self.region_kwarg(c) else {
                    return Ok(());
                };
                let Some(idx) = self.step_idx(create, c.line) else {
                    return Ok(());
                };
                for (dof, v) in forces {
                    self.steps[idx].cloads.push((region.clone(), dof, v));
                }
            }
            "pressure" => {
                self.warn_pos(c, "Pressure");
                match self.kwarg_num(c.kwargs, "magnitude") {
                    Ok(Some(mag)) => {
                        let region = self.region_kwarg(c);
                        let create = self
                            .kwarg_str(c.kwargs, "createStepName")
                            .unwrap_or_default();
                        if create.as_deref() == Some("Initial") {
                            self.warn(c.line, "loads cannot live in the Initial step — ignored");
                            return Ok(());
                        }
                        let Some(idx) = self.step_idx(create, c.line) else {
                            return Ok(());
                        };
                        if let Some(r) = region {
                            self.steps[idx].pressures.push((r, mag));
                        }
                    }
                    _ => self.warn(c.line, "Pressure without a numeric magnitude — ignored"),
                }
            }
            "gravityload" => {
                self.warn_pos(c, "GravityLoad");
                if self.kwarg(c.kwargs, "region").is_some() {
                    self.warn_once(
                        "grav-region",
                        c.line,
                        "GravityLoad region ignored (v1: gravity is model-wide)",
                    );
                }
                let g = match self.kwarg_val(c.kwargs, "gravitation") {
                    Ok(g) => g,
                    Err(e) => {
                        self.warn(c.line, e);
                        return Ok(());
                    }
                };
                let (Some(Value::Tuple(dir)) | Some(Value::List(dir))) = g else {
                    self.warn(c.line, "gravitation must be a 3-tuple — ignored");
                    return Ok(());
                };
                let v: Vec<f64> = dir.iter().filter_map(Value::as_f64).collect();
                if v.len() == 3 {
                    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                    let create = self
                        .kwarg_str(c.kwargs, "createStepName")
                        .unwrap_or_default();
                    let Some(idx) = self.step_idx(create, c.line) else {
                        return Ok(());
                    };
                    if mag > 0.0 {
                        let n = [v[0] / mag, v[1] / mag, v[2] / mag];
                        self.steps[idx].gravity = Some((mag, n));
                    }
                } else {
                    self.warn(c.line, "gravitation must have 3 components — ignored");
                }
            }
            _ => {
                self.unsupp(
                    c.line,
                    format!("mdb.models[..].{method} unsupported (v1: Material, HomogeneousSolidSection, StaticStep, DynamicStep, DisplacementBC, EncastreBC, ConcentratedForce, Pressure, GravityLoad)"),
                );
            }
        }
        Ok(())
    }

    fn displacement_bc(&mut self, c: &Call, encastre: bool) -> Result<(), TranslateError> {
        let mut dofs: Vec<usize> = Vec::new();
        if encastre {
            dofs = vec![1, 2];
        } else {
            for (i, key) in ["u1", "u2", "u3"].iter().enumerate() {
                match self.kwarg(c.kwargs, key) {
                    None => {}
                    Some(Expr::Num(_)) => dofs.push(i + 1),
                    Some(Expr::Name(n)) if n == "SET" => dofs.push(i + 1),
                    Some(Expr::Name(n)) if n == "FREE" || n == "UNSET" || n == "None" => {}
                    Some(other) => match eval_expr(other, &self.data) {
                        Ok(Value::Num(v)) => {
                            if v != 0.0 {
                                self.warn(
                                    c.line,
                                    format!("nonzero {key}={v} ignored (v1: homogeneous BCs only) — dof fixed"),
                                );
                            }
                            dofs.push(i + 1);
                        }
                        _ => self.warn(c.line, format!("cannot interpret {key} — ignored")),
                    },
                }
            }
        }
        for key in ["u4", "u5", "u6"] {
            if self.kwarg(c.kwargs, key).is_some() {
                self.unsupp(
                    c.line,
                    format!("{key} rotation dof not supported (plane elements)"),
                );
            }
        }
        if dofs.is_empty() {
            self.warn(c.line, "boundary condition fixes no dofs — ignored");
            return Ok(());
        }
        let create = self
            .kwarg_str(c.kwargs, "createStepName")
            .unwrap_or_default();
        if create.as_deref() == Some("Initial") {
            if let Some(r) = self.region_kwarg(c) {
                self.initial_bcs.push(PendingBc { region: r, dofs });
            }
            return Ok(());
        }
        let Some(region) = self.region_kwarg(c) else {
            return Ok(());
        };
        let Some(idx) = self.step_idx(create, c.line) else {
            return Ok(());
        };
        self.steps[idx].bcs.push(PendingBc { region, dofs });
        Ok(())
    }

    /// Index of the step named `create`, auto-creating it (with a
    /// warning) when the script references an unknown step.
    fn step_idx(&mut self, create: Option<String>, line: usize) -> Option<usize> {
        let name = create.unwrap_or_else(|| "Step-1".to_string());
        if let Some(i) = self.steps.iter().position(|s| s.name == name) {
            return Some(i);
        }
        self.warn(
            line,
            format!("unknown createStepName '{name}' — step auto-created"),
        );
        self.steps.push(StepAcc::new(&name));
        Some(self.steps.len() - 1)
    }

    /// The `region=` kwarg resolved to a [`Region`], if present and
    /// understandable.
    fn region_kwarg(&mut self, c: &Call) -> Option<Region> {
        let e = self.kwarg(c.kwargs, "region")?;
        self.region_of(e, c.line)
    }

    fn region_of(&mut self, e: &Expr, line: usize) -> Option<Region> {
        if let Expr::Name(n) = e {
            if let Some(Obj::Region(r)) = self.objs.get(n).map(clone_obj) {
                return Some(r);
            }
            self.warn(line, format!("unknown region variable '{n}'"));
            return None;
        }
        if let Some(c) = self.as_call(e, line) {
            if let Some(r) = self.region_from_call(&c) {
                return Some(r);
            }
        } else if let Some(path) = self.flatten(e) {
            if let Some(r) = self.region_from_path(&path) {
                return Some(r);
            }
        }
        self.unsupp(line, "unsupported region expression — record skipped");
        None
    }

    fn region_from_call(&mut self, c: &Call) -> Option<Region> {
        let p = &c.path;
        if p.len() >= 2 && p[p.len() - 1].eq_ignore_ascii_case("findAt") {
            let pts = self.points(c.pos, c.line);
            let head = &p[p.len() - 2];
            if head.eq_ignore_ascii_case("vertices") {
                return Some(Region::Coords(pts));
            }
            if head.eq_ignore_ascii_case("edges") {
                return Some(Region::EdgeMids(pts));
            }
        }
        None
    }

    fn region_from_path(&self, p: &[String]) -> Option<Region> {
        if p.len() >= 3 {
            let (head, key) = (&p[p.len() - 2], &p[p.len() - 1]);
            if head.eq_ignore_ascii_case("sets") || head.eq_ignore_ascii_case("surfaces") {
                return Some(Region::Named(key.clone()));
            }
        }
        if p.last().is_some_and(|s| {
            ["vertices", "edges", "faces", "cells"]
                .iter()
                .any(|x| s.eq_ignore_ascii_case(x))
        }) {
            return Some(Region::AllNodes);
        }
        if p.len() == 1 {
            if let Some(Obj::Region(r)) = self.objs.get(&p[0]).map(clone_obj) {
                return Some(r);
            }
        }
        None
    }

    /// `findAt(((x, y, z),), ...)` points from the positional args.
    fn points(&mut self, args: &[Expr], line: usize) -> Vec<[f64; 2]> {
        let mut pts = Vec::new();
        for a in args {
            match eval_expr(a, &self.data) {
                Ok(Value::Tuple(rows)) | Ok(Value::List(rows)) => {
                    for row in rows {
                        match row {
                            Value::Tuple(v) | Value::List(v) => {
                                let c: Vec<f64> = v.iter().filter_map(Value::as_f64).collect();
                                if c.len() >= 2 {
                                    pts.push([c[0], c[1]]);
                                }
                            }
                            Value::Num(x) => pts.push([x, 0.0]),
                            _ => {}
                        }
                    }
                }
                _ => self.warn(line, "findAt argument not understood — point skipped"),
            }
        }
        pts
    }

    // -- object method actions -----------------------------------------

    fn obj_action(&mut self, var: &str, c: &Call) -> Result<(), TranslateError> {
        let Some(obj) = self.objs.get(var).map(clone_obj) else {
            self.unsupp(c.line, format!("unknown variable '{var}'"));
            return Ok(());
        };
        match obj {
            Obj::Sketch(_) => self.sketch_action(var, c),
            Obj::Part(_) => self.part_action(var, c),
            Obj::Assembly => self.assembly_action(c),
            Obj::Region(_) | Obj::Opaque => {
                self.unsupp(c.line, "call on a non-actionable object — ignored");
                Ok(())
            }
        }
    }

    fn sketch_action(&mut self, var: &str, c: &Call) -> Result<(), TranslateError> {
        let method = c.path.last().map(|s| s.as_str()).unwrap_or_default();
        let m = method.to_ascii_lowercase();
        let line = c.line;
        let Some(Obj::Sketch(mut sk)) = self.objs.get(var).map(clone_obj) else {
            return Ok(());
        };
        match m.as_str() {
            "rectangle" => {
                let (Some(p1), Some(p2)) =
                    (self.kwarg_point(c, "point1"), self.kwarg_point(c, "point2"))
                else {
                    self.warn(line, "rectangle without point1/point2 — ignored");
                    return Ok(());
                };
                sk.segments.push((p1, [p2[0], p1[1]]));
                sk.segments.push(([p2[0], p1[1]], p2));
                sk.segments.push((p2, [p1[0], p2[1]]));
                sk.segments.push(([p1[0], p2[1]], p1));
            }
            "line" => {
                let (Some(p1), Some(p2)) =
                    (self.kwarg_point(c, "point1"), self.kwarg_point(c, "point2"))
                else {
                    self.warn(line, "Line without point1/point2 — ignored");
                    return Ok(());
                };
                sk.segments.push((p1, p2));
            }
            "circlebycenterradius" => {
                let (Some(center), radius) = (
                    self.kwarg_point(c, "center"),
                    self.kwarg_num(c.kwargs, "radius").ok().flatten(),
                ) else {
                    self.warn(line, "CircleByCenterRadius without center/radius — ignored");
                    return Ok(());
                };
                if let Some(r) = radius {
                    sk.circles.push((center[0], center[1], r));
                } else {
                    self.warn(line, "CircleByCenterRadius without radius — ignored");
                }
            }
            "setprimaryobject" | "unsetprimaryobject" | "unplaceprimaryobject" => {}
            _ => {
                self.unsupp(
                    line,
                    format!("sketch primitive '{method}' unsupported (v1: rectangle, Line, CircleByCenterRadius)"),
                );
                return Ok(());
            }
        }
        self.objs.insert(var.to_string(), Obj::Sketch(sk));
        Ok(())
    }

    fn part_action(&mut self, var: &str, c: &Call) -> Result<(), TranslateError> {
        let method = c.path.last().map(|s| s.as_str()).unwrap_or_default();
        let m = method.to_ascii_lowercase();
        let line = c.line;
        let Some(Obj::Part(idx)) = self.objs.get(var).map(clone_obj) else {
            return Ok(());
        };
        match m.as_str() {
            "baseshell" | "basesolid2d" => {
                let sk = self.kwarg(c.kwargs, "sketch").and_then(|e| match e {
                    Expr::Name(n) => Some(n.clone()),
                    _ => None,
                });
                let Some(skvar) = sk else {
                    self.warn(
                        line,
                        format!("{method} without sketch=<ConstrainedSketch> — ignored"),
                    );
                    return Ok(());
                };
                let Some(Obj::Sketch(data)) = self.objs.get(&skvar).map(clone_obj) else {
                    self.warn(line, format!("sketch '{skvar}' is not a ConstrainedSketch"));
                    return Ok(());
                };
                self.parts[idx].sketch = Some(data);
            }
            "seedpart" => match self.kwarg_num(c.kwargs, "size") {
                Ok(Some(s)) if s > 0.0 => self.parts[idx].seed = Some(s),
                _ => self.warn(line, "seedPart without a positive size — ignored"),
            },
            "generatemesh" => self.parts[idx].meshed = true,
            "sectionassignment" => match self.kwarg_str(c.kwargs, "sectionName") {
                Ok(Some(s)) => {
                    self.parts[idx].section = Some(s);
                    if self.kwarg(c.kwargs, "region").is_some() {
                        self.warn_once(
                            "section-region",
                            line,
                            "SectionAssignment region ignored (v1: single body)",
                        );
                    }
                }
                _ => self.warn(line, "SectionAssignment without sectionName — ignored"),
            },
            "setelementtype" => {
                self.unsupp(line, "setElementType ignored (v1: CPS4 structured grid)");
            }
            _ => {
                self.unsupp(line, format!("part method '{method}' unsupported"));
            }
        }
        Ok(())
    }

    fn assembly_action(&mut self, c: &Call) -> Result<(), TranslateError> {
        let method = c.path.last().map(|s| s.as_str()).unwrap_or_default();
        let m = method.to_ascii_lowercase();
        match m.as_str() {
            "instance" => self.assembly_instance(c),
            "set" | "surface" => self.assembly_set_or_surface(c),
            "regenerate" => {}
            _ => self.unsupp(c.line, format!("assembly method '{method}' unsupported")),
        }
        Ok(())
    }

    fn assembly_instance(&mut self, c: &Call) {
        let part = self.kwarg(c.kwargs, "part").and_then(|e| match e {
            Expr::Name(n) => Some(n.clone()),
            _ => None,
        });
        let Some(pvar) = part else {
            self.warn(c.line, "Instance without part=<Part> — ignored");
            return;
        };
        match self.objs.get(&pvar).map(clone_obj) {
            Some(Obj::Part(_)) => {}
            _ => self.warn(c.line, format!("Instance part '{pvar}' is not a Part")),
        }
    }

    fn assembly_set_or_surface(&mut self, c: &Call) {
        let line = c.line;
        let Some(name) = self.kwarg_str(c.kwargs, "name").unwrap_or_default() else {
            self.warn(line, "Set/Surface without a name — ignored");
            return;
        };
        for key in ["vertices", "edges", "side1Edges", "side2Edges", "nodes"] {
            if let Some(e) = self.kwarg(c.kwargs, key) {
                match self.region_of(e, line) {
                    Some(r) => self.sets.push((name.clone(), r)),
                    None => return,
                }
                return;
            }
        }
        self.warn(
            line,
            format!("Set/Surface '{name}' without a region keyword — ignored"),
        );
    }

    // -- finalization ----------------------------------------------------

    fn finish(mut self) -> Result<Translation, TranslateError> {
        let Some(idx) = self.parts.iter().position(|p| p.sketch.is_some()) else {
            return terr(
                self.parts.first().map(|p| p.line).unwrap_or(0),
                "no TWO_D_PLANAR part with a BaseShell/BaseSolid2D sketch found",
            );
        };
        let part_line = self.parts[idx].line;
        let part_name = self.parts[idx].name.clone();
        if !self.parts[idx].meshed {
            self.report
                .warnings
                .push("generateMesh() not called — mesh generated anyway".into());
        }
        let sk = self.parts[idx].sketch.clone().unwrap_or_default();

        let Some(p0) = sk.segments.first().map(|s| s.0) else {
            return terr(
                part_line,
                "sketch has no rectangle/Line geometry to bound a mesh",
            );
        };
        let (mut x0, mut y0, mut x1, mut y1) = (p0[0], p0[1], p0[0], p0[1]);
        for (a, b) in &sk.segments {
            for p in [*a, *b] {
                x0 = x0.min(p[0]);
                x1 = x1.max(p[0]);
                y0 = y0.min(p[1]);
                y1 = y1.max(p[1]);
            }
        }
        if (x1 - x0).abs() < 1e-12 || (y1 - y0).abs() < 1e-12 {
            return terr(
                part_line,
                "sketch bounding box is degenerate (zero width/height)",
            );
        }

        let mut holes = Vec::new();
        for &(cx, cy, r) in &sk.circles {
            if cx - r < x0 - 1e-9 || cx + r > x1 + 1e-9 || cy - r < y0 - 1e-9 || cy + r > y1 + 1e-9
            {
                self.report.warnings.push(format!(
                    "hole at ({cx}, {cy}) r={r} extends outside the rectangle — clipped"
                ));
            }
            holes.push(Hole { cx, cy, r });
        }

        let (w, h) = (x1 - x0, y1 - y0);
        let seed = match self.parts[idx].seed {
            Some(s) => s,
            None => {
                self.report
                    .warnings
                    .push("seedPart missing — defaulting to max-dimension/8".into());
                w.max(h) / 8.0
            }
        };
        let nx = mesh_count(w, seed);
        let ny = mesh_count(h, seed);

        let (mesh, notes) = meshgen::grid(x0, y0, w, h, nx, ny, &holes);
        self.report.warnings.extend(notes);

        let mut materials: Vec<Material> = Vec::new();
        for m in &self.mats {
            if let Some((e, nu)) = m.elastic {
                materials.push(Material {
                    name: m.name.clone(),
                    youngs_modulus: e,
                    poisson_ratio: nu,
                });
            }
        }
        if materials.is_empty() {
            return terr(
                part_line,
                "no Material with an Elastic table — nothing to solve",
            );
        }
        let density = self.mats.iter().find_map(|m| m.density);

        let (mat_name, thickness) = match &self.parts[idx].section {
            Some(sname) => match self.sections.iter().find(|(n, _, _)| n == sname) {
                Some((_, mat, th)) => {
                    if !materials.iter().any(|m| &m.name == mat) {
                        return terr(
                            part_line,
                            format!("section '{sname}' references material '{mat}' without an Elastic table"),
                        );
                    }
                    (mat.clone(), th.unwrap_or(1.0))
                }
                None => {
                    return terr(
                        part_line,
                        format!("SectionAssignment references unknown section '{sname}'"),
                    )
                }
            },
            None => {
                self.report.warnings.push(
                    "no SectionAssignment — first material applied with thickness 1.0".into(),
                );
                (materials[0].name.clone(), 1.0)
            }
        };
        if self.sections.len() > 1 || self.mats.len() > 1 {
            self.report.warnings.push(format!(
                "multiple materials/sections defined — v1 uses '{mat_name}' for the whole part '{part_name}'"
            ));
        }

        let mut model = Model {
            mesh,
            thickness,
            materials,
            elsets: vec![ElementSet {
                name: "ALL".into(),
                element_ids: Vec::new(),
            }],
            mat_assignments: vec![("ALL".into(), mat_name)],
            ..Default::default()
        };
        model.elsets[0].element_ids = model.mesh.elements.iter().map(|e| e.id).collect();
        if let Some(d) = density {
            model.density = d;
        }

        let geo = Geo {
            x0,
            y0,
            x1,
            y1,
            // Boundary classification must not swallow a whole cell.
            tol: 1e-7 * w.max(h).abs().max(1.0),
            cell: (w / nx as f64).max(h / ny as f64),
        };

        let steps = std::mem::take(&mut self.steps);
        let mut steps_out = Vec::with_capacity(steps.len());
        let mut carried_bcs = self.initial_bcs.clone();
        let mut carried_cloads: Vec<(Region, usize, f64)> = Vec::new();
        let mut carried_press: Vec<(Region, f64)> = Vec::new();
        let mut gravity = None;
        for st in &steps {
            carried_bcs.extend(st.bcs.iter().cloned());
            carried_cloads.extend(st.cloads.iter().cloned());
            carried_press.extend(st.pressures.iter().cloned());
            if st.gravity.is_some() {
                gravity = st.gravity;
            }
            let mut bcs = Vec::new();
            for pb in &carried_bcs {
                if let Some(nodes) = self.resolve_nodes(&pb.region, &geo, &model.mesh) {
                    if !nodes.is_empty() {
                        bcs.push(BoundaryCondition {
                            node_ids: nodes,
                            dofs: pb.dofs.clone(),
                        });
                    }
                }
            }
            let mut cloads = Vec::new();
            for (r, dof, val) in &carried_cloads {
                if let Some(nodes) = self.resolve_nodes(r, &geo, &model.mesh) {
                    for nid in nodes {
                        cloads.push(ConcentratedLoad {
                            node_id: nid,
                            dof: *dof,
                            value: *val,
                        });
                    }
                }
            }
            let mut dloads = Vec::new();
            for (r, mag) in &carried_press {
                dloads.extend(self.pressure_faces(r, *mag, &geo, &model.mesh));
            }
            steps_out.push(Step {
                name: st.name.clone(),
                boundary_conditions: bcs,
                cloads,
                dloads,
                gravity,
                dynamics: st.dynamics.clone(),
                ..Default::default()
            });
        }
        model.steps = steps_out;

        Ok(Translation {
            model,
            report: self.report,
        })
    }

    /// Resolve a region against the generated mesh.
    fn resolve_nodes(
        &mut self,
        r: &Region,
        geo: &Geo,
        mesh: &solvix_core::Mesh,
    ) -> Option<Vec<u32>> {
        let mut ids: Vec<u32> = match r {
            Region::Coords(pts) => {
                let mut v = Vec::new();
                for p in pts {
                    match nearest_node(mesh, *p) {
                        Some((id, d)) if d <= 2.0 * geo.cell => v.push(id),
                        _ => self.report.warnings.push(format!(
                            "no mesh node near ({}, {}) — point skipped",
                            p[0], p[1]
                        )),
                    }
                }
                v
            }
            Region::EdgeMids(mids) => {
                let mut v = Vec::new();
                for m in mids {
                    match geo.classify(*m) {
                        Some(side) => v.extend(nodes_on_side(mesh, geo, side)),
                        None => self.report.unsupported.push(format!(
                            "edge at ({}, {}) is not on the rectangle boundary — skipped",
                            m[0], m[1]
                        )),
                    }
                }
                v
            }
            Region::Named(n) => {
                let Some((_, r)) = self.sets.iter().find(|(name, _)| name == n) else {
                    self.report
                        .unsupported
                        .push(format!("unknown set/surface '{n}'"));
                    return None;
                };
                let r = r.clone();
                return self.resolve_nodes(&r, geo, mesh);
            }
            Region::AllNodes => {
                self.report
                    .warnings
                    .push("geometry region covers the whole part".into());
                mesh.nodes.iter().map(|n| n.id).collect()
            }
        };
        if ids.is_empty() {
            return Some(ids);
        }
        ids.sort_unstable();
        ids.dedup();
        Some(ids)
    }

    /// Best-effort edge pressure -> element P-face loads.
    fn pressure_faces(
        &mut self,
        r: &Region,
        mag: f64,
        geo: &Geo,
        mesh: &solvix_core::Mesh,
    ) -> Vec<DistributedLoad> {
        let mids: Vec<[f64; 2]> = match r {
            Region::EdgeMids(m) => m.clone(),
            Region::Named(n) => {
                let Some((_, inner)) = self.sets.iter().find(|(name, _)| name == n) else {
                    self.report
                        .unsupported
                        .push(format!("pressure references unknown surface '{n}'"));
                    return vec![];
                };
                match inner {
                    Region::EdgeMids(m) => m.clone(),
                    _ => {
                        self.report.unsupported.push(format!(
                            "pressure surface '{n}' is not edge-based — skipped"
                        ));
                        return vec![];
                    }
                }
            }
            _ => {
                self.report
                    .unsupported
                    .push("pressure region is not edge-based — skipped".into());
                return vec![];
            }
        };
        let mut sides = Vec::new();
        for m in &mids {
            match geo.classify(*m) {
                Some(s) => sides.push(s),
                None => self.report.unsupported.push(format!(
                    "pressure edge at ({}, {}) is not on the rectangle boundary — skipped",
                    m[0], m[1]
                )),
            }
        }
        let mut out = Vec::new();
        for e in &mesh.elements {
            for i in 0..4 {
                let (a, b) = (e.connectivity[i], e.connectivity[(i + 1) % 4]);
                let (na, nb) = match (mesh.node_index(a), mesh.node_index(b)) {
                    (Some(ia), Some(ib)) => (&mesh.nodes[ia], &mesh.nodes[ib]),
                    _ => continue,
                };
                for side in &sides {
                    if side.contains(na.coords, geo) && side.contains(nb.coords, geo) {
                        out.push(DistributedLoad {
                            element_id: e.id,
                            face: i + 1,
                            magnitude: mag,
                        });
                    }
                }
            }
        }
        if !out.is_empty() {
            self.report
                .warnings
                .push("Pressure: edge converted to element face loads (best effort)".into());
        } else {
            self.report
                .warnings
                .push("Pressure: no element faces found on the pressure edge — skipped".into());
        }
        out
    }
}

fn clone_obj(o: &Obj) -> Obj {
    match o {
        Obj::Sketch(s) => Obj::Sketch(s.clone()),
        Obj::Part(i) => Obj::Part(*i),
        Obj::Assembly => Obj::Assembly,
        Obj::Region(r) => Obj::Region(r.clone()),
        Obj::Opaque => Obj::Opaque,
    }
}

fn mesh_count(dim: f64, seed: f64) -> usize {
    let n = (dim / seed).round().clamp(1.0, 400.0);
    n as usize
}

fn nearest_node(mesh: &solvix_core::Mesh, p: [f64; 2]) -> Option<(u32, f64)> {
    mesh.nodes
        .iter()
        .map(|n| {
            (
                n.id,
                ((n.coords[0] - p[0]).powi(2) + (n.coords[1] - p[1]).powi(2)).sqrt(),
            )
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
}

struct Geo {
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
    /// Tight tolerance for on-boundary classification.
    tol: f64,
    /// Structured cell size (for nearest-node acceptance).
    cell: f64,
}

#[derive(Clone, Copy)]
enum Side {
    Left,
    Right,
    Bottom,
    Top,
}

impl Side {
    fn contains(&self, c: [f64; 3], geo: &Geo) -> bool {
        match self {
            Side::Left => (c[0] - geo.x0).abs() <= geo.tol,
            Side::Right => (c[0] - geo.x1).abs() <= geo.tol,
            Side::Bottom => (c[1] - geo.y0).abs() <= geo.tol,
            Side::Top => (c[1] - geo.y1).abs() <= geo.tol,
        }
    }
}

impl Geo {
    fn classify(&self, p: [f64; 2]) -> Option<Side> {
        if (p[0] - self.x0).abs() <= self.tol {
            Some(Side::Left)
        } else if (p[0] - self.x1).abs() <= self.tol {
            Some(Side::Right)
        } else if (p[1] - self.y0).abs() <= self.tol {
            Some(Side::Bottom)
        } else if (p[1] - self.y1).abs() <= self.tol {
            Some(Side::Top)
        } else {
            None
        }
    }
}

fn nodes_on_side(mesh: &solvix_core::Mesh, geo: &Geo, side: Side) -> Vec<u32> {
    mesh.nodes
        .iter()
        .filter(|n| side.contains(n.coords, geo))
        .map(|n| n.id)
        .collect()
}
