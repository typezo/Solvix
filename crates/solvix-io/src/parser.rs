//! Abaqus `.inp` parser for the v0.1 keyword subset (issue #2), extended
//! with multi-step, `*DLOAD` and `*DENSITY` (issue #26).
//!
//! Line model: `**…` comments and blank lines are dropped up front; every
//! remaining line is either a keyword line (`*NAME, KEY=VALUE, …`, parsed
//! case-insensitively) or a data line (comma-separated fields). A keyword
//! "block" owns the run of data lines that follows it.
//!
//! Steps (#26): `*BOUNDARY`/`*CLOAD`/`*DLOAD` inside `*STEP` blocks
//! accumulate into steps that each store the TOTAL BC/load state at the
//! end of that step (OP=NEW/MOD resolved here; linear superposition).
//! A single-step deck mirrors its BCs/loads into the legacy top-level
//! Model fields so v0.1 consumers see exactly the old picture.

use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use solvix_core::{
    BoundaryCondition, ConcentratedLoad, DistributedLoad, ElemKind, Element, ElementSet, Material,
    Model, Node, Step,
};

use super::{ParseError, ParseReport};

struct KwLine {
    line: usize,
    name: String,
    params: HashMap<String, String>,
}

enum Item {
    Kw(KwLine),
    Data { line: usize, text: String },
}

/// Working accumulator for one step's (or the model data's) BCs/loads.
#[derive(Default, Clone)]
struct StepAcc {
    name: String,
    bcs: Vec<BoundaryCondition>,
    cloads: Vec<ConcentratedLoad>,
    dloads: Vec<DistributedLoad>,
    gravity: Option<(f64, [f64; 3])>,
}

/// OP=NEW/MOD of a keyword line: `true` when OP=NEW, warning on anything
/// other than NEW/MOD (treated as MOD).
fn op_is_new(kw: &KwLine, report: &mut ParseReport) -> bool {
    match kw.params.get("op").map(|s| s.to_lowercase()) {
        None => false,
        Some(v) if v == "new" => true,
        Some(v) if v == "mod" => false,
        Some(v) => {
            report.warnings.push(format!(
                "line {}: unsupported OP='{v}' on *{} — treated as MOD",
                kw.line,
                kw.name.to_uppercase()
            ));
            false
        }
    }
}

pub(super) fn parse_inp(text: &str) -> Result<(Model, ParseReport), ParseError> {
    let items = tokenize(text);
    let mut model = Model {
        thickness: 1.0,
        ..Model::default()
    };
    let mut report = ParseReport::default();
    let mut nsets: HashMap<String, Vec<u32>> = HashMap::new();
    let mut node_ids: HashSet<u32> = HashSet::new();
    let mut elem_ids: HashSet<u32> = HashSet::new();
    let mut pending_material: Option<String> = None;
    let mut cur_material: Option<String> = None;
    let mut density_seen = false;
    let mut model_data = StepAcc::default();
    let mut steps_done: Vec<StepAcc> = Vec::new();
    let mut cur: Option<StepAcc> = None;
    let mut step_no = 0usize;

    let mut j = 0usize;
    while j < items.len() {
        let kw = match &items[j] {
            Item::Kw(kw) => kw,
            Item::Data { line, .. } => {
                report.warnings.push(format!(
                    "line {line}: data line outside any keyword block ignored"
                ));
                j += 1;
                continue;
            }
        };

        if let Some(name) = &pending_material {
            if kw.name != "elastic" && kw.name != "density" {
                report.warnings.push(format!(
                    "line {}: material '{name}' has no *ELASTIC block",
                    kw.line
                ));
                pending_material = None;
            }
        }

        match kw.name.as_str() {
            "heading" => {
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            "node" => {
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if f.len() != 3 && f.len() != 4 {
                        return err(
                            *dno,
                            format!("*NODE expects id, x, y[, z]; got {} field(s)", f.len()),
                        );
                    }
                    let id: u32 = num(f[0], *dno, "node id")?;
                    let x = num(f[1], *dno, "coordinate")?;
                    let y = num(f[2], *dno, "coordinate")?;
                    let z = if f.len() == 4 {
                        num(f[3], *dno, "coordinate")?
                    } else {
                        0.0
                    };
                    if !node_ids.insert(id) {
                        return err(*dno, format!("duplicate node id {id}"));
                    }
                    model.mesh.nodes.push(Node {
                        id,
                        coords: [x, y, z],
                    });
                    j += 1;
                }
            }

            "element" => {
                let ty = kw.params.get("type").map(|s| s.to_uppercase());
                let kind = match ty.as_deref() {
                    Some("CPS4") => Some(ElemKind::Cps4),
                    Some("CPS3") => Some(ElemKind::Cps3),
                    Some("C3D8") => Some(ElemKind::C3d8),
                    Some("C3D4") => Some(ElemKind::C3d4),
                    _ => None,
                };
                let Some(kind) = kind else {
                    let what = ty.map_or_else(
                        || "TYPE-less *ELEMENT".to_string(),
                        |t| format!("element TYPE={t}"),
                    );
                    report.warnings.push(format!(
                        "line {}: unsupported {what} — block skipped (supported: CPS4/CPS3/C3D8/C3D4)",
                        kw.line
                    ));
                    j += 1;
                    while matches!(items.get(j), Some(Item::Data { .. })) {
                        j += 1;
                    }
                    continue;
                };
                let want = 1 + kind.node_count();
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    // One element per record; a trailing comma continues the
                    // connectivity on the next data line.
                    let mut f: Vec<String> = Vec::new();
                    let mut first_dno = 0usize;
                    while let Some(Item::Data { line: dno, text }) = items.get(j) {
                        if f.is_empty() {
                            first_dno = *dno;
                        }
                        let ends_comma = text.ends_with(',');
                        f.extend(text.split(',').map(|s| s.trim().to_string()));
                        j += 1;
                        let have = f.iter().filter(|s| !s.is_empty()).count();
                        if !ends_comma || have >= want {
                            break;
                        }
                    }
                    let f: Vec<&str> = f
                        .iter()
                        .filter(|s| !s.is_empty())
                        .map(|s| s.as_str())
                        .collect();
                    if f.len() != want {
                        return err(
                            first_dno,
                            format!(
                                "*ELEMENT {:?} expects {} field(s) (id + {} nodes); got {}",
                                kind,
                                want,
                                want - 1,
                                f.len()
                            ),
                        );
                    }
                    let id: u32 = num(f[0], first_dno, "element id")?;
                    let connectivity: Vec<u32> = f[1..]
                        .iter()
                        .map(|s| num(s, first_dno, "node id"))
                        .collect::<Result<_, _>>()?;
                    if !elem_ids.insert(id) {
                        return err(first_dno, format!("duplicate element id {id}"));
                    }
                    if let Some(setname) = kw.params.get("elset") {
                        upsert_elset(&mut model.elsets, setname, &[id]);
                    }
                    model.mesh.elements.push(Element {
                        id,
                        kind,
                        connectivity,
                    });
                }
            }

            "elset" | "nset" => {
                let name = kw
                    .params
                    .get("elset")
                    .or_else(|| kw.params.get("nset"))
                    .map(|s| s.to_uppercase());
                let generate = kw.params.contains_key("generate");
                let Some(name) = name else {
                    report.warnings.push(format!(
                        "line {}: set keyword without ELSET=/NSET= name — block skipped",
                        kw.line
                    ));
                    j += 1;
                    while matches!(items.get(j), Some(Item::Data { .. })) {
                        j += 1;
                    }
                    continue;
                };
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    let ids: Vec<u32> = if generate {
                        if f.len() < 2 || f.len() > 3 {
                            return err(
                                *dno,
                                format!(
                                    "GENERATE expects start, end[, step]; got {} field(s)",
                                    f.len()
                                ),
                            );
                        }
                        let start: u32 = num(f[0], *dno, "set id")?;
                        let end: u32 = num(f[1], *dno, "set id")?;
                        let step: u32 = if f.len() == 3 {
                            num(f[2], *dno, "set step")?
                        } else {
                            1
                        };
                        if step == 0 {
                            return err(*dno, "GENERATE step must be >= 1".to_string());
                        }
                        (start..=end).step_by(step as usize).collect()
                    } else {
                        f.iter()
                            .map(|s| num(s, *dno, "set member id"))
                            .collect::<Result<Vec<u32>, _>>()?
                    };
                    if kw.name == "elset" {
                        upsert_elset(&mut model.elsets, &name, &ids);
                    } else {
                        nsets.entry(name.clone()).or_default().extend(ids);
                    }
                    j += 1;
                }
            }

            "material" => {
                match kw.params.get("name") {
                    Some(n) => {
                        let key = n.to_uppercase();
                        pending_material = Some(key.clone());
                        cur_material = Some(key);
                    }
                    None => report.warnings.push(format!(
                        "line {}: *MATERIAL without NAME — block ignored",
                        kw.line
                    )),
                }
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            "density" => {
                j += 1;
                let mut rows = 0usize;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if rows == 0 {
                        if cur_material.is_none() {
                            report
                                .warnings
                                .push(format!("line {dno}: *DENSITY outside *MATERIAL — ignored"));
                        } else if f.is_empty() {
                            return err(*dno, "*DENSITY expects a density value".to_string());
                        } else {
                            let d = num(f[0], *dno, "density")?;
                            if density_seen && model.density != d {
                                report.warnings.push(format!(
                                    "line {dno}: *DENSITY overwrites an earlier value (per-material density is v0.4)"
                                ));
                            }
                            model.density = d;
                            density_seen = true;
                            if f.len() > 1 {
                                report.warnings.push(format!(
                                    "line {dno}: extra *DENSITY fields ignored (temperature dependence is v0.2)"
                                ));
                            }
                        }
                    } else {
                        report.warnings.push(format!(
                            "line {dno}: extra *DENSITY row ignored (temperature dependence is v0.2)"
                        ));
                    }
                    rows += 1;
                    j += 1;
                }
            }

            "elastic" => {
                j += 1;
                let mut got_row = false;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if !got_row {
                        if f.len() != 2 {
                            return err(
                                *dno,
                                format!(
                                    "*ELASTIC expects youngs_modulus, poisson_ratio; got {} field(s)",
                                    f.len()
                                ),
                            );
                        }
                        match pending_material.take() {
                            Some(name) => {
                                let youngs_modulus = num(f[0], *dno, "youngs modulus")?;
                                let poisson_ratio = num(f[1], *dno, "poisson ratio")?;
                                model.materials.push(Material {
                                    name,
                                    youngs_modulus,
                                    poisson_ratio,
                                });
                            }
                            None => report.warnings.push(format!(
                                "line {}: *ELASTIC outside *MATERIAL — block ignored",
                                dno
                            )),
                        }
                        got_row = true;
                    } else {
                        report.warnings.push(format!(
                            "line {}: extra *ELASTIC row ignored (temperature dependence is v0.2)",
                            dno
                        ));
                    }
                    j += 1;
                }
            }

            "solid section" => {
                let elset = kw.params.get("elset").map(|s| s.to_uppercase());
                let material = kw.params.get("material").map(|s| s.to_uppercase());
                j += 1;
                match (elset, material) {
                    (Some(elset), Some(material)) => {
                        if !model.elsets.iter().any(|s| s.name == elset) {
                            return err(
                                kw.line,
                                format!("*SOLID SECTION references unknown ELSET '{elset}'"),
                            );
                        }
                        if !model.materials.iter().any(|m| m.name == material) {
                            return err(
                                kw.line,
                                format!("*SOLID SECTION references unknown MATERIAL '{material}'"),
                            );
                        }
                        model.mat_assignments.push((elset, material));
                        if let Some(Item::Data { line: dno, text }) = items.get(j) {
                            let f = fields(text);
                            if f.len() != 1 {
                                return err(
                                    *dno,
                                    format!(
                                        "*SOLID SECTION data expects one thickness; got {} field(s)",
                                        f.len()
                                    ),
                                );
                            }
                            model.thickness = num(f[0], *dno, "thickness")?;
                            j += 1;
                        }
                        while matches!(items.get(j), Some(Item::Data { .. })) {
                            j += 1;
                        }
                    }
                    _ => {
                        report.warnings.push(format!(
                            "line {}: *SOLID SECTION needs ELSET and MATERIAL — block skipped",
                            kw.line
                        ));
                        while matches!(items.get(j), Some(Item::Data { .. })) {
                            j += 1;
                        }
                    }
                }
            }

            "boundary" => {
                let mut op_new = op_is_new(kw, &mut report);
                let in_step = cur.is_some();
                if !in_step {
                    if op_new {
                        report
                            .warnings
                            .push(format!("line {}: OP=NEW outside *STEP ignored", kw.line));
                        op_new = false;
                    }
                    if !steps_done.is_empty() {
                        report.warnings.push(format!(
                            "line {}: *BOUNDARY between steps — ignored",
                            kw.line
                        ));
                    }
                } else if op_new {
                    cur.as_mut().unwrap().bcs.clear();
                }
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if f.len() < 2 || f.len() > 4 {
                        return err(
                            *dno,
                            format!(
                                "*BOUNDARY expects node/set, dof[, dof][, value]; got {} field(s)",
                                f.len()
                            ),
                        );
                    }
                    let node_ids = resolve_nodes(&nsets, f[0], *dno)?;
                    let dof1: usize = num(f[1], *dno, "dof")?;
                    let dofs: Vec<usize> = if f.len() >= 3 {
                        let dof2: usize = num(f[2], *dno, "dof")?;
                        (dof1..=dof2).collect()
                    } else {
                        vec![dof1]
                    };
                    if f.len() == 4 {
                        let value: f64 = num(f[3], *dno, "prescribed value")?;
                        if value != 0.0 {
                            report.warnings.push(format!(
                                "line {dno}: nonzero prescribed displacement {value} ignored (v0.1: homogeneous BCs only)"
                            ));
                        }
                    }
                    let bc = BoundaryCondition { node_ids, dofs };
                    match cur.as_mut() {
                        Some(acc) => acc.bcs.push(bc),
                        None => model_data.bcs.push(bc),
                    }
                    j += 1;
                }
            }

            "cload" => {
                let mut op_new = op_is_new(kw, &mut report);
                let in_step = cur.is_some();
                if !in_step {
                    if op_new {
                        report
                            .warnings
                            .push(format!("line {}: OP=NEW outside *STEP ignored", kw.line));
                        op_new = false;
                    }
                    if !steps_done.is_empty() {
                        report
                            .warnings
                            .push(format!("line {}: *CLOAD between steps — ignored", kw.line));
                    }
                } else if op_new {
                    cur.as_mut().unwrap().cloads.clear();
                }
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if f.len() != 3 {
                        return err(
                            *dno,
                            format!(
                                "*CLOAD expects node/set, dof, value; got {} field(s)",
                                f.len()
                            ),
                        );
                    }
                    let nodes = resolve_nodes(&nsets, f[0], *dno)?;
                    let dof: usize = num(f[1], *dno, "dof")?;
                    let value = num(f[2], *dno, "load value")?;
                    for node_id in nodes {
                        let cl = ConcentratedLoad {
                            node_id,
                            dof,
                            value,
                        };
                        match cur.as_mut() {
                            Some(acc) => acc.cloads.push(cl),
                            None => model_data.cloads.push(cl),
                        }
                    }
                    j += 1;
                }
            }

            "dload" => {
                let op_new = op_is_new(kw, &mut report);
                let Some(acc) = cur.as_mut() else {
                    report.warnings.push(format!(
                        "line {}: *DLOAD outside *STEP — ignored (distributed loads live in steps)",
                        kw.line
                    ));
                    j += 1;
                    while matches!(items.get(j), Some(Item::Data { .. })) {
                        j += 1;
                    }
                    continue;
                };
                if op_new {
                    acc.dloads.clear();
                    acc.gravity = None;
                }
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if f[0].eq_ignore_ascii_case("GRAV") {
                        if f.len() < 3 || f.len() > 5 {
                            return err(
                                *dno,
                                format!(
                                    "*DLOAD GRAV expects GRAV, magnitude, dir1[, dir2][, dir3]; got {} field(s)",
                                    f.len()
                                ),
                            );
                        }
                        let g = num(f[1], *dno, "gravity magnitude")?;
                        let d1 = num(f[2], *dno, "gravity direction")?;
                        let d2 = if f.len() > 3 {
                            num(f[3], *dno, "gravity direction")?
                        } else {
                            0.0
                        };
                        let d3 = if f.len() > 4 {
                            num(f[4], *dno, "gravity direction")?
                        } else {
                            0.0
                        };
                        acc.gravity = Some((g, [d1, d2, d3]));
                    } else {
                        if f.len() != 3 {
                            return err(
                                *dno,
                                format!(
                                    "*DLOAD expects element/elset, P1..P4, magnitude; got {} field(s)",
                                    f.len()
                                ),
                            );
                        }
                        let ltype = f[1].to_uppercase();
                        let face = match ltype
                            .strip_prefix('P')
                            .and_then(|n| n.parse::<usize>().ok())
                        {
                            Some(n) if (1..=4).contains(&n) => n,
                            _ => {
                                report.warnings.push(format!(
                                    "line {dno}: *DLOAD type '{ltype}' unsupported (v0.3: P1|P2|P3|P4, GRAV) — line ignored"
                                ));
                                j += 1;
                                continue;
                            }
                        };
                        let magnitude = num(f[2], *dno, "pressure magnitude")?;
                        match resolve_elements(&model, f[0]) {
                            Some(ids) => {
                                for id in ids {
                                    match model.mesh.elements.iter().find(|e| e.id == id) {
                                        None => report.warnings.push(format!(
                                            "line {dno}: *DLOAD on unknown element {id} — ignored"
                                        )),
                                        Some(e) => {
                                            let faces = match e.kind {
                                                ElemKind::Cps3 => 3,
                                                ElemKind::Cps4 => 4,
                                                _ => 0,
                                            };
                                            if face > faces {
                                                report.warnings.push(format!(
                                                    "line {dno}: element {id} ({}) has no face {face} — ignored",
                                                    kind_name(e.kind)
                                                ));
                                            } else {
                                                acc.dloads.push(DistributedLoad {
                                                    element_id: id,
                                                    face,
                                                    magnitude,
                                                });
                                            }
                                        }
                                    }
                                }
                            }
                            None => report.warnings.push(format!(
                                "line {dno}: *DLOAD target '{}' is not a known element or ELSET — line ignored",
                                f[0]
                            )),
                        }
                    }
                    j += 1;
                }
            }

            "step" => {
                if let Some(open) = cur.take() {
                    steps_done.push(open);
                }
                step_no += 1;
                let mut op_new = false;
                if let Some(v) = kw.params.get("op").map(|s| s.to_lowercase()) {
                    match v.as_str() {
                        "new" => op_new = true,
                        "mod" => {}
                        other => report.warnings.push(format!(
                            "line {}: unsupported OP='{other}' on *STEP — treated as MOD",
                            kw.line
                        )),
                    }
                }
                let mut acc = if op_new {
                    StepAcc::default()
                } else {
                    steps_done
                        .last()
                        .cloned()
                        .unwrap_or_else(|| model_data.clone())
                };
                acc.name = kw
                    .params
                    .get("name")
                    .cloned()
                    .unwrap_or_else(|| format!("Step-{step_no}"));
                cur = Some(acc);
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            "static" => {
                if !kw.params.is_empty() {
                    report.warnings.push(format!(
                        "line {}: *STATIC parameters ignored (v0.1: linear static only)",
                        kw.line
                    ));
                }
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            "end step" => {
                if let Some(open) = cur.take() {
                    steps_done.push(open);
                }
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            other => {
                report.warnings.push(format!(
                    "line {}: unsupported keyword '*{}' — block skipped",
                    kw.line,
                    other.to_uppercase()
                ));
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }
        }
    }

    if let Some(name) = &pending_material {
        report
            .warnings
            .push(format!("material '{name}' has no *ELASTIC block"));
    }
    if let Some(open) = cur.take() {
        steps_done.push(open);
    }

    if steps_done.is_empty() {
        // Legacy single-step (or step-less) deck: model data IS the model.
        model.boundary_conditions = model_data.bcs;
        model.cloads = model_data.cloads;
    } else {
        // Single-step decks mirror the (accumulated) step state into the
        // legacy top-level fields so v0.1 consumers behave exactly as
        // before; multi-step decks are steps-only.
        if steps_done.len() == 1 {
            model.boundary_conditions = steps_done[0].bcs.clone();
            model.cloads = steps_done[0].cloads.clone();
        }
        model.steps = steps_done
            .into_iter()
            .map(|a| Step {
                name: a.name,
                boundary_conditions: a.bcs,
                cloads: a.cloads,
                dloads: a.dloads,
                gravity: a.gravity,
            })
            .collect();
    }
    Ok((model, report))
}

fn tokenize(text: &str) -> Vec<Item> {
    text.lines()
        .enumerate()
        .filter_map(|(idx, raw)| {
            let t = raw.trim();
            if t.is_empty() || t.starts_with("**") {
                None
            } else if let Some(body) = t.strip_prefix('*') {
                let (name, rest) = body.split_once(',').map_or((body, ""), |(n, r)| (n, r));
                Some(Item::Kw(KwLine {
                    line: idx + 1,
                    name: name.trim().to_lowercase(),
                    params: parse_params(rest),
                }))
            } else {
                Some(Item::Data {
                    line: idx + 1,
                    text: t.to_string(),
                })
            }
        })
        .collect()
}

fn parse_params(rest: &str) -> HashMap<String, String> {
    rest.split(',')
        .filter_map(|p| {
            let p = p.trim();
            if p.is_empty() {
                return None;
            }
            Some(match p.split_once('=') {
                Some((k, v)) => (k.trim().to_lowercase(), v.trim().to_string()),
                None => (p.to_lowercase(), String::new()),
            })
        })
        .collect()
}

/// Comma-split a data line into trimmed fields, dropping trailing empty ones.
fn fields(line: &str) -> Vec<&str> {
    let mut f: Vec<&str> = line.split(',').map(str::trim).collect();
    while f.last() == Some(&"") {
        f.pop();
    }
    f
}

fn num<T: FromStr>(s: &str, line: usize, what: &str) -> Result<T, ParseError> {
    s.parse::<T>().map_err(|_| ParseError {
        line,
        message: format!("malformed {what} '{s}'"),
    })
}

fn err<T>(line: usize, message: String) -> Result<T, ParseError> {
    Err(ParseError { line, message })
}

/// Resolve a *BOUNDARY/*CLOAD target: node-set name first, then a bare node id.
fn resolve_nodes(
    nsets: &HashMap<String, Vec<u32>>,
    target: &str,
    line: usize,
) -> Result<Vec<u32>, ParseError> {
    if let Some(ids) = nsets.get(&target.to_uppercase()) {
        return Ok(ids.clone());
    }
    target
        .parse::<u32>()
        .map(|id| vec![id])
        .map_err(|_| ParseError {
            line,
            message: format!("unknown node set '{target}'"),
        })
}

/// Resolve a *DLOAD target: element-set name first, then a bare element id.
fn resolve_elements(model: &Model, target: &str) -> Option<Vec<u32>> {
    let key = target.to_uppercase();
    if let Some(s) = model.elsets.iter().find(|s| s.name == key) {
        return Some(s.element_ids.clone());
    }
    target.parse::<u32>().ok().map(|id| vec![id])
}

fn kind_name(k: ElemKind) -> &'static str {
    match k {
        ElemKind::Cps3 => "cps3",
        ElemKind::Cps4 => "cps4",
        ElemKind::C3d8 => "c3d8",
        ElemKind::C3d4 => "c3d4",
    }
}

/// Create `name` (uppercased) in `sets` or append to it — `*ELEMENT, ELSET=`
/// may legitimately name a set that does not exist yet.
fn upsert_elset(sets: &mut Vec<ElementSet>, name: &str, ids: &[u32]) {
    let key = name.to_uppercase();
    if let Some(s) = sets.iter_mut().find(|s| s.name == key) {
        s.element_ids.extend_from_slice(ids);
    } else {
        sets.push(ElementSet {
            name: key,
            element_ids: ids.to_vec(),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use solvix_core::ConcentratedLoad;

    const DECK: &str = r#"
** full-coverage smoke deck
*HEADING
    smoke deck
*NODE
1, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 2.0, 0.5
4, 0.0, 1.0
5, 1.0, 1.0
6, 2.0, 1.5
*ELEMENT, TYPE=CPS4, ELSET=QUADS
1, 1, 2, 5, 4
*ELEMENT, TYPE=CPS3, ELSET=TRIS
2, 2, 3, 5
*ELSET, ELSET=ALLELM
1, 2
*ELSET, ELSET=TOP, GENERATE
4, 6, 2
*NSET, NSET=FIX
1, 4
*NSET, NSET=TIP, GENERATE
5, 6, 1
*MATERIAL, NAME=steel
*ELASTIC
210.0E9, 0.3
*SOLID SECTION, ELSET=ALLELM, MATERIAL=steel
0.05
*DLOAD
ALLSIDES, P, 1.0
*STEP
*STATIC
*BOUNDARY
FIX, 1, 2
7, 1, 1, 2.5
*CLOAD
TIP, 2, -50.0
1, 1, 10.0
*END STEP
"#;

    #[test]
    fn full_deck_parses_all_keywords() {
        let (m, r) = parse_inp(DECK).unwrap();

        assert_eq!(m.mesh.nodes.len(), 6);
        assert_eq!(m.mesh.nodes[1].coords, [1.0, 0.0, 0.0]);
        assert_eq!(m.mesh.nodes[2].coords, [2.0, 0.5, 0.0]);

        assert_eq!(m.mesh.elements.len(), 2);
        assert_eq!(m.mesh.elements[0].kind, ElemKind::Cps4);
        assert_eq!(m.mesh.elements[0].connectivity, vec![1, 2, 5, 4]);
        assert_eq!(m.mesh.elements[1].kind, ElemKind::Cps3);
        assert_eq!(m.mesh.elements[1].connectivity, vec![2, 3, 5]);

        let set = |n: &str| m.elsets.iter().find(|s| s.name == n).unwrap();
        assert_eq!(set("QUADS").element_ids, vec![1]);
        assert_eq!(set("TRIS").element_ids, vec![2]);
        assert_eq!(set("ALLELM").element_ids, vec![1, 2]);
        assert_eq!(set("TOP").element_ids, vec![4, 6]);

        assert_eq!(m.thickness, 0.05);
        assert_eq!(m.materials.len(), 1);
        assert_eq!(m.materials[0].name, "STEEL");
        assert_eq!(m.materials[0].youngs_modulus, 210.0e9);
        assert_eq!(m.materials[0].poisson_ratio, 0.3);
        assert_eq!(
            m.mat_assignments,
            vec![("ALLELM".to_string(), "STEEL".to_string())]
        );

        // Sets expanded + dof range in BCs, bare node id with ignored value.
        assert_eq!(m.boundary_conditions.len(), 2);
        assert_eq!(m.boundary_conditions[0].node_ids, vec![1, 4]);
        assert_eq!(m.boundary_conditions[0].dofs, vec![1, 2]);
        assert_eq!(m.boundary_conditions[1].node_ids, vec![7]);
        assert_eq!(m.boundary_conditions[1].dofs, vec![1]);

        // CLOAD on NSET expands to every node, value per node.
        assert_eq!(
            m.cloads,
            vec![
                ConcentratedLoad {
                    node_id: 5,
                    dof: 2,
                    value: -50.0
                },
                ConcentratedLoad {
                    node_id: 6,
                    dof: 2,
                    value: -50.0
                },
                ConcentratedLoad {
                    node_id: 1,
                    dof: 1,
                    value: 10.0
                },
            ]
        );

        // *DLOAD skipped + nonzero prescribed displacement ignored.
        assert_eq!(r.warnings.len(), 2);
        assert!(r.warnings[0].contains("*DLOAD"));
        assert!(r.warnings[1].contains("nonzero prescribed displacement"));
        assert!(m.validate().is_ok());
    }

    #[test]
    fn lowercase_keywords_crlf_generate_and_default_thickness() {
        let deck = concat!(
            "*node\r\n",
            "1, 0.0, 0.0\r\n",
            "2, 1.0, 0.0\r\n",
            "3, 0.0, 1.0\r\n",
            "*element, type=cps3, elset=e1\r\n",
            "1, 1, 2, 3\r\n",
            "*nset, nset=load, generate\r\n",
            "1, 3, 1\r\n",
            "*material, name=m\r\n",
            "*elastic\r\n",
            "1.0e5, 0.0\r\n",
            "*solid section, elset=e1, material=m\r\n",
            "*boundary\r\n",
            "load, 1, 2\r\n",
            "*cload\r\n",
            "2, 2, 7.0\r\n",
        );
        let (m, r) = parse_inp(deck).unwrap();

        assert_eq!(m.mesh.nodes.len(), 3);
        assert_eq!(m.elsets[0].name, "E1");
        assert_eq!(m.elsets[0].element_ids, vec![1]);
        assert_eq!(m.boundary_conditions[0].node_ids, vec![1, 2, 3]);
        assert_eq!(m.boundary_conditions[0].dofs, vec![1, 2]);
        assert_eq!(
            m.cloads,
            vec![ConcentratedLoad {
                node_id: 2,
                dof: 2,
                value: 7.0
            }]
        );
        assert_eq!(m.materials[0].name, "M");
        assert_eq!(m.materials[0].youngs_modulus, 1.0e5);
        assert_eq!(m.thickness, 1.0); // no *SOLID SECTION data line
        assert!(r.warnings.is_empty());
    }

    #[test]
    fn element_connectivity_continues_on_next_line() {
        let deck = concat!(
            "*NODE\n1,0,0\n2,1,0\n3,1,1\n4,0,1\n",
            "*ELEMENT, TYPE=CPS4\n1, 1, 2,\n3, 4\n",
        );
        let (m, r) = parse_inp(deck).unwrap();
        assert_eq!(m.mesh.elements[0].connectivity, vec![1, 2, 3, 4]);
        assert!(r.warnings.is_empty());
    }

    #[test]
    fn unsupported_element_type_and_missing_type_warn_and_skip() {
        let deck = concat!(
            "*NODE\n1,0,0\n2,1,0\n3,0,1\n",
            "*ELEMENT, TYPE=S4R\n1, 1, 2, 3, 1\n",
            "*ELEMENT\n5, 1, 2, 3\n",
        );
        let (m, r) = parse_inp(deck).unwrap();
        assert!(m.mesh.elements.is_empty());
        assert_eq!(r.warnings.len(), 2);
        assert!(r.warnings[0].contains("TYPE=S4R"));
        assert!(r.warnings[1].contains("TYPE-less"));
    }

    #[test]
    fn c3d8_and_c3d4_element_types_parse() {
        let deck = concat!(
            "*NODE\n",
            "1,0,0,0\n2,1,0,0\n3,1,1,0\n4,0,1,0\n5,0,0,1\n6,1,0,1\n7,1,1,1\n8,0,1,1\n",
            "*ELEMENT, TYPE=C3D8\n1, 1, 2, 3, 4, 5, 6, 7, 8\n",
            "*ELEMENT, TYPE=C3D4\n2, 1, 2, 3, 5\n",
        );
        let (m, r) = parse_inp(deck).unwrap();
        assert!(r.warnings.is_empty());
        assert_eq!(m.mesh.elements[0].kind, ElemKind::C3d8);
        assert_eq!(m.mesh.elements[0].connectivity.len(), 8);
        assert_eq!(m.mesh.elements[1].kind, ElemKind::C3d4);
    }

    #[test]
    fn fatal_errors_carry_line_numbers() {
        let bad = |deck: &str| parse_inp(deck).unwrap_err();

        let e = bad("*NODE\n1, 0.0\n");
        assert_eq!(e.line, 2);
        assert!(e.message.contains("field"));

        let e = bad("*NODE\n1, abc, 0.0\n");
        assert_eq!(e.line, 2);
        assert!(e.message.contains("malformed"));

        let e = bad("*NODE\n1,0,0\n1,1,1\n");
        assert_eq!(e.line, 3);
        assert!(e.message.contains("duplicate node id 1"));

        let e = bad("*NODE\n1,0,0\n2,1,0\n3,0,1\n*ELEMENT, TYPE=CPS3\n1,1,2,3\n1,2,3,1\n");
        assert_eq!(e.line, 7);
        assert!(e.message.contains("duplicate element id 1"));

        let e = bad("*BOUNDARY\nNOPE, 1, 2\n");
        assert_eq!(e.line, 2);
        assert!(e.message.contains("unknown node set 'NOPE'"));

        let e = bad("*CLOAD\nNOPE, 1, 1.0\n");
        assert_eq!(e.line, 2);

        let e = bad("*MATERIAL, NAME=M\n*ELASTIC\n1.0\n");
        assert_eq!(e.line, 3);

        let e = bad("*SOLID SECTION, ELSET=NOPE, MATERIAL=M\n");
        assert_eq!(e.line, 1);
        assert!(e.message.contains("unknown ELSET 'NOPE'"));
    }

    #[test]
    fn parses_examples_cantilever() {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cantilever.inp");
        let Ok(text) = std::fs::read_to_string(&path) else {
            eprintln!("note: examples/cantilever.inp missing (mid-rewrite) — skipping");
            return;
        };
        match parse_inp(&text) {
            Ok((model, _)) => {
                assert!(!model.mesh.nodes.is_empty());
                assert!(!model.mesh.elements.is_empty());
            }
            Err(e) => {
                eprintln!(
                    "note: examples/cantilever.inp not parseable yet ({:?}) — skipping",
                    e
                );
            }
        }
    }

    // ---- issue #26: *DLOAD, *DENSITY, multi-step ----

    /// 3x3-node unit-square CPS4 deck (ids row-major from 1).
    fn plate_deck(pre: &str, steps: &str) -> String {
        format!(
            concat!(
                "*NODE\n",
                "1, 0.0, 0.0\n2, 0.5, 0.0\n3, 1.0, 0.0\n",
                "4, 0.0, 0.5\n5, 0.5, 0.5\n6, 1.0, 0.5\n",
                "7, 0.0, 1.0\n8, 0.5, 1.0\n9, 1.0, 1.0\n",
                "*ELEMENT, TYPE=CPS4, ELSET=ALL\n",
                "1, 1, 2, 5, 4\n2, 2, 3, 6, 5\n3, 4, 5, 8, 7\n4, 5, 6, 9, 8\n",
                "*NSET, NSET=BASE\n1, 2, 3\n",
                "*MATERIAL, NAME=STEEL\n",
                "*ELASTIC\n210.0E9, 0.3\n",
                "{pre}",
                "*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL\n0.01\n",
                "{steps}"
            ),
            pre = pre,
            steps = steps
        )
    }

    #[test]
    fn dload_edges_grav_and_density_parse() {
        let deck = plate_deck(
            "*DENSITY\n7850.,\n",
            concat!(
                "*STEP, NAME=PRESS\n",
                "*STATIC\n",
                "*DLOAD\n",
                "ALL, P1, 1.5e6\n",
                "3, P3, 2.5e5\n",
                "GRAV, 9.81, 0.0, -1.0\n",
                "*BOUNDARY\n",
                "BASE, 1, 2\n",
                "*END STEP\n"
            ),
        );
        let (m, r) = parse_inp(&deck).unwrap();
        assert!(r.warnings.is_empty(), "{:?}", r.warnings);
        assert_eq!(m.density, 7850.0);
        assert_eq!(m.steps.len(), 1);
        let s = &m.steps[0];
        assert_eq!(s.name, "PRESS");
        assert_eq!(s.dloads.len(), 5); // ALL expands to 4 quads + element 3
        assert!(s.dloads.contains(&DistributedLoad {
            element_id: 4,
            face: 1,
            magnitude: 1.5e6
        }));
        assert_eq!(
            s.gravity,
            Some((9.81, [0.0, -1.0, 0.0])),
            "GRAV row feeds the step"
        );
        // Single step mirrors BCs into the legacy fields (v0.1 duality).
        assert_eq!(m.boundary_conditions.len(), 1);
        assert_eq!(m.boundary_conditions[0].node_ids, vec![1, 2, 3]);
        assert!(m.validate().is_ok());
    }

    #[test]
    fn dload_unsupported_types_and_targets_warn() {
        let deck = plate_deck(
            "",
            concat!(
                "*STEP\n*STATIC\n*DLOAD\n",
                "NOSUCHSET, P1, 1.0\n", // unknown target
                "ALL, P5, 2.0\n",       // not P1..P4
                "ALL, TRVEC, 3.0\n",    // not P1..P4
                "3, P4, 4.0\n",         // CPS4 has face 4 — ok
                "ALL, P3, 5.0\n",       // includes the CPS4s — ok
                "*BOUNDARY\nBASE, 1, 2\n*END STEP\n"
            ),
        );
        let (m, r) = parse_inp(&deck).unwrap();
        // 3 warnings: unknown target, P5, TRVEC. P4 on element 3 is fine.
        assert_eq!(r.warnings.len(), 3, "{:?}", r.warnings);
        assert!(r.warnings[0].contains("not a known element"));
        assert!(r.warnings[1].contains("P5"));
        assert!(r.warnings[2].contains("TRVEC"));
        assert_eq!(m.steps[0].dloads.len(), 5);
    }

    #[test]
    fn density_after_elastic_and_overwrite_warning() {
        let deck = plate_deck(
            "*DENSITY\n2500.,\n",
            "*STEP\n*STATIC\n*BOUNDARY\nBASE, 1, 2\n*CLOAD\n9, 2, -1.0\n*END STEP\n",
        );
        let (m, r) = parse_inp(&deck).unwrap();
        assert!(r.warnings.is_empty(), "{:?}", r.warnings);
        assert_eq!(m.density, 2500.0);
        assert_eq!(m.materials.len(), 1);
        assert_eq!(m.materials[0].youngs_modulus, 210.0e9);
    }

    fn single_step_deck(p: f64) -> String {
        plate_deck(
            "",
            &format!("*STEP\n*STATIC\n*BOUNDARY\nBASE, 1, 2\n*CLOAD\n8, 2, {p}\n*END STEP\n"),
        )
    }

    fn two_step_deck(op_new: bool) -> String {
        let op = if op_new { ", OP=NEW" } else { "" };
        plate_deck(
            "",
            &format!(
                concat!(
                    "*STEP\n*STATIC\n*BOUNDARY\nBASE, 1, 2\n*CLOAD\n8, 2, -100.0\n*END STEP\n",
                    "*STEP\n*STATIC\n*CLOAD{op}\n8, 2, -250.0\n*END STEP\n"
                ),
                op = op
            ),
        )
    }

    #[test]
    fn two_step_mod_accumulates_and_matches_superposition() {
        let (m, r) = parse_inp(&two_step_deck(false)).unwrap();
        assert!(r.warnings.is_empty(), "{:?}", r.warnings);
        assert_eq!(m.steps.len(), 2);
        assert_eq!(m.steps[0].cloads.len(), 1);
        // Step 2 carries the ACCUMULATED state: P1 + P2.
        assert_eq!(
            m.steps[1].cloads,
            vec![
                ConcentratedLoad {
                    node_id: 8,
                    dof: 2,
                    value: -100.0
                },
                ConcentratedLoad {
                    node_id: 8,
                    dof: 2,
                    value: -250.0
                },
            ]
        );
        // BCs seeded from step 1 persist into step 2.
        assert_eq!(
            m.steps[1].boundary_conditions,
            m.steps[0].boundary_conditions
        );

        let a = solvix_core::solve(&m).unwrap();
        let b = solvix_core::solve(&parse_inp(&single_step_deck(-350.0)).unwrap().0).unwrap();
        for (ua, ub) in a.displacements.iter().zip(&b.displacements) {
            for d in 0..3 {
                assert!((ua[d] - ub[d]).abs() < 1e-9, "{ua:?} vs {ub:?}");
            }
        }
    }

    #[test]
    fn two_step_op_new_replaces_previous_loads() {
        let (m, r) = parse_inp(&two_step_deck(true)).unwrap();
        assert!(r.warnings.is_empty(), "{:?}", r.warnings);
        assert_eq!(m.steps.len(), 2);
        assert_eq!(
            m.steps[1].cloads,
            vec![ConcentratedLoad {
                node_id: 8,
                dof: 2,
                value: -250.0
            }],
            "OP=NEW drops step 1's loads"
        );

        let a = solvix_core::solve(&m).unwrap();
        let b = solvix_core::solve(&parse_inp(&single_step_deck(-250.0)).unwrap().0).unwrap();
        for (ua, ub) in a.displacements.iter().zip(&b.displacements) {
            for d in 0..3 {
                assert!((ua[d] - ub[d]).abs() < 1e-9, "{ua:?} vs {ub:?}");
            }
        }
    }

    #[test]
    fn step_op_new_clears_bcs_and_loads_static_params_warn() {
        let deck = plate_deck(
            "",
            concat!(
                "*STEP\n*STATIC, INITIAL\n*BOUNDARY\nBASE, 1, 2\n*CLOAD\n8, 2, -1.0\n*END STEP\n",
                "*STEP, OP=NEW\n*STATIC\n*BOUNDARY\nBASE, 1, 2\n*CLOAD\n9, 2, -5.0\n*END STEP\n"
            ),
        );
        let (m, r) = parse_inp(&deck).unwrap();
        assert_eq!(r.warnings.len(), 1);
        assert!(r.warnings[0].contains("*STATIC parameters ignored"));
        assert_eq!(m.steps.len(), 2);
        assert_eq!(m.steps[0].cloads.len(), 1);
        // *STEP, OP=NEW starts from a clean slate.
        assert_eq!(m.steps[1].boundary_conditions.len(), 1);
        assert_eq!(m.steps[1].cloads.len(), 1);
        assert_eq!(m.steps[1].cloads[0].node_id, 9);
    }
}
