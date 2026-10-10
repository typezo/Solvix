//! # solvix-ui — the Abaqus-CAE-style web app for Solvix
//!
//! `solvix serve` starts an HTTP server: static UI under `/` and `/static/`,
//! JSON REST API under `/api/*` (issue #13). The same API serves the browser
//! UI (issue #14) and AI agents (issue #15).
//!
//! API v1 contract:
//! - `GET  /api/model?job=NAME` → mesh+model JSON (see issue #13)
//! - `POST /api/run` body `{"path": "examples/cantilever.inp"}` (optional
//!   `"modes": N` also extracts the N lowest natural frequencies, issue #40)
//!   → run summary JSON (CLI schema; gains `modes`/`frequencies_hz` when run
//!   with modes, `thermal`/`max_temperature_c` for auto-detected thermal
//!   decks (issue #46) and `dynamic`/`steps`/`dt`/`duration`/`final_disp_mm`
//!   for auto-detected `*DYNAMIC` decks (issue #47))
//! - `GET  /api/results/{job}` → displacements + von Mises arrays; after a
//!   modal run also `frequencies_hz` and the first 3 `modes` as per-node
//!   vectors scaled to max |component| = 1.0; thermal runs add
//!   `thermal.temperatures` (issue #46); dynamic runs add
//!   `dynamic.time_history` (t, max |u| mm, <= 201 points) and
//!   `dynamic.frames` (per-node displacement triples at the same times,
//!   driving the frame slider) (issue #47)
//! - `GET  /api/agent` → machine-readable API description
//!
//! std-only: TcpListener + thread per connection, Connection: close.

use std::collections::HashMap;
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use solvix_core::{
    solve_modal, DynSolution, ElemKind, ModalOptions, ModalSolution, Model, Solution, SolveError,
};

const USAGE: &str = "\
solvix-ui 0.1.0 — Abaqus-CAE-style web interface for Solvix

Usage:
  solvix-ui serve [--port 8877] [--root .]   start HTTP server + JSON REST API
  solvix-ui help | --help | -h               show this text

API: GET /api/model?job=NAME · POST /api/run · GET /api/results/{job} · GET /api/agent
";

// ---------------------------------------------------------------- arguments

#[derive(Debug, PartialEq, Eq)]
struct ServeConfig {
    port: u16,
    root: PathBuf,
}

fn parse_args(args: &[String]) -> Result<ServeConfig, String> {
    let mut cfg = ServeConfig {
        port: 8877,
        root: PathBuf::from("."),
    };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "serve" => {}
            "--port" => {
                i += 1;
                let Some(v) = args.get(i) else {
                    return Err("--port requires a value".into());
                };
                cfg.port = v.parse().map_err(|_| format!("invalid port: {v}"))?;
            }
            "--root" => {
                i += 1;
                let Some(v) = args.get(i) else {
                    return Err("--root requires a value".into());
                };
                cfg.root = PathBuf::from(v);
            }
            "help" | "--help" | "-h" => return Err(String::new()),
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    Ok(cfg)
}

// ------------------------------------------------------------------- router

/// Where a request routes to. Pure function of (method, path) — no I/O.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Route {
    /// Sanitized relative path under the static dir.
    Static(String),
    Model,
    Run,
    Results(String),
    Agent,
}

fn route(method: &str, raw_path: &str) -> Option<Route> {
    let path = raw_path.split('?').next().unwrap_or(raw_path);
    match (method, path) {
        ("GET", "/") => Some(Route::Static("index.html".into())),
        ("GET", "/api/model") => Some(Route::Model),
        ("GET", "/api/agent") => Some(Route::Agent),
        ("POST", "/api/run") => Some(Route::Run),
        _ => {
            if method != "GET" {
                return None;
            }
            if let Some(rest) = path.strip_prefix("/static/") {
                return safe_static_rel(rest).map(Route::Static);
            }
            // Root-relative assets referenced by index.html (href="style.css").
            if let Some(file) = path.strip_prefix('/') {
                if matches!(file, "app.js" | "style.css" | "index.html" | "favicon.svg") {
                    return Some(Route::Static(file.into()));
                }
            }
            if let Some(job) = path.strip_prefix("/api/results/") {
                if valid_name(job) {
                    return Some(Route::Results(job.into()));
                }
            }
            None
        }
    }
}

/// Path-traversal guard: a relative static path is safe only if no component
/// is `..` (absolute paths, `\`, and NUL are rejected outright).
fn safe_static_rel(rel: &str) -> Option<String> {
    if rel.is_empty() || rel.starts_with('/') || rel.contains('\\') || rel.contains('\0') {
        return None;
    }
    for seg in rel.split('/') {
        if seg == ".." {
            return None;
        }
    }
    Some(rel.to_string())
}

/// Job/set names used in paths or file joins: no separators, no `..`.
fn valid_name(name: &str) -> bool {
    !name.is_empty() && !name.contains('/') && !name.contains('\\') && !name.contains("..")
}

// ----------------------------------------------------------------- response

struct Response {
    status: u16,
    content_type: &'static str,
    body: Vec<u8>,
}

impl Response {
    fn json(status: u16, body: String) -> Response {
        Response {
            status,
            content_type: "application/json",
            body: body.into_bytes(),
        }
    }

    fn error(status: u16, stage: &str, msg: &str) -> Response {
        Response::json(status, error_json(stage, msg))
    }
}

fn status_text(status: u16) -> &'static str {
    match status {
        200 => "OK",
        400 => "Bad Request",
        404 => "Not Found",
        405 => "Method Not Allowed",
        500 => "Internal Server Error",
        _ => "Unknown",
    }
}

fn content_type(path: &Path) -> &'static str {
    match path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase()
        .as_str()
    {
        "html" | "htm" => "text/html; charset=utf-8",
        "js" => "text/javascript; charset=utf-8",
        "css" => "text/css; charset=utf-8",
        "json" => "application/json",
        "svg" => "image/svg+xml",
        "png" => "image/png",
        _ => "application/octet-stream",
    }
}

// ------------------------------------------------------------- JSON helpers

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// f64 as a valid JSON number (never NaN/inf).
fn num(x: f64) -> String {
    if x.is_finite() {
        format!("{x:?}")
    } else {
        "0.0".into()
    }
}

fn fixed3(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.3}")
    } else {
        "0.000".into()
    }
}

fn exp3(x: f64) -> String {
    if x.is_finite() {
        format!("{x:.3e}")
    } else {
        "0.000e0".into()
    }
}

fn error_json(stage: &str, msg: &str) -> String {
    format!(
        "{{\"ok\":false,\"stage\":\"{}\",\"error\":\"{}\"}}",
        stage,
        json_escape(msg)
    )
}

fn ids(v: &[u32]) -> String {
    v.iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn kind_name(k: ElemKind) -> &'static str {
    match k {
        ElemKind::Cps3 => "cps3",
        ElemKind::Cps4 => "cps4",
        ElemKind::Cps8 => "cps8",
        ElemKind::Cpe4 => "cpe4",
        ElemKind::Cpe8 => "cpe8",
        ElemKind::C3d8 => "c3d8",
        ElemKind::C3d4 => "c3d4",
        ElemKind::B31 => "b31",
    }
}

/// Pull one top-level-ish string field out of a small JSON object body.
/// Handles `\"` escapes; enough for `{"path": "..."}` without a JSON lib.
fn json_string_field(body: &[u8], key: &str) -> Option<String> {
    let s = std::str::from_utf8(body).ok()?;
    let pat = format!("\"{key}\"");
    let mut from = 0usize;
    while let Some(pos) = s[from..].find(&pat) {
        let rest = s[from + pos + pat.len()..].trim_start();
        if let Some(rest) = rest.strip_prefix(':').map(str::trim_start) {
            let mut chars = rest.strip_prefix('"')?.chars();
            let mut out = String::new();
            while let Some(c) = chars.next() {
                match c {
                    '\\' => match chars.next() {
                        Some('n') => out.push('\n'),
                        Some('t') => out.push('\t'),
                        Some('r') => out.push('\r'),
                        Some(o) => out.push(o),
                        None => return None,
                    },
                    '"' => return Some(out),
                    c => out.push(c),
                }
            }
            return None;
        }
        from += pos + pat.len();
    }
    None
}

/// Pull one non-negative integer field out of a small JSON object body
/// (`{"modes": 3}`). `None` when absent or not a plain integer — a
/// non-integer `modes` value keeps the run static-only (the modal count
/// itself is validated by `solve_modal`).
fn json_uint_field(body: &[u8], key: &str) -> Option<usize> {
    let s = std::str::from_utf8(body).ok()?;
    let pat = format!("\"{key}\"");
    let mut from = 0usize;
    while let Some(pos) = s[from..].find(&pat) {
        let rest = s[from + pos + pat.len()..].trim_start();
        if let Some(rest) = rest.strip_prefix(':').map(str::trim_start) {
            let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
            return digits.parse().ok();
        }
        from += pos + pat.len();
    }
    None
}

/// Value of `name` in the query string of a request target (no URL decoding).
fn query_param(target: &str, name: &str) -> Option<String> {
    let q = target.split_once('?')?.1;
    for seg in q.split('&') {
        if let Some((k, v)) = seg.split_once('=') {
            if k == name {
                return Some(v.to_string());
            }
        }
    }
    None
}

// --------------------------------------------------------------- app state

/// Auto-detected transient run data (issue #47): summary scalars plus the
/// time history and per-frame displacement triples the results API serves
/// to the frontend's sparkline + frame slider.
struct UiDynamic {
    steps: usize,
    dt: f64,
    duration: f64,
    /// Max over ALL time points of the largest nodal |u| (m).
    max_disp_m: f64,
    /// Largest nodal |u| at the final time (m).
    final_disp_m: f64,
    /// (t in s, max nodal |u| in mm), subsampled to <= MAX_HISTORY_POINTS.
    time_history: Vec<(f64, f64)>,
    /// Per-node displacement triples (m) at the SAME subsample indices —
    /// the frame slider's deformed shapes.
    frames: Vec<Vec<[f64; 3]>>,
}

/// One cached solved job: model, solution, modal run if requested, plus
/// auto-detected thermal temperatures (issue #46) and dynamics (issue #47).
struct JobResult {
    model: Model,
    /// Static/thermomech solution; for dynamic jobs the final-time frame
    /// packed as a Solution (v1 dynamics recovers no stresses).
    sol: Solution,
    modal: Option<ModalSolution>,
    /// Nodal temperatures (the deck's temperature units) of a thermal run.
    thermal: Option<Vec<f64>>,
    dynamic: Option<UiDynamic>,
}

/// Time-history/frames subsample ceiling served by /api/results.
const MAX_HISTORY_POINTS: usize = 201;

/// Evenly spaced indices over 0..n-1 capped at `cap` (first and last
/// always included, strictly increasing when n > cap).
fn subsample_indices(n: usize, cap: usize) -> Vec<usize> {
    if n == 0 {
        Vec::new()
    } else if n <= cap {
        (0..n).collect()
    } else {
        (0..cap).map(|i| i * (n - 1) / (cap - 1)).collect()
    }
}

/// The run pipeline's analysis-type detection (issues #46/#47): any
/// `*DYNAMIC` step wins; temperature data without dynamics is thermal.
fn has_dynamics(model: &Model) -> bool {
    model.steps.iter().any(|s| s.dynamics.is_some())
}

fn has_thermal(model: &Model) -> bool {
    model
        .steps
        .iter()
        .any(|s| !s.temp_bcs.is_empty() || !s.heat_fluxes.is_empty() || s.heat_source.is_some())
}

/// Largest nodal displacement magnitude in a flat dof vector.
fn max_abs_u(flat: &[f64], ndof: usize) -> f64 {
    let mut max: f64 = 0.0;
    for c in flat.chunks(ndof.max(1)) {
        let mut s = 0.0;
        for &v in c {
            s += v * v;
        }
        max = max.max(s.sqrt());
    }
    max
}

/// Per-node displacement triples of one flat dof vector.
fn flat_triples(flat: &[f64], nodes: usize) -> Vec<[f64; 3]> {
    let ndof = flat.len().checked_div(nodes).unwrap_or(2);
    (0..nodes)
        .map(|i| {
            let b = i * ndof;
            [
                flat.get(b).copied().unwrap_or(0.0),
                flat.get(b + 1).copied().unwrap_or(0.0),
                if ndof == 3 {
                    flat.get(b + 2).copied().unwrap_or(0.0)
                } else {
                    0.0
                },
            ]
        })
        .collect()
}

/// Summarize a transient solution for the cache: scalars + subsampled
/// history and frames.
fn ui_dynamic(model: &Model, dyn_sol: &DynSolution, dt: f64, duration: f64) -> UiDynamic {
    let nodes = model.mesh.nodes.len();
    let last = dyn_sol.points.last().expect("solve_dynamic emits t=0");
    let ndof = last.displacements.len() / nodes.max(1);
    let max_disp_m = dyn_sol
        .points
        .iter()
        .map(|p| max_abs_u(&p.displacements, ndof))
        .fold(0.0, f64::max);
    let final_disp_m = max_abs_u(&last.displacements, ndof);
    let mut time_history = Vec::new();
    let mut frames = Vec::new();
    for i in subsample_indices(dyn_sol.points.len(), MAX_HISTORY_POINTS) {
        let p = &dyn_sol.points[i];
        time_history.push((p.time, max_abs_u(&p.displacements, ndof) * 1000.0));
        frames.push(flat_triples(&p.displacements, nodes));
    }
    UiDynamic {
        steps: dyn_sol.points.len().saturating_sub(1),
        dt,
        duration,
        max_disp_m,
        final_disp_m,
        time_history,
        frames,
    }
}

struct AppState {
    /// Workspace root: `examples/` decks resolve under it, `/api/run` paths
    /// resolve relative to it.
    root: PathBuf,
    static_dir: PathBuf,
    /// job name → last solved JobResult.
    results: Mutex<HashMap<String, JobResult>>,
}

/// `<root>/crates/solvix-ui/static` when run from the workspace root,
/// `<root>/static` when run from the crate dir.
fn find_static_dir(root: &Path) -> PathBuf {
    for cand in ["crates/solvix-ui/static", "static"] {
        let p = root.join(cand);
        if p.is_dir() {
            return p;
        }
    }
    root.join("crates/solvix-ui/static")
}

// ------------------------------------------------------------- API payload

/// Model JSON — exact shape of `static/mock/model.json`.
fn model_json(model: &Model, job: &str) -> String {
    let nodes = model
        .mesh
        .nodes
        .iter()
        .map(|n| format!("[{},{},{}]", n.id, num(n.coords[0]), num(n.coords[1])))
        .collect::<Vec<_>>()
        .join(",");
    let elements = model
        .mesh
        .elements
        .iter()
        .map(|e| {
            format!(
                "[{},\"{}\",[{}]]",
                e.id,
                kind_name(e.kind),
                ids(&e.connectivity)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    // ponytail: Model carries no named nsets (the parser expands *NSET into
    // BCs); emit {} — revisit if the parser ever keeps set names.
    let elsets = model
        .elsets
        .iter()
        .map(|s| format!("\"{}\":[{}]", json_escape(&s.name), ids(&s.element_ids)))
        .collect::<Vec<_>>()
        .join(",");
    let materials = model
        .materials
        .iter()
        .map(|m| {
            format!(
                "{{\"name\":\"{}\",\"young\":{},\"poisson\":{}}}",
                json_escape(&m.name),
                num(m.youngs_modulus),
                num(m.poisson_ratio)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let bcs = model
        .boundary_conditions
        .iter()
        .map(|b| {
            format!(
                "{{\"nodes\":[{}],\"dofs\":[{}]}}",
                ids(&b.node_ids),
                b.dofs
                    .iter()
                    .map(|d| d.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let cloads = model
        .cloads
        .iter()
        .map(|c| format!("[{},{},{}]", c.node_id, c.dof, num(c.value)))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"ok\":true,\"job\":\"{}\",\"nodes\":[{}],\"elements\":[{}],\"nsets\":{{}},\"elsets\":{{{}}},\"materials\":[{}],\"bcs\":[{}],\"cloads\":[{}]}}",
        json_escape(job),
        nodes,
        elements,
        elsets,
        materials,
        bcs,
        cloads,
    )
}

/// Run summary JSON — same schema as `solvix run --json` (issue #3); a
/// non-empty `frequencies_hz` appends the modal fields (issue #40), a
/// thermal run appends `thermal`/`max_temperature_c` (issue #46) and a
/// transient run `dynamic`/`steps`/`dt`/`duration`/`final_disp_mm`
/// (issue #47).
fn run_summary_json(jr: &JobResult, job: &str, solve_ms: u64, outputs: &[String]) -> String {
    let max_disp_m = match &jr.dynamic {
        // Dynamic runs report the max over ALL time points.
        Some(d) => d.max_disp_m,
        None => jr
            .sol
            .displacements
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max),
    };
    let max_mises = jr.sol.von_mises.iter().copied().fold(0.0, f64::max);
    let outs = outputs
        .iter()
        .map(|o| format!("\"{}\"", json_escape(o)))
        .collect::<Vec<_>>()
        .join(",");
    let modal = match &jr.modal {
        Some(m) if !m.frequencies_hz.is_empty() => {
            let freqs = m
                .frequencies_hz
                .iter()
                .map(|f| fixed3(*f))
                .collect::<Vec<_>>()
                .join(",");
            format!(
                ",\"modes\":{},\"frequencies_hz\":[{}]",
                m.frequencies_hz.len(),
                freqs
            )
        }
        _ => String::new(),
    };
    let thermal = match &jr.thermal {
        Some(temps) => format!(
            ",\"thermal\":true,\"max_temperature_c\":{}",
            fixed3(temps.iter().copied().fold(f64::NEG_INFINITY, f64::max))
        ),
        None => String::new(),
    };
    let dynamic = match &jr.dynamic {
        Some(d) => format!(
            ",\"dynamic\":true,\"steps\":{},\"dt\":{},\"duration\":{},\"final_disp_mm\":{}",
            d.steps,
            num(d.dt),
            num(d.duration),
            fixed3(d.final_disp_m * 1000.0)
        ),
        None => String::new(),
    };
    let mut j = format!(
        "{{\"ok\":true,\"job\":\"{}\",\"nodes\":{},\"elements\":{},\"iterations\":{},\"residual\":{},\"solve_ms\":{},\"max_disp_mm\":{},\"max_von_mises_mpa\":{},\"outputs\":[{}]",
        json_escape(job),
        jr.model.mesh.nodes.len(),
        jr.model.mesh.elements.len(),
        jr.sol.iterations,
        exp3(jr.sol.final_residual),
        solve_ms,
        fixed3(max_disp_m * 1000.0),
        fixed3(max_mises / 1e6),
        outs,
    );
    j.push_str(&modal);
    j.push_str(&thermal);
    j.push_str(&dynamic);
    j.push('}');
    j
}

/// Dofs per node for an already-validated model: 3 when any solid or beam
/// element is present, else 2 (mixing is rejected by `solve` before we get
/// here). ponytail: for B31 the third dof is rotz, rendered as z — a
/// display-only approximation.
fn ndof_for_display(model: &Model) -> usize {
    if model
        .mesh
        .elements
        .iter()
        .any(|e| matches!(e.kind, ElemKind::C3d8 | ElemKind::C3d4 | ElemKind::B31))
    {
        3
    } else {
        2
    }
}

/// First [`MAX_DISPLAY_MODES`] mode shapes as per-node `[x,y,z]` triples,
/// each scaled to max |component| = 1.0 for viewport display (issue #40).
fn mode_display_triples(model: &Model, modal: &ModalSolution) -> Vec<Vec<[f64; 3]>> {
    const MAX_DISPLAY_MODES: usize = 3;
    let ndof = ndof_for_display(model);
    let nodes = model.mesh.nodes.len();
    modal
        .modes
        .iter()
        .take(MAX_DISPLAY_MODES)
        .map(|flat| {
            let scale = flat.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let s = if scale > 0.0 { scale } else { 1.0 };
            (0..nodes)
                .map(|n| {
                    let dx = flat.get(n * ndof).copied().unwrap_or(0.0) / s;
                    let dy = flat.get(n * ndof + 1).copied().unwrap_or(0.0) / s;
                    let dz = if ndof == 3 {
                        flat.get(n * ndof + 2).copied().unwrap_or(0.0) / s
                    } else {
                        0.0
                    };
                    [dx, dy, dz]
                })
                .collect()
        })
        .collect()
}

/// Results JSON — exact shape of `static/mock/results.json`; after a modal
/// run also `frequencies_hz` plus display-scaled `modes` triples (issue #40);
/// thermal runs add `thermal.temperatures` (issue #46); dynamic runs add
/// `dynamic.time_history` + `dynamic.frames` (issue #47).
fn results_json(jr: &JobResult, job: &str) -> String {
    let sol = &jr.sol;
    let disp = sol
        .displacements
        .iter()
        .map(|d| format!("[{},{},{}]", num(d[0]), num(d[1]), num(d[2])))
        .collect::<Vec<_>>()
        .join(",");
    let mises = sol
        .von_mises
        .iter()
        .map(|v| num(*v))
        .collect::<Vec<_>>()
        .join(",");
    let nodal = sol
        .nodal_von_mises
        .iter()
        .map(|v| num(*v))
        .collect::<Vec<_>>()
        .join(",");
    let max_disp_m = match &jr.dynamic {
        Some(d) => d.max_disp_m,
        None => sol
            .displacements
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max),
    };
    let max_mises = sol.von_mises.iter().copied().fold(0.0, f64::max);
    let modal_json = match (&jr.modal, &jr.model) {
        (Some(m), model) => {
            let freqs = m
                .frequencies_hz
                .iter()
                .map(|f| fixed3(*f))
                .collect::<Vec<_>>()
                .join(",");
            let modes = mode_display_triples(model, m)
                .iter()
                .map(|mode| {
                    mode.iter()
                        .map(|d| format!("[{},{},{}]", num(d[0]), num(d[1]), num(d[2])))
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .collect::<Vec<_>>()
                .join("],[");
            format!(",\"frequencies_hz\":[{}],\"modes\":[[{}]]", freqs, modes)
        }
        (None, _) => String::new(),
    };
    let thermal_json = match &jr.thermal {
        Some(temps) => {
            let values = temps.iter().map(|t| num(*t)).collect::<Vec<_>>().join(",");
            format!(
                ",\"thermal\":{{\"temperatures\":[{}],\"thermal\":true}}",
                values
            )
        }
        None => String::new(),
    };
    let dynamic_json = match &jr.dynamic {
        Some(d) => {
            let history = d
                .time_history
                .iter()
                .map(|(t, u)| format!("{{\"t\":{},\"max_abs_u_mm\":{}}}", num(*t), fixed3(*u)))
                .collect::<Vec<_>>()
                .join(",");
            let frames = d
                .frames
                .iter()
                .map(|frame| {
                    frame
                        .iter()
                        .map(|v| format!("[{},{},{}]", num(v[0]), num(v[1]), num(v[2])))
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .collect::<Vec<_>>()
                .join("],[");
            format!(
                ",\"dynamic\":{{\"dynamic\":true,\"time_history\":[{}],\"frames\":[[{}]]}}",
                history, frames
            )
        }
        None => String::new(),
    };
    let mut j = format!(
        "{{\"ok\":true,\"job\":\"{}\",\"displacements\":[{}],\"von_mises\":[{}],\"nodal_von_mises\":[{}],\"max_disp_mm\":{},\"max_von_mises_mpa\":{},\"iterations\":{},\"residual\":{}",
        json_escape(job),
        disp,
        mises,
        nodal,
        fixed3(max_disp_m * 1000.0),
        fixed3(max_mises / 1e6),
        sol.iterations,
        exp3(sol.final_residual),
    );
    j.push_str(&thermal_json);
    j.push_str(&dynamic_json);
    j.push_str(&modal_json);
    j.push('}');
    j
}

/// Self-describing API index for humans and AI agents (issue #15).
fn agent_json() -> String {
    let endpoints = [
        (
            "GET",
            "/api/model?job=NAME",
            "null",
            "shape of /static/mock/model.json with ok=true",
            "parse <root>/examples/NAME.inp and return mesh + model JSON",
        ),
        (
            "POST",
            "/api/run",
            "{\"path\": \"examples/NAME.inp\"} (+ optional \"modes\": N for modal analysis)",
            "{\"ok\":true,\"job\":..,\"nodes\":..,\"elements\":..,\"iterations\":..,\"residual\":..,\"solve_ms\":..,\"max_disp_mm\":..,\"max_von_mises_mpa\":..,\"outputs\":[\"NAME.vtu\"]} (CLI run --json schema; with modes also \"modes\":N,\"frequencies_hz\":[Hz]; thermal decks gain \"thermal\":true,\"max_temperature_c\":C; *DYNAMIC decks gain \"dynamic\":true,\"steps\":N,\"dt\":s,\"duration\":s,\"final_disp_mm\":mm with max_disp_mm = max over time)",
            "parse, solve and write NAME.vtu next to the deck; caches results for /api/results/{job}; \"modes\": N also extracts natural frequencies (cached for the mode viewer); the analysis type is auto-detected from the deck (no flags): *DYNAMIC runs transient Newmark dynamics, temperature data without dynamics runs a thermal-stress pass, else linear statics",
        ),
        (
            "GET",
            "/api/results/{job}",
            "null",
            "shape of /static/mock/results.json with ok=true; thermal jobs add \"thermal\":{\"temperatures\":[C per node],\"thermal\":true}; dynamic jobs add \"dynamic\":{\"dynamic\":true,\"time_history\":[{\"t\":s,\"max_abs_u_mm\":mm}..<=201 points],\"frames\":[[[x,y,z] per node]..same times]}",
            "displacements and von Mises arrays of the last POST /api/run of {job} (dynamic jobs: final-time displacements + time history and frames for the player)",
        ),
        (
            "GET",
            "/api/agent",
            "null",
            "this document",
            "machine-readable API self-description",
        ),
        (
            "GET",
            "/ and /static/*",
            "null",
            "static files (html/js/css/json/svg/png)",
            "browser UI (issue #14); mocks live under /static/mock/",
        ),
    ];
    let list = endpoints
        .iter()
        .map(|(m, p, b, r, d)| {
            format!(
                "{{\"method\":\"{m}\",\"path\":\"{}\",\"body\":{b},\"response\":\"{}\",\"description\":\"{}\"}}",
                json_escape(p),
                json_escape(r),
                json_escape(d)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("{{\"ok\":true,\"endpoints\":[{list}]}}")
}

// ---------------------------------------------------------------- handlers

fn serve_static(static_dir: &Path, rel: &str) -> Response {
    let full = static_dir.join(rel);
    match std::fs::read(&full) {
        Ok(bytes) => Response {
            status: 200,
            content_type: content_type(&full),
            body: bytes,
        },
        Err(_) => Response::error(404, "io", &format!("no such file: /static/{rel}")),
    }
}

fn api_model(state: &AppState, job: Option<String>) -> Response {
    let Some(job) = job else {
        return Response::error(400, "usage", "missing ?job=NAME");
    };
    if !valid_name(&job) {
        return Response::error(400, "usage", "invalid job name");
    }
    let deck = state.root.join("examples").join(format!("{job}.inp"));
    let text = match std::fs::read_to_string(&deck) {
        Ok(t) => t,
        Err(e) => {
            return Response::error(404, "io", &format!("cannot read {}: {e}", deck.display()))
        }
    };
    match solvix_io::parse_inp(&text) {
        Ok((model, _)) => Response::json(200, model_json(&model, &job)),
        Err(p) => Response::error(400, "parse", &format!("line {}: {}", p.line, p.message)),
    }
}

fn solve_error_message(e: SolveError) -> String {
    match e {
        SolveError::EmptyMesh => "empty mesh".into(),
        SolveError::NoConstraints => "no constraints".into(),
        SolveError::InvalidModel(m) => format!("invalid model: {m}"),
        SolveError::NotConverged {
            iterations,
            residual,
        } => {
            format!(
                "not converged after {iterations} iterations (residual {})",
                exp3(residual)
            )
        }
    }
}

fn api_run(state: &AppState, body: &[u8]) -> Response {
    let Some(path) = json_string_field(body, "path") else {
        return Response::error(400, "usage", "body must be {\"path\": \"job.inp\"}");
    };
    let deck = if Path::new(&path).is_absolute() {
        PathBuf::from(&path)
    } else {
        state.root.join(&path)
    };
    let text = match std::fs::read_to_string(&deck) {
        Ok(t) => t,
        Err(e) => {
            return Response::error(404, "io", &format!("cannot read {}: {e}", deck.display()))
        }
    };
    let (model, _) = match solvix_io::parse_inp(&text) {
        Ok(m) => m,
        Err(p) => return Response::error(400, "parse", &format!("line {}: {}", p.line, p.message)),
    };
    let t0 = Instant::now();
    // Analysis-type dispatch (issues #46/#47): *DYNAMIC wins; temperature
    // data without dynamics runs the thermal-stress pass; else statics.
    let (sol, thermal, dynamic) = if has_dynamics(&model) {
        let dyn_sol = match solvix_core::solve_dynamic(&model) {
            Ok(d) => d,
            Err(e) => return Response::error(400, "solve", &solve_error_message(e)),
        };
        let spec = model
            .steps
            .iter()
            .rev()
            .find(|s| s.dynamics.is_some())
            .and_then(|s| s.dynamics.as_ref())
            .expect("has_dynamics checked is_some");
        let info = ui_dynamic(&model, &dyn_sol, spec.dt, spec.duration);
        let nodes = model.mesh.nodes.len();
        let flat = dyn_sol
            .points
            .last()
            .map(|p| p.displacements.as_slice())
            .unwrap_or(&[]);
        let final_sol = Solution {
            displacements: flat_triples(flat, nodes),
            von_mises: Vec::new(),
            nodal_von_mises: Vec::new(),
            reactions: Vec::new(),
            iterations: dyn_sol.iterations_total,
            final_residual: 0.0,
            solver_method: "newmark".into(),
            cond_estimate: 0.0,
        };
        (final_sol, None, Some(info))
    } else if has_thermal(&model) {
        // solve_thermomech runs the conduction pass internally; the extra
        // solve_thermal recovers the temperature field for the API
        // (core offers no temps-out variant — the double conduction solve
        // is the v0.4 price).
        let temps = match solvix_core::solve_thermal(&model) {
            Ok(t) => t,
            Err(e) => return Response::error(400, "solve", &solve_error_message(e)),
        };
        let sol = match solvix_core::solve_thermomech(&model) {
            Ok(s) => s,
            Err(e) => return Response::error(400, "solve", &solve_error_message(e)),
        };
        (sol, Some(temps.temperatures), None)
    } else {
        let sol = match solvix_core::solve(&model) {
            Ok(s) => s,
            Err(e) => return Response::error(400, "solve", &solve_error_message(e)),
        };
        (sol, None, None)
    };
    // Optional modal pass (issue #40); count validated by solve_modal.
    let modal = match json_uint_field(body, "modes") {
        Some(n) => match solve_modal(&model, &ModalOptions { num_modes: n }) {
            Ok(m) => Some(m),
            Err(e) => return Response::error(400, "solve", &solve_error_message(e)),
        },
        None => None,
    };
    let solve_ms = t0.elapsed().as_millis() as u64;
    let job = deck
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("job")
        .to_string();
    let vtu = deck.with_extension("vtu");
    if let Err(e) = solvix_io::write_vtu(&vtu, &model, &sol) {
        return Response::error(500, "io", &format!("cannot write {}: {e}", vtu.display()));
    }
    let jr = JobResult {
        model,
        sol,
        modal,
        thermal,
        dynamic,
    };
    let summary = run_summary_json(&jr, &job, solve_ms, &[vtu.display().to_string()]);
    if let Ok(mut cache) = state.results.lock() {
        cache.insert(job, jr);
    }
    Response::json(200, summary)
}

fn api_results(state: &AppState, job: &str) -> Response {
    let cache = state.results.lock().unwrap_or_else(|p| p.into_inner());
    match cache.get(job) {
        Some(jr) => Response::json(200, results_json(jr, job)),
        None => Response::error(
            404,
            "io",
            &format!("unknown job '{job}'; POST /api/run {{\"path\":\"{job}.inp\"}} first"),
        ),
    }
}

/// Dispatch one request to a response. Pure given `state` + inputs — the
/// unit-testable core of the server (no sockets).
fn handle(state: &AppState, method: &str, target: &str, body: &[u8]) -> Response {
    match route(method, target) {
        Some(Route::Static(rel)) => serve_static(&state.static_dir, &rel),
        Some(Route::Model) => api_model(state, query_param(target, "job")),
        Some(Route::Run) => api_run(state, body),
        Some(Route::Results(job)) => api_results(state, &job),
        Some(Route::Agent) => Response::json(200, agent_json()),
        None => Response::error(404, "usage", &format!("no route: {method} {target}")),
    }
}

// ------------------------------------------------------------------- server

fn handle_conn(state: &AppState, mut stream: TcpStream) {
    let mut buf = Vec::new();
    let mut chunk = [0u8; 4096];
    // Read until end of headers (or a sane cap).
    while buf.len() <= 64 * 1024 {
        let n = match stream.read(&mut chunk) {
            Ok(0) | Err(_) => return,
            Ok(n) => n,
        };
        buf.extend_from_slice(&chunk[..n]);
        if find_subslice(&buf, b"\r\n\r\n").is_some() {
            break;
        }
    }
    let Some(hdr_end) = find_subslice(&buf, b"\r\n\r\n") else {
        let resp = Response::error(400, "usage", "malformed HTTP request");
        let _ = write_response(&mut stream, &resp);
        return;
    };
    let head = String::from_utf8_lossy(&buf[..hdr_end]).into_owned();
    let mut lines = head.split("\r\n");
    let Some(request_line) = lines.next() else {
        let resp = Response::error(400, "usage", "empty request");
        let _ = write_response(&mut stream, &resp);
        return;
    };
    let mut parts = request_line.split_whitespace();
    let method = parts.next().unwrap_or("").to_string();
    let target = parts.next().unwrap_or("").to_string();
    if method.is_empty() || target.is_empty() {
        let resp = Response::error(400, "usage", "malformed request line");
        let _ = write_response(&mut stream, &resp);
        return;
    }
    let mut content_length = 0usize;
    for line in lines {
        if let Some((k, v)) = line.split_once(':') {
            if k.trim().eq_ignore_ascii_case("content-length") {
                content_length = v.trim().parse().unwrap_or(0);
            }
        }
    }
    if content_length > 1024 * 1024 {
        let resp = Response::error(400, "usage", "body too large");
        let _ = write_response(&mut stream, &resp);
        return;
    }
    let mut body: Vec<u8> = buf[hdr_end + 4..].to_vec();
    while body.len() < content_length {
        let n = match stream.read(&mut chunk) {
            Ok(0) | Err(_) => break,
            Ok(n) => n,
        };
        body.extend_from_slice(&chunk[..n]);
    }
    body.truncate(content_length);
    let resp = handle(state, &method, &target, &body);
    let _ = write_response(&mut stream, &resp);
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|w| w == needle)
}

fn write_response(stream: &mut TcpStream, resp: &Response) -> std::io::Result<()> {
    let head = format!(
        "HTTP/1.1 {} {}\r\nContent-Type: {}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
        resp.status,
        status_text(resp.status),
        resp.content_type,
        resp.body.len()
    );
    stream.write_all(head.as_bytes())?;
    stream.write_all(&resp.body)?;
    stream.flush()
}

fn serve(cfg: ServeConfig) -> std::io::Result<()> {
    let static_dir = find_static_dir(&cfg.root);
    let state = Arc::new(AppState {
        root: cfg.root.clone(),
        static_dir,
        results: Mutex::new(HashMap::new()),
    });
    let listener = TcpListener::bind(("127.0.0.1", cfg.port))?;
    eprintln!(
        "solvix-ui: serving http://127.0.0.1:{} (root: {}, static: {})",
        cfg.port,
        cfg.root.display(),
        state.static_dir.display()
    );
    for stream in listener.incoming() {
        let Ok(stream) = stream else { continue };
        let state = Arc::clone(&state);
        std::thread::spawn(move || handle_conn(&state, stream));
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&args) {
        Ok(cfg) => {
            if let Err(e) = serve(cfg) {
                eprintln!("error (io): {e}");
                std::process::exit(1);
            }
        }
        Err(msg) => {
            eprint!("{USAGE}");
            if !msg.is_empty() {
                eprintln!("error (usage): {msg}");
                std::process::exit(1);
            }
        }
    }
}

// -------------------------------------------------------------------- tests

#[cfg(test)]
mod tests {
    use super::*;

    const TINY_DECK: &str = "\
*HEADING
    tiny fixture
*NODE
1, 0.0, 0.0
2, 1.0, 0.0
3, 0.0, 1.0
*ELEMENT, TYPE=CPS3, ELSET=ALL
1, 1, 2, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210e9, 0.3
*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL
1.0
*STEP
*STATIC
*BOUNDARY
1, 1, 2
*CLOAD
3, 2, -100.0
*END STEP
";

    fn sv(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    /// Two-triangle square, left edge clamped: 4 free dofs, stable, so the
    /// modal pencil is SPD (the single pinned triangle of TINY_DECK is a
    /// rotational mechanism and legitimately fails the modal factorization).
    const TINY_MODAL_DECK: &str = "\
*HEADING
    modal fixture
*NODE
1, 0.0, 0.0
2, 1.0, 0.0
3, 1.0, 1.0
4, 0.0, 1.0
*ELEMENT, TYPE=CPS3, ELSET=ALL
1, 1, 2, 3
2, 1, 3, 4
*MATERIAL, NAME=STEEL
*ELASTIC
210e9, 0.3
*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL
1.0
*STEP
*STATIC
*BOUNDARY
1, 1, 2
4, 1, 2
*CLOAD
2, 2, -100.0
3, 2, -100.0
*END STEP
";

    fn fixture_root(tag: &str) -> PathBuf {
        let root =
            std::env::temp_dir().join(format!("solvix-ui-test-{}-{tag}", std::process::id()));
        let ex = root.join("examples");
        std::fs::create_dir_all(&ex).unwrap();
        std::fs::write(ex.join("tiny.inp"), TINY_DECK).unwrap();
        std::fs::write(ex.join("tiny_modal.inp"), TINY_MODAL_DECK).unwrap();
        root
    }

    fn test_state(root: &Path) -> AppState {
        AppState {
            root: root.to_path_buf(),
            static_dir: PathBuf::from("static"),
            results: Mutex::new(HashMap::new()),
        }
    }

    fn job_result(sol: Solution, modal: Option<ModalSolution>) -> JobResult {
        let (model, _) = solvix_io::parse_inp(TINY_DECK).unwrap();
        JobResult {
            model,
            sol,
            modal,
            thermal: None,
            dynamic: None,
        }
    }

    // -- router dispatch ---------------------------------------------------

    #[test]
    fn router_dispatch_table() {
        assert_eq!(route("GET", "/"), Some(Route::Static("index.html".into())));
        assert_eq!(
            route("GET", "/style.css"),
            Some(Route::Static("style.css".into()))
        );
        assert_eq!(
            route("GET", "/app.js"),
            Some(Route::Static("app.js".into()))
        );
        assert_eq!(
            route("GET", "/static/mock/model.json"),
            Some(Route::Static("mock/model.json".into()))
        );
        assert_eq!(
            route("GET", "/static/a/b.js"),
            Some(Route::Static("a/b.js".into()))
        );
        assert_eq!(
            route("GET", "/api/model?job=cantilever"),
            Some(Route::Model)
        );
        assert_eq!(route("GET", "/api/model"), Some(Route::Model));
        assert_eq!(route("POST", "/api/run"), Some(Route::Run));
        assert_eq!(
            route("GET", "/api/results/cantilever"),
            Some(Route::Results("cantilever".into()))
        );
        assert_eq!(route("GET", "/api/agent"), Some(Route::Agent));
        assert_eq!(route("GET", "/nope"), None);
        assert_eq!(route("GET", "/api/run"), None, "POST-only route");
        assert_eq!(route("POST", "/api/model?job=x"), None, "GET-only route");
        assert_eq!(route("DELETE", "/api/results/x"), None);
    }

    #[test]
    fn router_rejects_traversal() {
        assert_eq!(route("GET", "/static/../Cargo.toml"), None);
        assert_eq!(route("GET", "/static/a/../../Cargo.toml"), None);
        assert!(route("GET", "/static//etc/passwd").is_none());
        assert_eq!(
            route("GET", "/api/results/../run"),
            None,
            "job segment with .."
        );
        assert_eq!(route("GET", "/api/results/a/b"), None, "job segment with /");
        assert_eq!(safe_static_rel("../x"), None);
        assert_eq!(safe_static_rel("a/../x"), None);
        assert_eq!(safe_static_rel("/abs"), None);
        assert_eq!(safe_static_rel("a\\b"), None);
        assert_eq!(
            safe_static_rel("mock/model.json").as_deref(),
            Some("mock/model.json")
        );
    }

    #[test]
    fn traversal_through_handle_is_404() {
        let state = test_state(Path::new("."));
        let resp = handle(&state, "GET", "/static/../Cargo.toml", b"");
        assert_eq!(resp.status, 404);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.starts_with(r#"{"ok":false,"stage":"usage""#), "{body}");
    }

    // -- static ------------------------------------------------------------

    #[test]
    fn static_mock_served_with_json_type() {
        let state = test_state(Path::new("."));
        let resp = handle(&state, "GET", "/static/mock/model.json", b"");
        assert_eq!(resp.status, 200);
        assert_eq!(resp.content_type, "application/json");
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.contains(r#""nsets": { "FIX""#), "{body}");
    }

    #[test]
    fn static_unknown_file_404() {
        let state = test_state(Path::new("."));
        let resp = handle(&state, "GET", "/static/nope.js", b"");
        assert_eq!(resp.status, 404);
        assert_eq!(resp.content_type, "application/json");
    }

    // -- /api/model ----------------------------------------------------------

    #[test]
    fn model_json_exact_shape_from_real_parser() {
        let (model, _) = solvix_io::parse_inp(TINY_DECK).unwrap();
        assert_eq!(
            model_json(&model, "tiny"),
            r#"{"ok":true,"job":"tiny","nodes":[[1,0.0,0.0],[2,1.0,0.0],[3,0.0,1.0]],"elements":[[1,"cps3",[1,2,3]]],"nsets":{},"elsets":{"ALL":[1]},"materials":[{"name":"STEEL","young":210000000000.0,"poisson":0.3}],"bcs":[{"nodes":[1],"dofs":[1,2]}],"cloads":[[3,2,-100.0]]}"#
        );
    }

    #[test]
    fn model_endpoint_end_to_end() {
        let root = fixture_root("modelerr");
        let state = test_state(&root);
        let resp = handle(&state, "GET", "/api/model?job=tiny", b"");
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.starts_with(r#"{"ok":true,"job":"tiny""#), "{body}");
        assert!(
            body.contains(r#""elements":[[1,"cps3",[1,2,3]]]"#),
            "{body}"
        );
    }

    #[test]
    fn model_errors() {
        let root = fixture_root("modele2e");
        let state = test_state(&root);
        let resp = handle(&state, "GET", "/api/model", b"");
        assert_eq!(resp.status, 400);
        assert_eq!(
            String::from_utf8(resp.body).unwrap(),
            r#"{"ok":false,"stage":"usage","error":"missing ?job=NAME"}"#
        );

        let resp = handle(&state, "GET", "/api/model?job=missing", b"");
        assert_eq!(resp.status, 404);
        assert_eq!(
            String::from_utf8(resp.body).unwrap(),
            error_json(
                "io",
                &format!(
                    "cannot read {}/examples/missing.inp: No such file or directory (os error 2)",
                    root.display()
                )
            )
        );

        let resp = handle(&state, "GET", "/api/model?job=../secrets", b"");
        assert_eq!(resp.status, 400);
        assert_eq!(
            String::from_utf8(resp.body).unwrap(),
            r#"{"ok":false,"stage":"usage","error":"invalid job name"}"#
        );

        std::fs::write(root.join("examples/bad.inp"), "*NODE\n1, x.y, 0.0\n").unwrap();
        let resp = handle(&state, "GET", "/api/model?job=bad", b"");
        assert_eq!(resp.status, 400);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(
            body.starts_with(r#"{"ok":false,"stage":"parse","error":"line "#),
            "{body}"
        );
    }

    // -- /api/run + /api/results ---------------------------------------------

    #[test]
    fn run_then_results_end_to_end() {
        let root = fixture_root("rune2e");
        let state = test_state(&root);

        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"  { "path" : "examples/tiny.inp" } "#,
        );
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(
            body.starts_with(r#"{"ok":true,"job":"tiny","nodes":3,"elements":1,"iterations":"#),
            "{body}"
        );
        assert!(body.contains(r#""max_von_mises_mpa":"#), "{body}");
        assert!(body.contains(r#""outputs":["#), "{body}");
        let vtu = root.join("examples/tiny.vtu");
        assert!(vtu.is_file(), "vtu written next to deck");

        let resp = handle(&state, "GET", "/api/results/tiny", b"");
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(
            body.starts_with(r#"{"ok":true,"job":"tiny","displacements":[[0.0"#),
            "{body}"
        );
        assert!(body.contains(r#""von_mises":["#), "{body}");
        assert!(body.contains(r#""max_disp_mm":"#), "{body}");
        assert!(body.contains(r#""iterations":"#), "{body}");
        assert!(
            body.ends_with(&format!(
                r#""residual":{}}}"#,
                exp3(
                    solvix_core::solve(&solvix_io::parse_inp(TINY_DECK).unwrap().0)
                        .unwrap()
                        .final_residual
                )
            )),
            "{body}"
        );
    }

    #[test]
    fn run_errors() {
        let root = fixture_root("runerr");
        let state = test_state(&root);

        let resp = handle(&state, "POST", "/api/run", b"not json");
        assert_eq!(resp.status, 400);
        assert_eq!(
            String::from_utf8(resp.body).unwrap(),
            r#"{"ok":false,"stage":"usage","error":"body must be {\"path\": \"job.inp\"}"}"#
        );

        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/nope.inp"}"#,
        );
        assert_eq!(resp.status, 404);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.starts_with(r#"{"ok":false,"stage":"io""#), "{body}");

        let resp = handle(&state, "GET", "/api/results/never-ran", b"");
        assert_eq!(resp.status, 404);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(
            body.starts_with(r#"{"ok":false,"stage":"io","error":"unknown job 'never-ran'"#),
            "{body}"
        );
    }

    #[test]
    fn results_json_exact_shape() {
        let sol = Solution {
            displacements: vec![[0.0, 0.0, 0.0], [4.0e-6, -8.0e-6, 0.0]],
            von_mises: vec![120.0e6, 119.0e6],
            nodal_von_mises: vec![120.0e6, 119.0e6],
            reactions: vec![],
            iterations: 390,
            final_residual: 6.933e-11,
            ..Solution::default()
        };
        let jr = job_result(sol, None);
        assert_eq!(
            results_json(&jr, "cantilever"),
            r#"{"ok":true,"job":"cantilever","displacements":[[0.0,0.0,0.0],[4e-6,-8e-6,0.0]],"von_mises":[120000000.0,119000000.0],"nodal_von_mises":[120000000.0,119000000.0],"max_disp_mm":0.009,"max_von_mises_mpa":120.000,"iterations":390,"residual":6.933e-11}"#
        );
    }

    #[test]
    fn run_summary_json_appends_modal_fields() {
        let (model, _) = solvix_io::parse_inp(TINY_DECK).unwrap();
        let sol = solvix_core::solve(&model).unwrap();
        let modal = ModalSolution {
            frequencies_hz: vec![41.5, 260.25],
            modes: Vec::new(),
            iterations: 3,
        };
        let j = run_summary_json(
            &job_result(sol.clone(), Some(modal)),
            "tiny",
            3,
            &["tiny.vtu".into()],
        );
        assert!(
            j.ends_with(r#""outputs":["tiny.vtu"],"modes":2,"frequencies_hz":[41.500,260.250]}"#),
            "{j}"
        );
        let static_only = run_summary_json(&job_result(sol, None), "tiny", 3, &["tiny.vtu".into()]);
        assert!(
            static_only.ends_with(r#""outputs":["tiny.vtu"]}"#),
            "{static_only}"
        );
    }

    #[test]
    fn mode_display_triples_scale_and_cap() {
        let (model, _) = solvix_io::parse_inp(TINY_DECK).unwrap();
        // 3 nodes x 2 dofs; 4 modes asked, only 3 may come back for display.
        let modal = ModalSolution {
            frequencies_hz: vec![10.0, 40.0, 90.0, 160.0],
            modes: vec![
                vec![0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, 0.0, 0.5, -0.5],
                vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
                vec![9.0, 9.0, 9.0, 9.0, 9.0, 9.0],
            ],
            iterations: 7,
        };
        let triples = mode_display_triples(&model, &modal);
        assert_eq!(triples.len(), 3, "first 3 modes only");
        for mode in &triples {
            let max = mode
                .iter()
                .flat_map(|d| [d[0].abs(), d[1].abs(), d[2].abs()])
                .fold(0.0f64, f64::max);
            assert!((max - 1.0).abs() < 1e-12, "max {max}");
        }
        assert_eq!(triples[0][0], [0.0, 0.0, 0.0]); // node 1 constrained
        assert_eq!(triples[0][1], [1.0, 0.0, 0.0]); // 2.0 / 2.0
        assert_eq!(triples[1][2], [1.0, -1.0, 0.0]); // 0.5 / 0.5, -0.5 / 0.5
    }

    #[test]
    fn results_json_with_modal_exact_tail() {
        let (model, _) = solvix_io::parse_inp(TINY_DECK).unwrap();
        let sol = solvix_core::solve(&model).unwrap();
        let modal = ModalSolution {
            frequencies_hz: vec![12.5, 99.125],
            modes: vec![
                vec![0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.0],
            ],
            iterations: 5,
        };
        let mut jr = job_result(sol, Some(modal));
        jr.model = model;
        let j = results_json(&jr, "tiny");
        assert!(
            j.ends_with(
                r#""frequencies_hz":[12.500,99.125],"modes":[[[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.0],[1.0,0.0,0.0]]]}"#
            ),
            "{j}"
        );
    }

    #[test]
    fn json_uint_field_parses() {
        assert_eq!(
            json_uint_field(br#"{"path":"a.inp","modes":3}"#, "modes"),
            Some(3)
        );
        assert_eq!(
            json_uint_field(br#"  { "modes" : 12 , "x":1}"#, "modes"),
            Some(12)
        );
        assert_eq!(json_uint_field(br#"{"path":"a.inp"}"#, "modes"), None);
        assert_eq!(json_uint_field(br#"{"modes":"many"}"#, "modes"), None);
        assert_eq!(json_uint_field(br#"{"modes":-2}"#, "modes"), None);
    }

    #[test]
    fn run_with_modes_end_to_end() {
        let root = fixture_root("modes");
        let state = test_state(&root);
        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/tiny_modal.inp","modes":2}"#,
        );
        assert_eq!(resp.status, 200, "{}", String::from_utf8_lossy(&resp.body));
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.contains(r#""modes":2,"frequencies_hz":["#), "{body}");
        assert!(body.ends_with("]}"), "{body}");

        let resp = handle(&state, "GET", "/api/results/tiny_modal", b"");
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.contains(r#""frequencies_hz":["#), "{body}");
        assert!(body.contains(r#""modes":[[["#), "{body}");
        assert!(body.ends_with("]]}"), "{body}");

        // static-only run for the same job drops the modal fields again
        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/tiny_modal.inp"}"#,
        );
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(!body.contains("frequencies_hz"), "{body}");
        let resp = handle(&state, "GET", "/api/results/tiny_modal", b"");
        let body = String::from_utf8(resp.body).unwrap();
        assert!(!body.contains("frequencies_hz"), "{body}");
    }

    #[test]
    fn run_with_invalid_modes_count_is_solve_error() {
        let root = fixture_root("modes0");
        let state = test_state(&root);
        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/tiny.inp","modes":0}"#,
        );
        assert_eq!(resp.status, 400);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.contains(r#""stage":"solve""#), "{body}");
        assert!(body.contains("num_modes"), "{body}");
    }

    // -- issues #46/#47: auto-detected thermal and dynamic runs ----------

    /// 2x1 CPS4 bar: clamped left edge, T=100 at x=0, T=0 at x=0.1.
    const TINY_THERMAL_DECK: &str = "\
*HEADING
    thermal fixture
*NODE
1, 0.0, 0.0
2, 0.05, 0.0
3, 0.1, 0.0
4, 0.0, 0.025
5, 0.05, 0.025
6, 0.1, 0.025
*ELEMENT, TYPE=CPS4, ELSET=ALL
1, 1, 2, 5, 4
2, 2, 3, 6, 5
*NSET, NSET=HOT
1, 4
*NSET, NSET=COLD
3, 6
*MATERIAL, NAME=STEEL
*ELASTIC
210e9, 0.3
*CONDUCTIVITY
45.0
*EXPANSION
1.2e-5
*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL
0.01
*STEP
*HEAT TRANSFER, STEADY STATE
*TEMPERATURE
HOT, 100.0
COLD, 0.0
*BOUNDARY
1, 1, 2
4, 1, 2
*END STEP
";

    /// Same mesh clamped left, tip load, transient: dt x duration.
    fn tiny_dynamic_deck(dt: &str, duration: &str) -> String {
        format!(
            "\
*HEADING
    dynamics fixture
*NODE
1, 0.0, 0.0
2, 0.05, 0.0
3, 0.1, 0.0
4, 0.0, 0.025
5, 0.05, 0.025
6, 0.1, 0.025
*ELEMENT, TYPE=CPS4, ELSET=ALL
1, 1, 2, 5, 4
2, 2, 3, 6, 5
*MATERIAL, NAME=STEEL
*ELASTIC
210e9, 0.3
*DENSITY
7850.
*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL
0.01
*STEP
*DYNAMIC
{dt}, {duration}
*BOUNDARY
1, 1, 2
4, 1, 2
*CLOAD
3, 2, -100.0
6, 2, -100.0
*END STEP
"
        )
    }

    #[test]
    fn subsample_indices_caps_and_spans() {
        assert_eq!(subsample_indices(0, 201), Vec::<usize>::new());
        assert_eq!(subsample_indices(5, 201), vec![0, 1, 2, 3, 4]);
        let idx = subsample_indices(501, 201);
        assert_eq!(idx.len(), 201);
        assert_eq!(idx[0], 0);
        assert_eq!(idx[200], 500);
        assert!(idx.windows(2).all(|w| w[0] < w[1]), "strictly increasing");
        let idx = subsample_indices(202, 201);
        assert_eq!(idx.len(), 201);
        assert_eq!(idx[0], 0);
        assert_eq!(idx[200], 201);
    }

    #[test]
    fn run_thermal_end_to_end() {
        let root = fixture_root("thermal");
        std::fs::write(root.join("examples/tiny_thermal.inp"), TINY_THERMAL_DECK).unwrap();
        let state = test_state(&root);
        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/tiny_thermal.inp"}"#,
        );
        assert_eq!(resp.status, 200, "{}", String::from_utf8_lossy(&resp.body));
        let body = String::from_utf8(resp.body).unwrap();
        assert!(
            body.contains(r#""thermal":true,"max_temperature_c":100.000"#),
            "{body}"
        );
        assert!(!body.contains(r#""dynamic""#), "{body}");

        let resp = handle(&state, "GET", "/api/results/tiny_thermal", b"");
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.contains(r#""thermal":{"temperatures":["#), "{body}");
        assert!(body.contains(r#""thermal":true}}"#), "{body}");
        // Linear gradient: node 1 (x=0) is 100, node 3 (x=0.1) is 0.
        assert!(body.contains(r#""temperatures":[100.0,"#), "{body}");
        assert!(body.contains(r#",0.0],"thermal":true}}"#), "{body}");

        // static deck: no thermal fields
        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/tiny.inp"}"#,
        );
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(!body.contains(r#""thermal""#), "{body}");
    }

    #[test]
    fn run_dynamic_end_to_end() {
        let root = fixture_root("dyn");
        std::fs::write(
            root.join("examples/tiny_dyn.inp"),
            tiny_dynamic_deck("1.0E-3", "0.01"),
        )
        .unwrap();
        let state = test_state(&root);
        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/tiny_dyn.inp"}"#,
        );
        assert_eq!(resp.status, 200, "{}", String::from_utf8_lossy(&resp.body));
        let body = String::from_utf8(resp.body).unwrap();
        assert!(
            body.contains(
                r#""dynamic":true,"steps":10,"dt":0.001,"duration":0.01,"final_disp_mm":"#
            ),
            "{body}"
        );
        assert!(!body.contains(r#""thermal""#), "{body}");

        let resp = handle(&state, "GET", "/api/results/tiny_dyn", b"");
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(
            body.contains(r#""dynamic":{"dynamic":true,"time_history":["#),
            "{body}"
        );
        assert_eq!(
            body.matches("\"t\":").count(),
            11,
            "11 points, no subsample"
        );
        assert!(body.contains(r#""t":0.0,"max_abs_u_mm":0.000"#), "{body}");
        // frames: 11 frames of 6 nodes; frame 0 is all-zero (rest state).
        assert!(body.contains(r#""frames":[[[0.0,0.0,0.0]"#), "{body}");
    }

    #[test]
    fn dynamic_history_is_subsampled_to_cap() {
        let root = fixture_root("dynsub");
        // 501 output points -> capped at MAX_HISTORY_POINTS (201).
        std::fs::write(
            root.join("examples/tiny_dynsub.inp"),
            tiny_dynamic_deck("1.0E-3", "0.5"),
        )
        .unwrap();
        let state = test_state(&root);
        let resp = handle(
            &state,
            "POST",
            "/api/run",
            br#"{"path":"examples/tiny_dynsub.inp"}"#,
        );
        assert_eq!(resp.status, 200, "{}", String::from_utf8_lossy(&resp.body));
        let resp = handle(&state, "GET", "/api/results/tiny_dynsub", b"");
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert_eq!(body.matches("\"t\":").count(), MAX_HISTORY_POINTS, "{body}");
        assert!(body.contains(r#""t":0.0,"max_abs_u_mm":0.000"#), "{body}");
        // last subsample index is the final point: t = 0.5.
        assert!(body.contains(r#""t":0.5,"#), "{body}");
    }

    // -- /api/agent ------------------------------------------------------------

    #[test]
    fn agent_self_describes_all_routes() {
        let state = test_state(Path::new("."));
        let resp = handle(&state, "GET", "/api/agent", b"");
        assert_eq!(resp.status, 200);
        let body = String::from_utf8(resp.body).unwrap();
        assert!(body.starts_with(r#"{"ok":true,"endpoints":["#), "{body}");
        for frag in [
            r#""method":"GET","path":"/api/model?job=NAME""#,
            r#""method":"POST","path":"/api/run""#,
            r#""path":"/api/results/{job}""#,
            r#""path":"/api/agent""#,
            r#""body":{"path": "examples/NAME.inp"}"#,
            r#""description":""#,
            r#""response":""#,
        ] {
            assert!(body.contains(frag), "missing {frag} in {body}");
        }
    }

    // -- JSON escaping & helpers -------------------------------------------------

    #[test]
    fn json_escaping_in_payloads() {
        assert_eq!(json_escape("a\"b\\c\nd"), "a\\\"b\\\\c\\nd");
        let (mut model, _) = solvix_io::parse_inp(TINY_DECK).unwrap();
        model.materials[0].name = "MA\"T".into();
        model.elsets[0].name = "SE\\T".into();
        let j = model_json(&model, "jo\"b");
        assert!(j.contains(r#""job":"jo\"b""#), "{j}");
        assert!(j.contains(r#""name":"MA\"T""#), "{j}");
        assert!(j.contains(r#""SE\\T":[1]"#), "{j}");
    }

    #[test]
    fn json_string_field_parses_and_unescapes() {
        assert_eq!(
            json_string_field(br#"{"path": "a\"b.inp", "x": 1}"#, "path").as_deref(),
            Some("a\"b.inp")
        );
        assert_eq!(
            json_string_field(br#"  {"other":1, "path":"ex/t.inp"}"#, "path").as_deref(),
            Some("ex/t.inp")
        );
        assert_eq!(json_string_field(br#"{"nope": 1}"#, "path"), None);
        assert_eq!(json_string_field(b"\xff\xfe", "path"), None);
    }

    #[test]
    fn query_param_extracts_last_value() {
        assert_eq!(query_param("/api/model?job=a", "job").as_deref(), Some("a"));
        assert_eq!(
            query_param("/api/model?x=1&job=b&y=2", "job").as_deref(),
            Some("b")
        );
        assert_eq!(query_param("/api/model", "job"), None);
        assert_eq!(query_param("/api/model?jobs=a", "job"), None);
    }

    #[test]
    fn content_types_by_extension() {
        assert_eq!(
            content_type(Path::new("a/index.html")),
            "text/html; charset=utf-8"
        );
        assert_eq!(
            content_type(Path::new("a/app.JS")),
            "text/javascript; charset=utf-8"
        );
        assert_eq!(
            content_type(Path::new("a/x.css")),
            "text/css; charset=utf-8"
        );
        assert_eq!(content_type(Path::new("a/x.json")), "application/json");
        assert_eq!(content_type(Path::new("a/x.svg")), "image/svg+xml");
        assert_eq!(content_type(Path::new("a/x.png")), "image/png");
        assert_eq!(
            content_type(Path::new("a/x.bin")),
            "application/octet-stream"
        );
    }

    #[test]
    fn args_serve_config() {
        assert_eq!(
            parse_args(&sv(&["serve"])),
            Ok(ServeConfig {
                port: 8877,
                root: PathBuf::from(".")
            })
        );
        assert_eq!(
            parse_args(&sv(&["serve", "--port", "9000", "--root", "/tmp/r"])),
            Ok(ServeConfig {
                port: 9000,
                root: PathBuf::from("/tmp/r")
            })
        );
        assert!(parse_args(&sv(&["--port"])).is_err());
        assert!(parse_args(&sv(&["--port", "abc"])).is_err());
        assert!(parse_args(&sv(&["--root"])).is_err());
        assert!(parse_args(&sv(&["frobnicate"])).is_err());
    }

    #[test]
    fn num_and_fixed_formatters() {
        assert_eq!(num(0.05), "0.05");
        assert_eq!(num(210e9), "210000000000.0");
        assert_eq!(num(f64::NAN), "0.0");
        assert_eq!(fixed3(0.0138319 * 1000.0), "13.832");
        assert_eq!(exp3(6.933e-11), "6.933e-11");
    }
}
