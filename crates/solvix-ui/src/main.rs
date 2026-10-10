//! # solvix-ui — the Abaqus-CAE-style web app for Solvix
//!
//! `solvix serve` starts an HTTP server: static UI under `/` and `/static/`,
//! JSON REST API under `/api/*` (issue #13). The same API serves the browser
//! UI (issue #14) and AI agents (issue #15).
//!
//! API v1 contract:
//! - `GET  /api/model?job=NAME` → mesh+model JSON (see issue #13)
//! - `POST /api/run` body `{"path": "examples/cantilever.inp"}` → run
//!   summary JSON (CLI schema)
//! - `GET  /api/results/{job}` → displacements + von Mises arrays
//! - `GET  /api/agent` → machine-readable API description
//!
//! std-only: TcpListener + thread per connection, Connection: close.

use std::collections::HashMap;
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use solvix_core::{ElemKind, Model, Solution, SolveError};

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
        ElemKind::C3d8 => "c3d8",
        ElemKind::C3d4 => "c3d4",
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

struct AppState {
    /// Workspace root: `examples/` decks resolve under it, `/api/run` paths
    /// resolve relative to it.
    root: PathBuf,
    static_dir: PathBuf,
    /// job name → last solved (model, solution).
    results: Mutex<HashMap<String, (Model, Solution)>>,
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

/// Run summary JSON — same schema as `solvix run --json` (issue #3).
fn run_summary_json(
    model: &Model,
    sol: &Solution,
    job: &str,
    solve_ms: u64,
    outputs: &[String],
) -> String {
    let max_disp_m = sol
        .displacements
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
        .fold(0.0, f64::max);
    let max_mises = sol.von_mises.iter().copied().fold(0.0, f64::max);
    let outs = outputs
        .iter()
        .map(|o| format!("\"{}\"", json_escape(o)))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"ok\":true,\"job\":\"{}\",\"nodes\":{},\"elements\":{},\"iterations\":{},\"residual\":{},\"solve_ms\":{},\"max_disp_mm\":{},\"max_von_mises_mpa\":{},\"outputs\":[{}]}}",
        json_escape(job),
        model.mesh.nodes.len(),
        model.mesh.elements.len(),
        sol.iterations,
        exp3(sol.final_residual),
        solve_ms,
        fixed3(max_disp_m * 1000.0),
        fixed3(max_mises / 1e6),
        outs,
    )
}

/// Results JSON — exact shape of `static/mock/results.json`.
fn results_json(job: &str, sol: &Solution) -> String {
    let disp = sol
        .displacements
        .iter()
        .map(|d| format!("[{},{}]", num(d[0]), num(d[1])))
        .collect::<Vec<_>>()
        .join(",");
    let mises = sol
        .von_mises
        .iter()
        .map(|v| num(*v))
        .collect::<Vec<_>>()
        .join(",");
    let max_disp_m = sol
        .displacements
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
        .fold(0.0, f64::max);
    let max_mises = sol.von_mises.iter().copied().fold(0.0, f64::max);
    format!(
        "{{\"ok\":true,\"job\":\"{}\",\"displacements\":[{}],\"von_mises\":[{}],\"max_disp_mm\":{},\"max_von_mises_mpa\":{},\"iterations\":{},\"residual\":{}}}",
        json_escape(job),
        disp,
        mises,
        fixed3(max_disp_m * 1000.0),
        fixed3(max_mises / 1e6),
        sol.iterations,
        exp3(sol.final_residual),
    )
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
            "{\"path\": \"examples/NAME.inp\"}",
            "{\"ok\":true,\"job\":..,\"nodes\":..,\"elements\":..,\"iterations\":..,\"residual\":..,\"solve_ms\":..,\"max_disp_mm\":..,\"max_von_mises_mpa\":..,\"outputs\":[\"NAME.vtu\"]} (CLI run --json schema)",
            "parse, solve and write NAME.vtu next to the deck; caches results for /api/results/{job}",
        ),
        (
            "GET",
            "/api/results/{job}",
            "null",
            "shape of /static/mock/results.json with ok=true",
            "displacements and von Mises arrays of the last POST /api/run of {job}",
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
    let sol = match solvix_core::solve(&model) {
        Ok(s) => s,
        Err(e) => return Response::error(400, "solve", &solve_error_message(e)),
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
    let summary = run_summary_json(&model, &sol, &job, solve_ms, &[vtu.display().to_string()]);
    if let Ok(mut cache) = state.results.lock() {
        cache.insert(job, (model, sol));
    }
    Response::json(200, summary)
}

fn api_results(state: &AppState, job: &str) -> Response {
    let cache = state.results.lock().unwrap_or_else(|p| p.into_inner());
    match cache.get(job) {
        Some((_, sol)) => Response::json(200, results_json(job, sol)),
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

    fn fixture_root(tag: &str) -> PathBuf {
        let root =
            std::env::temp_dir().join(format!("solvix-ui-test-{}-{tag}", std::process::id()));
        let ex = root.join("examples");
        std::fs::create_dir_all(&ex).unwrap();
        std::fs::write(ex.join("tiny.inp"), TINY_DECK).unwrap();
        root
    }

    fn test_state(root: &Path) -> AppState {
        AppState {
            root: root.to_path_buf(),
            static_dir: PathBuf::from("static"),
            results: Mutex::new(HashMap::new()),
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
        assert_eq!(
            results_json("cantilever", &sol),
            r#"{"ok":true,"job":"cantilever","displacements":[[0.0,0.0],[4e-6,-8e-6]],"von_mises":[120000000.0,119000000.0],"max_disp_mm":0.009,"max_von_mises_mpa":120.000,"iterations":390,"residual":6.933e-11}"#
        );
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
