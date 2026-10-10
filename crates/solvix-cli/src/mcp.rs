//! `solvix mcp` — Model Context Protocol server over stdio (issue #8).
//!
//! Speaks the 2024-11-05 MCP subset over newline-delimited JSON-RPC 2.0:
//! one JSON object per request line, one JSON object per response line
//! (no Content-Length framing — matches common MCP stdio client behavior;
//! batch arrays are not handled). stdout carries ONLY responses; any
//! server-side diagnostics go to stderr. Parse/solve failures inside a
//! tool call become `isError:true` tool results — the server never
//! panics and keeps serving. EOF exits 0.
//!
//! Exposes the CLI JSON layer as four tools: solvix_info, solvix_run,
//! solvix_results (solve cached by path), solvix_bench.

use std::io::{BufRead, Write};

use solvix_core::Solution;

use crate::{
    bench, error_json, fin, fmt_mm, fmt_mpa, info_json, info_report, job_name, json_escape,
    parse_deck, read_deck, run_summary, solve_auto, solve_error_message, success_json, Stage,
};

// ---------------------------------------------------------------------------
// Minimal JSON value + parser (zero-dep workspace policy).
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(fields) => fields.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    #[cfg(test)]
    fn at(&self, i: usize) -> Option<&Json> {
        match self {
            Json::Arr(a) => a.get(i),
            _ => None,
        }
    }

    fn parse(s: &str) -> Result<Json, &'static str> {
        let mut p = Parser {
            b: s.as_bytes(),
            i: 0,
        };
        let v = p.value(0)?;
        p.ws();
        if p.i != p.b.len() {
            return Err("trailing characters after JSON value");
        }
        Ok(v)
    }

    fn to_json_string(&self) -> String {
        match self {
            Json::Null => "null".into(),
            Json::Bool(b) => b.to_string(),
            Json::Num(n) => format!("{}", fin(*n)),
            Json::Str(s) => format!("\"{}\"", json_escape(s)),
            Json::Arr(a) => format!(
                "[{}]",
                a.iter()
                    .map(|v| v.to_json_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            Json::Obj(o) => format!(
                "{{{}}}",
                o.iter()
                    .map(|(k, v)| format!("\"{}\":{}", json_escape(k), v.to_json_string()))
                    .collect::<Vec<_>>()
                    .join(",")
            ),
        }
    }
}

struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Parser<'a> {
    fn ws(&mut self) {
        while self.i < self.b.len() && matches!(self.b[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }

    fn value(&mut self, depth: usize) -> Result<Json, &'static str> {
        if depth > 64 {
            return Err("JSON nesting too deep");
        }
        self.ws();
        match self.b.get(self.i) {
            None => Err("unexpected end of input"),
            Some(b'{') => self.object(depth),
            Some(b'[') => self.array(depth),
            Some(b'"') => Ok(Json::Str(self.string()?)),
            Some(b't') => self.literal("true", Json::Bool(true)),
            Some(b'f') => self.literal("false", Json::Bool(false)),
            Some(b'n') => self.literal("null", Json::Null),
            Some(c) if c.is_ascii_digit() || *c == b'-' => self.number(),
            Some(_) => Err("unexpected character"),
        }
    }

    fn literal(&mut self, lit: &str, v: Json) -> Result<Json, &'static str> {
        if self.b[self.i..].starts_with(lit.as_bytes()) {
            self.i += lit.len();
            Ok(v)
        } else {
            Err("invalid literal")
        }
    }

    fn object(&mut self, depth: usize) -> Result<Json, &'static str> {
        self.i += 1; // '{'
        let mut fields = Vec::new();
        self.ws();
        if self.b.get(self.i) == Some(&b'}') {
            self.i += 1;
            return Ok(Json::Obj(fields));
        }
        loop {
            self.ws();
            if self.b.get(self.i) != Some(&b'"') {
                return Err("expected object key string");
            }
            let key = self.string()?;
            self.ws();
            if self.b.get(self.i) != Some(&b':') {
                return Err("expected ':' in object");
            }
            self.i += 1;
            let val = self.value(depth + 1)?;
            fields.push((key, val));
            self.ws();
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b'}') => {
                    self.i += 1;
                    return Ok(Json::Obj(fields));
                }
                _ => return Err("expected ',' or '}' in object"),
            }
        }
    }

    fn array(&mut self, depth: usize) -> Result<Json, &'static str> {
        self.i += 1; // '['
        let mut items = Vec::new();
        self.ws();
        if self.b.get(self.i) == Some(&b']') {
            self.i += 1;
            return Ok(Json::Arr(items));
        }
        loop {
            items.push(self.value(depth + 1)?);
            self.ws();
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b']') => {
                    self.i += 1;
                    return Ok(Json::Arr(items));
                }
                _ => return Err("expected ',' or ']' in array"),
            }
        }
    }

    fn string(&mut self) -> Result<String, &'static str> {
        self.i += 1; // '"'
        let mut out = Vec::new();
        loop {
            let Some(&c) = self.b.get(self.i) else {
                return Err("unterminated string");
            };
            match c {
                b'"' => {
                    self.i += 1;
                    return String::from_utf8(out).map_err(|_| "invalid utf-8 in string");
                }
                b'\\' => {
                    self.i += 1;
                    let Some(&e) = self.b.get(self.i) else {
                        return Err("unterminated escape");
                    };
                    self.i += 1;
                    match e {
                        b'"' => out.push(b'"'),
                        b'\\' => out.push(b'\\'),
                        b'/' => out.push(b'/'),
                        b'b' => out.push(0x08),
                        b'f' => out.push(0x0C),
                        b'n' => out.push(b'\n'),
                        b'r' => out.push(b'\r'),
                        b't' => out.push(b'\t'),
                        b'u' => {
                            let hi = self.hex4()?;
                            let cp = if (0xD800..=0xDBFF).contains(&hi) {
                                if self.b.get(self.i) != Some(&b'\\')
                                    || self.b.get(self.i + 1) != Some(&b'u')
                                {
                                    return Err("lone high surrogate");
                                }
                                self.i += 2;
                                let lo = self.hex4()?;
                                if !(0xDC00..=0xDFFF).contains(&lo) {
                                    return Err("invalid low surrogate");
                                }
                                0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00)
                            } else if (0xDC00..=0xDFFF).contains(&hi) {
                                return Err("lone low surrogate");
                            } else {
                                hi
                            };
                            let ch = char::from_u32(cp).ok_or("invalid unicode escape")?;
                            let mut buf = [0u8; 4];
                            out.extend_from_slice(ch.encode_utf8(&mut buf).as_bytes());
                        }
                        _ => return Err("invalid escape"),
                    }
                }
                c if c < 0x20 => return Err("unescaped control character in string"),
                c => {
                    out.push(c);
                    self.i += 1;
                }
            }
        }
    }

    fn hex4(&mut self) -> Result<u32, &'static str> {
        let mut v = 0u32;
        for _ in 0..4 {
            let Some(&c) = self.b.get(self.i) else {
                return Err("truncated \\u escape");
            };
            let d = (c as char)
                .to_digit(16)
                .ok_or("invalid hex in \\u escape")?;
            v = v * 16 + d;
            self.i += 1;
        }
        Ok(v)
    }

    fn number(&mut self) -> Result<Json, &'static str> {
        let start = self.i;
        if self.b.get(self.i) == Some(&b'-') {
            self.i += 1;
        }
        match self.b.get(self.i) {
            Some(b'0') => self.i += 1,
            Some(c) if c.is_ascii_digit() => {
                while self.b.get(self.i).is_some_and(|c| c.is_ascii_digit()) {
                    self.i += 1;
                }
            }
            _ => return Err("invalid number"),
        }
        if self.b.get(self.i) == Some(&b'.') {
            self.i += 1;
            if !self.b.get(self.i).is_some_and(|c| c.is_ascii_digit()) {
                return Err("invalid number fraction");
            }
            while self.b.get(self.i).is_some_and(|c| c.is_ascii_digit()) {
                self.i += 1;
            }
        }
        if matches!(self.b.get(self.i), Some(b'e') | Some(b'E')) {
            self.i += 1;
            if matches!(self.b.get(self.i), Some(b'+') | Some(b'-')) {
                self.i += 1;
            }
            if !self.b.get(self.i).is_some_and(|c| c.is_ascii_digit()) {
                return Err("invalid number exponent");
            }
            while self.b.get(self.i).is_some_and(|c| c.is_ascii_digit()) {
                self.i += 1;
            }
        }
        let s = std::str::from_utf8(&self.b[start..self.i]).map_err(|_| "invalid number")?;
        s.parse::<f64>()
            .map(Json::Num)
            .map_err(|_| "invalid number")
    }
}

// ---------------------------------------------------------------------------
// MCP server
// ---------------------------------------------------------------------------

struct CachedSolve {
    node_ids: Vec<u32>,
    solution: Solution,
}

pub struct McpState {
    /// ponytail: unbounded growth guard — small FIFO of solved decks by path.
    cache: Vec<(String, CachedSolve)>,
}

impl McpState {
    pub fn new() -> McpState {
        McpState { cache: Vec::new() }
    }

    fn cached(&mut self, path: &str) -> Option<&CachedSolve> {
        self.cache.iter().find(|(p, _)| p == path).map(|(_, c)| c)
    }

    fn insert(&mut self, path: &str, c: CachedSolve) {
        if self.cache.len() >= 4 {
            self.cache.remove(0);
        }
        self.cache.push((path.to_string(), c));
    }
}

impl Default for McpState {
    fn default() -> Self {
        Self::new()
    }
}

const TOOL_INFO: &str = r#"{"name":"solvix_info","description":"Parse a Solvix .inp deck and report mesh, BC and load statistics without solving. Units: SI in (m, N, Pa).","inputSchema":{"type":"object","properties":{"path":{"type":"string","description":"Path to the .inp deck (SI units: m, N, Pa)"}},"required":["path"]}}"#;
const TOOL_RUN: &str = r#"{"name":"solvix_run","description":"Parse, solve and write <jobname>.vtu for a .inp deck; returns the run summary with solver_method (pcg or ldlt) and cond_estimate (a diagonal conditioning indicator). The analysis type is AUTO-DETECTED from the deck (no flags): a *DYNAMIC step runs transient Newmark dynamics — JSON gains \"dynamic\":true, steps, dt, duration, final_disp_mm, and max_disp_mm is the max over all time points; step temperature data (*HEAT TRANSFER, *TEMPERATURE, *DFLUX) without dynamics runs a thermal-stress pass — JSON gains \"thermal\":true and max_temperature_c; anything else is linear statics. Optional \"modes\": N also extracts the N lowest natural frequencies (modal analysis; JSON gains modes and frequencies_hz). Optional \"odb\": true (the CLI's `run --odb` flag) also writes <jobname>.solvix, a binary output database with result frames and per-node |u| history — inspect it with `solvix odb <file> --inspect` / `--frame N`. Units: SI in (m, N, Pa); mm and MPa out.","inputSchema":{"type":"object","properties":{"path":{"type":"string","description":"Path to the .inp deck (SI units: m, N, Pa)"},"outdir":{"type":"string","description":"Output directory for the .vtu file (default: working directory)"},"modes":{"type":"integer","minimum":1,"description":"Also extract the N lowest natural frequencies (modal analysis)"},"odb":{"type":"boolean","description":"Also write <jobname>.solvix, a binary output database (run --odb)"}},"required":["path"]}}"#;
const TOOL_RESULTS: &str = r#"{"name":"solvix_results","description":"Solve a .inp deck (reuses a cached solution per path) and return per-node displacements plus nodal and element von Mises arrays. Units: SI in (m, N, Pa); mm and MPa out.","inputSchema":{"type":"object","properties":{"path":{"type":"string","description":"Path to the .inp deck (SI units: m, N, Pa)"}},"required":["path"]}}"#;
const TOOL_BENCH: &str = r#"{"name":"solvix_bench","description":"Run the built-in verification suite (cantilever40x4, kirsch, patch_cps4) and report values vs analytical references with pass/fail gates. Metrics in mm / MPa.","inputSchema":{"type":"object","properties":{},"required":[]}}"#;

fn rpc_result(id: &Json, result: &str) -> String {
    format!(
        r#"{{"jsonrpc":"2.0","id":{},"result":{}}}"#,
        id.to_json_string(),
        result
    )
}

fn rpc_error(id: Json, code: i32, msg: &str) -> String {
    format!(
        r#"{{"jsonrpc":"2.0","id":{},"error":{{"code":{},"message":"{}"}}}}"#,
        id.to_json_string(),
        code,
        json_escape(msg)
    )
}

fn text_envelope(text: &str, is_error: bool) -> String {
    format!(
        r#"{{"content":[{{"type":"text","text":"{}"}}],"isError":{}}}"#,
        json_escape(text),
        is_error
    )
}

fn arg_str<'a>(args: &'a Json, key: &str) -> Result<&'a str, (Stage, String)> {
    args.get(key)
        .and_then(|v| v.as_str())
        .ok_or_else(|| (Stage::Usage, format!("missing required argument: {key}")))
}

fn tool_info(args: &Json) -> Result<String, (Stage, String)> {
    let path = arg_str(args, "path")?;
    let text = read_deck(path)?;
    let (model, report) = parse_deck(&text)?;
    Ok(info_json(&info_report(&model, &job_name(path), &report)))
}

/// Optional non-negative integer argument (validated further by
/// `solve_modal`, which rejects `num_modes` outside `1..free_dof_count`);
/// `None` when absent or not a JSON number — the run stays static-only.
fn arg_modes(args: &Json) -> Option<usize> {
    let n = match args.get("modes") {
        Some(Json::Num(n)) => *n,
        _ => return None,
    };
    if n.fract() == 0.0 && n >= 0.0 && n <= u32::MAX as f64 {
        Some(n as usize)
    } else {
        None
    }
}

fn tool_run(args: &Json) -> Result<String, (Stage, String)> {
    let path = arg_str(args, "path")?;
    let outdir = args
        .get("outdir")
        .and_then(|v| v.as_str())
        .unwrap_or(".")
        .to_string();
    let text = read_deck(path)?;
    let (model, _) = parse_deck(&text)?;
    let t0 = std::time::Instant::now();
    let outcome = solve_auto(&model).map_err(|e| (Stage::Solve, solve_error_message(e)))?;
    let frequencies_hz = match arg_modes(args) {
        Some(n) => {
            solvix_core::solve_modal(&model, &solvix_core::ModalOptions { num_modes: n })
                .map_err(|e| (Stage::Solve, solve_error_message(e)))?
                .frequencies_hz
        }
        None => Vec::new(),
    };
    let solve_ms = t0.elapsed().as_millis() as u64;
    std::fs::create_dir_all(&outdir)
        .map_err(|e| (Stage::Io, format!("cannot create directory {outdir}: {e}")))?;
    let vtu = std::path::Path::new(&outdir).join(format!("{}.vtu", job_name(path)));
    let sol_for_vtu = match &outcome {
        crate::AutoOutcome::Static { sol, .. } => sol,
        crate::AutoOutcome::Dynamic { final_sol, .. } => final_sol,
    };
    solvix_io::write_vtu(&vtu, &model, sol_for_vtu)
        .map_err(|e| (Stage::Io, format!("cannot write {}: {e}", vtu.display())))?;
    let mut outputs = vec![vtu.display().to_string()];
    if matches!(args.get("odb"), Some(Json::Bool(true))) {
        // Issue #32: binary output database alongside the VTU.
        let odb = std::path::Path::new(&outdir).join(format!("{}.solvix", job_name(path)));
        crate::write_job_odb(&odb, &model, &outcome)?;
        outputs.push(odb.display().to_string());
    }
    let s = run_summary(
        &model,
        &outcome,
        &job_name(path),
        solve_ms,
        outputs,
        frequencies_hz,
    );
    Ok(success_json(&s))
}

fn results_json(job: &str, node_ids: &[u32], sol: &Solution) -> String {
    let ids = node_ids
        .iter()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let disps = sol
        .displacements
        .iter()
        .map(|d| format!("[{},{},{}]", fmt_mm(d[0]), fmt_mm(d[1]), fmt_mm(d[2])))
        .collect::<Vec<_>>()
        .join(",");
    let nodal = sol
        .nodal_von_mises
        .iter()
        .map(|v| fmt_mpa(*v))
        .collect::<Vec<_>>()
        .join(",");
    let elems = sol
        .von_mises
        .iter()
        .map(|v| fmt_mpa(*v))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        r#"{{"ok":true,"job":"{}","nodes":{},"node_ids":[{}],"displacements_mm":[{}],"nodal_von_mises_mpa":[{}],"element_von_mises_mpa":[{}]}}"#,
        json_escape(job),
        node_ids.len(),
        ids,
        disps,
        nodal,
        elems,
    )
}

fn tool_results(state: &mut McpState, args: &Json) -> Result<String, (Stage, String)> {
    let path = arg_str(args, "path")?;
    if let Some(hit) = state.cached(path) {
        return Ok(results_json(&job_name(path), &hit.node_ids, &hit.solution));
    }
    let text = read_deck(path)?;
    let (model, _) = parse_deck(&text)?;
    let sol = solvix_core::solve(&model).map_err(|e| (Stage::Solve, solve_error_message(e)))?;
    let node_ids = model.mesh.nodes.iter().map(|n| n.id).collect::<Vec<_>>();
    let out = results_json(&job_name(path), &node_ids, &sol);
    state.insert(
        path,
        CachedSolve {
            node_ids,
            solution: sol,
        },
    );
    Ok(out)
}

fn tool_bench() -> Result<String, (Stage, String)> {
    let mut outcomes = Vec::with_capacity(bench::CASES.len());
    for def in bench::CASES {
        outcomes.push(
            bench::run_case(def)
                .map_err(|m| (Stage::Solve, format!("bench case '{}': {m}", def.name)))?,
        );
    }
    Ok(bench::json_line(&outcomes))
}

fn tools_call(state: &mut McpState, params: &Json) -> String {
    let Some(name) = params.get("name").and_then(|v| v.as_str()) else {
        return rpc_error(Json::Null, -32602, "invalid params: missing tool name");
    };
    let args = params
        .get("arguments")
        .cloned()
        .unwrap_or(Json::Obj(vec![]));
    let outcome = match name {
        "solvix_info" => tool_info(&args),
        "solvix_run" => tool_run(&args),
        "solvix_results" => tool_results(state, &args),
        "solvix_bench" => tool_bench(),
        other => Err((Stage::Usage, format!("unknown tool: {other}"))),
    };
    match outcome {
        Ok(payload) => text_envelope(&payload, false),
        Err((stage, msg)) => text_envelope(&error_json(stage, &msg), true),
    }
}

/// Pure message handler: one request line in, one response line out.
/// `None` means "produce no output" (notification, empty line).
pub fn handle_message(state: &mut McpState, line: &str) -> Option<String> {
    let line = line.trim();
    if line.is_empty() {
        return None;
    }
    let msg = match Json::parse(line) {
        Ok(v) => v,
        Err(e) => return Some(rpc_error(Json::Null, -32700, &format!("parse error: {e}"))),
    };
    if !matches!(msg, Json::Obj(_)) {
        return Some(rpc_error(
            Json::Null,
            -32600,
            "invalid request: expected object",
        ));
    }
    let id = msg.get("id").cloned();
    let Some(id) = id else {
        return None; // notification (no id): never respond
    };
    let Some(method) = msg.get("method").and_then(|m| m.as_str()) else {
        return Some(rpc_error(id, -32600, "invalid request: missing method"));
    };
    let params = msg.get("params").cloned().unwrap_or(Json::Obj(vec![]));
    let result = match method {
        "initialize" => format!(
            r#"{{"protocolVersion":"2024-11-05","capabilities":{{"tools":{{}}}},"serverInfo":{{"name":"solvix","version":"{}"}}}}"#,
            env!("CARGO_PKG_VERSION")
        ),
        "initialized" => "null".to_string(),
        "ping" => "{}".to_string(),
        "tools/list" => {
            format!(r#"{{"tools":[{TOOL_INFO},{TOOL_RUN},{TOOL_RESULTS},{TOOL_BENCH}]}}"#)
        }
        "tools/call" => tools_call(state, &params),
        other => {
            return Some(rpc_error(id, -32601, &format!("method not found: {other}")));
        }
    };
    Some(rpc_result(&id, &result))
}

/// stdio main loop: line in, response line out, flush, EOF exits 0.
pub fn serve() -> i32 {
    let mut state = McpState::new();
    let stdout = std::io::stdout();
    let mut w = stdout.lock();
    for line in std::io::stdin().lock().lines() {
        let Ok(line) = line else { break };
        if let Some(resp) = handle_message(&mut state, &line) {
            if writeln!(w, "{resp}").is_err() {
                break; // stdout closed
            }
            let _ = w.flush();
        }
    }
    0
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Tiny 2-element CPS4 cantilever: 6 nodes, left edge clamped,
    /// -100 N on each right-edge node (dof 2). SI units.
    fn tiny_deck() -> String {
        "\
*HEADING
    mcp test tiny cantilever
*NODE
1, 0.0, 0.0
2, 0.05, 0.0
3, 0.1, 0.0
4, 0.0, 0.05
5, 0.05, 0.05
6, 0.1, 0.05
*ELEMENT, TYPE=CPS4, ELSET=ALL
1, 1, 2, 5, 4
2, 2, 3, 6, 5
*ELSET, ELSET=ALL, GENERATE
1, 2, 1
*MATERIAL, NAME=STEEL
*ELASTIC
210e9, 0.3
*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL
0.01
*STEP
*STATIC
*BOUNDARY
1, 1, 2
4, 1, 2
*CLOAD
3, 2, -100.0
6, 2, -100.0
*END STEP
"
        .to_string()
    }

    fn tmp_path(name: &str) -> String {
        std::env::temp_dir()
            .join(format!("solvix-mcp-test-{}-{name}.inp", std::process::id()))
            .display()
            .to_string()
    }

    fn tmp_deck(name: &str) -> String {
        let p = tmp_path(name);
        std::fs::write(&p, tiny_deck()).unwrap();
        p
    }

    /// Tiny thermal deck: 2x1 CPS4 bar, T=100 at x=0, T=0 at x=0.1,
    /// clamped left edge (issue #46 auto-detection).
    fn thermal_deck() -> String {
        "\
*HEADING
    mcp thermal fixture
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
"
        .to_string()
    }

    /// Tiny dynamics deck: same mesh, *DYNAMIC 1e-3 over 0.01 s
    /// (issue #47 auto-detection).
    fn dynamic_deck() -> String {
        "\
*HEADING
    mcp dynamics fixture
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
1.0E-3, 0.01
*BOUNDARY
1, 1, 2
4, 1, 2
*CLOAD
3, 2, -100.0
6, 2, -100.0
*END STEP
"
        .to_string()
    }

    fn tmp_deck_with(name: &str, deck: &str) -> String {
        let p = tmp_path(name);
        std::fs::write(&p, deck).unwrap();
        p
    }

    fn call_tool(state: &mut McpState, name: &str, args_json: &str) -> String {
        handle_message(
            state,
            &format!(
                r#"{{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{{"name":"{name}","arguments":{args_json}}}}}"#
            ),
        )
        .unwrap()
    }

    fn text_of(response: &str) -> String {
        let v = Json::parse(response).unwrap();
        let text = v
            .get("result")
            .unwrap()
            .get("content")
            .unwrap()
            .at(0)
            .and_then(|c| c.get("text"))
            .and_then(|t| t.as_str())
            .unwrap()
            .to_string();
        Json::parse(&text).unwrap().to_json_string()
    }

    // JSON parser ---------------------------------------------------------

    #[test]
    fn json_parse_roundtrip() {
        let v = Json::parse(r#"{"a":[1,-2.5e3,true,false,null],"s":"he said \"hi\"\n"}"#).unwrap();
        assert_eq!(
            v.to_json_string(),
            r#"{"a":[1,-2500,true,false,null],"s":"he said \"hi\"\n"}"#
        );
        assert_eq!(v.get("s").unwrap().as_str().unwrap(), "he said \"hi\"\n");
    }

    #[test]
    fn json_parse_unicode_escapes() {
        let v = Json::parse(r#""A😀""#).unwrap();
        assert_eq!(v.as_str().unwrap(), "A😀");
    }

    #[test]
    fn json_parse_rejects_garbage() {
        for bad in [
            "",
            "{",
            "{\"a\":}",
            "[1,]",
            "{\"a\" 1}",
            "tru",
            "\"\\q\"",
            "01",
            "1.",
            "{\"a\":1} x",
            "\"\u{1}\"",
        ] {
            assert!(Json::parse(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn json_parse_rejects_deep_nesting() {
        let deep = "[".repeat(100) + &"]".repeat(100);
        assert_eq!(Json::parse(&deep), Err("JSON nesting too deep"));
    }

    // initialize / notifications / errors ----------------------------------

    #[test]
    fn initialize_round_trip() {
        let mut s = McpState::new();
        let resp = handle_message(
            &mut s,
            r#"{"jsonrpc":"2.0","id":7,"method":"initialize","params":{}}"#,
        )
        .unwrap();
        assert!(
            resp.starts_with(r#"{"jsonrpc":"2.0","id":7,"result":{"protocolVersion":"2024-11-05""#),
            "{resp}"
        );
        assert!(resp.contains(r#""capabilities":{"tools":{}}"#), "{resp}");
        assert!(
            resp.contains(&format!(
                r#""serverInfo":{{"name":"solvix","version":"{}"}}"#,
                env!("CARGO_PKG_VERSION")
            )),
            "{resp}"
        );
    }

    #[test]
    fn string_id_is_escaped_and_echoed() {
        let mut s = McpState::new();
        let resp = handle_message(
            &mut s,
            r#"{"jsonrpc":"2.0","id":"a\"b","method":"initialize"}"#,
        )
        .unwrap();
        assert!(resp.contains(r#""id":"a\"b""#), "{resp}");
    }

    #[test]
    fn initialized_notification_no_response() {
        let mut s = McpState::new();
        assert_eq!(
            handle_message(&mut s, r#"{"jsonrpc":"2.0","method":"initialized"}"#),
            None
        );
    }

    #[test]
    fn any_notification_gets_no_response() {
        let mut s = McpState::new();
        assert_eq!(
            handle_message(
                &mut s,
                r#"{"jsonrpc":"2.0","method":"notifications/whatever"}"#
            ),
            None
        );
        assert_eq!(handle_message(&mut s, "   "), None);
    }

    #[test]
    fn malformed_json_is_32700_with_null_id() {
        let mut s = McpState::new();
        let resp = handle_message(&mut s, r#"{"jsonrpc":"2.0","id":5,"method":"#).unwrap();
        assert!(
            resp.contains(r#""id":null,"error":{"code":-32700"#),
            "{resp}"
        );
    }

    #[test]
    fn unknown_method_is_32601() {
        let mut s = McpState::new();
        let resp =
            handle_message(&mut s, r#"{"jsonrpc":"2.0","id":2,"method":"frobnicate"}"#).unwrap();
        assert!(
            resp.contains(r#""error":{"code":-32601,"message":"method not found: frobnicate"}"#),
            "{resp}"
        );
    }

    #[test]
    fn non_object_request_is_32600() {
        let mut s = McpState::new();
        let resp = handle_message(&mut s, r#"[1,2]"#).unwrap();
        assert!(resp.contains("-32600"), "{resp}");
    }

    // tools/list -----------------------------------------------------------

    #[test]
    fn tools_list_shape() {
        let mut s = McpState::new();
        let resp =
            handle_message(&mut s, r#"{"jsonrpc":"2.0","id":1,"method":"tools/list"}"#).unwrap();
        for tool in [
            "solvix_info",
            "solvix_run",
            "solvix_results",
            "solvix_bench",
        ] {
            assert!(
                resp.contains(&format!(r#""name":"{tool}""#)),
                "{tool}: {resp}"
            );
        }
        assert_eq!(resp.matches("\"required\":").count(), 4, "{resp}");
        assert_eq!(resp.matches("\"inputSchema\":").count(), 4, "{resp}");
        assert!(resp.contains(r#""required":["path"]"#), "{resp}");
        assert!(resp.contains("mm and MPa out"), "{resp}");
        assert!(resp.contains("solver_method"), "{resp}");
        assert!(resp.contains("cond_estimate"), "{resp}");
    }

    // tools/call -----------------------------------------------------------

    #[test]
    fn call_info_on_fixture_deck() {
        let mut s = McpState::new();
        let p = tmp_deck("info");
        let resp = call_tool(&mut s, "solvix_info", &format!(r#"{{"path":"{p}"}}"#));
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains(r#""ok":true"#), "{payload}");
        assert!(payload.contains(r#""nodes":6"#), "{payload}");
        assert!(payload.contains(r#""elements":2"#), "{payload}");
        assert!(payload.contains(r#""warnings":"#), "{payload}");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn call_run_returns_real_numbers() {
        let mut s = McpState::new();
        let p = tmp_deck("run");
        let outdir = std::env::temp_dir()
            .join(format!("solvix-mcp-test-{}-run", std::process::id()))
            .display()
            .to_string();
        let resp = call_tool(
            &mut s,
            "solvix_run",
            &format!(r#"{{"path":"{p}","outdir":"{outdir}"}}"#),
        );
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains(r#""ok":true"#), "{payload}");
        assert!(payload.contains(r#""solver_method":"#), "{payload}");
        assert!(payload.contains(r#""cond_estimate":"#), "{payload}");
        assert!(payload.contains(r#""max_disp_mm":"#), "{payload}");
        assert!(!payload.contains(r#""max_disp_mm":0.000"#), "{payload}");
        assert!(
            !payload.contains(r#""max_von_mises_mpa":0.000"#),
            "{payload}"
        );
        assert!(payload.contains(".vtu"), "{payload}");
        let stem = std::path::Path::new(&p)
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();
        assert!(
            std::path::Path::new(&outdir)
                .join(format!("{stem}.vtu"))
                .exists(),
            "{payload}"
        );
        std::fs::remove_file(p).ok();
        std::fs::remove_dir_all(outdir).ok();
    }

    #[test]
    fn call_run_with_modes_returns_frequencies() {
        let mut s = McpState::new();
        let p = tmp_deck("runmodes");
        let outdir = std::env::temp_dir()
            .join(format!("solvix-mcp-test-{}-runmodes", std::process::id()))
            .display()
            .to_string();
        let resp = call_tool(
            &mut s,
            "solvix_run",
            &format!(r#"{{"path":"{p}","outdir":"{outdir}","modes":2}}"#),
        );
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains(r#""modes":2,"#), "{payload}");
        let start = payload
            .find(r#""frequencies_hz":["#)
            .expect("frequencies present");
        let freq_json = &payload[start..];
        let freqs = freq_json
            .trim_start_matches(r#""frequencies_hz":["#)
            .split([',', ']'])
            .filter_map(|t| t.parse::<f64>().ok())
            .collect::<Vec<_>>();
        assert_eq!(freqs.len(), 2, "{payload}");
        assert!(freqs[0] > 0.0 && freqs[1] > freqs[0], "{payload}");
        assert!(payload.ends_with(r#"}"#), "{payload}");
        std::fs::remove_file(p).ok();
        std::fs::remove_dir_all(outdir).ok();
    }

    #[test]
    fn call_run_modes_garbage_is_static_only() {
        let mut s = McpState::new();
        let p = tmp_deck("runmodes2");
        let resp = call_tool(
            &mut s,
            "solvix_run",
            &format!(r#"{{"path":"{p}","outdir":".","modes":"many"}}"#),
        );
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(!payload.contains("frequencies_hz"), "{payload}");
        std::fs::remove_file(&p).ok();
        std::fs::remove_file(format!("{}.vtu", p.trim_end_matches(".inp"))).ok();
    }

    #[test]
    fn call_run_modes_zero_is_solve_error() {
        let mut s = McpState::new();
        let p = tmp_deck("runmodes3");
        let resp = call_tool(
            &mut s,
            "solvix_run",
            &format!(r#"{{"path":"{p}","outdir":".","modes":0}}"#),
        );
        assert!(resp.contains(r#""isError":true"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains("num_modes"), "{payload}");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn call_results_arrays_and_cache() {
        let mut s = McpState::new();
        let p = tmp_deck("results");
        let resp = call_tool(&mut s, "solvix_results", &format!(r#"{{"path":"{p}"}}"#));
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains(r#""node_ids":[1,2,3,4,5,6]"#), "{payload}");
        assert!(payload.contains(r#""displacements_mm":["#), "{payload}");
        assert!(payload.contains(r#""nodal_von_mises_mpa":["#), "{payload}");
        assert!(
            payload.contains(r#""element_von_mises_mpa":["#),
            "{payload}"
        );
        // second call hits the cache and returns identical payload
        let resp2 = call_tool(&mut s, "solvix_results", &format!(r#"{{"path":"{p}"}}"#));
        assert_eq!(text_of(&resp), text_of(&resp2));
        assert_eq!(s.cache.len(), 1);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn call_run_thermal_deck_reports_temperature() {
        let mut s = McpState::new();
        let p = tmp_deck_with("runthermal", &thermal_deck());
        let outdir = std::env::temp_dir()
            .join(format!("solvix-mcp-test-{}-runthermal", std::process::id()))
            .display()
            .to_string();
        let resp = call_tool(
            &mut s,
            "solvix_run",
            &format!(r#"{{"path":"{p}","outdir":"{outdir}"}}"#),
        );
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(
            payload.contains(r#""thermal":true,"max_temperature_c":100"#),
            "{payload}"
        );
        assert!(!payload.contains(r#""dynamic""#), "{payload}");
        std::fs::remove_file(p).ok();
        std::fs::remove_dir_all(outdir).ok();
    }

    #[test]
    fn call_run_dynamic_deck_reports_transient_fields() {
        let mut s = McpState::new();
        let p = tmp_deck_with("rundyn", &dynamic_deck());
        let outdir = std::env::temp_dir()
            .join(format!("solvix-mcp-test-{}-rundyn", std::process::id()))
            .display()
            .to_string();
        let resp = call_tool(
            &mut s,
            "solvix_run",
            &format!(r#"{{"path":"{p}","outdir":"{outdir}"}}"#),
        );
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(
            payload.contains(
                r#""dynamic":true,"steps":10,"dt":0.001,"duration":0.01,"final_disp_mm":"#
            ),
            "{payload}"
        );
        assert!(!payload.contains(r#""thermal""#), "{payload}");
        std::fs::remove_file(p).ok();
        std::fs::remove_dir_all(outdir).ok();
    }

    #[test]
    fn call_bench_ok() {
        let mut s = McpState::new();
        let resp = call_tool(&mut s, "solvix_bench", "{}");
        assert!(resp.contains(r#""isError":false"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains(r#""ok":true"#), "{payload}");
        assert!(payload.contains("cantilever40x4"), "{payload}");
        assert!(payload.contains("kirsch"), "{payload}");
        assert!(payload.contains("patch_cps4"), "{payload}");
    }

    // error paths ----------------------------------------------------------

    #[test]
    fn call_missing_file_is_error_not_panic() {
        let mut s = McpState::new();
        let resp = call_tool(&mut s, "solvix_info", r#"{"path":"/no/such/deck.inp"}"#);
        assert!(resp.contains(r#""isError":true"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains(r#""ok":false"#), "{payload}");
        assert!(payload.contains(r#""stage":"io""#), "{payload}");
    }

    #[test]
    fn call_missing_argument_is_error() {
        let mut s = McpState::new();
        let resp = call_tool(&mut s, "solvix_info", "{}");
        assert!(resp.contains(r#""isError":true"#), "{resp}");
        assert!(resp.contains("missing required argument: path"), "{resp}");
    }

    #[test]
    fn call_parse_error_is_error() {
        let mut s = McpState::new();
        let p = tmp_path("broken");
        std::fs::write(&p, "*NODE\n1, 0.0\n*GARBAGE_KEYWORD\n").unwrap();
        let resp = call_tool(&mut s, "solvix_info", &format!(r#"{{"path":"{p}"}}"#));
        assert!(resp.contains(r#""isError":true"#), "{resp}");
        let payload = text_of(&resp);
        assert!(payload.contains(r#""stage":"parse""#), "{payload}");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn call_unknown_tool_is_error() {
        let mut s = McpState::new();
        let resp = call_tool(&mut s, "solvix_nope", "{}");
        assert!(resp.contains(r#""isError":true"#), "{resp}");
        assert!(resp.contains("unknown tool: solvix_nope"), "{resp}");
    }

    #[test]
    fn call_missing_name_is_32602() {
        let mut s = McpState::new();
        let resp = handle_message(
            &mut s,
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"arguments":{}}}"#,
        )
        .unwrap();
        assert!(resp.contains("-32602"), "{resp}");
    }

    #[test]
    fn server_keeps_serving_after_bad_line() {
        let mut s = McpState::new();
        let bad = handle_message(&mut s, "not json").unwrap();
        assert!(bad.contains("-32700"), "{bad}");
        let good = handle_message(&mut s, r#"{"jsonrpc":"2.0","id":9,"method":"ping"}"#).unwrap();
        assert!(good.contains(r#""id":9,"result":{}"#), "{good}");
    }
}
