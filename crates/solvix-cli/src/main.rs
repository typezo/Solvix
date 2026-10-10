//! # solvix — the agent-era FEA CLI
//!
//! Implementation owned by issue #3; `bench` by issue #17. Subcommands:
//!   solvix run <job.inp> [-o outdir] [--json]
//!   solvix info <job.inp> [--json]
//!   solvix bench [--json]
//!
//! Exit codes: 0 ok · 1 usage/IO · 2 parse error · 3 solve error.
//! stdout carries the payload (human table or single-line JSON), stderr
//! carries warnings — designed for AI agents.

use std::path::{Path, PathBuf};
use std::time::Instant;

use solvix_core::{ElemKind, Model, Solution, SolveError};
use solvix_io::ParseReport;

mod bench;

const USAGE: &str = "\
solvix 0.1.0 — the agent-era FEA CLI

solvix is designed for AI agents: all subcommands accept --json and exit 0/1/2/3.

Usage:
  solvix run <job.inp> [-o outdir] [--json]   parse, solve, write <jobname>.vtu
  solvix info <job.inp> [--json]              parse-only model report, never solves
  solvix bench [--json]                       built-in verification suite (all-pass exit 0, else 3)
  solvix help | --help | -h                   show this text

Exit codes: 0 ok · 1 usage/IO · 2 parse error · 3 solve error.
Units: SI in (m, N, Pa) · mm and MPa out for display.
Errors with --json: {\"ok\":false,\"stage\":\"parse|solve|io|usage\",\"error\":\"...\"} on stdout.
";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Stage {
    Usage,
    Io,
    Parse,
    Solve,
}

impl Stage {
    fn as_str(self) -> &'static str {
        match self {
            Stage::Usage => "usage",
            Stage::Io => "io",
            Stage::Parse => "parse",
            Stage::Solve => "solve",
        }
    }

    fn exit_code(self) -> i32 {
        match self {
            Stage::Usage | Stage::Io => 1,
            Stage::Parse => 2,
            Stage::Solve => 3,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Command {
    Run {
        inp: String,
        outdir: String,
        json: bool,
    },
    Info {
        inp: String,
        json: bool,
    },
    Bench {
        json: bool,
    },
    Help,
}

fn parse_args(args: &[String]) -> Result<Command, String> {
    let mut it = args.iter();
    let Some(first) = it.next() else {
        return Ok(Command::Help);
    };
    match first.as_str() {
        "help" | "--help" | "-h" => Ok(Command::Help),
        "run" | "info" => {
            let run = first == "run";
            let mut inp: Option<String> = None;
            let mut outdir = ".".to_string();
            let mut json = false;
            let rest: Vec<&String> = it.collect();
            let mut i = 0;
            while i < rest.len() {
                let a = rest[i].as_str();
                match a {
                    "--json" => json = true,
                    "-o" if run => {
                        i += 1;
                        let Some(v) = rest.get(i) else {
                            return Err("-o requires a value".into());
                        };
                        outdir = v.to_string();
                    }
                    _ if a.starts_with('-') => return Err(format!("unknown flag: {a}")),
                    _ if inp.is_none() => inp = Some(a.to_string()),
                    _ => return Err(format!("unexpected argument: {a}")),
                }
                i += 1;
            }
            let Some(inp) = inp else {
                return Err(format!(
                    "usage: solvix {first} <job.inp>{}",
                    if run {
                        " [-o outdir] [--json]"
                    } else {
                        " [--json]"
                    }
                ));
            };
            if run {
                Ok(Command::Run { inp, outdir, json })
            } else {
                Ok(Command::Info { inp, json })
            }
        }
        "bench" => {
            let mut json = false;
            for a in it {
                match a.as_str() {
                    "--json" => json = true,
                    other => return Err(format!("unexpected argument for bench: {other}")),
                }
            }
            Ok(Command::Bench { json })
        }
        other => Err(format!("unknown subcommand: {other}")),
    }
}

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

fn fin(x: f64) -> f64 {
    if x.is_finite() {
        x
    } else {
        0.0
    }
}

fn fmt_mm(m: f64) -> String {
    format!("{:.3}", fin(m) * 1000.0)
}

fn fmt_mpa(pa: f64) -> String {
    format!("{:.3}", fin(pa) / 1e6)
}

fn fmt_e(x: f64) -> String {
    format!("{:.3e}", fin(x))
}

fn names_json(names: &[String]) -> String {
    names
        .iter()
        .map(|n| format!("\"{}\"", json_escape(n)))
        .collect::<Vec<_>>()
        .join(",")
}

fn error_json(stage: Stage, msg: &str) -> String {
    format!(
        "{{\"ok\":false,\"stage\":\"{}\",\"error\":\"{}\"}}",
        stage.as_str(),
        json_escape(msg)
    )
}

fn warnings_json(warnings: &[String]) -> String {
    format!(
        "{{\"ok\":false,\"stage\":\"parse\",\"warnings\":[{}]}}",
        names_json(warnings)
    )
}

struct RunSummary {
    job: String,
    nodes: usize,
    elements: usize,
    iterations: usize,
    residual: f64,
    solve_ms: u64,
    max_disp_m: f64,
    max_von_mises_pa: f64,
    outputs: Vec<String>,
}

fn run_summary(
    model: &Model,
    sol: &Solution,
    job: &str,
    solve_ms: u64,
    outputs: Vec<String>,
) -> RunSummary {
    let max_disp_m = sol
        .displacements
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
        .fold(0.0, f64::max);
    let max_von_mises_pa = sol.von_mises.iter().copied().fold(0.0, f64::max);
    RunSummary {
        job: job.to_string(),
        nodes: model.mesh.nodes.len(),
        elements: model.mesh.elements.len(),
        iterations: sol.iterations,
        residual: sol.final_residual,
        solve_ms,
        max_disp_m,
        max_von_mises_pa,
        outputs,
    }
}

fn success_json(s: &RunSummary) -> String {
    format!(
        "{{\"ok\":true,\"job\":\"{}\",\"nodes\":{},\"elements\":{},\"iterations\":{},\"residual\":{},\"solve_ms\":{},\"max_disp_mm\":{},\"max_von_mises_mpa\":{},\"outputs\":[{}]}}",
        json_escape(&s.job),
        s.nodes,
        s.elements,
        s.iterations,
        fmt_e(s.residual),
        s.solve_ms,
        fmt_mm(s.max_disp_m),
        fmt_mpa(s.max_von_mises_pa),
        names_json(&s.outputs),
    )
}

fn row(label: &str, value: &str) -> String {
    format!("{label:<14}{value}\n")
}

fn run_human(s: &RunSummary) -> String {
    let mut out = String::new();
    out.push_str(&row("job", &s.job));
    out.push_str(&row("nodes", &s.nodes.to_string()));
    out.push_str(&row("elements", &s.elements.to_string()));
    out.push_str(&row("iterations", &s.iterations.to_string()));
    out.push_str(&row("residual", &fmt_e(s.residual)));
    out.push_str(&row("solve time", &format!("{} ms", s.solve_ms)));
    out.push_str(&row("max |disp|", &format!("{} mm", fmt_mm(s.max_disp_m))));
    out.push_str(&row(
        "max von Mises",
        &format!("{} MPa", fmt_mpa(s.max_von_mises_pa)),
    ));
    for o in &s.outputs {
        out.push_str(&row("output", o));
    }
    out
}

struct InfoReport {
    job: String,
    nodes: usize,
    elements: usize,
    by_kind: Vec<(&'static str, usize)>,
    elsets: Vec<String>,
    materials: Vec<String>,
    bc_nodes: usize,
    cload_sum: Vec<(usize, f64)>,
}

fn kind_name(k: ElemKind) -> &'static str {
    match k {
        ElemKind::Cps3 => "cps3",
        ElemKind::Cps4 => "cps4",
    }
}

fn info_report(model: &Model, job: &str) -> InfoReport {
    let mut by_kind: Vec<(&'static str, usize)> = [ElemKind::Cps3, ElemKind::Cps4]
        .iter()
        .map(|k| (kind_name(*k), 0))
        .collect();
    for e in &model.mesh.elements {
        let name = kind_name(e.kind);
        if let Some(c) = by_kind.iter_mut().find(|(n, _)| *n == name) {
            c.1 += 1;
        }
    }
    let mut bc: Vec<u32> = model
        .boundary_conditions
        .iter()
        .flat_map(|b| b.node_ids.iter().copied())
        .collect();
    bc.sort_unstable();
    bc.dedup();
    let mut dofs: Vec<usize> = model.cloads.iter().map(|c| c.dof).collect();
    dofs.sort_unstable();
    dofs.dedup();
    let cload_sum = dofs
        .into_iter()
        .map(|d| {
            let sum = model
                .cloads
                .iter()
                .filter(|c| c.dof == d)
                .map(|c| c.value)
                .sum::<f64>();
            (d, sum)
        })
        .collect();
    InfoReport {
        job: job.to_string(),
        nodes: model.mesh.nodes.len(),
        elements: model.mesh.elements.len(),
        by_kind,
        elsets: model.elsets.iter().map(|s| s.name.clone()).collect(),
        materials: model.materials.iter().map(|m| m.name.clone()).collect(),
        bc_nodes: bc.len(),
        cload_sum,
    }
}

fn info_json(r: &InfoReport) -> String {
    let kinds = r
        .by_kind
        .iter()
        .map(|(k, n)| format!("\"{k}\":{n}"))
        .collect::<Vec<_>>()
        .join(",");
    let loads = r
        .cload_sum
        .iter()
        .map(|(d, v)| format!("\"{d}\":{:.3}", fin(*v)))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"ok\":true,\"job\":\"{}\",\"nodes\":{},\"elements\":{},\"elements_by_kind\":{{{}}},\"elsets\":{},\"elset_names\":[{}],\"materials\":{},\"material_names\":[{}],\"bc_nodes\":{},\"cload_sum_per_dof\":{{{}}}}}",
        json_escape(&r.job),
        r.nodes,
        r.elements,
        kinds,
        r.elsets.len(),
        names_json(&r.elsets),
        r.materials.len(),
        names_json(&r.materials),
        r.bc_nodes,
        loads,
    )
}

fn list_or_count(names: &[String]) -> String {
    if names.is_empty() {
        "0".to_string()
    } else {
        format!("{} ({})", names.len(), names.join(", "))
    }
}

fn info_human(r: &InfoReport) -> String {
    let kinds = r
        .by_kind
        .iter()
        .map(|(k, n)| format!("{k}: {n}"))
        .collect::<Vec<_>>()
        .join(", ");
    let loads = if r.cload_sum.is_empty() {
        "(none)".to_string()
    } else {
        r.cload_sum
            .iter()
            .map(|(d, v)| format!("dof {d}: {:.3} N", fin(*v)))
            .collect::<Vec<_>>()
            .join(" | ")
    };
    let mut out = String::new();
    out.push_str(&row("job", &r.job));
    out.push_str(&row("nodes", &r.nodes.to_string()));
    out.push_str(&row("elements", &format!("{} ({})", r.elements, kinds)));
    out.push_str(&row("elsets", &list_or_count(&r.elsets)));
    out.push_str(&row("materials", &list_or_count(&r.materials)));
    out.push_str(&row("bc nodes", &r.bc_nodes.to_string()));
    out.push_str(&row("cload sum", &loads));
    out
}

fn job_name(inp: &str) -> String {
    Path::new(inp)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("job")
        .to_string()
}

fn read_deck(inp: &str) -> Result<String, (Stage, String)> {
    std::fs::read_to_string(inp).map_err(|e| (Stage::Io, format!("cannot read {inp}: {e}")))
}

fn parse_deck(text: &str) -> Result<(Model, ParseReport), (Stage, String)> {
    solvix_io::parse_inp(text)
        .map_err(|p| (Stage::Parse, format!("line {}: {}", p.line, p.message)))
}

fn emit_warnings(report: &ParseReport) {
    if !report.warnings.is_empty() {
        eprintln!("{}", warnings_json(&report.warnings));
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
        } => format!(
            "not converged after {iterations} iterations (residual {})",
            fmt_e(residual)
        ),
    }
}

fn do_run(inp: &str, outdir: &str, json: bool) -> Result<(), (Stage, String)> {
    let text = read_deck(inp)?;
    let (model, report) = parse_deck(&text)?;
    emit_warnings(&report);
    let t0 = Instant::now();
    let sol = solvix_core::solve(&model).map_err(|e| (Stage::Solve, solve_error_message(e)))?;
    let solve_ms = t0.elapsed().as_millis() as u64;
    let job = job_name(inp);
    std::fs::create_dir_all(outdir)
        .map_err(|e| (Stage::Io, format!("cannot create directory {outdir}: {e}")))?;
    let vtu: PathBuf = Path::new(outdir).join(format!("{job}.vtu"));
    solvix_io::write_vtu(&vtu, &model, &sol)
        .map_err(|e| (Stage::Io, format!("cannot write {}: {e}", vtu.display())))?;
    let s = run_summary(
        &model,
        &sol,
        &job,
        solve_ms,
        vec![vtu.display().to_string()],
    );
    if json {
        println!("{}", success_json(&s));
    } else {
        print!("{}", run_human(&s));
    }
    Ok(())
}

fn do_info(inp: &str, json: bool) -> Result<(), (Stage, String)> {
    let text = read_deck(inp)?;
    let (model, report) = parse_deck(&text)?;
    emit_warnings(&report);
    let r = info_report(&model, &job_name(inp));
    if json {
        println!("{}", info_json(&r));
    } else {
        print!("{}", info_human(&r));
    }
    Ok(())
}

fn finish(r: Result<(), (Stage, String)>, json: bool) -> i32 {
    match r {
        Ok(()) => 0,
        Err((stage, msg)) => {
            if json {
                println!("{}", error_json(stage, &msg));
            } else {
                eprintln!("error ({}): {}", stage.as_str(), msg);
            }
            stage.exit_code()
        }
    }
}

fn cli(args: &[String]) -> i32 {
    let wants_json = args.iter().any(|a| a == "--json");
    let cmd = match parse_args(args) {
        Ok(c) => c,
        Err(msg) => {
            if wants_json {
                println!("{}", error_json(Stage::Usage, &msg));
            } else {
                eprintln!("error (usage): {msg}");
                eprint!("{USAGE}");
            }
            return Stage::Usage.exit_code();
        }
    };
    match cmd {
        Command::Help => {
            print!("{USAGE}");
            0
        }
        Command::Info { inp, json } => finish(do_info(&inp, json), json),
        Command::Run { inp, outdir, json } => finish(do_run(&inp, &outdir, json), json),
        Command::Bench { json } => bench::run(json),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    std::process::exit(cli(&args));
}

#[cfg(test)]
mod tests {
    use super::*;
    use solvix_core::{
        BoundaryCondition, ConcentratedLoad, Element, ElementSet, Material, Mesh, Node,
    };

    fn sv(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    fn fake_model() -> Model {
        Model {
            mesh: Mesh {
                nodes: vec![
                    Node {
                        id: 1,
                        coords: [0.0, 0.0, 0.0],
                    },
                    Node {
                        id: 2,
                        coords: [1.0, 0.0, 0.0],
                    },
                ],
                elements: vec![Element {
                    id: 1,
                    kind: ElemKind::Cps4,
                    connectivity: vec![1, 2, 1, 2],
                }],
            },
            thickness: 1.0,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: 210e9,
                poisson_ratio: 0.3,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: vec![1],
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            boundary_conditions: vec![BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1, 2],
            }],
            cloads: vec![ConcentratedLoad {
                node_id: 2,
                dof: 2,
                value: 1000.0,
            }],
        }
    }

    fn fake_solution() -> Solution {
        Solution {
            displacements: vec![[0.0, 0.0, 0.0], [0.0047619, 0.0, 0.0]],
            von_mises: vec![38.2e6],
            nodal_von_mises: vec![38.2e6, 38.2e6],
            reactions: vec![],
            iterations: 42,
            final_residual: 1.5e-12,
        }
    }

    #[test]
    fn unit_fmt_mm() {
        assert_eq!(fmt_mm(0.00476), "4.760");
        assert_eq!(fmt_mm(0.0047619), "4.762");
        assert_eq!(fmt_mm(0.0), "0.000");
        assert_eq!(fmt_mm(f64::NAN), "0.000");
        assert_eq!(fmt_mm(f64::INFINITY), "0.000");
    }

    #[test]
    fn unit_fmt_mpa() {
        assert_eq!(fmt_mpa(38.2e6), "38.200");
        assert_eq!(fmt_mpa(0.0), "0.000");
        assert_eq!(fmt_mpa(f64::NEG_INFINITY), "0.000");
    }

    #[test]
    fn args_help_variants() {
        assert_eq!(parse_args(&[]), Ok(Command::Help));
        for h in ["help", "--help", "-h"] {
            assert_eq!(parse_args(&sv(&[h])), Ok(Command::Help));
        }
    }

    #[test]
    fn args_run() {
        assert_eq!(
            parse_args(&sv(&["run", "job.inp"])),
            Ok(Command::Run {
                inp: "job.inp".into(),
                outdir: ".".into(),
                json: false
            })
        );
        assert_eq!(
            parse_args(&sv(&["run", "job.inp", "-o", "out", "--json"])),
            Ok(Command::Run {
                inp: "job.inp".into(),
                outdir: "out".into(),
                json: true
            })
        );
        assert_eq!(
            parse_args(&sv(&["run", "--json", "a.inp", "-o", "b"])),
            Ok(Command::Run {
                inp: "a.inp".into(),
                outdir: "b".into(),
                json: true
            })
        );
    }

    #[test]
    fn args_info() {
        assert_eq!(
            parse_args(&sv(&["info", "job.inp", "--json"])),
            Ok(Command::Info {
                inp: "job.inp".into(),
                json: true
            })
        );
    }

    #[test]
    fn args_bench() {
        assert_eq!(
            parse_args(&sv(&["bench"])),
            Ok(Command::Bench { json: false })
        );
        assert_eq!(
            parse_args(&sv(&["bench", "--json"])),
            Ok(Command::Bench { json: true })
        );
        assert!(parse_args(&sv(&["bench", "extra"])).is_err());
        assert!(parse_args(&sv(&["bench", "-o", "x"])).is_err());
    }

    #[test]
    fn args_errors() {
        assert!(parse_args(&sv(&["run"])).is_err());
        assert!(parse_args(&sv(&["info"])).is_err());
        assert!(parse_args(&sv(&["run", "a.inp", "b.inp"])).is_err());
        assert!(parse_args(&sv(&["run", "a.inp", "--pretty"])).is_err());
        assert!(parse_args(&sv(&["run", "a.inp", "-o"])).is_err());
        assert!(parse_args(&sv(&["info", "a.inp", "-o", "x"])).is_err());
        let err = parse_args(&sv(&["frobnicate"])).unwrap_err();
        assert!(err.contains("unknown subcommand"), "{err}");
        let err = parse_args(&sv(&["run", "a.inp", "-x"])).unwrap_err();
        assert!(err.contains("unknown flag"), "{err}");
    }

    #[test]
    fn exit_code_mapping() {
        assert_eq!(Stage::Usage.exit_code(), 1);
        assert_eq!(Stage::Io.exit_code(), 1);
        assert_eq!(Stage::Parse.exit_code(), 2);
        assert_eq!(Stage::Solve.exit_code(), 3);
        assert_eq!(Stage::Parse.as_str(), "parse");
    }

    #[test]
    fn success_json_exact() {
        let s = run_summary(
            &fake_model(),
            &fake_solution(),
            "demo",
            7,
            vec!["demo.vtu".into()],
        );
        assert_eq!(
            success_json(&s),
            r#"{"ok":true,"job":"demo","nodes":2,"elements":1,"iterations":42,"residual":1.500e-12,"solve_ms":7,"max_disp_mm":4.762,"max_von_mises_mpa":38.200,"outputs":["demo.vtu"]}"#
        );
    }

    #[test]
    fn success_json_escapes_job() {
        let s = run_summary(
            &fake_model(),
            &fake_solution(),
            "my \"job\"",
            0,
            vec!["a b.vtu".into()],
        );
        let j = success_json(&s);
        assert!(j.starts_with(r#"{"ok":true,"job":"my \"job\"""#), "{j}");
        assert!(j.ends_with(r#""outputs":["a b.vtu"]}"#), "{j}");
        assert!(j.contains("\"solve_ms\":0,"));
    }

    #[test]
    fn error_json_exact() {
        assert_eq!(
            error_json(Stage::Parse, "line 3: bad keyword"),
            r#"{"ok":false,"stage":"parse","error":"line 3: bad keyword"}"#
        );
        assert_eq!(
            error_json(Stage::Io, "cannot read x: \"nope\""),
            r#"{"ok":false,"stage":"io","error":"cannot read x: \"nope\""}"#
        );
        assert_eq!(
            error_json(Stage::Usage, "a\nb"),
            r#"{"ok":false,"stage":"usage","error":"a\nb"}"#
        );
    }

    #[test]
    fn warnings_json_exact() {
        assert_eq!(
            warnings_json(&["skipped *HEADING".into()]),
            r#"{"ok":false,"stage":"parse","warnings":["skipped *HEADING"]}"#
        );
        assert_eq!(
            warnings_json(&["a".into(), "b".into()]),
            r#"{"ok":false,"stage":"parse","warnings":["a","b"]}"#
        );
        assert_eq!(
            warnings_json(&[]),
            r#"{"ok":false,"stage":"parse","warnings":[]}"#
        );
    }

    #[test]
    fn info_json_exact() {
        let r = info_report(&fake_model(), "demo");
        assert_eq!(
            info_json(&r),
            r#"{"ok":true,"job":"demo","nodes":2,"elements":1,"elements_by_kind":{"cps3":0,"cps4":1},"elsets":1,"elset_names":["all"],"materials":1,"material_names":["steel"],"bc_nodes":1,"cload_sum_per_dof":{"2":1000.000}}"#
        );
    }

    #[test]
    fn info_json_no_loads() {
        let mut m = fake_model();
        m.cloads.clear();
        let r = info_report(&m, "demo");
        assert!(info_json(&r).ends_with(r#""cload_sum_per_dof":{}}"#));
    }

    #[test]
    fn info_report_sums_per_dof() {
        let mut m = fake_model();
        m.cloads.push(ConcentratedLoad {
            node_id: 1,
            dof: 2,
            value: 250.5,
        });
        m.cloads.push(ConcentratedLoad {
            node_id: 2,
            dof: 1,
            value: -10.0,
        });
        let r = info_report(&m, "demo");
        assert_eq!(r.cload_sum, vec![(1, -10.0), (2, 1250.5)]);
    }

    #[test]
    fn run_human_has_units() {
        let s = run_summary(
            &fake_model(),
            &fake_solution(),
            "demo",
            7,
            vec!["demo.vtu".into()],
        );
        let h = run_human(&s);
        assert!(h.contains("4.762 mm"), "{h}");
        assert!(h.contains("38.200 MPa"), "{h}");
        assert!(h.contains("demo.vtu"), "{h}");
        assert!(h.contains("42"), "{h}");
    }

    #[test]
    fn usage_states_agent_contract() {
        assert!(USAGE
            .contains("designed for AI agents: all subcommands accept --json and exit 0/1/2/3"));
    }

    #[test]
    fn job_name_from_path() {
        assert_eq!(job_name("examples/cantilever.inp"), "cantilever");
        assert_eq!(job_name("weird"), "weird");
        assert_eq!(job_name("/a/b/c.job.inp"), "c.job");
    }
}
