//! # abaqus2solvix — translate an Abaqus Python script to a Solvix deck
//!
//! Usage:
//!   abaqus2solvix <input.py> [-o <out.inp>] [--json] [--validate]
//!
//! `--validate` re-parses the generated deck with `solvix_io` and, when
//! solvable, solves it with `solvix_core` — the migration confidence
//! check.
//!
//! JSON-first UX: `--json` prints `{"ok":true,"output":...,
//! "warnings":[...],"unsupported":[...]}` (plus a `validate` object
//! with `--validate`) on stdout; human diagnostics go to stderr.
//! Exit codes: 0 ok · 1 usage/IO/validate failure · 2 script parse or
//! translation error.

use std::path::PathBuf;
use std::process::ExitCode;

use solvix_tools::deck::write_deck;
use solvix_tools::translate::translate;

const USAGE: &str = "\
abaqus2solvix — translate a bounded subset of Abaqus Python scripts to Solvix decks

Usage:
  abaqus2solvix <input.py> [-o <out.inp>] [--json] [--validate]
                                             translate the script to a .inp deck.
                                             Default output: <input stem>.inp next
                                             to the input. --validate re-parses the
                                             deck with solvix_io and solves it when
                                             solvable (migration confidence check).
  abaqus2solvix help | --help | -h           show this text

Exit codes: 0 ok · 1 usage/IO/validate failure · 2 script parse/translate error.
--json: {\"ok\":true,\"output\":path,\"warnings\":[...],\"unsupported\":[...]} on stdout.
";

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

fn jarr(items: &[String]) -> String {
    let inner: Vec<String> = items
        .iter()
        .map(|s| format!("\"{}\"", json_escape(s)))
        .collect();
    format!("[{}]", inner.join(","))
}

struct Args {
    input: PathBuf,
    output: PathBuf,
    json: bool,
    validate: bool,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut json = false;
    let mut validate = false;
    let mut it = argv.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "-h" | "--help" | "help" => {
                print!("{USAGE}");
                std::process::exit(0);
            }
            "--json" => json = true,
            "--validate" => validate = true,
            "-o" | "--output" => match it.next() {
                Some(o) => output = Some(PathBuf::from(o)),
                None => return Err(format!("{a} needs a value")),
            },
            other if other.starts_with('-') => {
                return Err(format!("unknown flag '{other}'"));
            }
            other => {
                if input.is_some() {
                    return Err(format!("unexpected extra argument '{other}'"));
                }
                input = Some(PathBuf::from(other));
            }
        }
    }
    let Some(input) = input else {
        return Err("missing <input.py>".to_string());
    };
    let output = output.unwrap_or_else(|| {
        let stem = input.file_stem().map(|s| s.to_string_lossy().into_owned());
        match (input.parent(), stem) {
            (Some(dir), Some(stem)) if !dir.as_os_str().is_empty() => {
                dir.join(format!("{stem}.inp"))
            }
            (_, Some(stem)) => PathBuf::from(format!("{stem}.inp")),
            _ => PathBuf::from("out.inp"),
        }
    });
    Ok(Args {
        input,
        output,
        json,
        validate,
    })
}

fn error_json(stage: &str, msg: &str) -> String {
    format!(
        "{{\"ok\":false,\"stage\":\"{stage}\",\"error\":\"{}\"}}",
        json_escape(msg)
    )
}

fn fail(stage: &str, msg: &str, code: u8, json: bool) -> ExitCode {
    if json {
        println!("{}", error_json(stage, msg));
    } else {
        eprintln!("abaqus2solvix: {msg}");
    }
    ExitCode::from(code)
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = match parse_args(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("abaqus2solvix: {e}\n");
            eprint!("{USAGE}");
            return ExitCode::from(1);
        }
    };

    let src = match std::fs::read_to_string(&args.input) {
        Ok(s) => s,
        Err(e) => {
            return fail(
                "io",
                &format!("cannot read {}: {e}", args.input.display()),
                1,
                args.json,
            )
        }
    };

    let tr = match translate(&src) {
        Ok(t) => t,
        Err(e) => return fail("parse", &format!("{e}"), 2, args.json),
    };

    let deck = write_deck(&tr.model);
    if let Err(e) = std::fs::write(&args.output, &deck) {
        return fail(
            "io",
            &format!("cannot write {}: {e}", args.output.display()),
            1,
            args.json,
        );
    }

    if !args.json {
        println!(
            "abaqus2solvix: {} -> {} ({} nodes, {} elements, {} steps)",
            args.input.display(),
            args.output.display(),
            tr.model.mesh.nodes.len(),
            tr.model.mesh.elements.len(),
            tr.model.steps.len()
        );
        for w in &tr.report.warnings {
            eprintln!("warning: {w}");
        }
        for u in &tr.report.unsupported {
            eprintln!("unsupported: {u}");
        }
    }

    // --validate: migration confidence check.
    let mut validate_json = String::new();
    let mut validate_failed = false;
    if args.validate {
        match solvix_io::parse_inp(&deck) {
            Err(e) => {
                validate_failed = true;
                let msg = format!(
                    "generated deck failed to re-parse: line {}: {}",
                    e.line, e.message
                );
                if args.json {
                    validate_json =
                        format!("{{\"parsed\":false,\"error\":\"{}\"}}", json_escape(&msg));
                } else {
                    eprintln!("validate: {msg}");
                }
            }
            Ok((m2, _)) => match solvix_core::solve(&m2) {
                Ok(sol) => {
                    let max_disp = sol
                        .displacements
                        .iter()
                        .map(|d| (d[0] * d[0] + d[1] * d[1]).sqrt())
                        .fold(0.0f64, f64::max);
                    let max_vm = sol.von_mises.iter().cloned().fold(0.0f64, f64::max);
                    if args.json {
                        validate_json = format!(
                            "{{\"parsed\":true,\"solved\":true,\"nodes\":{},\"elements\":{},\"iterations\":{},\"solver\":\"{}\",\"max_disp_mm\":{:.3},\"max_von_mises_mpa\":{:.3}}}",
                            m2.mesh.nodes.len(),
                            m2.mesh.elements.len(),
                            sol.iterations,
                            sol.solver_method,
                            max_disp * 1000.0,
                            max_vm / 1e6
                        );
                    } else {
                        println!(
                            "validate: deck re-parsed ok; solved via {} in {} iterations; max |u| = {:.3} mm, max von Mises = {:.3} MPa",
                            sol.solver_method,
                            sol.iterations,
                            max_disp * 1000.0,
                            max_vm / 1e6
                        );
                    }
                }
                Err(e) => {
                    let reason = format!("{e:?}");
                    if args.json {
                        validate_json = format!(
                            "{{\"parsed\":true,\"solved\":false,\"reason\":\"{}\"}}",
                            json_escape(&reason)
                        );
                    } else {
                        eprintln!("validate: deck re-parsed ok but is not solvable: {reason}");
                    }
                }
            },
        }
    }

    if args.json {
        let mut out = format!(
            "{{\"ok\":true,\"output\":\"{}\",\"nodes\":{},\"elements\":{},\"steps\":{},\"warnings\":{},\"unsupported\":{}",
            json_escape(&args.output.display().to_string()),
            tr.model.mesh.nodes.len(),
            tr.model.mesh.elements.len(),
            tr.model.steps.len(),
            jarr(&tr.report.warnings),
            jarr(&tr.report.unsupported)
        );
        if args.validate {
            out.push_str(&format!(",\"validate\":{validate_json}"));
        }
        out.push('}');
        println!("{out}");
    }

    if validate_failed {
        ExitCode::from(1)
    } else {
        ExitCode::SUCCESS
    }
}
