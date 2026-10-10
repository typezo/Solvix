//! # solvix — the agent-era FEA CLI
//!
//! Implementation owned by issue #3; `bench` by issue #17. Subcommands:
//!   solvix run <job.inp> [-o outdir] [--odb] [--json]
//!   solvix info <job.inp> [--json]
//!   solvix odb <file.solvix> [--inspect | --frame N] [--json]
//!   solvix bench [--json]
//!
//! `run` auto-detects the analysis type from the deck (issues #46/#47):
//! a `*DYNAMIC` step runs transient Newmark dynamics, temperature data
//! without dynamics runs the thermal-stress pass, anything else is
//! linear statics. No flags — see [`solve_auto`].
//!
//! `--odb` additionally writes a binary `<job>.solvix` output database
//! (issue #32): static jobs store one frame, NLGEOM/plasticity/contact
//! jobs one frame per load increment, dynamics up to 200 subsampled
//! frames. `solvix odb` inspects them.
//!
//! Exit codes: 0 ok · 1 usage/IO · 2 parse error · 3 solve error.
//! stdout carries the payload (human table or single-line JSON), stderr
//! carries warnings — designed for AI agents.

use std::path::{Path, PathBuf};
use std::time::Instant;

use solvix_core::{
    DynSolution, ElemKind, Model, NlFrame, NlOptions, Solution, SolveError, TimePoint,
};
use solvix_io::{Odb, OdbFrame, OdbWriter, ParseReport};

mod bench;
mod mcp;

const USAGE: &str = "\
solvix 0.1.0 — the agent-era FEA CLI

solvix is designed for AI agents: all subcommands accept --json and exit 0/1/2/3.

Usage:
  solvix run <job.inp> [-o outdir] [--modes N] [--odb] [--json]
                                               parse, solve, write <jobname>.vtu;
                                               --modes N also extracts the N lowest
                                               natural frequencies (modal analysis).
                                               --odb also writes <jobname>.solvix, a
                                               binary output database (frames with
                                               displacements, von Mises, temperatures
                                               or velocities, plus per-node |u|
                                               history; inspect with `solvix odb`).
                                               Analysis type is AUTO-DETECTED from
                                               the deck (no flags): a *DYNAMIC step
                                               runs transient Newmark dynamics (JSON
                                               gains \"dynamic\":true, steps, dt,
                                               duration, final_disp_mm; max_disp_mm
                                               is the max over all time points); step
                                               temperature data (*HEAT TRANSFER,
                                               *TEMPERATURE, *DFLUX) without dynamics
                                               runs a thermal-stress pass (JSON gains
                                               \"thermal\":true, max_temperature_c);
                                               anything else is linear statics
  solvix info <job.inp> [--json]              parse-only model report, never solves
  solvix odb <file.solvix> [--inspect | --frame N] [--json]
                                               binary output database reader: default
                                               (or --inspect) prints the file schema
                                               (frames, fields, model size, bytes);
                                               --frame N prints that frame's summary
                                               (time, load fraction, max displacement
                                               / stress / temperature / velocity)
  solvix bench [--json]                       built-in verification suite (all-pass exit 0, else 3)
  solvix mcp                                  MCP server on stdio (JSON-RPC 2.0, one object per line)
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
        modes: Option<usize>,
        odb: bool,
        json: bool,
    },
    Info {
        inp: String,
        json: bool,
    },
    Odb {
        file: String,
        /// `None` means inspect (the schema); `Some(n)` reads frame n.
        frame: Option<usize>,
        json: bool,
    },
    Bench {
        json: bool,
    },
    Mcp,
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
            let mut modes: Option<usize> = None;
            let mut odb = false;
            let mut json = false;
            let rest: Vec<&String> = it.collect();
            let mut i = 0;
            while i < rest.len() {
                let a = rest[i].as_str();
                match a {
                    "--json" => json = true,
                    "--odb" if run => odb = true,
                    "--modes" if run => {
                        i += 1;
                        let Some(v) = rest.get(i) else {
                            return Err("--modes requires a value".into());
                        };
                        modes = Some(
                            v.parse()
                                .map_err(|_| format!("invalid --modes value: {v}"))?,
                        );
                    }
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
                        " [-o outdir] [--modes N] [--odb] [--json]"
                    } else {
                        " [--json]"
                    }
                ));
            };
            if run {
                Ok(Command::Run {
                    inp,
                    outdir,
                    modes,
                    odb,
                    json,
                })
            } else {
                Ok(Command::Info { inp, json })
            }
        }
        "odb" => {
            let mut file: Option<String> = None;
            let mut inspect = false;
            let mut frame: Option<usize> = None;
            let mut json = false;
            let rest: Vec<&String> = it.collect();
            let mut i = 0;
            while i < rest.len() {
                let a = rest[i].as_str();
                match a {
                    "--json" => json = true,
                    "--inspect" => inspect = true,
                    "--frame" => {
                        i += 1;
                        let Some(v) = rest.get(i) else {
                            return Err("--frame requires a value".into());
                        };
                        frame = Some(
                            v.parse()
                                .map_err(|_| format!("invalid --frame value: {v}"))?,
                        );
                    }
                    _ if a.starts_with('-') => return Err(format!("unknown flag: {a}")),
                    _ if file.is_none() => file = Some(a.to_string()),
                    _ => return Err(format!("unexpected argument: {a}")),
                }
                i += 1;
            }
            if inspect && frame.is_some() {
                return Err("--inspect and --frame are mutually exclusive".into());
            }
            let Some(file) = file else {
                return Err(
                    "usage: solvix odb <file.solvix> [--inspect | --frame N] [--json]".into(),
                );
            };
            Ok(Command::Odb { file, frame, json })
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
        "mcp" => {
            if let Some(extra) = it.next() {
                return Err(format!("unexpected argument for mcp: {extra}"));
            }
            Ok(Command::Mcp)
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
    solver_method: String,
    cond_estimate: f64,
    residual: f64,
    solve_ms: u64,
    max_disp_m: f64,
    max_von_mises_pa: f64,
    outputs: Vec<String>,
    /// Natural frequencies (Hz) when `--modes N` ran a modal pass; empty
    /// otherwise.
    frequencies_hz: Vec<f64>,
    /// Auto-detected thermal pass (issue #46): nodal temperature range.
    thermal: Option<ThermalInfo>,
    /// Auto-detected transient dynamics (issue #47).
    dynamic: Option<DynamicInfo>,
}

fn run_summary(
    model: &Model,
    outcome: &AutoOutcome,
    job: &str,
    solve_ms: u64,
    outputs: Vec<String>,
    frequencies_hz: Vec<f64>,
) -> RunSummary {
    let (sol, thermal, dynamic) = match outcome {
        AutoOutcome::Static { sol, thermal, .. } => (sol, *thermal, None),
        AutoOutcome::Dynamic {
            info, final_sol, ..
        } => (final_sol, None, Some(*info)),
    };
    // Dynamic runs report the max over ALL time points, not just the
    // final frame packed into `sol`.
    let max_disp_m = match &dynamic {
        Some(d) => d.max_disp_m,
        None => sol
            .displacements
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max),
    };
    let max_von_mises_pa = sol.von_mises.iter().copied().fold(0.0, f64::max);
    RunSummary {
        job: job.to_string(),
        nodes: model.mesh.nodes.len(),
        elements: model.mesh.elements.len(),
        iterations: sol.iterations,
        solver_method: sol.solver_method.clone(),
        cond_estimate: sol.cond_estimate,
        residual: sol.final_residual,
        solve_ms,
        max_disp_m,
        max_von_mises_pa,
        outputs,
        frequencies_hz,
        thermal,
        dynamic,
    }
}

fn fmt_hz(x: f64) -> String {
    format!("{:.3}", fin(x))
}

fn fmt_c(c: f64) -> String {
    format!("{:.3}", fin(c))
}

/// Seconds in plain decimals where reasonable, scientific otherwise.
fn fmt_secs(x: f64) -> String {
    let x = fin(x);
    if x != 0.0 && (x.abs() < 1e-3 || x.abs() >= 1e6) {
        fmt_e(x)
    } else {
        let s = format!("{x:.6}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

/// Nodal temperature range of an auto-detected thermal run (the deck's
/// temperature units — only differences enter the stress pass).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ThermalInfo {
    pub(crate) min_c: f64,
    pub(crate) max_c: f64,
}

/// Summary of an auto-detected transient run (`*DYNAMIC`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct DynamicInfo {
    pub(crate) steps: usize,
    pub(crate) dt: f64,
    pub(crate) duration: f64,
    /// Max over ALL time points of the largest nodal |u| (m).
    pub(crate) max_disp_m: f64,
    /// Largest nodal |u| at the final time (m).
    pub(crate) final_disp_m: f64,
}

/// What a run actually solved (issues #46/#47): `*DYNAMIC` in any step
/// wins (transient Newmark; thermal data in the same deck is ignored in
/// v1 with a stderr warning), temperature data without dynamics runs the
/// thermal-stress pass, anything else is plain linear statics (NLGEOM,
/// plasticity and contact ride the NL machinery with one frame per load
/// increment, kept in [`AutoOutcome::Static::nl_frames`]).
pub(crate) enum AutoOutcome {
    Static {
        sol: Solution,
        thermal: Option<ThermalInfo>,
        /// Nodal temperatures of a thermal run (parallel to
        /// `mesh.nodes`); kept for the `.solvix` output database.
        temps: Option<Vec<f64>>,
        /// Converged NL increments (NLGEOM/plasticity/contact decks);
        /// empty for linear statics. `sol` is the last frame's solution.
        nl_frames: Vec<NlFrame>,
    },
    Dynamic {
        info: DynamicInfo,
        /// Final-time displacements packed as a [`Solution`] for the VTU
        /// writer and summary fields (v1 dynamics recovers no stresses).
        final_sol: Solution,
        /// All output points (kept for the `.solvix` output database;
        /// v1 dynamics recovers no stresses).
        dyn_points: Vec<TimePoint>,
    },
}

pub(crate) fn has_dynamics(model: &Model) -> bool {
    model.steps.iter().any(|s| s.dynamics.is_some())
}

pub(crate) fn has_thermal(model: &Model) -> bool {
    model
        .steps
        .iter()
        .any(|s| !s.temp_bcs.is_empty() || !s.heat_fluxes.is_empty() || s.heat_source.is_some())
}

/// Whether `solvix_core::solve` would route through the NL machinery
/// (contact pairs, plasticity or NLGEOM on the last step).
fn is_nl(model: &Model) -> bool {
    !model.contact_pairs.is_empty()
        || model.plasticity.is_some()
        || model.steps.last().is_some_and(|s| s.nlgeom)
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

/// Final-time displacements of a transient run as a [`Solution`]: the VTU
/// writer and the summary consume per-node triples; stress recovery is
/// not part of v1 dynamics (empty arrays render as 0 in the VTU).
fn final_frame_solution(model: &Model, dyn_sol: &DynSolution) -> Solution {
    let nodes = model.mesh.nodes.len();
    let flat = dyn_sol
        .points
        .last()
        .map(|p| p.displacements.as_slice())
        .unwrap_or(&[]);
    Solution {
        displacements: flat_triples(flat, nodes),
        von_mises: Vec::new(),
        nodal_von_mises: Vec::new(),
        reactions: Vec::new(),
        iterations: dyn_sol.iterations_total,
        final_residual: 0.0,
        solver_method: "newmark".into(),
        contact_forces: Vec::new(),
        cond_estimate: 0.0,
    }
}

/// The `run` pipeline's analysis-type dispatch (issues #46/#47): detect
/// from step content, solve with the matching core entry point.
pub(crate) fn solve_auto(model: &Model) -> Result<AutoOutcome, SolveError> {
    if has_dynamics(model) {
        let dyn_sol = solvix_core::solve_dynamic(model)?;
        let spec = model
            .steps
            .iter()
            .rev()
            .find(|s| s.dynamics.is_some())
            .and_then(|s| s.dynamics.as_ref())
            .expect("has_dynamics checked is_some");
        let nodes = model.mesh.nodes.len();
        let ndof = dyn_sol
            .points
            .last()
            .map(|p| p.displacements.len() / nodes.max(1))
            .unwrap_or(2);
        let max_disp_m = dyn_sol
            .points
            .iter()
            .map(|p| max_abs_u(&p.displacements, ndof))
            .fold(0.0, f64::max);
        let final_disp_m = dyn_sol
            .points
            .last()
            .map(|p| max_abs_u(&p.displacements, ndof))
            .unwrap_or(0.0);
        let info = DynamicInfo {
            steps: dyn_sol.points.len().saturating_sub(1),
            dt: spec.dt,
            duration: spec.duration,
            max_disp_m,
            final_disp_m,
        };
        let final_sol = final_frame_solution(model, &dyn_sol);
        return Ok(AutoOutcome::Dynamic {
            info,
            final_sol,
            dyn_points: dyn_sol.points,
        });
    }
    if has_thermal(model) {
        // solve_thermomech runs the conduction pass internally; the extra
        // solve_thermal recovers the temperature field for the summary
        // and the output database (core offers no temps-out variant —
        // the double conduction solve is the v0.4 price).
        let temps = solvix_core::solve_thermal(model)?;
        let sol = solvix_core::solve_thermomech(model)?;
        let info = ThermalInfo {
            min_c: temps
                .temperatures
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min),
            max_c: temps
                .temperatures
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
        };
        return Ok(AutoOutcome::Static {
            sol,
            thermal: Some(info),
            temps: Some(temps.temperatures),
            nl_frames: Vec::new(),
        });
    }
    if is_nl(model) {
        // Same routing solvix_core::solve applies, but the frames are
        // kept for the output database.
        let nl_frames = solvix_core::solve_nl(model, &NlOptions::default())?;
        let sol = nl_frames
            .last()
            .map(|f| f.solution.clone())
            .unwrap_or_default();
        return Ok(AutoOutcome::Static {
            sol,
            thermal: None,
            temps: None,
            nl_frames,
        });
    }
    Ok(AutoOutcome::Static {
        sol: solvix_core::solve(model)?,
        thermal: None,
        temps: None,
        nl_frames: Vec::new(),
    })
}

// ---- issue #32: .solvix output database ----

/// Frame ceiling for dynamic jobs written to the output database.
const ODB_MAX_FRAMES: usize = 200;

/// Indices of the output points kept when subsampling to
/// [`ODB_MAX_FRAMES`]: all of them when there are fewer, else 200
/// evenly spaced from first to last (inclusive).
fn dyn_frame_indices(n: usize) -> Vec<usize> {
    if n <= ODB_MAX_FRAMES {
        (0..n).collect()
    } else {
        (0..ODB_MAX_FRAMES)
            .map(|k| k * (n - 1) / (ODB_MAX_FRAMES - 1))
            .collect()
    }
}

/// The output-database frames of a run: static jobs one frame, NL jobs
/// one per converged increment, dynamic jobs the subsampled time points
/// with velocities.
fn odb_frames(model: &Model, outcome: &AutoOutcome) -> Vec<OdbFrame> {
    let nodes = model.mesh.nodes.len();
    match outcome {
        AutoOutcome::Static {
            sol,
            temps,
            nl_frames,
            ..
        } => {
            if nl_frames.is_empty() {
                vec![OdbFrame {
                    time: 0.0,
                    load_fraction: 1.0,
                    displacements: sol.displacements.clone(),
                    von_mises: sol.von_mises.clone(),
                    nodal_von_mises: sol.nodal_von_mises.clone(),
                    temperatures: temps.clone(),
                    velocities: None,
                }]
            } else {
                nl_frames
                    .iter()
                    .map(|f| OdbFrame {
                        time: 0.0,
                        load_fraction: f.load_fraction,
                        displacements: f.solution.displacements.clone(),
                        von_mises: f.solution.von_mises.clone(),
                        nodal_von_mises: f.solution.nodal_von_mises.clone(),
                        temperatures: None,
                        velocities: None,
                    })
                    .collect()
            }
        }
        AutoOutcome::Dynamic {
            info, dyn_points, ..
        } => dyn_frame_indices(dyn_points.len())
            .into_iter()
            .map(|i| {
                let p = &dyn_points[i];
                let load_fraction = if info.duration > 0.0 {
                    p.time / info.duration
                } else {
                    0.0
                };
                OdbFrame {
                    time: p.time,
                    load_fraction,
                    displacements: flat_triples(&p.displacements, nodes),
                    von_mises: Vec::new(),
                    nodal_von_mises: Vec::new(),
                    temperatures: None,
                    velocities: Some(flat_triples(&p.velocities, nodes)),
                }
            })
            .collect(),
    }
}

/// Write a run's output database to `path` (model once, then frames).
fn write_job_odb(path: &Path, model: &Model, outcome: &AutoOutcome) -> Result<(), (Stage, String)> {
    let wrap = |e: std::io::Error| (Stage::Io, format!("cannot write {}: {e}", path.display()));
    let mut w = OdbWriter::create(path, &model.mesh).map_err(wrap)?;
    for frame in odb_frames(model, outcome) {
        w.append_frame(&frame).map_err(wrap)?;
    }
    w.finish().map_err(wrap)?;
    Ok(())
}

/// Schema report for `solvix odb --inspect` (issue #32).
struct OdbInspect {
    file: String,
    frames: usize,
    fields: Vec<&'static str>,
    nodes: usize,
    elements: usize,
    bytes: u64,
}

/// One frame's summary for `solvix odb --frame N`.
struct FrameSummary {
    frame: usize,
    time: f64,
    load_fraction: f64,
    max_disp_m: f64,
    max_von_mises_pa: f64,
    max_nodal_von_mises_pa: f64,
    max_temperature: Option<f64>,
    max_velocity_m: Option<f64>,
}

fn frame_summary(frame: usize, f: &OdbFrame) -> FrameSummary {
    let max_disp_m = f
        .displacements
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
        .fold(0.0, f64::max);
    let max_vel_m = f.velocities.as_ref().map(|vs| {
        vs.iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0, f64::max)
    });
    FrameSummary {
        frame,
        time: f.time,
        load_fraction: f.load_fraction,
        max_disp_m,
        max_von_mises_pa: f.von_mises.iter().copied().fold(0.0, f64::max),
        max_nodal_von_mises_pa: f.nodal_von_mises.iter().copied().fold(0.0, f64::max),
        max_temperature: f
            .temperatures
            .as_ref()
            .map(|t| t.iter().copied().fold(f64::NEG_INFINITY, f64::max))
            .filter(|t| t.is_finite()),
        max_velocity_m: max_vel_m,
    }
}

fn inspect_json(r: &OdbInspect) -> String {
    format!(
        "{{\"ok\":true,\"file\":\"{}\",\"version\":1,\"nodes\":{},\"elements\":{},\"frames\":{},\"fields\":[{}],\"bytes\":{}}}",
        json_escape(&r.file),
        r.nodes,
        r.elements,
        r.frames,
        r.fields
            .iter()
            .map(|f| format!("\"{f}\""))
            .collect::<Vec<_>>()
            .join(","),
        r.bytes,
    )
}

fn inspect_human(r: &OdbInspect) -> String {
    let mut out = String::new();
    out.push_str(&row("file", &r.file));
    out.push_str(&row("format", "SLVXODB1 v1"));
    out.push_str(&row("nodes", &r.nodes.to_string()));
    out.push_str(&row("elements", &r.elements.to_string()));
    out.push_str(&row("frames", &r.frames.to_string()));
    out.push_str(&row("fields", &r.fields.join(", ")));
    out.push_str(&row("bytes", &r.bytes.to_string()));
    out
}

fn frame_json(file: &str, s: &FrameSummary) -> String {
    let mut j = format!(
        "{{\"ok\":true,\"file\":\"{}\",\"frame\":{},\"time\":{},\"load_fraction\":{},\"max_disp_mm\":{},\"max_von_mises_mpa\":{},\"max_nodal_von_mises_mpa\":{}",
        json_escape(file),
        s.frame,
        fmt_secs(s.time),
        fmt_secs(s.load_fraction),
        fmt_mm(s.max_disp_m),
        fmt_mpa(s.max_von_mises_pa),
        fmt_mpa(s.max_nodal_von_mises_pa),
    );
    if let Some(t) = s.max_temperature {
        j.push_str(&format!(",\"max_temperature_c\":{}", fmt_c(t)));
    }
    if let Some(v) = s.max_velocity_m {
        j.push_str(&format!(",\"max_velocity_mm_s\":{}", fmt_mm(v)));
    }
    j.push('}');
    j
}

fn frame_human(file: &str, s: &FrameSummary) -> String {
    let mut out = String::new();
    out.push_str(&row("file", file));
    out.push_str(&row("frame", &s.frame.to_string()));
    out.push_str(&row("time", &format!("{} s", fmt_secs(s.time))));
    out.push_str(&row("load fraction", &fmt_secs(s.load_fraction)));
    out.push_str(&row("max |disp|", &format!("{} mm", fmt_mm(s.max_disp_m))));
    out.push_str(&row(
        "max von Mises",
        &format!("{} MPa", fmt_mpa(s.max_von_mises_pa)),
    ));
    out.push_str(&row(
        "max nodal vM",
        &format!("{} MPa", fmt_mpa(s.max_nodal_von_mises_pa)),
    ));
    if let Some(t) = s.max_temperature {
        out.push_str(&row("max temp", &format!("{} C", fmt_c(t))));
    }
    if let Some(v) = s.max_velocity_m {
        out.push_str(&row("max |vel|", &format!("{} mm/s", fmt_mm(v))));
    }
    out
}

fn success_json(s: &RunSummary) -> String {
    let modal = if s.frequencies_hz.is_empty() {
        String::new()
    } else {
        let freqs = s
            .frequencies_hz
            .iter()
            .map(|f| fmt_hz(*f))
            .collect::<Vec<_>>()
            .join(",");
        format!(
            ",\"modes\":{},\"frequencies_hz\":[{}]",
            s.frequencies_hz.len(),
            freqs
        )
    };
    let thermal = match s.thermal {
        Some(t) => format!(",\"thermal\":true,\"max_temperature_c\":{}", fmt_c(t.max_c)),
        None => String::new(),
    };
    let dynamic = match s.dynamic {
        Some(d) => format!(
            ",\"dynamic\":true,\"steps\":{},\"dt\":{},\"duration\":{},\"final_disp_mm\":{}",
            d.steps,
            fmt_secs(d.dt),
            fmt_secs(d.duration),
            fmt_mm(d.final_disp_m)
        ),
        None => String::new(),
    };
    let mut j = format!(
        "{{\"ok\":true,\"job\":\"{}\",\"nodes\":{},\"elements\":{},\"iterations\":{},\"solver_method\":\"{}\",\"cond_estimate\":{},\"residual\":{},\"solve_ms\":{},\"max_disp_mm\":{},\"max_von_mises_mpa\":{},\"outputs\":[{}]",
        json_escape(&s.job),
        s.nodes,
        s.elements,
        s.iterations,
        json_escape(&s.solver_method),
        fmt_e(s.cond_estimate),
        fmt_e(s.residual),
        s.solve_ms,
        fmt_mm(s.max_disp_m),
        fmt_mpa(s.max_von_mises_pa),
        names_json(&s.outputs),
    );
    j.push_str(&modal);
    j.push_str(&thermal);
    j.push_str(&dynamic);
    j.push('}');
    j
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
    out.push_str(&row("solver method", &s.solver_method));
    out.push_str(&row("cond estimate", &fmt_e(s.cond_estimate)));
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
    if !s.frequencies_hz.is_empty() {
        out.push_str(&row("modes", &s.frequencies_hz.len().to_string()));
        for (i, f) in s.frequencies_hz.iter().enumerate() {
            out.push_str(&row(
                &format!("mode {}", i + 1),
                &format!("{} Hz", fmt_hz(*f)),
            ));
        }
    }
    if let Some(t) = s.thermal {
        out.push_str(&row("analysis", "thermal-stress (auto)"));
        out.push_str(&row("max temp", &format!("{} C", fmt_c(t.max_c))));
        out.push_str(&row("min temp", &format!("{} C", fmt_c(t.min_c))));
    }
    if let Some(d) = s.dynamic {
        out.push_str(&row("analysis", "dynamics, Newmark (auto)"));
        out.push_str(&row("time steps", &d.steps.to_string()));
        out.push_str(&row("dt", &format!("{} s", fmt_secs(d.dt))));
        out.push_str(&row("duration", &format!("{} s", fmt_secs(d.duration))));
        out.push_str(&row(
            "final |disp|",
            &format!("{} mm", fmt_mm(d.final_disp_m)),
        ));
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
    warnings: usize,
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

fn info_report(model: &Model, job: &str, report: &ParseReport) -> InfoReport {
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
        warnings: report.warnings.len(),
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
        "{{\"ok\":true,\"job\":\"{}\",\"nodes\":{},\"elements\":{},\"elements_by_kind\":{{{}}},\"elsets\":{},\"elset_names\":[{}],\"materials\":{},\"material_names\":[{}],\"bc_nodes\":{},\"cload_sum_per_dof\":{{{}}},\"warnings\":{}}}",
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
        r.warnings,
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
    if r.warnings > 0 {
        out.push_str(&row("warnings", &r.warnings.to_string()));
    }
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

fn do_run(
    inp: &str,
    outdir: &str,
    modes: Option<usize>,
    odb: bool,
    json: bool,
) -> Result<(), (Stage, String)> {
    let text = read_deck(inp)?;
    let (model, report) = parse_deck(&text)?;
    emit_warnings(&report);
    if has_dynamics(&model) && has_thermal(&model) {
        eprintln!("warning: *DYNAMIC step present — thermal data ignored (v1: dynamics wins)");
    }
    let t0 = Instant::now();
    let outcome = solve_auto(&model).map_err(|e| (Stage::Solve, solve_error_message(e)))?;
    let solve_ms = t0.elapsed().as_millis() as u64;
    // Modal pass (issue #40): after the solve, extract the N lowest
    // natural frequencies. VTU mode-shape output is deferred: solvix-io is
    // owned by issue #25 (write_vtu extension is a #25-followup item), as
    // is temperature PointData for thermal jobs (issue #46 follow-up).
    let frequencies_hz = match modes {
        Some(n) => {
            solvix_core::solve_modal(&model, &solvix_core::ModalOptions { num_modes: n })
                .map_err(|e| (Stage::Solve, solve_error_message(e)))?
                .frequencies_hz
        }
        None => Vec::new(),
    };
    let job = job_name(inp);
    std::fs::create_dir_all(outdir)
        .map_err(|e| (Stage::Io, format!("cannot create directory {outdir}: {e}")))?;
    let vtu: PathBuf = Path::new(outdir).join(format!("{job}.vtu"));
    let sol_for_vtu = match &outcome {
        AutoOutcome::Static { sol, .. } => sol,
        AutoOutcome::Dynamic { final_sol, .. } => final_sol,
    };
    solvix_io::write_vtu(&vtu, &model, sol_for_vtu)
        .map_err(|e| (Stage::Io, format!("cannot write {}: {e}", vtu.display())))?;
    let mut outputs = vec![vtu.display().to_string()];
    if odb {
        // Issue #32: binary output database alongside the VTU.
        let odb_path: PathBuf = Path::new(outdir).join(format!("{job}.solvix"));
        write_job_odb(&odb_path, &model, &outcome)?;
        outputs.push(odb_path.display().to_string());
    }
    let s = run_summary(&model, &outcome, &job, solve_ms, outputs, frequencies_hz);
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
    let r = info_report(&model, &job_name(inp), &report);
    if json {
        println!("{}", info_json(&r));
    } else {
        print!("{}", info_human(&r));
    }
    Ok(())
}

/// `solvix odb <file>`: `frame == None` inspects the schema, `Some(n)`
/// prints frame n's summary. A malformed database (bad magic/version)
/// maps to the parse stage; other IO failures to the io stage.
fn do_odb(file: &str, frame: Option<usize>, json: bool) -> Result<(), (Stage, String)> {
    let odb = Odb::open(Path::new(file)).map_err(|e| {
        let stage = if e.kind() == std::io::ErrorKind::InvalidData {
            Stage::Parse
        } else {
            Stage::Io
        };
        (stage, format!("cannot open {file}: {e}"))
    })?;
    match frame {
        None => {
            let mut fields = vec!["displacements", "von_mises", "nodal_von_mises"];
            if odb.has_temperatures() {
                fields.push("temperatures");
            }
            if odb.has_velocities() {
                fields.push("velocities");
            }
            let r = OdbInspect {
                file: file.to_string(),
                frames: odb.frame_count(),
                fields,
                nodes: odb.model().nodes.len(),
                elements: odb.model().elements.len(),
                bytes: odb.file_len(),
            };
            if json {
                println!("{}", inspect_json(&r));
            } else {
                print!("{}", inspect_human(&r));
            }
        }
        Some(n) => {
            if n >= odb.frame_count() {
                return Err((
                    Stage::Usage,
                    format!("frame {n} out of range ({} frames)", odb.frame_count()),
                ));
            }
            let f = odb
                .frame(n)
                .map_err(|e| (Stage::Parse, format!("cannot read frame {n}: {e}")))?;
            let s = frame_summary(n, &f);
            if json {
                println!("{}", frame_json(file, &s));
            } else {
                print!("{}", frame_human(file, &s));
            }
        }
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
        Command::Odb { file, frame, json } => finish(do_odb(&file, frame, json), json),
        Command::Run {
            inp,
            outdir,
            modes,
            odb,
            json,
        } => finish(do_run(&inp, &outdir, modes, odb, json), json),
        Command::Bench { json } => bench::run(json),
        Command::Mcp => mcp::serve(),
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
            density: 7850.0,
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
            ..Default::default()
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
            solver_method: "pcg".into(),
            contact_forces: Vec::new(),
            cond_estimate: 12.5,
        }
    }

    fn static_outcome() -> AutoOutcome {
        AutoOutcome::Static {
            sol: fake_solution(),
            thermal: None,
            temps: None,
            nl_frames: Vec::new(),
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
                modes: None,
                odb: false,
                json: false
            })
        );
        assert_eq!(
            parse_args(&sv(&["run", "job.inp", "-o", "out", "--json"])),
            Ok(Command::Run {
                inp: "job.inp".into(),
                outdir: "out".into(),
                modes: None,
                odb: false,
                json: true
            })
        );
        assert_eq!(
            parse_args(&sv(&["run", "--json", "a.inp", "-o", "b"])),
            Ok(Command::Run {
                inp: "a.inp".into(),
                outdir: "b".into(),
                modes: None,
                odb: false,
                json: true
            })
        );
        assert_eq!(
            parse_args(&sv(&["run", "job.inp", "--modes", "3", "--json"])),
            Ok(Command::Run {
                inp: "job.inp".into(),
                outdir: ".".into(),
                modes: Some(3),
                odb: false,
                json: true
            })
        );
        assert_eq!(
            parse_args(&sv(&["run", "--modes", "6", "job.inp"])),
            Ok(Command::Run {
                inp: "job.inp".into(),
                outdir: ".".into(),
                modes: Some(6),
                odb: false,
                json: false
            })
        );
    }

    #[test]
    fn args_modes_errors() {
        assert!(parse_args(&sv(&["run", "a.inp", "--modes"])).is_err());
        let err = parse_args(&sv(&["run", "a.inp", "--modes", "x"])).unwrap_err();
        assert!(err.contains("invalid --modes value"), "{err}");
        let err = parse_args(&sv(&["run", "a.inp", "--modes", "-1"])).unwrap_err();
        assert!(err.contains("invalid --modes value"), "{err}");
        let err = parse_args(&sv(&["info", "a.inp", "--modes", "3"])).unwrap_err();
        assert!(err.contains("unknown flag"), "{err}");
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
            &static_outcome(),
            "demo",
            7,
            vec!["demo.vtu".into()],
            Vec::new(),
        );
        assert_eq!(
            success_json(&s),
            r#"{"ok":true,"job":"demo","nodes":2,"elements":1,"iterations":42,"solver_method":"pcg","cond_estimate":1.250e1,"residual":1.500e-12,"solve_ms":7,"max_disp_mm":4.762,"max_von_mises_mpa":38.200,"outputs":["demo.vtu"]}"#
        );
    }

    #[test]
    fn success_json_with_modes_exact() {
        let s = run_summary(
            &fake_model(),
            &static_outcome(),
            "demo",
            7,
            vec!["demo.vtu".into()],
            vec![12.5, 73.21, 198.0],
        );
        assert_eq!(
            success_json(&s),
            r#"{"ok":true,"job":"demo","nodes":2,"elements":1,"iterations":42,"solver_method":"pcg","cond_estimate":1.250e1,"residual":1.500e-12,"solve_ms":7,"max_disp_mm":4.762,"max_von_mises_mpa":38.200,"outputs":["demo.vtu"],"modes":3,"frequencies_hz":[12.500,73.210,198.000]}"#
        );
    }

    #[test]
    fn run_human_modal_table() {
        let s = run_summary(
            &fake_model(),
            &static_outcome(),
            "demo",
            7,
            vec!["demo.vtu".into()],
            vec![12.5, 73.21],
        );
        let h = run_human(&s);
        assert!(h.contains("modes         2"), "{h}");
        assert!(h.contains("mode 1        12.500 Hz"), "{h}");
        assert!(h.contains("mode 2        73.210 Hz"), "{h}");
    }

    #[test]
    fn run_human_static_has_no_mode_rows() {
        let s = run_summary(
            &fake_model(),
            &static_outcome(),
            "demo",
            7,
            vec!["demo.vtu".into()],
            Vec::new(),
        );
        let h = run_human(&s);
        assert!(!h.contains("mode 1"), "{h}");
    }

    #[test]
    fn success_json_escapes_job() {
        let s = run_summary(
            &fake_model(),
            &static_outcome(),
            "my \"job\"",
            0,
            vec!["a b.vtu".into()],
            Vec::new(),
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
        let report = ParseReport {
            warnings: vec!["skipped *HEADING".into()],
        };
        let r = info_report(&fake_model(), "demo", &report);
        assert_eq!(
            info_json(&r),
            r#"{"ok":true,"job":"demo","nodes":2,"elements":1,"elements_by_kind":{"cps3":0,"cps4":1},"elsets":1,"elset_names":["all"],"materials":1,"material_names":["steel"],"bc_nodes":1,"cload_sum_per_dof":{"2":1000.000},"warnings":1}"#
        );
    }

    #[test]
    fn info_json_no_loads() {
        let mut m = fake_model();
        m.cloads.clear();
        let r = info_report(&m, "demo", &ParseReport::default());
        assert!(info_json(&r).ends_with(r#""cload_sum_per_dof":{},"warnings":0}"#));
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
        let r = info_report(&m, "demo", &ParseReport::default());
        assert_eq!(r.cload_sum, vec![(1, -10.0), (2, 1250.5)]);
    }

    #[test]
    fn run_human_has_units() {
        let s = run_summary(
            &fake_model(),
            &static_outcome(),
            "demo",
            7,
            vec!["demo.vtu".into()],
            Vec::new(),
        );
        let h = run_human(&s);
        assert!(h.contains("4.762 mm"), "{h}");
        assert!(h.contains("38.200 MPa"), "{h}");
        assert!(h.contains("demo.vtu"), "{h}");
        assert!(h.contains("42"), "{h}");
        assert!(h.contains("solver method pcg"), "{h}");
        assert!(h.contains("cond estimate 1.250e1"), "{h}");
    }

    #[test]
    fn info_human_prints_warning_count_when_present() {
        let report = ParseReport {
            warnings: vec!["skipped *HEADING".into()],
        };
        let r = info_report(&fake_model(), "demo", &report);
        let h = info_human(&r);

        assert!(h.contains("warnings      1"), "{h}");
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

    // ---- issues #46/#47: auto-detected thermal and dynamics runs ----

    /// 2x1 CPS4 bar, 0.1 x 0.05: clamped left edge, HOT=100 at x=0,
    /// COLD=0 at x=0.1. T(x) = 100(1 - x/0.1).
    const THERMAL_DECK: &str = "\
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
*STEP, NAME=HEAT
*HEAT TRANSFER, STEADY STATE
*TEMPERATURE
HOT, 100.0
COLD, 0.0
*BOUNDARY
1, 1, 2
4, 1, 2
*END STEP
";

    /// Same mesh clamped at the left, tip load, transient Newmark:
    /// dt = 1e-3 over 0.01 s -> 10 steps, 11 output points.
    const DYNAMIC_DECK: &str = "\
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
1.0E-3, 0.01
*BOUNDARY
1, 1, 2
4, 1, 2
*CLOAD
3, 2, -100.0
6, 2, -100.0
*END STEP
";

    #[test]
    fn auto_detection_helpers() {
        let (model, _) = solvix_io::parse_inp(THERMAL_DECK).unwrap();
        assert!(has_thermal(&model));
        assert!(!has_dynamics(&model));
        let (model, _) = solvix_io::parse_inp(DYNAMIC_DECK).unwrap();
        assert!(has_dynamics(&model));
        assert!(!has_thermal(&model));
        let (model, _) = solvix_io::parse_inp(
            "*NODE\n1, 0.0, 0.0\n2, 1.0, 0.0\n3, 0.0, 1.0\n*ELEMENT, TYPE=CPS3, ELSET=A\n1, 1, 2, 3\n*MATERIAL, NAME=M\n*ELASTIC\n1.0, 0.0\n*SOLID SECTION, ELSET=A, MATERIAL=M\n1.0\n*STEP\n*STATIC\n*BOUNDARY\n1, 1, 2\n*END STEP\n",
        )
        .unwrap();
        assert!(!has_thermal(&model));
        assert!(!has_dynamics(&model));
    }

    #[test]
    fn auto_thermal_run_json_and_human() {
        let (model, _) = solvix_io::parse_inp(THERMAL_DECK).unwrap();
        let outcome = solve_auto(&model).unwrap();
        let AutoOutcome::Static { thermal, .. } = &outcome else {
            panic!("thermal deck must take the thermomech path");
        };
        let t = thermal.expect("thermal info present");
        assert!((t.max_c - 100.0).abs() < 1e-6, "max {t:?}");
        assert!(t.min_c.abs() < 1e-6, "min {t:?}");
        let s = run_summary(
            &model,
            &outcome,
            "demo",
            7,
            vec!["demo.vtu".into()],
            Vec::new(),
        );
        let j = success_json(&s);
        assert!(
            j.contains(r#""thermal":true,"max_temperature_c":100.000"#),
            "{j}"
        );
        assert!(!j.contains("dynamic"), "{j}");
        let h = run_human(&s);
        assert!(h.contains("analysis      thermal-stress (auto)"), "{h}");
        assert!(h.contains("max temp      100.000 C"), "{h}");
        assert!(h.contains("min temp      0.000 C"), "{h}");
    }

    #[test]
    fn auto_dynamic_run_json_and_human() {
        let (model, _) = solvix_io::parse_inp(DYNAMIC_DECK).unwrap();
        let outcome = solve_auto(&model).unwrap();
        let AutoOutcome::Dynamic {
            info, final_sol, ..
        } = &outcome
        else {
            panic!("dynamic deck must take the transient path");
        };
        assert_eq!(info.steps, 10);
        assert_eq!(final_sol.displacements.len(), model.mesh.nodes.len());
        // Undamped step load: the oscillation must exceed the static
        // equilibrium (max over time beats the final value somewhere in
        // the window) and be non-zero.
        assert!(info.max_disp_m > 0.0);
        assert!(info.max_disp_m >= info.final_disp_m, "{info:?}");
        let s = run_summary(
            &model,
            &outcome,
            "demo",
            7,
            vec!["demo.vtu".into()],
            Vec::new(),
        );
        assert_eq!(s.max_disp_m, info.max_disp_m, "summary max is over time");
        assert_eq!(s.solver_method, "newmark");
        let j = success_json(&s);
        assert!(
            j.contains(r#""dynamic":true,"steps":10,"dt":0.001,"duration":0.01,"final_disp_mm":"#),
            "{j}"
        );
        assert!(!j.contains("thermal"), "{j}");
        let h = run_human(&s);
        assert!(h.contains("analysis      dynamics, Newmark (auto)"), "{h}");
        assert!(h.contains("time steps    10"), "{h}");
        assert!(h.contains("dt            0.001 s"), "{h}");
        assert!(h.contains("duration      0.01 s"), "{h}");
        assert!(h.contains("final |disp|  "), "{h}");
    }

    #[test]
    fn auto_static_run_has_no_new_fields() {
        let (model, _) = solvix_io::parse_inp(
            "*NODE\n1, 0.0, 0.0\n2, 1.0, 0.0\n3, 0.0, 1.0\n*ELEMENT, TYPE=CPS3, ELSET=A\n1, 1, 2, 3\n*MATERIAL, NAME=M\n*ELASTIC\n1.0, 0.0\n*SOLID SECTION, ELSET=A, MATERIAL=M\n1.0\n*STEP\n*STATIC\n*BOUNDARY\n1, 1, 2\n*CLOAD\n2, 1, 1.0\n*END STEP\n",
        )
        .unwrap();
        let outcome = solve_auto(&model).unwrap();
        assert!(matches!(outcome, AutoOutcome::Static { thermal: None, .. }));
        let s = run_summary(
            &model,
            &outcome,
            "demo",
            1,
            vec!["demo.vtu".into()],
            Vec::new(),
        );
        let j = success_json(&s);
        assert!(!j.contains("thermal"), "{j}");
        assert!(!j.contains("dynamic"), "{j}");
    }

    #[test]
    fn fmt_secs_plain_and_scientific() {
        assert_eq!(fmt_secs(0.1), "0.1");
        assert_eq!(fmt_secs(0.001), "0.001");
        assert_eq!(fmt_secs(1e-4), "1.000e-4");
        assert_eq!(fmt_secs(1.5), "1.5");
        assert_eq!(fmt_secs(0.0), "0");
        assert_eq!(fmt_secs(f64::NAN), "0");
    }

    // ---- issue #32: .solvix output database ----

    /// Same mesh as THERMAL_DECK, NLGEOM static step.
    const NLGEOM_DECK: &str = "\
*HEADING
    nl fixture
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
*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL
0.01
*STEP, NLGEOM=YES
*STATIC
*BOUNDARY
1, 1, 2
4, 1, 2
*CLOAD
3, 2, -500.0
6, 2, -500.0
*END STEP
";

    fn tmp_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "solvix_cli_odb_{name}_{}.solvix",
            std::process::id()
        ))
    }

    #[test]
    fn args_odb() {
        // Bare odb defaults to inspect (frame: None).
        assert_eq!(
            parse_args(&sv(&["odb", "job.solvix"])),
            Ok(Command::Odb {
                file: "job.solvix".into(),
                frame: None,
                json: false
            })
        );
        assert_eq!(
            parse_args(&sv(&["odb", "job.solvix", "--inspect", "--json"])),
            Ok(Command::Odb {
                file: "job.solvix".into(),
                frame: None,
                json: true
            })
        );
        assert_eq!(
            parse_args(&sv(&["odb", "job.solvix", "--frame", "3", "--json"])),
            Ok(Command::Odb {
                file: "job.solvix".into(),
                frame: Some(3),
                json: true
            })
        );
        assert!(parse_args(&sv(&["odb"])).is_err());
        assert!(parse_args(&sv(&["odb", "f.solvix", "--frame"])).is_err());
        let err = parse_args(&sv(&["odb", "f.solvix", "--frame", "x"])).unwrap_err();
        assert!(err.contains("invalid --frame value"), "{err}");
        let err = parse_args(&sv(&["odb", "f.solvix", "--inspect", "--frame", "1"])).unwrap_err();
        assert!(err.contains("mutually exclusive"), "{err}");
        let err = parse_args(&sv(&["odb", "f.solvix", "--wat"])).unwrap_err();
        assert!(err.contains("unknown flag"), "{err}");
        let err = parse_args(&sv(&["odb", "a.solvix", "b.solvix"])).unwrap_err();
        assert!(err.contains("unexpected argument"), "{err}");
    }

    #[test]
    fn args_run_odb_flag() {
        assert_eq!(
            parse_args(&sv(&["run", "job.inp", "--odb", "--json"])),
            Ok(Command::Run {
                inp: "job.inp".into(),
                outdir: ".".into(),
                modes: None,
                odb: true,
                json: true
            })
        );
        let err = parse_args(&sv(&["info", "a.inp", "--odb"])).unwrap_err();
        assert!(err.contains("unknown flag"), "{err}");
    }

    #[test]
    fn dyn_frame_indices_subsamples() {
        assert_eq!(dyn_frame_indices(0), Vec::<usize>::new());
        assert_eq!(dyn_frame_indices(11), (0..11).collect::<Vec<_>>());
        let idx = dyn_frame_indices(1001);
        assert_eq!(idx.len(), ODB_MAX_FRAMES);
        assert_eq!(idx.first(), Some(&0));
        assert_eq!(idx.last(), Some(&1000));
        assert!(idx.windows(2).all(|w| w[0] < w[1]));
        let idx = dyn_frame_indices(201);
        assert_eq!(idx.len(), ODB_MAX_FRAMES);
        assert_eq!(idx.last(), Some(&200));
    }

    #[test]
    fn odb_static_one_frame() {
        let model = fake_model();
        let frames = odb_frames(&model, &static_outcome());
        assert_eq!(frames.len(), 1);
        let f = &frames[0];
        assert_eq!(f.load_fraction, 1.0);
        assert_eq!(f.displacements, fake_solution().displacements);
        assert_eq!(f.von_mises.len(), 1);
        assert!(f.temperatures.is_none() && f.velocities.is_none());
    }

    #[test]
    fn odb_run_pipeline_static_round_trip() {
        let path = tmp_path("static");
        let model = fake_model();
        write_job_odb(&path, &model, &static_outcome()).unwrap();
        let odb = Odb::open(&path).unwrap();
        assert_eq!(odb.frame_count(), 1);
        assert_eq!(odb.model().nodes.len(), 2);
        let f = odb.frame(0).unwrap();
        assert_eq!(f.displacements[1][0].to_bits(), 0.0047619f64.to_bits());
        assert!(!odb.has_temperatures() && !odb.has_velocities());
        assert_eq!(odb.history(1).unwrap().len(), 1);
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn odb_thermal_stores_temperatures() {
        let (model, _) = solvix_io::parse_inp(THERMAL_DECK).unwrap();
        let outcome = solve_auto(&model).unwrap();
        let frames = odb_frames(&model, &outcome);
        assert_eq!(frames.len(), 1);
        let temps = frames[0].temperatures.as_ref().unwrap();
        assert_eq!(temps.len(), model.mesh.nodes.len());
        assert!((temps.iter().copied().fold(f64::NEG_INFINITY, f64::max) - 100.0).abs() < 1e-6);
    }

    #[test]
    fn odb_nl_one_frame_per_increment() {
        let (model, _) = solvix_io::parse_inp(NLGEOM_DECK).unwrap();
        let outcome = solve_auto(&model).unwrap();
        let AutoOutcome::Static { nl_frames, .. } = &outcome else {
            panic!("NLGEOM deck must keep its frames");
        };
        assert_eq!(nl_frames.len(), 20, "default max_increments = 20");
        assert!((nl_frames.last().unwrap().load_fraction - 1.0).abs() < 1e-12);
        let path = tmp_path("nl");
        write_job_odb(&path, &model, &outcome).unwrap();
        let odb = Odb::open(&path).unwrap();
        assert_eq!(odb.frame_count(), 20);
        let f = odb.frame(19).unwrap();
        assert!((f.load_fraction - 1.0).abs() < 1e-12);
        let hist = odb.history(2).unwrap();
        assert_eq!(hist.len(), 20);
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn odb_dynamic_keeps_small_runs_and_velocities() {
        let (model, _) = solvix_io::parse_inp(DYNAMIC_DECK).unwrap();
        let outcome = solve_auto(&model).unwrap();
        // 11 output points stay 11 frames (no sub-sampling needed).
        let frames = odb_frames(&model, &outcome);
        assert_eq!(frames.len(), 11);
        let last = &frames[10];
        assert!((last.time - 0.01).abs() < 1e-12);
        assert!((last.load_fraction - 1.0).abs() < 1e-12);
        let v = last.velocities.as_ref().unwrap();
        assert_eq!(v.len(), model.mesh.nodes.len());
        assert!(last.von_mises.is_empty(), "v1 dynamics has no stresses");
    }

    #[test]
    fn inspect_json_exact() {
        let r = OdbInspect {
            file: "job.solvix".into(),
            frames: 11,
            fields: vec![
                "displacements",
                "von_mises",
                "nodal_von_mises",
                "temperatures",
            ],
            nodes: 6,
            elements: 2,
            bytes: 4096,
        };
        assert_eq!(
            inspect_json(&r),
            r#"{"ok":true,"file":"job.solvix","version":1,"nodes":6,"elements":2,"frames":11,"fields":["displacements","von_mises","nodal_von_mises","temperatures"],"bytes":4096}"#
        );
        let h = inspect_human(&r);
        assert!(h.contains("frames        11"), "{h}");
        assert!(
            h.contains("fields        displacements, von_mises, nodal_von_mises, temperatures"),
            "{h}"
        );
    }

    #[test]
    fn frame_json_exact() {
        let s = FrameSummary {
            frame: 3,
            time: 0.003,
            load_fraction: 0.3,
            max_disp_m: 0.0047619,
            max_von_mises_pa: 38.2e6,
            max_nodal_von_mises_pa: 40.0e6,
            max_temperature: Some(100.0),
            max_velocity_m: None,
        };
        assert_eq!(
            frame_json("job.solvix", &s),
            r#"{"ok":true,"file":"job.solvix","frame":3,"time":0.003,"load_fraction":0.3,"max_disp_mm":4.762,"max_von_mises_mpa":38.200,"max_nodal_von_mises_mpa":40.000,"max_temperature_c":100.000}"#
        );
        let h = frame_human("job.solvix", &s);
        assert!(h.contains("max |disp|    4.762 mm"), "{h}");
        assert!(h.contains("max temp      100.000 C"), "{h}");
        // Empty optional fields are omitted entirely.
        let bare = FrameSummary {
            max_temperature: None,
            max_velocity_m: Some(0.5),
            ..s
        };
        let j = frame_json("job.solvix", &bare);
        assert!(!j.contains("max_temperature_c"), "{j}");
        assert!(j.contains(r#""max_velocity_mm_s":500.000"#), "{j}");
    }

    #[test]
    fn frame_summary_maxes() {
        let f = OdbFrame {
            time: 1.0,
            load_fraction: 0.5,
            displacements: vec![[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]],
            von_mises: vec![1.0, 9.0, 4.0],
            nodal_von_mises: vec![2.0],
            temperatures: Some(vec![-5.0, 80.0]),
            velocities: Some(vec![[0.0, 0.0, 2.0]]),
        };
        let s = frame_summary(0, &f);
        assert_eq!(s.max_disp_m, 5.0);
        assert_eq!(s.max_von_mises_pa, 9.0);
        assert_eq!(s.max_nodal_von_mises_pa, 2.0);
        assert_eq!(s.max_temperature, Some(80.0));
        assert_eq!(s.max_velocity_m, Some(2.0));
    }

    #[test]
    fn do_odb_error_stages() {
        // Missing file: io.
        let (stage, msg) = do_odb("/nonexistent/x.solvix", None, true).unwrap_err();
        assert_eq!(stage, Stage::Io);
        assert!(msg.contains("cannot open"), "{msg}");
        // Garbage magic: parse.
        let path = tmp_path("garbage");
        std::fs::write(&path, b"definitely not an odb").unwrap();
        let (stage, _) = do_odb(&path.display().to_string(), None, true).unwrap_err();
        assert_eq!(stage, Stage::Parse);
        // Valid file, frame out of range: usage.
        write_job_odb(&path, &fake_model(), &static_outcome()).unwrap();
        let (stage, msg) = do_odb(&path.display().to_string(), Some(5), true).unwrap_err();
        assert_eq!(stage, Stage::Usage);
        assert!(msg.contains("out of range"), "{msg}");
        // Success paths print and return Ok.
        assert!(do_odb(&path.display().to_string(), None, true).is_ok());
        assert!(do_odb(&path.display().to_string(), Some(0), true).is_ok());
        std::fs::remove_file(&path).unwrap();
    }
}
