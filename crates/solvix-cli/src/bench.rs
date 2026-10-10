//! Built-in verification suite: `solvix bench` (issue #17).
//!
//! Cases are embedded decks built by pure functions (loops, not literals, and
//! no dependency on `examples/`), so unit tests exercise parse+solve+gate
//! in-process without spawning the binary. Each case gates one key metric
//! against an analytical reference:
//!
//! * `cantilever40x4` — Euler-Bernoulli tip deflection P L³/(3 E I) = 15.238 mm.
//! * `kirsch` — quarter plate with a hole, far-field 1 MPa, SCF 3 → 3.0 MPa.
//!   Same radial 8×8 pattern as `examples/plate_hole.inp`; coarse element
//!   centroid sampling underpredicts the hole-surface peak (documented in
//!   docs/benchmarks/kirsch.md).
//! * `patch_cps4` — distorted 9-node patch under uniform sigma_x; every
//!   element centroid von Mises must equal the applied traction.

use solvix_core::{Model, Solution};

use crate::{fin, json_escape, Stage};

/// Metric extractor: value in display unit (mm / MPa).
type Extract = fn(&Model, &Solution) -> f64;

/// One benchmark case as data: deck builder + gate definition.
pub struct CaseDef {
    pub name: &'static str,
    pub metric: &'static str,
    pub unit: &'static str,
    pub expected: f64,
    pub tol_pct: f64,
    pub deck: fn() -> String,
    pub extract: Extract,
}

pub const CASES: &[CaseDef] = &[
    CaseDef {
        name: "cantilever40x4",
        metric: "max_disp_mm",
        unit: "mm",
        expected: 15.238,
        tol_pct: 12.0,
        deck: cantilever_deck,
        extract: max_disp_mm,
    },
    CaseDef {
        name: "kirsch",
        metric: "max_nodal_von_mises_mpa",
        unit: "MPa",
        expected: 3.0,
        tol_pct: 15.0,
        deck: kirsch_deck,
        extract: max_nodal_von_mises_mpa,
    },
    CaseDef {
        name: "patch_cps4",
        metric: "max_centroid_von_mises_mpa",
        unit: "MPa",
        expected: 1.0,
        tol_pct: 2.0,
        deck: patch_deck,
        extract: max_centroid_von_mises_mpa,
    },
];

#[derive(Clone)]
pub struct CaseOutcome {
    pub name: &'static str,
    pub metric: &'static str,
    pub unit: &'static str,
    pub nodes: usize,
    pub elements: usize,
    pub iterations: usize,
    pub value: f64,
    pub expected: f64,
    pub tol_pct: f64,
    pub pct_error: f64,
    pub pass: bool,
}

fn max_disp_mm(_model: &Model, sol: &Solution) -> f64 {
    sol.displacements
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
        .fold(0.0, f64::max)
        * 1000.0
}

fn max_centroid_von_mises_mpa(_model: &Model, sol: &Solution) -> f64 {
    sol.von_mises.iter().copied().fold(0.0, f64::max) / 1e6
}

fn max_nodal_von_mises_mpa(_model: &Model, sol: &Solution) -> f64 {
    sol.nodal_von_mises.iter().copied().fold(0.0, f64::max) / 1e6
}

/// Relative error in percent and the pass/fail verdict; the boundary value
/// (|pct| == tol) passes.
pub fn gate(value: f64, expected: f64, tol_pct: f64) -> (f64, bool) {
    let pct = if expected != 0.0 {
        (value - expected) / expected * 100.0
    } else {
        f64::NAN
    };
    (pct, pct.is_finite() && pct.abs() <= tol_pct + 1e-9)
}

pub fn run_case(def: &CaseDef) -> Result<CaseOutcome, String> {
    let (model, _report) = solvix_io::parse_inp(&(def.deck)())
        .map_err(|p| format!("line {}: {}", p.line, p.message))?;
    let sol = solvix_core::solve(&model).map_err(crate::solve_error_message)?;
    let value = (def.extract)(&model, &sol);
    let (pct_error, pass) = gate(value, def.expected, def.tol_pct);
    Ok(CaseOutcome {
        name: def.name,
        metric: def.metric,
        unit: def.unit,
        nodes: model.mesh.nodes.len(),
        elements: model.mesh.elements.len(),
        iterations: sol.iterations,
        value,
        expected: def.expected,
        tol_pct: def.tol_pct,
        pct_error,
        pass,
    })
}

pub fn run(json: bool) -> i32 {
    let mut outcomes = Vec::with_capacity(CASES.len());
    for def in CASES {
        match run_case(def) {
            Ok(o) => outcomes.push(o),
            Err(msg) => {
                let msg = format!("bench case '{}': {msg}", def.name);
                if json {
                    println!("{}", crate::error_json(Stage::Solve, &msg));
                } else {
                    eprintln!("error (solve): {msg}");
                }
                return Stage::Solve.exit_code();
            }
        }
    }
    let all_pass = outcomes.iter().all(|o| o.pass);
    if json {
        println!("{}", json_line(&outcomes));
    } else {
        print!("{}", human(&outcomes));
    }
    if all_pass {
        0
    } else {
        Stage::Solve.exit_code()
    }
}

pub fn human(outcomes: &[CaseOutcome]) -> String {
    let mut out = format!(
        "{:<14}{:>6}{:>7}{:>7}  {:<26}{:>11}{:>11}{:>9}  {}\n",
        "case", "nodes", "elems", "iter", "metric", "value", "expected", "error", "verdict"
    );
    for o in outcomes {
        out.push_str(&format!(
            "{:<14}{:>6}{:>7}{:>7}  {:<26}{:>11}{:>11}{:>9}  {}\n",
            o.name,
            o.nodes,
            o.elements,
            o.iterations,
            o.metric,
            format!("{:.3} {}", fin(o.value), o.unit),
            format!("{:.3}", fin(o.expected)),
            format!("{:+.2}%", fin(o.pct_error)),
            if o.pass { "PASS" } else { "FAIL" },
        ));
    }
    let passed = outcomes.iter().filter(|o| o.pass).count();
    out.push_str(&format!("{}/{} passed\n", passed, outcomes.len()));
    out
}

pub fn json_line(outcomes: &[CaseOutcome]) -> String {
    let cases = outcomes
        .iter()
        .map(|o| {
            format!(
                "{{\"name\":\"{}\",\"metric\":\"{}\",\"value\":{:.3},\"expected\":{:.3},\"tolerance\":{:.1},\"pct_error\":{:.2},\"pass\":{}}}",
                json_escape(o.name),
                json_escape(o.metric),
                fin(o.value),
                fin(o.expected),
                o.tol_pct,
                fin(o.pct_error),
                o.pass,
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"ok\":{},\"cases\":[{}]}}",
        outcomes.iter().all(|o| o.pass),
        cases
    )
}

fn material_and_step(s: &mut String, thickness: f64) {
    s.push_str("*MATERIAL, NAME=STEEL\n*ELASTIC\n210e9, 0.3\n");
    s.push_str(&format!(
        "*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL\n{thickness}\n"
    ));
    s.push_str("*STEP\n*STATIC\n");
}

fn close_step(s: &mut String) {
    s.push_str("*END STEP\n");
}

/// 40×4 CPS4 cantilever, L=2.0 m, h=0.1 m, t=0.01 m; left edge clamped,
/// five tip nodes carry −200 N each (1000 N total).
fn cantilever_deck() -> String {
    let (nx, ny) = (40usize, 4usize);
    let (l, h) = (2.0_f64, 0.1_f64);
    let (wx, wy) = (l / nx as f64, h / ny as f64);
    let id = |i: usize, j: usize| (j * (nx + 1) + i + 1) as u32;
    let mut s = String::from(
        "** 40x4 CPS4 cantilever, SI units (m, Pa, N).\n\
         ** Euler-Bernoulli reference: P*L^3/(3*E*I) = 15.238 mm.\n\
         *HEADING\n    bench cantilever40x4\n*NODE\n",
    );
    for j in 0..=ny {
        for i in 0..=nx {
            s.push_str(&format!(
                "{}, {:.6}, {:.6}\n",
                id(i, j),
                i as f64 * wx,
                j as f64 * wy
            ));
        }
    }
    s.push_str("*ELEMENT, TYPE=CPS4, ELSET=ALL\n");
    for j in 0..ny {
        for i in 0..nx {
            s.push_str(&format!(
                "{}, {}, {}, {}, {}\n",
                j * nx + i + 1,
                id(i, j),
                id(i + 1, j),
                id(i + 1, j + 1),
                id(i, j + 1)
            ));
        }
    }
    let fix: Vec<String> = (0..=ny).map(|j| id(0, j).to_string()).collect();
    let tip: Vec<String> = (0..=ny).map(|j| id(nx, j).to_string()).collect();
    s.push_str(&format!("*NSET, NSET=FIX\n{}\n", fix.join(", ")));
    s.push_str(&format!("*NSET, NSET=TIP\n{}\n", tip.join(", ")));
    s.push_str(&format!("*ELSET, ELSET=ALL, GENERATE\n1, {}, 1\n", nx * ny));
    material_and_step(&mut s, 0.01);
    s.push_str("*BOUNDARY\nFIX, 1, 2\n*CLOAD\n");
    for id in tip {
        s.push_str(&format!("{id}, 2, -200.0\n"));
    }
    close_step(&mut s);
    s
}

/// Quarter plate with a circular hole: the radial 8×8 pattern of
/// examples/plate_hole.inp — rings at r = 0.2…1.0 m step 0.1, sectors
/// θ = 0…90° step π/16. Symmetry BCs on the axes; consistent uniaxial
/// far-field traction (total 1e4 N) gives nominal 1 MPa far field.
fn kirsch_deck() -> String {
    let (nr, nt) = (8usize, 8usize);
    let (r0, r1) = (0.2_f64, 1.0_f64);
    let id = |i: usize, j: usize| (j * (nt + 1) + i + 1) as u32;
    let mut s = String::from(
        "** Quarter plate with a circular hole, SI units (m, Pa, N).\n\
         ** Kirsch reference: 3.0 MPa at the hole for 1 MPa far field.\n\
         *HEADING\n    bench kirsch\n*NODE\n",
    );
    for j in 0..=nr {
        let r = r0 + (r1 - r0) * j as f64 / nr as f64;
        for i in 0..=nt {
            let th = i as f64 * (std::f64::consts::FRAC_PI_2) / nt as f64;
            s.push_str(&format!(
                "{}, {:.6}, {:.6}\n",
                id(i, j),
                r * th.cos(),
                r * th.sin()
            ));
        }
    }
    s.push_str("*ELEMENT, TYPE=CPS4, ELSET=ALL\n");
    for j in 0..nr {
        for i in 0..nt {
            s.push_str(&format!(
                "{}, {}, {}, {}, {}\n",
                j * nt + i + 1,
                id(i, j),
                id(i, j + 1),
                id(i + 1, j + 1),
                id(i + 1, j)
            ));
        }
    }
    let ysym: Vec<String> = (0..=nr).map(|j| id(0, j).to_string()).collect();
    let xsym: Vec<String> = (0..=nr).map(|j| id(nt, j).to_string()).collect();
    let far: Vec<String> = (0..=nt).map(|i| id(i, nr).to_string()).collect();
    s.push_str(&format!("*NSET, NSET=YSYM\n{}\n", ysym.join(", ")));
    s.push_str(&format!("*NSET, NSET=XSYM\n{}\n", xsym.join(", ")));
    s.push_str(&format!("*NSET, NSET=FAR_EDGE\n{}\n", far.join(", ")));
    s.push_str(&format!("*ELSET, ELSET=ALL, GENERATE\n1, {}, 1\n", nr * nt));
    material_and_step(&mut s, 0.01);
    s.push_str("*BOUNDARY\nYSYM, 2, 2\nXSYM, 1, 1\n*CLOAD\n");
    // Consistent far-field uniaxial traction: t = (sigma*cos(theta), 0) on
    // the arc (uniaxial sigma_x state: no y-component), per-node
    // F = sigma*t*∫N_i·cos(theta) ds, 2-point Gauss per arc segment.
    // Total Fx = sigma*t*r = 1e4 N.
    let (sigma, t) = (1.0e6_f64, 0.01_f64);
    let mut fx = vec![0.0_f64; nt + 1];
    let gauss = [(-1.0 / 3.0_f64.sqrt(), 1.0), (1.0 / 3.0_f64.sqrt(), 1.0)];
    for i in 0..nt {
        let (ta, tb) = (
            i as f64 * std::f64::consts::FRAC_PI_2 / nt as f64,
            (i + 1) as f64 * std::f64::consts::FRAC_PI_2 / nt as f64,
        );
        for (s_gp, w) in gauss {
            let theta = 0.5 * (ta + tb) + 0.5 * (tb - ta) * s_gp;
            let ds = r1 * 0.5 * (tb - ta) * w;
            let nx = theta.cos();
            fx[i] += sigma * t * 0.5 * (1.0 - s_gp) * nx * ds;
            fx[i + 1] += sigma * t * 0.5 * (1.0 + s_gp) * nx * ds;
        }
    }
    for i in 0..=nt {
        s.push_str(&format!("{}, 1, {:.6}\n", id(i, nr), fx[i]));
    }
    close_step(&mut s);
    s
}

/// 9-node distorted unit-square patch: four CPS4 elements around an
/// off-center interior node. Consistent traction loads (1/4, 1/2, 1/4 of
/// σ·t·L = 10⁴ N — bilinear edges split linearly, the shared edge node
/// accumulates half from each element) on the right edge give an exact
/// uniform σx = 1 MPa; left edge ux = 0 with one uy pin.
fn patch_deck() -> String {
    let nodes: [(f64, f64); 9] = [
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (0.5, 0.0),
        (1.0, 0.5),
        (0.5, 1.0),
        (0.0, 0.5),
        (0.35, 0.55),
    ];
    let elems: [[u32; 4]; 4] = [[1, 5, 9, 8], [5, 2, 6, 9], [9, 6, 3, 7], [8, 9, 7, 4]];
    let mut s = String::from(
        "** Distorted unit-square CPS4 patch test, SI units (m, Pa, N).\n\
         ** Exact uniform sigma_x = 1 MPa; centroid von Mises must match.\n\
         *HEADING\n    bench patch_cps4\n*NODE\n",
    );
    for (k, (x, y)) in nodes.iter().enumerate() {
        s.push_str(&format!("{}, {:.6}, {:.6}\n", k + 1, x, y));
    }
    s.push_str("*ELEMENT, TYPE=CPS4, ELSET=ALL\n");
    for (k, e) in elems.iter().enumerate() {
        s.push_str(&format!(
            "{}, {}, {}, {}, {}\n",
            k + 1,
            e[0],
            e[1],
            e[2],
            e[3]
        ));
    }
    s.push_str("*ELSET, ELSET=ALL, GENERATE\n1, 4, 1\n");
    material_and_step(&mut s, 0.01);
    s.push_str("*BOUNDARY\n1, 1, 1\n8, 1, 1\n4, 1, 1\n1, 2, 2\n*CLOAD\n");
    s.push_str("2, 1, 2500.0\n6, 1, 5000.0\n3, 1, 2500.0\n");
    close_step(&mut s);
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn solved(name: &str) -> CaseOutcome {
        let def = CASES.iter().find(|d| d.name == name).unwrap();
        run_case(def).unwrap()
    }

    #[test]
    fn decks_parse_with_expected_sizes() {
        for (name, nodes, elements) in [
            ("cantilever40x4", 205, 160),
            ("kirsch", 81, 64),
            ("patch_cps4", 9, 4),
        ] {
            let def = CASES.iter().find(|d| d.name == name).unwrap();
            let (model, report) = solvix_io::parse_inp(&(def.deck)()).unwrap();
            assert_eq!(model.mesh.nodes.len(), nodes, "{name}");
            assert_eq!(model.mesh.elements.len(), elements, "{name}");
            assert!(report.warnings.is_empty(), "{name}: {:?}", report.warnings);
        }
    }

    #[test]
    fn all_cases_solve_to_finite_metrics() {
        for def in CASES {
            let o = run_case(def).unwrap();
            assert!(
                o.value.is_finite() && o.value > 0.0,
                "{}: {}",
                def.name,
                o.value
            );
            assert!(o.pct_error.is_finite(), "{}", def.name);
        }
    }

    #[test]
    fn cantilever_and_patch_pass_gates() {
        for name in ["cantilever40x4", "patch_cps4"] {
            let o = solved(name);
            assert!(o.pass, "{name}: value {} vs {}", o.value, o.expected);
        }
    }

    #[test]
    fn kirsch_deck_matches_example_recipe() {
        let d = kirsch_deck();
        assert!(
            d.contains("2, 0.196157, 0.039018\n"),
            "ring 0 sector 1 node"
        );
        assert!(
            d.contains("77, 0.707107, 0.707107\n"),
            "r=1.0 45-degree node"
        );
        assert!(
            d.contains("1, 1, 10, 11, 2\n"),
            "first element connectivity"
        );
    }

    #[test]
    fn gate_inclusive_boundary() {
        for (v, ok) in [
            (1.12, true),
            (0.88, true),
            (1.0, true),
            (1.1201, false),
            (0.8799, false),
        ] {
            let (pct, pass) = gate(v, 1.0, 12.0);
            assert_eq!(pass, ok, "value {v} -> pct {pct}");
        }
        let (pct, pass) = gate(1.0, 1.0, 12.0);
        assert_eq!((pct, pass), (0.0, true));
        let (_, pass) = gate(f64::NAN, 1.0, 12.0);
        assert!(!pass);
        let (_, pass) = gate(3.0, 0.0, 12.0);
        assert!(!pass, "zero expected must fail, not divide by zero");
    }

    #[test]
    fn human_table_and_summary() {
        let o = solved("patch_cps4");
        let h = human(&[o.clone(), o]);
        assert!(h.contains("patch_cps4"), "{h}");
        assert!(h.contains("PASS"), "{h}");
        assert!(h.contains("2/2 passed"), "{h}");
        assert!(h.ends_with("2/2 passed\n"), "{h}");
    }

    #[test]
    fn json_single_line_schema() {
        let o = solved("patch_cps4");
        let j = json_line(&[o]);
        assert!(!j.contains('\n'), "{j}");
        assert!(
            j.starts_with(r#"{"ok":true,"cases":[{"name":"patch_cps4""#),
            "{j}"
        );
        assert!(j.ends_with(r#""pass":true}]}"#), "{j}");
        assert!(
            j.contains(r#""metric":"max_centroid_von_mises_mpa""#),
            "{j}"
        );
        assert!(j.contains(r#""tolerance":2.0"#), "{j}");
        let failing = CaseOutcome {
            name: "x",
            metric: "m",
            unit: "mm",
            nodes: 1,
            elements: 1,
            iterations: 1,
            value: 1.0,
            expected: 2.0,
            tol_pct: 5.0,
            pct_error: -50.0,
            pass: false,
        };
        assert!(json_line(&[failing]).starts_with(r#"{"ok":false,"cases":[{"name":"x""#));
    }
}
