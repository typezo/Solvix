//! Issue #35 acceptance tests: three canonical Abaqus tutorial scripts
//! translated, round-tripped through `solvix_io::parse_inp` and solved
//! against benchmark tolerances.

use solvix_core::solve;
use solvix_io::parse_inp;
use solvix_tools::deck::write_deck;
use solvix_tools::translate::translate;

const E: f64 = 210e9;

fn tip_deflection(model: &solvix_core::Model, x: f64, y: f64) -> f64 {
    let sol = solve(model).unwrap();
    let n = model
        .mesh
        .nodes
        .iter()
        .find(|n| (n.coords[0] - x).abs() < 1e-9 && (n.coords[1] - y).abs() < 1e-9)
        .unwrap_or_else(|| panic!("no node at ({x}, {y})"));
    sol.displacements[model.mesh.node_index(n.id).unwrap()][1]
}

/// An independently written 40x4 CPS4 cantilever deck (same geometry,
/// load and BCs as fixtures/cantilever.py at seed 0.025): 41x5 nodes,
/// clamped left edge, -1000 N at the middle tip node.
fn handwritten_cantilever_deck() -> String {
    let (nx, ny, l, h) = (40usize, 4usize, 1.0f64, 0.1f64);
    let id = |i: usize, j: usize| (j * (nx + 1) + i + 1) as u32;
    let mut d = String::from("*HEADING\n handwritten cantilever\n*NODE\n");
    for j in 0..=ny {
        for i in 0..=nx {
            d.push_str(&format!(
                "{}, {}, {}\n",
                id(i, j),
                i as f64 * l / nx as f64,
                j as f64 * h / ny as f64
            ));
        }
    }
    d.push_str("*ELEMENT, TYPE=CPS4, ELSET=ALL\n");
    let mut eid = 0u32;
    for j in 0..ny {
        for i in 0..nx {
            eid += 1;
            d.push_str(&format!(
                "{}, {}, {}, {}, {}\n",
                eid,
                id(i, j),
                id(i + 1, j),
                id(i + 1, j + 1),
                id(i, j + 1)
            ));
        }
    }
    let fix: Vec<String> = (0..=ny).map(|j| id(0, j).to_string()).collect();
    d.push_str(&format!("*NSET, NSET=FIX\n{}\n", fix.join(", ")));
    d.push_str("*MATERIAL, NAME=STEEL\n*ELASTIC\n210e9, 0.3\n*DENSITY\n7850.0\n");
    d.push_str("*SOLID SECTION, ELSET=ALL, MATERIAL=STEEL\n0.01\n");
    d.push_str("*STEP\n*STATIC\n*BOUNDARY\nFIX, 1, 2\n*CLOAD\n");
    d.push_str(&format!("{}, 2, -1000.0\n", id(nx, 2)));
    d.push_str("*END STEP\n");
    d
}

#[test]
fn cantilever_matches_handwritten_deck_and_beam_theory() {
    let src = include_str!("fixtures/cantilever.py");
    let tr = translate(src).unwrap();
    assert!(
        tr.report.unsupported.is_empty(),
        "canonical cantilever should be fully supported: {:?}",
        tr.report.unsupported
    );

    // Round-trip: the generated deck must re-parse to the same model.
    let deck = write_deck(&tr.model);
    let (m2, report) = parse_inp(&deck).unwrap();
    assert!(report.warnings.is_empty(), "{:?}", report.warnings);
    assert_eq!(m2.mesh.nodes.len(), tr.model.mesh.nodes.len());
    assert_eq!(m2.mesh.elements.len(), tr.model.mesh.elements.len());
    assert_eq!(m2.steps.len(), 1);
    assert_eq!(m2.thickness, 0.01);
    assert_eq!(m2.density, 7850.0);
    assert_eq!(m2.steps[0].cloads.len(), 1);
    assert_eq!(m2.steps[0].cloads[0].value, -1000.0);
    assert_eq!(m2.steps[0].boundary_conditions.len(), 1);
    assert_eq!(m2.steps[0].boundary_conditions[0].dofs, vec![1, 2]);

    // Solve the translated model: tip deflection within 5% of the same
    // handwritten deck.
    let translated = tip_deflection(&tr.model, 1.0, 0.05);
    let (hw, _) = parse_inp(&handwritten_cantilever_deck()).unwrap();
    let handwritten = tip_deflection(&hw, 1.0, 0.05);
    assert!(
        ((translated - handwritten) / handwritten).abs() < 0.05,
        "translated {translated:.6e} vs handwritten {handwritten:.6e}"
    );

    // Physical sanity: Euler-Bernoulli within 15% (shear locking on a
    // 40x4 CPS4 grid, same tolerance as the solvix-core benchmark).
    let inertia = 0.01 * 0.1f64.powi(3) / 12.0;
    let beam = 1000.0 / (3.0 * E * inertia);
    assert!(
        (translated.abs() - beam).abs() / beam < 0.15,
        "tip {:.6e} vs beam theory {beam:.6e}",
        translated.abs()
    );
}

#[test]
fn plate_with_hole_kirsch() {
    let src = include_str!("fixtures/plate_hole.py");
    let tr = translate(src).unwrap();
    assert!(
        tr.report.unsupported.is_empty(),
        "{:?}",
        tr.report.unsupported
    );

    // Round-trip through the deck parser, then solve the re-parsed deck
    // (full pipeline, exactly what --validate does).
    let deck = write_deck(&tr.model);
    let (m2, report) = parse_inp(&deck).unwrap();
    assert!(report.warnings.is_empty(), "{:?}", report.warnings);
    let sol = solve(&m2).unwrap();

    // Kirsch: far-field 1 MPa uniaxial -> 3 MPa von Mises at the hole
    // edge (coarse faceted hole, 20% tolerance).
    let peak = sol.nodal_von_mises.iter().cloned().fold(0.0f64, f64::max);
    assert!(
        (peak - 3.0e6).abs() / 3.0e6 < 0.20,
        "Kirsch peak {peak:.4e} Pa, want 3.0 MPa +- 20%"
    );
    // Pressure became element face loads on the right edge.
    assert!(!m2.steps[0].dloads.is_empty());
    assert!(m2.steps[0].dloads.iter().all(|d| d.face == 2));
}

#[test]
fn bracket_multi_feature_reports_unsupported_and_solves() {
    let src = include_str!("fixtures/bracket.py");
    let tr = translate(src).unwrap();

    // Fillet and element-type override must be surfaced, not silently
    // dropped.
    let joined = tr.report.unsupported.join("\n");
    assert!(joined.contains("Fillet"), "{joined}");
    assert!(joined.contains("ElemType"), "{joined}");

    let deck = write_deck(&tr.model);
    let (m2, report) = parse_inp(&deck).unwrap();
    assert!(report.warnings.is_empty(), "{:?}", report.warnings);

    // Two steps; the second carries the first step's BCs and loads
    // (accumulated) plus gravity and the top pressure.
    assert_eq!(m2.steps.len(), 2);
    let s2 = &m2.steps[1];
    assert_eq!(s2.name, "Gravity");
    assert_eq!(s2.boundary_conditions.len(), 1); // the EncastreBC, carried
    assert_eq!(s2.cloads.len(), 2); // cf1 + cf2, carried
    let (g, dir) = s2.gravity.unwrap();
    assert!((g - 9.81).abs() < 1e-9);
    assert!((dir[1] + 1.0).abs() < 1e-9);
    // Top-edge pressure became P3 face loads; gravity adds GRAV.
    assert!(!s2.dloads.is_empty());
    assert!(s2.dloads.iter().all(|d| d.face == 3));

    let sol = solve(&m2).unwrap();
    assert!(sol.final_residual <= 1e-10);
}

#[test]
fn cli_smoke_validate_json() {
    let bin = env!("CARGO_BIN_EXE_abaqus2solvix");
    let out = std::process::Command::new(bin)
        .arg("tests/fixtures/cantilever.py")
        .arg("-o")
        .arg(concat!(env!("CARGO_TARGET_TMPDIR"), "/cant_cli.inp"))
        .arg("--json")
        .arg("--validate")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("run abaqus2solvix");
    assert!(out.status.success(), "status: {out:?}");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.starts_with(r#"{"ok":true,"output":"#), "{stdout}");
    assert!(
        stdout.contains(r#""validate":{"parsed":true,"solved":true"#),
        "{stdout}"
    );
}
