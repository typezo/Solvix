//! Issue #53 acceptance tests: the v0.8 geometry kernel (sketch
//! booleans + conformal Delaunay meshing) — mesh validity, hole
//! boundary recovery, L-bracket quality, the plate-with-hole Kirsch
//! benchmark, and the `--mesh2d` translation path.

use std::collections::HashSet;

use solvix_core::{
    solve, BoundaryCondition, ConcentratedLoad, ElemKind, ElementSet, Material, Model,
};
use solvix_io::parse_inp;
use solvix_tools::deck::write_deck;
use solvix_tools::geom::Primitive;
use solvix_tools::meshgen2::{mesh_sketch, to_cps3_mesh, MeshResult};
use solvix_tools::translate::{translate_with, Options};

fn line(x1: f64, y1: f64, x2: f64, y2: f64) -> Primitive {
    Primitive::Line { x1, y1, x2, y2 }
}

fn square_sketch() -> Vec<Primitive> {
    vec![
        line(0.0, 0.0, 1.0, 0.0),
        line(1.0, 0.0, 1.0, 1.0),
        line(1.0, 1.0, 0.0, 1.0),
        line(0.0, 1.0, 0.0, 0.0),
    ]
}

fn signed_area(m: &MeshResult) -> f64 {
    m.tris
        .iter()
        .map(|t| {
            0.5 * ((m.nodes[t[1]][0] - m.nodes[t[0]][0]) * (m.nodes[t[2]][1] - m.nodes[t[0]][1])
                - (m.nodes[t[1]][1] - m.nodes[t[0]][1]) * (m.nodes[t[2]][0] - m.nodes[t[0]][0]))
        })
        .sum()
}

/// Mesh validity: CCW triangles, every edge shared by at most two
/// triangles, Euler characteristic V - E + T = 1 - #holes.
fn check_valid(m: &MeshResult, holes: usize) {
    let mut incidences: Vec<(usize, usize)> = Vec::new();
    for t in &m.tris {
        let (a, b, c) = (m.nodes[t[0]], m.nodes[t[1]], m.nodes[t[2]]);
        let cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        assert!(cross > 0.0, "triangle not CCW: {t:?}");
        for i in 0..3 {
            let (p, q) = (t[i], t[(i + 1) % 3]);
            assert_ne!(p, q);
            incidences.push((p.min(q), p.max(q)));
        }
    }
    incidences.sort_unstable();
    incidences.dedup();
    for w in incidences.windows(3) {
        assert!(
            !(w[0] == w[2]),
            "edge {:?} shared by more than two triangles",
            w[0]
        );
    }
    let (n_tris, n_edges) = (m.tris.len(), incidences.len());
    let chi = m.nodes.len() as i64 - n_edges as i64 + n_tris as i64;
    assert_eq!(
        chi,
        (1 - holes as i64),
        "Euler characteristic: V={} E={} T={}",
        m.nodes.len(),
        n_edges,
        n_tris
    );
    assert_eq!(m.boundary.len(), 1 + holes);
}

#[test]
fn a_unit_square_mesh_valid_and_exact_area() {
    let m = mesh_sketch(&square_sketch(), 0.15).unwrap();
    check_valid(&m, 0);
    assert!(
        (signed_area(&m) - 1.0).abs() / 1.0 < 1e-3,
        "area {}",
        signed_area(&m)
    );
    assert!(m.quality.min_angle_deg > 15.0, "{:?}", m.quality);
    // CPS3 deck round-trip via the existing writer.
    let mesh = to_cps3_mesh(&m);
    assert_eq!(mesh.elements.len(), m.tris.len());
    assert!(mesh.elements.iter().all(|e| e.kind == ElemKind::Cps3));
}

#[test]
fn b_square_with_hole_recovers_boundary_exactly() {
    let mut sk = square_sketch();
    sk.push(Primitive::Circle {
        cx: 0.5,
        cy: 0.5,
        r: 0.2,
    });
    let m = mesh_sketch(&sk, 0.05).unwrap();
    check_valid(&m, 1);

    // Every hole-loop edge is present as a mesh edge.
    let mesh_edges: HashSet<(usize, usize)> = m
        .tris
        .iter()
        .flat_map(|t| {
            (0..3).map(move |i| {
                let (p, q) = (t[i], t[(i + 1) % 3]);
                (p.min(q), p.max(q))
            })
        })
        .collect();
    let hole = &m.boundary[1];
    let n = hole.len();
    for i in 0..n {
        let (a, b) = (hole[i], hole[(i + 1) % n]);
        assert!(
            mesh_edges.contains(&(a.min(b), a.max(b))),
            "hole edge {i} missing from the mesh"
        );
    }

    // Every hole node lies on the circle within 2% of r.
    for &v in hole {
        let d = ((m.nodes[v][0] - 0.5).hypot(m.nodes[v][1] - 0.5) - 0.2).abs();
        assert!(d < 0.02 * 0.2, "hole node off radius by {d}");
    }
    assert!(m.quality.min_angle_deg > 15.0, "{:?}", m.quality);
}

#[test]
fn c_l_bracket_no_slivers() {
    let sk = vec![
        line(0.0, 0.0, 1.0, 0.0),
        line(1.0, 0.0, 1.0, 0.4),
        line(1.0, 0.4, 0.4, 0.4),
        line(0.4, 0.4, 0.4, 1.0),
        line(0.4, 1.0, 0.0, 1.0),
        line(0.0, 1.0, 0.0, 0.0),
    ];
    let m = mesh_sketch(&sk, 0.08).unwrap();
    check_valid(&m, 0);
    assert!(
        (signed_area(&m) - 0.64).abs() / 0.64 < 1e-2,
        "area {}",
        signed_area(&m)
    );
    assert!(m.quality.min_angle_deg > 15.0, "{:?}", m.quality);
    assert!(m.quality.max_aspect < 20.0, "{:?}", m.quality);
}

#[test]
fn d_plate_hole_kirsch_benchmark() {
    // Quarter-symmetric tension strip with a hole: x, y in [0, 1]^2,
    // r = 0.2 at the origin, sigma = 1 MPa far field, t = 0.01.
    // Kirsch hole-edge peak von Mises = 3 MPa at (0, r).
    let sk = vec![
        line(0.2, 0.0, 1.0, 0.0),
        line(1.0, 0.0, 1.0, 1.0),
        line(1.0, 1.0, 0.0, 1.0),
        line(0.0, 1.0, 0.0, 0.2),
        Primitive::Arc {
            cx: 0.0,
            cy: 0.0,
            r: 0.2,
            a0: std::f64::consts::FRAC_PI_2,
            a1: 0.0,
        },
    ];
    let m = mesh_sketch(&sk, 0.02).unwrap();
    check_valid(&m, 0); // the hole is a boundary notch on the quarter model
    assert!(m.quality.min_angle_deg > 15.0, "{:?}", m.quality);

    let mut model = Model {
        mesh: to_cps3_mesh(&m),
        thickness: 0.01,
        materials: vec![Material {
            name: "Steel".into(),
            youngs_modulus: 210e9,
            poisson_ratio: 0.3,
        }],
        elsets: vec![ElementSet {
            name: "all".into(),
            element_ids: (1..=m.tris.len() as u32).collect(),
        }],
        mat_assignments: vec![("all".into(), "Steel".into())],
        ..Default::default()
    };

    // Symmetry: u2 = 0 on y = 0, u1 = 0 on x = 0 (hole mouths
    // excluded — no nodes there anyway).
    let (mut xs, mut ys): (Vec<_>, Vec<_>) = (Vec::new(), Vec::new());
    for n in &model.mesh.nodes {
        if n.coords[1].abs() < 1e-7 {
            xs.push(n.id);
        }
        if n.coords[0].abs() < 1e-7 {
            ys.push(n.id);
        }
    }
    model.boundary_conditions = vec![
        BoundaryCondition {
            node_ids: xs,
            dofs: vec![2],
        },
        BoundaryCondition {
            node_ids: ys,
            dofs: vec![1],
        },
    ];

    // Consistent x-traction on x = 1: F_i = sigma t (h_prev + h_next)/2.
    let mut right: Vec<_> = model
        .mesh
        .nodes
        .iter()
        .filter(|n| (n.coords[0] - 1.0).abs() < 1e-7)
        .map(|n| (n.id, n.coords[1]))
        .collect();
    right.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let sigma = 1e6f64;
    for w in right.windows(2) {
        let h = w[1].1 - w[0].1;
        for (id, _) in w {
            model.cloads.push(ConcentratedLoad {
                node_id: *id,
                dof: 1,
                value: sigma * 0.01 * h / 2.0,
            });
        }
    }

    let sol = solve(&model).unwrap();
    assert!(sol.final_residual <= 1e-10);
    let peak = sol.nodal_von_mises.iter().cloned().fold(0.0f64, f64::max);
    eprintln!(
        "KIRSCH mesh2d: peak = {peak:.6e} Pa ({:.2}% vs Kirsch 3.0 MPa); {} nodes / {} tris, min angle {:.1} deg, max aspect {:.2}, iterations {}",
        100.0 * (peak - 3.0e6) / 3.0e6,
        m.nodes.len(),
        m.tris.len(),
        m.quality.min_angle_deg,
        m.quality.max_aspect,
        sol.iterations
    );
    let rel = (peak - 3.0e6) / 3.0e6;
    assert!(
        rel.abs() < 0.12,
        "Kirsch peak {peak:.4e} Pa = {:.2}% off Kirsch 3.0 MPa (mesh: {} nodes / {} tris, min angle {:.1} deg)",
        100.0 * rel,
        m.nodes.len(),
        m.tris.len(),
        m.quality.min_angle_deg
    );
}

#[test]
fn e_translated_script_uses_conformal_mesher_with_mesh2d() {
    let src = include_str!("fixtures/plate_hole.py");
    let tr = translate_with(src, &Options { mesh2d: true }).unwrap();
    let joined = tr.report.warnings.join("\n");
    assert!(joined.contains("conformal mesh"), "{joined}");

    // CPS3 triangles, not the structured grid.
    assert!(!tr.model.mesh.elements.is_empty());
    assert!(tr
        .model
        .mesh
        .elements
        .iter()
        .all(|e| e.kind == ElemKind::Cps3));

    // Hole boundary recovered: nodes on the r = 0.2 circle, and no
    // element centroid inside it.
    let mut on_circle = 0;
    for n in &tr.model.mesh.nodes {
        let d = (n.coords[0].hypot(n.coords[1]) - 0.2).abs();
        if d < 1e-6 {
            on_circle += 1;
        }
    }
    assert!(on_circle >= 16, "only {on_circle} nodes on the hole circle");
    for e in &tr.model.mesh.elements {
        let (mut cx, mut cy) = (0.0, 0.0);
        for &id in &e.connectivity {
            let n = &tr.model.mesh.nodes[tr.model.mesh.node_index(id).unwrap()];
            cx += n.coords[0];
            cy += n.coords[1];
        }
        assert!(
            (cx / 3.0).hypot(cy / 3.0) > 0.2,
            "element {} centroid inside the hole",
            e.id
        );
    }

    // Full pipeline: deck round-trip, solve, Kirsch sanity.
    let deck = write_deck(&tr.model);
    let (m2, report) = parse_inp(&deck).unwrap();
    assert!(report.warnings.is_empty(), "{:?}", report.warnings);
    let sol = solve(&m2).unwrap();
    let peak = sol.nodal_von_mises.iter().cloned().fold(0.0f64, f64::max);
    assert!(
        (2.2e6..3.6e6).contains(&peak),
        "Kirsch sanity peak {peak:.4e} Pa"
    );
}

#[test]
fn e_cli_mesh2d_flag() {
    let bin = env!("CARGO_BIN_EXE_abaqus2solvix");
    let out = std::process::Command::new(bin)
        .arg("tests/fixtures/plate_hole.py")
        .arg("--mesh2d")
        .arg("-o")
        .arg(concat!(env!("CARGO_TARGET_TMPDIR"), "/plate_mesh2d.inp"))
        .arg("--json")
        .arg("--validate")
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .expect("run abaqus2solvix --mesh2d");
    assert!(out.status.success(), "status: {out:?}");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.starts_with(r#"{"ok":true,"output":"#), "{stdout}");
    assert!(stdout.contains("conformal mesh"), "{stdout}");
    assert!(
        stdout.contains(r#""validate":{"parsed":true,"solved":true"#),
        "{stdout}"
    );
}
