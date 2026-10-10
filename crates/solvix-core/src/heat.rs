//! # Steady-state heat conduction and sequential thermal stress
//!
//! [`solve_thermal`]: scalar temperature dofs (one per node), Fourier
//! conduction `C_e = k ∫ ∇N^T ∇N dV`. The scalar gradients reuse the
//! mechanical B-matrix machinery — first rows of `quad_b`/`quad8_b`/
//! `hex_b`, constant closed-form gradients for CPS3 and C3D4 — with the
//! same Gauss rules as the stiffness (2x2 for quads, 3x3 for
//! serendipity quads, 2x2x2 for hexes). Assembled on the shared CSR
//! machinery at node-count size and solved with the Jacobi-PCG.
//! Prescribed temperatures (`*TEMPERATURE`) are nonhomogeneous Dirichlet
//! rows; `*DFLUX` surface fluxes (positive INTO the element) and
//! volumetric generation (BN) feed the RHS.
//!
//! [`solve_thermomech`]: sequential coupling — thermal solve, then the
//! mechanical solve of the last step's state plus equivalent thermal
//! loads `f_th = ∫ B^T sigma_th dV` with
//! `sigma_th = E alpha dT / (1 - nu)` (plane stress) or
//! `E alpha dT / (1 - 2 nu)` (plane strain, solids) on the normal
//! components. Stress recovery subtracts `sigma_th` before von Mises,
//! so returned stresses include the thermal strain.
//!
//! v1 ceilings: uniform dT per element (element-average nodal
//! temperature) against T_ref = 0; one conductivity/expansion per
//! model; no beam thermal expansion; PCG only for the thermal solve
//! (conduction matrices are Laplacian-like and benign — wire the LDLT
//! fallback in here if a mesh ever stalls it).

use crate::assembly::{coords2, csr_from_triplets, Csr};
use crate::element::{
    hex_b, hex_shapes, invalid, quad8_b, quad8_shapes, quad_b, quad_shapes, tet_gradients, G3, W3,
};
use crate::solver;
use crate::{
    ElemKind, Element, HeatFlux, Method, Model, Solution, SolveError, TempBc, ThermalLoads,
};

/// Steady-state thermal solve output. `temperatures` is parallel to
/// `mesh.nodes` (units as prescribed by the BCs).
#[derive(Debug, Clone, Default)]
pub struct ThermalSolution {
    pub temperatures: Vec<f64>,
    pub iterations: usize,
    pub residual: f64,
}

/// Solve the steady-state conduction problem for the LAST step's
/// thermal state (temperature BCs and fluxes live on steps — a step-less
/// model has none and fails with `NoConstraints`, as does a step without
/// any `*TEMPERATURE`: a pure-Neumann problem is singular). One
/// conductivity covers the whole model. Relative tolerance 1e-10,
/// `max_iter = max(1000, 10 * nodes)`.
pub fn solve_thermal(model: &Model) -> Result<ThermalSolution, SolveError> {
    if model.mesh.nodes.is_empty() || model.mesh.elements.is_empty() {
        return Err(SolveError::EmptyMesh);
    }
    if model.mesh.elements.iter().any(|e| e.kind == ElemKind::B31) {
        return Err(SolveError::InvalidModel(
            "beam elements have no conduction formulation (v1: plane and solid only)".into(),
        ));
    }
    let Some(step) = model.steps.last() else {
        return Err(SolveError::NoConstraints);
    };
    if step.temp_bcs.is_empty() {
        return Err(SolveError::NoConstraints);
    }
    if model.thickness <= 0.0
        || model.thickness.is_nan()
        || model.conductivity <= 0.0
        || model.conductivity.is_nan()
    {
        return Err(SolveError::InvalidModel(
            "thickness and conductivity must be positive".into(),
        ));
    }

    let n = model.mesh.nodes.len();
    let mut trips: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 0.0)).collect();
    let mut f = vec![0.0; n];
    for e in &model.mesh.elements {
        conduction_element(model, e, step.heat_source, &mut trips, &mut f)?;
    }
    heat_flux_loads(model, &step.heat_fluxes, &mut f)?;
    let mut c = csr_from_triplets(n, trips);
    apply_temperatures(model, &step.temp_bcs, &mut c, &mut f)?;

    let (temperatures, iterations, residual) = solver::pcg(&c, &f, 1e-10, 1000.max(10 * n))?;
    Ok(ThermalSolution {
        temperatures,
        iterations,
        residual,
    })
}

/// Sequential thermal-stress solve: [`solve_thermal`] on the last step,
/// then the linear-static solve of the accumulated mechanical state
/// plus equivalent thermal loads from the temperature field. Returns
/// the mechanical [`Solution`]; von Mises includes thermal strain
/// (sigma = D (B u - eps_th)). Uniform dT per element (element-average
/// nodal temperature), T_ref = 0; no beam thermal expansion (v1).
pub fn solve_thermomech(model: &Model) -> Result<Solution, SolveError> {
    let thermal = solve_thermal(model)?;
    let alpha = model.thermal_expansion;
    let ndof = crate::ndof_per_node(model)?;
    let mut f = vec![0.0; model.mesh.nodes.len() * ndof];
    let mut sigma_th = vec![[0.0f64; 6]; model.mesh.elements.len()];
    for (ei, e) in model.mesh.elements.iter().enumerate() {
        if e.kind == ElemKind::B31 {
            continue; // v1: no beam thermal expansion
        }
        let mat = model
            .material_for_element(ei)
            .ok_or_else(|| invalid(e.id, "no material assigned"))?;
        let idx: Vec<usize> = e
            .connectivity
            .iter()
            .map(|&nid| model.mesh.node_index(nid).unwrap())
            .collect();
        let dt = idx.iter().map(|&i| thermal.temperatures[i]).sum::<f64>() / idx.len() as f64;
        // Equivalent thermal stress on the normal components: for plane
        // stress E alpha dT / (1 - nu) with sigma_zz = 0; for plane
        // strain and solids E alpha dT / (1 - 2 nu) including sigma_zz
        // (the plane-strain out-of-plane constraint carries stress).
        let constrained_zz = !matches!(e.kind, ElemKind::Cps3 | ElemKind::Cps4 | ElemKind::Cps8);
        let c = mat.youngs_modulus * alpha * dt
            / if constrained_zz {
                1.0 - 2.0 * mat.poisson_ratio
            } else {
                1.0 - mat.poisson_ratio
            };
        let s = [c, c, if constrained_zz { c } else { 0.0 }, 0.0, 0.0, 0.0];
        sigma_th[ei] = s;
        thermal_load_element(model, &idx, e, s, ndof, &mut f)?;
    }
    crate::solve_with_thermal(model, Method::Auto, Some(ThermalLoads { f, sigma_th }))
}

/// CPS3 constant-gradient ingredients (b_i, c_i and the signed 2A),
/// same definitions as the CPS3 stiffness.
fn tri_constants(p: &[[f64; 2]; 3]) -> ([f64; 3], [f64; 3], f64) {
    let b = [p[1][1] - p[2][1], p[2][1] - p[0][1], p[0][1] - p[1][1]];
    let c = [p[2][0] - p[1][0], p[0][0] - p[2][0], p[1][0] - p[0][0]];
    let two_a = b[0] * c[1] - b[1] * c[0];
    (b, c, two_a)
}

/// Element conduction matrix and consistent heat-source load, scattered
/// straight into the global triplets/RHS.
fn conduction_element(
    model: &Model,
    e: &Element,
    q: Option<f64>,
    trips: &mut Vec<(usize, usize, f64)>,
    f: &mut [f64],
) -> Result<(), SolveError> {
    let idx: Vec<usize> = e
        .connectivity
        .iter()
        .map(|&nid| model.mesh.node_index(nid).unwrap())
        .collect();
    let (k, t) = (model.conductivity, model.thickness);
    let jac = || invalid(e.id, "non-positive Jacobian determinant at a Gauss point");
    match e.kind {
        ElemKind::Cps3 => {
            let p = [
                coords2(model, &idx, 0),
                coords2(model, &idx, 1),
                coords2(model, &idx, 2),
            ];
            let (b, c, two_a) = tri_constants(&p);
            if two_a == 0.0 {
                return Err(invalid(e.id, "degenerate (zero-area) triangle"));
            }
            // ∫∇N·∇N dA = (b_i b_j + c_i c_j) / (4A), times k*t.
            let w = k * t / (2.0 * two_a.abs());
            for i in 0..3 {
                if let Some(q) = q {
                    f[idx[i]] += q * t * two_a.abs() / 6.0; // q t A / 3
                }
                for j in 0..3 {
                    trips.push((idx[i], idx[j], w * (b[i] * b[j] + c[i] * c[j])));
                }
            }
        }
        ElemKind::Cps4 | ElemKind::Cpe4 => {
            let coords: [[f64; 2]; 4] = core::array::from_fn(|m| coords2(model, &idx, m));
            let g = 1.0f64 / (3.0f64).sqrt();
            let mut ke = [[0.0f64; 4]; 4];
            let mut fe = [0.0f64; 4];
            for &xi in &[-g, g] {
                for &eta in &[-g, g] {
                    let (bm, detj) = quad_b(&coords, xi, eta).ok_or_else(jac)?;
                    let n = quad_shapes(xi, eta);
                    for i in 0..4 {
                        fe[i] += q.unwrap_or(0.0) * t * detj * n[i];
                        for j in 0..4 {
                            ke[i][j] += k
                                * t
                                * detj
                                * (bm[0][2 * i] * bm[0][2 * j]
                                    + bm[1][2 * i + 1] * bm[1][2 * j + 1]);
                        }
                    }
                }
            }
            for i in 0..4 {
                f[idx[i]] += fe[i];
                for j in 0..4 {
                    trips.push((idx[i], idx[j], ke[i][j]));
                }
            }
        }
        ElemKind::Cps8 | ElemKind::Cpe8 => {
            let coords: [[f64; 2]; 8] = core::array::from_fn(|m| coords2(model, &idx, m));
            let mut ke = [[0.0f64; 8]; 8];
            let mut fe = [0.0f64; 8];
            for (a, &xi) in [-G3, 0.0, G3].iter().enumerate() {
                for (b, &eta) in [-G3, 0.0, G3].iter().enumerate() {
                    let (bm, detj) = quad8_b(&coords, xi, eta).ok_or_else(jac)?;
                    let n = quad8_shapes(xi, eta);
                    let w = k * t * W3[a] * W3[b];
                    for i in 0..8 {
                        fe[i] += q.unwrap_or(0.0) * t * detj * W3[a] * W3[b] * n[i];
                        for j in 0..8 {
                            ke[i][j] += w
                                * detj
                                * (bm[0][2 * i] * bm[0][2 * j]
                                    + bm[1][2 * i + 1] * bm[1][2 * j + 1]);
                        }
                    }
                }
            }
            for i in 0..8 {
                f[idx[i]] += fe[i];
                for j in 0..8 {
                    trips.push((idx[i], idx[j], ke[i][j]));
                }
            }
        }
        ElemKind::C3d8 => {
            let coords: [[f64; 3]; 8] = core::array::from_fn(|m| model.mesh.nodes[idx[m]].coords);
            let g = 1.0f64 / (3.0f64).sqrt();
            let mut ke = [[0.0f64; 8]; 8];
            let mut fe = [0.0f64; 8];
            for &xi in &[-g, g] {
                for &eta in &[-g, g] {
                    for &zeta in &[-g, g] {
                        let (bm, detj) = hex_b(&coords, xi, eta, zeta).ok_or_else(jac)?;
                        let n = hex_shapes(xi, eta, zeta);
                        for i in 0..8 {
                            fe[i] += q.unwrap_or(0.0) * detj * n[i];
                            for j in 0..8 {
                                let mut v = 0.0;
                                for r in 0..3 {
                                    v += bm[r][3 * i + r] * bm[r][3 * j + r];
                                }
                                ke[i][j] += k * detj * v;
                            }
                        }
                    }
                }
            }
            for i in 0..8 {
                f[idx[i]] += fe[i];
                for j in 0..8 {
                    trips.push((idx[i], idx[j], ke[i][j]));
                }
            }
        }
        ElemKind::C3d4 => {
            let coords: [[f64; 3]; 4] = core::array::from_fn(|m| model.mesh.nodes[idx[m]].coords);
            let (gs, six_v) = tet_gradients(&coords)
                .ok_or_else(|| invalid(e.id, "degenerate (zero-volume) tetrahedron"))?;
            let v = six_v.abs() / 6.0;
            for i in 0..4 {
                if let Some(q) = q {
                    f[idx[i]] += q * v / 4.0;
                }
                for j in 0..4 {
                    let dot = gs[i][0] * gs[j][0] + gs[i][1] * gs[j][1] + gs[i][2] * gs[j][2];
                    trips.push((idx[i], idx[j], k * v * dot));
                }
            }
        }
        ElemKind::B31 => {
            return Err(invalid(
                e.id,
                "beam elements have no conduction formulation (v1: plane and solid only)",
            ));
        }
    }
    Ok(())
}

/// Consistent nodal heat loads for `*DFLUX` S-faces: q * t * len split
/// 1/2 - 1/2 over a linear edge (corner/midside 1/6 - 2/3 - 1/6 on the
/// quadratic quads). Positive magnitude flows INTO the element.
fn heat_flux_loads(model: &Model, fluxes: &[HeatFlux], f: &mut [f64]) -> Result<(), SolveError> {
    for hf in fluxes {
        let e = model
            .mesh
            .elements
            .iter()
            .find(|e| e.id == hf.element_id)
            .ok_or_else(|| {
                SolveError::InvalidModel(format!("heat flux on unknown element {}", hf.element_id))
            })?;
        let quadratic = matches!(e.kind, ElemKind::Cps8 | ElemKind::Cpe8);
        let max_face = match e.kind {
            ElemKind::Cps3 => 3,
            ElemKind::Cps4 | ElemKind::Cpe4 | ElemKind::Cps8 | ElemKind::Cpe8 => 4,
            _ => {
                return Err(SolveError::InvalidModel(format!(
                    "heat flux on element {} ({:?}) unsupported (plane edges only)",
                    e.id, e.kind
                )))
            }
        };
        if hf.face == 0 || hf.face > max_face {
            return Err(SolveError::InvalidModel(format!(
                "element {} ({:?}) has no face {}",
                e.id, e.kind, hf.face
            )));
        }
        let corners = if quadratic { 4 } else { e.connectivity.len() };
        let na = e.connectivity[hf.face - 1];
        let nb = e.connectivity[hf.face % corners];
        let pos = |nid: u32| {
            let i = model.mesh.node_index(nid).unwrap();
            model.mesh.nodes[i].coords
        };
        let (pa, pb) = (pos(na), pos(nb));
        let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2)).sqrt();
        if len == 0.0 {
            return Err(SolveError::InvalidModel(format!(
                "zero-length flux edge on element {}",
                e.id
            )));
        }
        let total = hf.magnitude * model.thickness * len;
        let ia = model.mesh.node_index(na).unwrap();
        let ib = model.mesh.node_index(nb).unwrap();
        if quadratic {
            let mid = e.connectivity[hf.face + 3]; // midside 4 + face - 1
            let im = model.mesh.node_index(mid).unwrap();
            f[im] += total * 2.0 / 3.0;
            f[ia] += total / 6.0;
            f[ib] += total / 6.0;
        } else {
            f[ia] += total / 2.0;
            f[ib] += total / 2.0;
        }
    }
    Ok(())
}

/// Nonhomogeneous Dirichlet temperatures: zero the row and column of
/// every prescribed node, 1 on the diagonal, the temperature in the RHS
/// (later `TempBc`s win on shared nodes). Same single-pass structure as
/// the mechanical `apply_dirichlet`.
fn apply_temperatures(
    model: &Model,
    temp_bcs: &[TempBc],
    c: &mut Csr,
    f: &mut [f64],
) -> Result<(), SolveError> {
    let n = f.len();
    let mut prescribed = vec![f64::NAN; n];
    for tb in temp_bcs {
        for &nid in &tb.node_ids {
            let i = model.mesh.node_index(nid).ok_or_else(|| {
                SolveError::InvalidModel(format!("temperature bc on unknown node {nid}"))
            })?;
            prescribed[i] = tb.temperature;
        }
    }
    let fixed = |i: usize| !prescribed[i].is_nan();
    // Free rows: move the prescribed-value coupling K[a,r]*T_r to the
    // RHS, then zero the column entry (for homogeneous Dirichlet the
    // compensation is zero — that is why the mechanical version can
    // skip it).
    for (r, fr) in f.iter_mut().enumerate() {
        if fixed(r) {
            continue;
        }
        for ii in c.row_ptr[r]..c.row_ptr[r + 1] {
            let col = c.col_idx[ii];
            if fixed(col) {
                *fr -= c.vals[ii] * prescribed[col];
                c.vals[ii] = 0.0;
            }
        }
    }
    // Prescribed rows: identity row, RHS = temperature.
    for (r, &t) in prescribed.iter().enumerate() {
        if t.is_nan() {
            continue;
        }
        for ii in c.row_ptr[r]..c.row_ptr[r + 1] {
            let col = c.col_idx[ii];
            c.vals[ii] = if r == col { 1.0 } else { 0.0 };
        }
        f[r] = t;
    }
    Ok(())
}

/// Scatter `f_th = ∫ B^T sigma_th dV` (thickness-weighted in 2D) into
/// the global load vector, model dof layout.
fn thermal_load_element(
    model: &Model,
    idx: &[usize],
    e: &Element,
    s: [f64; 6],
    ndof: usize,
    f: &mut [f64],
) -> Result<(), SolveError> {
    let t = model.thickness;
    let jac = || invalid(e.id, "non-positive Jacobian determinant at a Gauss point");
    match e.kind {
        ElemKind::Cps3 => {
            let p = [
                coords2(model, idx, 0),
                coords2(model, idx, 1),
                coords2(model, idx, 2),
            ];
            let (b, c, two_a) = tri_constants(&p);
            if two_a == 0.0 {
                return Err(invalid(e.id, "degenerate (zero-area) triangle"));
            }
            for i in 0..3 {
                let col = idx[i] * ndof;
                f[col] += s[0] * t * b[i] / 2.0;
                f[col + 1] += s[1] * t * c[i] / 2.0;
            }
        }
        ElemKind::Cps4 | ElemKind::Cpe4 => {
            let coords: [[f64; 2]; 4] = core::array::from_fn(|m| coords2(model, idx, m));
            let g = 1.0f64 / (3.0f64).sqrt();
            for &xi in &[-g, g] {
                for &eta in &[-g, g] {
                    let (bm, detj) = quad_b(&coords, xi, eta).ok_or_else(jac)?;
                    for i in 0..4 {
                        let col = idx[i] * ndof;
                        f[col] += s[0] * t * detj * bm[0][2 * i];
                        f[col + 1] += s[1] * t * detj * bm[1][2 * i + 1];
                    }
                }
            }
        }
        ElemKind::Cps8 | ElemKind::Cpe8 => {
            let coords: [[f64; 2]; 8] = core::array::from_fn(|m| coords2(model, idx, m));
            for (a, &xi) in [-G3, 0.0, G3].iter().enumerate() {
                for (b, &eta) in [-G3, 0.0, G3].iter().enumerate() {
                    let (bm, detj) = quad8_b(&coords, xi, eta).ok_or_else(jac)?;
                    let w = t * detj * W3[a] * W3[b];
                    for i in 0..8 {
                        let col = idx[i] * ndof;
                        f[col] += s[0] * w * bm[0][2 * i];
                        f[col + 1] += s[1] * w * bm[1][2 * i + 1];
                    }
                }
            }
        }
        ElemKind::C3d8 => {
            let coords: [[f64; 3]; 8] = core::array::from_fn(|m| model.mesh.nodes[idx[m]].coords);
            let g = 1.0f64 / (3.0f64).sqrt();
            for &xi in &[-g, g] {
                for &eta in &[-g, g] {
                    for &zeta in &[-g, g] {
                        let (bm, detj) = hex_b(&coords, xi, eta, zeta).ok_or_else(jac)?;
                        for i in 0..8 {
                            let col = idx[i] * ndof;
                            for r in 0..3 {
                                f[col + r] += s[r] * detj * bm[r][3 * i + r];
                            }
                        }
                    }
                }
            }
        }
        ElemKind::C3d4 => {
            let coords: [[f64; 3]; 4] = core::array::from_fn(|m| model.mesh.nodes[idx[m]].coords);
            let (gs, _) = tet_gradients(&coords)
                .ok_or_else(|| invalid(e.id, "degenerate (zero-volume) tetrahedron"))?;
            // ∫B dV = B * V with B entries g_i/(6V): the integral of a
            // gradient row is g_i/6 regardless of the volume.
            for i in 0..4 {
                let col = idx[i] * ndof;
                for r in 0..3 {
                    f[col + r] += s[r] * gs[i][r] / 6.0;
                }
            }
        }
        ElemKind::B31 => {} // caller skips beams
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::tests::{hex_mesh, quad_mesh};
    use crate::tests::solid_model;
    use crate::{BoundaryCondition, ConcentratedLoad, ElementSet, Material, Step};

    const E: f64 = 210e9;
    const NU: f64 = 0.3;

    fn strip(nx: usize, ny: usize) -> Model {
        let n = (nx * ny) as u32;
        Model {
            mesh: quad_mesh(nx, ny, 1.0, 0.1),
            thickness: 1.0,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: E,
                poisson_ratio: NU,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=n).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            ..Model::default()
        }
    }

    #[test]
    fn bar_fixed_temps_linear_gradient_exact() {
        let (nx, ny) = (8, 1);
        let mut m = strip(nx, ny);
        m.conductivity = 1.0;
        let left: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        let right: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + nx + 1) as u32).collect();
        m.steps = vec![Step {
            name: "heat".into(),
            temp_bcs: vec![
                TempBc {
                    node_ids: left,
                    temperature: 100.0,
                },
                TempBc {
                    node_ids: right,
                    temperature: 0.0,
                },
            ],
            ..Step::default()
        }];
        let sol = solve_thermal(&m).unwrap();
        assert!(sol.residual <= 1e-10);
        assert!(sol.iterations > 0);
        for node in &m.mesh.nodes {
            let exact = 100.0 * (1.0 - node.coords[0]);
            let got = sol.temperatures[m.mesh.node_index(node.id).unwrap()];
            assert!((got - exact).abs() < 1e-8, "node: {got} vs {exact}");
        }
    }

    #[test]
    fn hex_bar_fixed_temps_linear_gradient_exact() {
        let (nx, ny, nz) = (6, 1, 1);
        let mut m = solid_model(hex_mesh(nx, ny, nz, 1.0, 0.1, 0.1));
        m.conductivity = 5.0;
        let nid =
            |i: usize, j: usize, k: usize| (k * (ny + 1) * (nx + 1) + j * (nx + 1) + i + 1) as u32;
        let cold: Vec<u32> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| nid(0, j, k)))
            .collect();
        let hot: Vec<u32> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| nid(nx, j, k)))
            .collect();
        m.steps = vec![Step {
            name: "heat".into(),
            temp_bcs: vec![
                TempBc {
                    node_ids: cold,
                    temperature: 20.0,
                },
                TempBc {
                    node_ids: hot,
                    temperature: 120.0,
                },
            ],
            ..Step::default()
        }];
        let sol = solve_thermal(&m).unwrap();
        for node in &m.mesh.nodes {
            let exact = 20.0 + 100.0 * node.coords[0];
            let got = sol.temperatures[m.mesh.node_index(node.id).unwrap()];
            assert!((got - exact).abs() < 1e-8, "node: {got} vs {exact}");
        }
    }

    #[test]
    fn surface_flux_matches_equivalent_dirichlet() {
        // Heat flux q into the right edge of a k=2 bar is equivalent to
        // fixing T_right = q*L/k: both give T(x) = q(L-x)/k.
        let (nx, ny) = (6, 1);
        let (k, q) = (2.0, 50.0);
        let mut fluxed = strip(nx, ny);
        fluxed.conductivity = k;
        let left: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        let right: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + nx + 1) as u32).collect();
        fluxed.steps = vec![Step {
            name: "flux".into(),
            temp_bcs: vec![TempBc {
                node_ids: left.clone(),
                temperature: 0.0,
            }],
            heat_fluxes: (0..ny)
                .map(|j| HeatFlux {
                    element_id: (j * nx + nx) as u32,
                    face: 2,
                    magnitude: q,
                })
                .collect(),
            ..Step::default()
        }];
        let mut fixed = fluxed.clone();
        fixed.steps = vec![Step {
            name: "fixed".into(),
            temp_bcs: vec![
                TempBc {
                    node_ids: left.clone(),
                    temperature: 0.0,
                },
                TempBc {
                    node_ids: right,
                    temperature: q / k, // L = 1
                },
            ],
            ..Step::default()
        }];
        let a = solve_thermal(&fluxed).unwrap();
        let b = solve_thermal(&fixed).unwrap();
        for (i, (ta, tb)) in a.temperatures.iter().zip(&b.temperatures).enumerate() {
            let x = fluxed.mesh.nodes[i].coords[0];
            let exact = q * x / k; // T(0)=0 fixed, q flows in at x=L=1
            assert!(
                (ta - exact).abs() < 1e-8 && (tb - exact).abs() < 1e-8,
                "node {i}: flux={ta} fixed={tb} exact={exact}"
            );
        }
    }

    #[test]
    fn uniform_dt_fully_constrained_bar_sigma_is_e_alpha_dt() {
        // Plane-stress bar with ux fixed on both ends and a uniform dT:
        // sigma_x = -E alpha dT (sigma_y = 0), vM = E alpha dT exactly.
        let (nx, ny) = (8, 2);
        let mut m = strip(nx, ny);
        let dt = 100.0;
        let left: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        let right: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + nx + 1) as u32).collect();
        m.steps = vec![Step {
            name: "thermal".into(),
            boundary_conditions: vec![
                BoundaryCondition {
                    node_ids: left,
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: right,
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: vec![1],
                    dofs: vec![2],
                },
            ],
            temp_bcs: vec![TempBc {
                node_ids: (1..=((nx + 1) * (ny + 1)) as u32).collect(),
                temperature: dt,
            }],
            ..Step::default()
        }];
        let sol = solve_thermomech(&m).unwrap();
        let expect = E * m.thermal_expansion * dt;
        for (i, &vm) in sol.von_mises.iter().enumerate() {
            assert!(
                (vm - expect).abs() / expect < 1e-6,
                "element {i}: vM={vm} vs {expect}"
            );
        }
        for (i, &nvm) in sol.nodal_von_mises.iter().enumerate() {
            assert!(
                (nvm - expect).abs() / expect < 1e-6,
                "node {i}: vM={nvm} vs {expect}"
            );
        }
        // Fully constrained in x: no x displacement anywhere.
        for d in &sol.displacements {
            assert!(d[0].abs() < 1e-12);
        }
    }

    #[test]
    fn thermal_and_mechanical_superpose_linearly() {
        let (nx, ny) = (8, 2);
        let bcs = vec![
            BoundaryCondition {
                node_ids: (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect(),
                dofs: vec![1],
            },
            BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![2],
            },
        ];
        let loads = vec![ConcentratedLoad {
            node_id: (ny * (nx + 1) + nx + 1) as u32,
            dof: 1,
            value: 5.0e5,
        }];
        let hot = |bcs: Vec<BoundaryCondition>, loads: Vec<ConcentratedLoad>| {
            let mut m = strip(nx, ny);
            m.steps = vec![Step {
                name: "s".into(),
                boundary_conditions: bcs,
                cloads: loads,
                temp_bcs: vec![TempBc {
                    node_ids: (1..=((nx + 1) * (ny + 1)) as u32).collect(),
                    temperature: 50.0,
                }],
                ..Step::default()
            }];
            m
        };
        let full = solve_thermomech(&hot(bcs.clone(), loads.clone())).unwrap();
        let th = solve_thermomech(&hot(bcs.clone(), vec![])).unwrap();
        let mut mech = hot(bcs, loads);
        mech.steps[0].temp_bcs.clear();
        let mech = crate::solve(&mech).unwrap();
        let scale = full
            .displacements
            .iter()
            .flat_map(|d| d.iter())
            .fold(0.0f64, |a, &v| a.max(v.abs()));
        for (uf, (um, ut)) in full
            .displacements
            .iter()
            .zip(mech.displacements.iter().zip(&th.displacements))
        {
            for d in 0..3 {
                let sum = um[d] + ut[d];
                assert!(
                    (uf[d] - sum).abs() < 1e-9 * scale + 1e-15,
                    "dof {d}: full={} mech+th={sum}",
                    uf[d]
                );
            }
        }
    }

    #[test]
    fn thermal_errors() {
        // No step at all.
        let m = strip(2, 1);
        assert!(matches!(solve_thermal(&m), Err(SolveError::NoConstraints)));
        // Step without temperature BCs (pure Neumann = singular).
        let mut m = strip(2, 1);
        m.steps = vec![Step::default()];
        assert!(matches!(solve_thermal(&m), Err(SolveError::NoConstraints)));
        // Beams have no conduction.
        let mut m = strip(2, 1);
        m.mesh = crate::assembly::tests::beam_mesh(2, 1.0);
        m.steps = vec![Step {
            temp_bcs: vec![TempBc {
                node_ids: vec![1],
                temperature: 1.0,
            }],
            ..Step::default()
        }];
        match solve_thermal(&m) {
            Err(SolveError::InvalidModel(msg)) => assert!(msg.contains("conduction")),
            other => panic!("expected InvalidModel(conduction), got {other:?}"),
        }
    }
}
