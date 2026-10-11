//! Small-strain J2 plasticity with linear isotropic hardening.
//!
//! The constitutive update is the classical closed-form radial return at
//! each integration point.  Stress and equivalent plastic strain are
//! history variables: every Newton trial starts from the last committed
//! increment, and only a line-search-accepted trial is committed.  The
//! returned stiffness is the closed-form consistent elastoplastic tangent
//! (not the elastic secant), so full Newton retains its local convergence
//! rate after yield.
//!
//! Element support is deliberate.  C3D8/C3D4 use the direct six-component
//! return.  Plane elements carry three B rows `(eps_xx, eps_yy, gamma_xy)`
//! mapped into the six-component order's in-plane slots `[0, 1, 3]`:
//! plane-strain CPE4/CPE8 (eps_zz = 0) reuse the six-component return,
//! plane-stress CPS3/CPS4/CPS8 use the dedicated three-component
//! fixed-point return below that enforces `sigma_zz = 0`.  B31 remains
//! rejected.  Plasticity also takes precedence over NLGEOM and uses the
//! small-strain branch documented here.

use crate::assembly::{self, Csr};
use crate::element::{
    cps3_b_matrix, d3_matrix, d_matrix, hex_b, quad8_b, quad_b, tet_gradients, von_mises, G3, W3,
};
use crate::ldl::Ldlt;
use crate::solver::{apply_dirichlet, constrained_dofs};
use crate::{ElemKind, Model, NlFrame, NlOptions, PlasticitySpec, Reaction, Solution, SolveError};

/// Strain slots of the six-component order `(xx, yy, zz, xy, yz, xz)`
/// that plane elements' three B rows occupy.
const IN_PLANE_SLOTS: [usize; 3] = [0, 1, 3];
/// All six slots, in order (solid elements).
const SOLID_SLOTS: [usize; 6] = [0, 1, 2, 3, 4, 5];
/// Plane-stress fixed-point convergence tolerance, relative to the
/// current yield stress.
const PS_TOL: f64 = 1.0e-10;
/// Plane-stress fixed-point iteration cap (3-5 passes is typical).
const PS_MAX_ITER: usize = 50;

/// One integration point's flattened B matrix (row-major) and integration
/// weight (including thickness for plane elements).
struct IntegrationPoint {
    b: Vec<Vec<f64>>,
    weight: f64,
}

/// Per-Gauss-point history carried between Newton increments.
#[derive(Clone, Copy)]
struct HistoryPoint {
    sigma: [f64; 6],
    equivalent_plastic_strain: f64,
}

struct PlasticElem {
    node_idx: Vec<usize>,
    points: Vec<IntegrationPoint>,
    /// Full 3D elastic matrix; plane-strain elements feed only their
    /// in-plane strain components into the six-component return.
    /// Plane-stress elements instead carry the plane-stress 3x3
    /// embedded at the in-plane slots (their elastic-fallback tangent).
    d: [[f64; 6]; 6],
    youngs_modulus: f64,
    poisson_ratio: f64,
    /// CPS elements route to the plane-stress fixed-point return.
    plane_stress: bool,
    history: Vec<HistoryPoint>,
}

struct Evaluation {
    kt: Csr,
    f_int: Vec<f64>,
    history: Vec<Vec<HistoryPoint>>,
}

/// An internal inhomogeneous Dirichlet value.  The public model still has
/// homogeneous BCs only; this hook lets the core's displacement-controlled
/// tests exercise the material path without changing that public contract.
#[derive(Clone, Copy)]
pub(crate) struct PrescribedDof {
    node_id: u32,
    dof: usize,
    value: f64,
}

#[derive(Clone, Copy)]
struct ReturnResult {
    sigma: [f64; 6],
    equivalent_plastic_strain: f64,
    tangent: [[f64; 6]; 6],
}

fn mat_vec6(a: &[[f64; 6]; 6], x: &[f64; 6]) -> [f64; 6] {
    let mut out = [0.0; 6];
    for (row, value) in a.iter().zip(out.iter_mut()) {
        *value = row.iter().zip(x).map(|(a, b)| a * b).sum();
    }
    out
}

fn deviator(sigma: &[f64; 6]) -> [f64; 6] {
    let mean = (sigma[0] + sigma[1] + sigma[2]) / 3.0;
    let mut out = *sigma;
    out[0] -= mean;
    out[1] -= mean;
    out[2] -= mean;
    out
}

/// Von Mises equivalent stress.  Engineering shear components are weighted
/// by the tensor contraction convention used throughout the crate.
fn equivalent_stress(dev: &[f64; 6]) -> f64 {
    (1.5 * (dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2])
        + 3.0 * (dev[3] * dev[3] + dev[4] * dev[4] + dev[5] * dev[5]))
        .sqrt()
}

/// Closed-form J2 radial return and consistent tangent for one history
/// point.  Component order is `(xx, yy, zz, xy, yz, xz)` with engineering
/// shear strain in the incremental strain vector.
fn radial_return(
    spec: PlasticitySpec,
    elastic: &[[f64; 6]; 6],
    youngs_modulus: f64,
    poisson_ratio: f64,
    old: HistoryPoint,
    strain_increment: &[f64; 6],
) -> ReturnResult {
    let trial = {
        let elastic_stress = mat_vec6(elastic, strain_increment);
        let mut trial = old.sigma;
        for (a, b) in trial.iter_mut().zip(elastic_stress) {
            *a += b;
        }
        trial
    };
    let s_trial = deviator(&trial);
    let q_trial = equivalent_stress(&s_trial);
    let current_yield = spec.yield_stress + spec.hardening_modulus * old.equivalent_plastic_strain;

    if q_trial <= current_yield {
        return ReturnResult {
            sigma: trial,
            equivalent_plastic_strain: old.equivalent_plastic_strain,
            tangent: *elastic,
        };
    }

    let shear = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    let denominator = 3.0 * shear + spec.hardening_modulus;
    let delta_lambda = (q_trial - current_yield) / denominator;
    let scale = 3.0 * shear * delta_lambda / q_trial;
    let mut sigma = trial;
    for (sigma_i, s_i) in sigma.iter_mut().zip(s_trial) {
        *sigma_i -= scale * s_i;
    }

    // Algorithmic tangent.  `r` is normalized with the tensor inner
    // product, whose engineering-shear metric is (1, 1, 1, 2, 2, 2).  The
    // radial term includes that metric because q's derivative has twice
    // the engineering-shear weight; the perpendicular term is the derivative
    // of the return direction at finite plastic strain.
    let r_norm = q_trial / (1.5f64).sqrt();
    let metric = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];
    let mut radial = [[0.0; 6]; 6];
    for (row, s_i) in radial.iter_mut().zip(s_trial) {
        for (cell, (s_j, m_j)) in row.iter_mut().zip(s_trial.iter().zip(metric)) {
            *cell = s_i * m_j * s_j / (r_norm * r_norm);
        }
    }
    let mut projection = [[0.0; 6]; 6];
    for (i, row) in projection.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for row in projection.iter_mut().take(3) {
        for cell in row.iter_mut().take(3) {
            *cell -= 1.0 / 3.0;
        }
    }
    for (row, radial_row) in projection.iter_mut().zip(radial) {
        for (cell, radial_cell) in row.iter_mut().zip(radial_row) {
            *cell -= radial_cell;
        }
    }
    let radial_coefficient = 3.0 * shear / denominator;
    let perpendicular_coefficient = (6.0f64).sqrt() * shear * delta_lambda / r_norm;
    let mut algorithmic = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            algorithmic[i][j] = if i == j { 1.0 } else { 0.0 }
                - radial_coefficient * radial[i][j]
                - perpendicular_coefficient * projection[i][j];
        }
    }
    let mut tangent = [[0.0; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut value = 0.0;
            for k in 0..6 {
                value += algorithmic[i][k] * elastic[k][j];
            }
            tangent[i][j] = value;
        }
    }

    ReturnResult {
        sigma,
        equivalent_plastic_strain: old.equivalent_plastic_strain + delta_lambda,
        tangent,
    }
}

/// One pass of the classical plane-stress fixed-point return: the
/// three-component J2 update enforcing `sigma_zz = 0` by iterating the
/// out-of-plane strain increment (de Souza Neto, Peric & Owen style).
///
/// Each pass builds the six-component trial with the current `d_eps_zz`
/// (seeded with the elastic elimination
/// `-lambda/(lambda+2G) * (d_eps_xx + d_eps_yy)`), performs the standard
/// radial return, and corrects `d_eps_zz` by the `sigma_zz` the return
/// leaves behind, scaled by `1/(lambda + 2G)`.  The plastic flow's zz
/// component is what moves the constraint around, so the iteration
/// converges when the returned `sigma_zz` does; 3-5 passes is typical
/// for hardening plastics (`PS_TOL` relative to the current yield
/// stress, `PS_MAX_ITER` cap).  In-plane strain components are read
/// from slots `[0, 1, 3]` of `strain_increment`; the returned stress
/// carries `sigma_yz = sigma_xz = 0`.  Returns the stress, the new
/// equivalent plastic strain, and whether the step stayed elastic.
fn plane_stress_core(
    spec: PlasticitySpec,
    youngs_modulus: f64,
    poisson_ratio: f64,
    old: HistoryPoint,
    strain_increment: [f64; 6],
) -> ([f64; 6], f64, bool) {
    let shear = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    let lambda =
        youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    let zz_stiffness = lambda + 2.0 * shear;
    let current_yield = spec.yield_stress + spec.hardening_modulus * old.equivalent_plastic_strain;
    // Elastic d_eps_zz eliminating sigma_zz; the old.sigma[2] term keeps
    // any converged fixed-point residual from accumulating across steps.
    let elastic_dzz =
        -(lambda * (strain_increment[0] + strain_increment[1]) + old.sigma[2]) / zz_stiffness;
    let trial_at = |d_eps_zz: f64| {
        let mut trial = old.sigma;
        trial[0] += zz_stiffness * strain_increment[0] + lambda * (strain_increment[1] + d_eps_zz);
        trial[1] += zz_stiffness * strain_increment[1] + lambda * (strain_increment[0] + d_eps_zz);
        trial[2] += lambda * (strain_increment[0] + strain_increment[1]) + zz_stiffness * d_eps_zz;
        trial[3] += shear * strain_increment[3];
        trial
    };

    let mut d_eps_zz = elastic_dzz;
    let mut sigma = old.sigma;
    let mut equivalent_plastic_strain = old.equivalent_plastic_strain;
    let mut elastic = false;
    for _ in 0..PS_MAX_ITER {
        let trial = trial_at(d_eps_zz);
        let s_trial = deviator(&trial);
        let q_trial = equivalent_stress(&s_trial);
        if q_trial <= current_yield {
            sigma = trial_at(elastic_dzz);
            elastic = true;
            break;
        }
        let delta_lambda = (q_trial - current_yield) / (3.0 * shear + spec.hardening_modulus);
        let scale = 3.0 * shear * delta_lambda / q_trial;
        let mut returned = trial;
        for (sigma_i, s_i) in returned.iter_mut().zip(s_trial) {
            *sigma_i -= scale * s_i;
        }
        sigma = returned;
        equivalent_plastic_strain = old.equivalent_plastic_strain + delta_lambda;
        if returned[2].abs() <= PS_TOL * current_yield {
            break;
        }
        d_eps_zz -= returned[2] / zz_stiffness;
    }
    // ponytail: a non-converged fixed point (perfect plasticity at
    // extreme nu is the known slow case) accepts the last iterate; the
    // line-searched outer Newton still enforces equilibrium.
    (sigma, equivalent_plastic_strain, elastic)
}

/// Embed a 3x3 in-plane `(xx, yy, xy)` matrix into the six-component
/// order at slots `[0, 1, 3]`.
fn embed_in_plane(d: &[[f64; 3]; 3]) -> [[f64; 6]; 6] {
    let mut out = [[0.0; 6]; 6];
    for (i, &si) in IN_PLANE_SLOTS.iter().enumerate() {
        for (j, &sj) in IN_PLANE_SLOTS.iter().enumerate() {
            out[si][sj] = d[i][j];
        }
    }
    out
}

/// Plane-stress J2 return with the numerical consistent tangent.
///
/// The tangent is formed by three extra fixed-point returns per Gauss
/// point: forward differences (h = 1e-9) of the returned in-plane stress
/// with respect to the three in-plane strain components, symmetrized.
/// Documented v1 choice — it preserves Newton's quadratic convergence in
/// practice at a fraction of the derivation risk of the closed-form
/// plane-stress algorithmic tangent (with its delta-lambda correction
/// terms), which is the upgrade path.
fn plane_stress_return(
    spec: PlasticitySpec,
    youngs_modulus: f64,
    poisson_ratio: f64,
    old: HistoryPoint,
    strain_increment: &[f64; 6],
) -> ReturnResult {
    let (sigma, equivalent_plastic_strain, elastic) =
        plane_stress_core(spec, youngs_modulus, poisson_ratio, old, *strain_increment);
    let mut c = [[0.0f64; 3]; 3];
    if !elastic {
        let h = 1.0e-9;
        for (j, &slot) in IN_PLANE_SLOTS.iter().enumerate() {
            let mut perturbed = *strain_increment;
            perturbed[slot] += h;
            let (sigma_p, _, _) =
                plane_stress_core(spec, youngs_modulus, poisson_ratio, old, perturbed);
            for (i, &row_slot) in IN_PLANE_SLOTS.iter().enumerate() {
                c[i][j] = (sigma_p[row_slot] - sigma[row_slot]) / h;
            }
        }
        for &(i, j) in &[(0, 1), (0, 2), (1, 2)] {
            let average = 0.5 * (c[i][j] + c[j][i]);
            c[i][j] = average;
            c[j][i] = average;
        }
    } else {
        c = d_matrix(youngs_modulus, poisson_ratio);
    }
    ReturnResult {
        sigma,
        equivalent_plastic_strain,
        tangent: embed_in_plane(&c),
    }
}

fn element_points(
    model: &Model,
    element: &crate::Element,
    node_idx: &[usize],
) -> Result<Vec<IntegrationPoint>, SolveError> {
    let invalid =
        |message: &str| SolveError::InvalidModel(format!("element {}: {message}", element.id));
    match element.kind {
        ElemKind::Cpe4 | ElemKind::Cps4 => {
            let coords: [[f64; 2]; 4] = core::array::from_fn(|i| {
                let coordinates = model.mesh.nodes[node_idx[i]].coords;
                [coordinates[0], coordinates[1]]
            });
            let g = 1.0f64 / (3.0f64).sqrt();
            let mut out = Vec::with_capacity(4);
            for xi in [-g, g] {
                for eta in [-g, g] {
                    let (b, detj) = quad_b(&coords, xi, eta)
                        .ok_or_else(|| invalid("non-positive Jacobian determinant"))?;
                    out.push(IntegrationPoint {
                        b: b.iter().map(|r| r.to_vec()).collect(),
                        weight: model.thickness * detj,
                    });
                }
            }
            Ok(out)
        }
        ElemKind::Cpe8 | ElemKind::Cps8 => {
            let coords: [[f64; 2]; 8] = core::array::from_fn(|i| {
                let coordinates = model.mesh.nodes[node_idx[i]].coords;
                [coordinates[0], coordinates[1]]
            });
            let mut out = Vec::with_capacity(9);
            for (i, xi) in [-G3, 0.0, G3].iter().enumerate() {
                for (j, eta) in [-G3, 0.0, G3].iter().enumerate() {
                    let (b, detj) = quad8_b(&coords, *xi, *eta)
                        .ok_or_else(|| invalid("non-positive Jacobian determinant"))?;
                    out.push(IntegrationPoint {
                        b: b.iter().map(|r| r.to_vec()).collect(),
                        weight: model.thickness * detj * W3[i] * W3[j],
                    });
                }
            }
            Ok(out)
        }
        ElemKind::C3d8 => {
            let coords: [[f64; 3]; 8] =
                core::array::from_fn(|i| model.mesh.nodes[node_idx[i]].coords);
            let g = 1.0f64 / (3.0f64).sqrt();
            let mut out = Vec::with_capacity(8);
            for xi in [-g, g] {
                for eta in [-g, g] {
                    for zeta in [-g, g] {
                        let (b, detj) = hex_b(&coords, xi, eta, zeta)
                            .ok_or_else(|| invalid("non-positive Jacobian determinant"))?;
                        out.push(IntegrationPoint {
                            b: b.iter().map(|r| r.to_vec()).collect(),
                            weight: detj,
                        });
                    }
                }
            }
            Ok(out)
        }
        ElemKind::C3d4 => {
            let coords: [[f64; 3]; 4] =
                core::array::from_fn(|i| model.mesh.nodes[node_idx[i]].coords);
            let (gradients, determinant) = tet_gradients(&coords)
                .ok_or_else(|| invalid("degenerate (zero-volume) tetrahedron"))?;
            let inverse = 1.0 / determinant;
            let mut b = vec![[0.0; 12]; 6];
            for (i, gradient) in gradients.iter().enumerate() {
                let column = 3 * i;
                b[0][column] = gradient[0] * inverse;
                b[1][column + 1] = gradient[1] * inverse;
                b[2][column + 2] = gradient[2] * inverse;
                b[3][column] = gradient[1] * inverse;
                b[3][column + 1] = gradient[0] * inverse;
                b[4][column + 1] = gradient[2] * inverse;
                b[4][column + 2] = gradient[1] * inverse;
                b[5][column] = gradient[2] * inverse;
                b[5][column + 2] = gradient[0] * inverse;
            }
            Ok(vec![IntegrationPoint {
                b: b.iter().map(|r| r.to_vec()).collect(),
                weight: determinant.abs() / 6.0,
            }])
        }
        ElemKind::Cps3 => {
            let coords: [[f64; 2]; 3] = core::array::from_fn(|i| {
                let coordinates = model.mesh.nodes[node_idx[i]].coords;
                [coordinates[0], coordinates[1]]
            });
            let (b, two_a) =
                cps3_b_matrix(&coords).ok_or_else(|| invalid("degenerate (zero-area) triangle"))?;
            Ok(vec![IntegrationPoint {
                b: b.iter().map(|r| r.to_vec()).collect(),
                weight: model.thickness * two_a.abs() / 2.0,
            }])
        }
        ElemKind::B31 => Err(SolveError::InvalidModel(format!(
            "element {}: J2 plasticity is not supported for beam kind {:?}",
            element.id, element.kind
        ))),
    }
}

fn build_elements(model: &Model, ndof: usize) -> Result<Vec<PlasticElem>, SolveError> {
    let loads = model.load_state();
    let (_, _, infos) = assembly::assemble(model, ndof, &loads)?;
    let mut out = Vec::with_capacity(infos.len());
    for (index, info) in infos.iter().enumerate() {
        let element = &model.mesh.elements[index];
        let material = model.material_for_element(index).ok_or_else(|| {
            SolveError::InvalidModel(format!("no material for element {}", element.id))
        })?;
        let points = element_points(model, element, &info.node_idx)?;
        let plane_stress = matches!(
            element.kind,
            ElemKind::Cps3 | ElemKind::Cps4 | ElemKind::Cps8
        );
        let history = points
            .iter()
            .map(|_| HistoryPoint {
                sigma: [0.0; 6],
                equivalent_plastic_strain: 0.0,
            })
            .collect();
        out.push(PlasticElem {
            node_idx: info.node_idx.clone(),
            points,
            d: if plane_stress {
                embed_in_plane(&d_matrix(material.youngs_modulus, material.poisson_ratio))
            } else {
                d3_matrix(material.youngs_modulus, material.poisson_ratio)
            },
            youngs_modulus: material.youngs_modulus,
            poisson_ratio: material.poisson_ratio,
            plane_stress,
            history,
        });
    }
    Ok(out)
}

fn evaluate(
    model: &Model,
    ndof: usize,
    elements: &[PlasticElem],
    step: &[f64],
) -> Result<Evaluation, SolveError> {
    evaluate_with_tangent(model, ndof, elements, step, false)
}

fn evaluate_with_tangent(
    model: &Model,
    ndof: usize,
    elements: &[PlasticElem],
    step: &[f64],
    elastic_tangent: bool,
) -> Result<Evaluation, SolveError> {
    let n = model.mesh.nodes.len() * ndof;
    let mut triplets = Vec::new();
    for row in 0..n {
        triplets.push((row, row, 0.0));
    }
    let mut f_int = vec![0.0; n];
    let mut history = Vec::with_capacity(elements.len());

    for element in elements {
        let spec = model.plasticity.unwrap();
        let mut new_history = Vec::with_capacity(element.points.len());
        let mut element_force = vec![0.0; element.node_idx.len() * ndof];

        for (point, old) in element.points.iter().zip(&element.history) {
            // Plane elements carry three B rows (eps_xx, eps_yy,
            // gamma_xy) in the in-plane slots; solids carry six in order.
            let slots: &[usize] = if point.b.len() == 3 {
                &IN_PLANE_SLOTS
            } else {
                &SOLID_SLOTS
            };
            let mut strain = [0.0; 6];
            for (r, b_row) in point.b.iter().enumerate() {
                strain[slots[r]] = b_row
                    .iter()
                    .enumerate()
                    .map(|(column, b)| {
                        let node = element.node_idx[column / ndof];
                        b * step[node * ndof + column % ndof]
                    })
                    .sum();
            }
            let returned = if element.plane_stress {
                plane_stress_return(
                    spec,
                    element.youngs_modulus,
                    element.poisson_ratio,
                    *old,
                    &strain,
                )
            } else {
                radial_return(
                    spec,
                    &element.d,
                    element.youngs_modulus,
                    element.poisson_ratio,
                    *old,
                    &strain,
                )
            };
            new_history.push(HistoryPoint {
                sigma: returned.sigma,
                equivalent_plastic_strain: returned.equivalent_plastic_strain,
            });
            let tangent = if elastic_tangent {
                element.d
            } else {
                returned.tangent
            };

            for (r, b_row) in point.b.iter().enumerate() {
                let sigma_r = returned.sigma[slots[r]];
                for (a, force) in element_force.iter_mut().enumerate() {
                    *force += point.weight * b_row[a] * sigma_r;
                }
            }
            for a in 0..element_force.len() {
                for c in a..element_force.len() {
                    let mut value = 0.0;
                    for (r, b_row) in point.b.iter().enumerate() {
                        let mut db = 0.0;
                        for (k, b_k) in point.b.iter().enumerate() {
                            db += tangent[slots[r]][slots[k]] * b_k[c];
                        }
                        value += b_row[a] * db;
                    }
                    value *= point.weight;
                    let global_a = element.node_idx[a / ndof] * ndof + a % ndof;
                    let global_c = element.node_idx[c / ndof] * ndof + c % ndof;
                    triplets.push((global_a, global_c, value));
                    if a != c {
                        triplets.push((global_c, global_a, value));
                    }
                }
            }
        }
        for (local, force) in element_force.iter().enumerate() {
            let global = element.node_idx[local / ndof] * ndof + local % ndof;
            f_int[global] += force;
        }
        history.push(new_history);
    }

    Ok(Evaluation {
        kt: assembly::csr_from_triplets(n, triplets),
        f_int,
        history,
    })
}

#[allow(clippy::too_many_arguments)]
fn stress_solution(
    model: &Model,
    ndof: usize,
    elements: &[PlasticElem],
    u: &[f64],
    ev: &Evaluation,
    f_ext: &[f64],
    constrained: &[crate::solver::ConstrainedDof],
    iterations: usize,
    residual: f64,
    cond: f64,
) -> Solution {
    let element_stress: Vec<[f64; 6]> = elements
        .iter()
        .map(|element| {
            let count = element.history.len() as f64;
            let mut sum = [0.0; 6];
            for point in &element.history {
                for (a, b) in sum.iter_mut().zip(point.sigma) {
                    *a += b;
                }
            }
            for a in sum.iter_mut() {
                *a /= count;
            }
            sum
        })
        .collect();
    let mut nodal_sum = vec![[0.0; 6]; model.mesh.nodes.len()];
    let mut nodal_count = vec![0usize; model.mesh.nodes.len()];
    for (element, stress) in elements.iter().zip(&element_stress) {
        for &node in &element.node_idx {
            for (a, b) in nodal_sum[node].iter_mut().zip(stress) {
                *a += b;
            }
            nodal_count[node] += 1;
        }
    }
    Solution {
        displacements: model
            .mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let base = i * ndof;
                [
                    u[base],
                    u.get(base + 1).copied().unwrap_or(0.0),
                    u.get(base + 2).copied().unwrap_or(0.0),
                ]
            })
            .collect(),
        von_mises: element_stress.iter().map(von_mises).collect(),
        nodal_von_mises: nodal_sum
            .iter()
            .zip(nodal_count)
            .map(|(stress, count)| {
                if count == 0 {
                    0.0
                } else {
                    let mut average = *stress;
                    for value in average.iter_mut() {
                        *value /= count as f64;
                    }
                    von_mises(&average)
                }
            })
            .collect(),
        reactions: constrained
            .iter()
            .map(|c| Reaction {
                node_id: c.node_id,
                dof: c.dof,
                value: ev.f_int[c.row] - f_ext[c.row],
            })
            .collect(),
        contact_forces: Vec::new(),
        iterations,
        final_residual: residual,
        solver_method: "newton-j2".to_string(),
        cond_estimate: cond,
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_targets(
    model: &Model,
    opts: &NlOptions,
    targets: &[f64],
    prescribed: &[PrescribedDof],
) -> Result<Vec<NlFrame>, SolveError> {
    let loads = model.load_state();
    if loads.bcs.is_empty() {
        return Err(SolveError::NoConstraints);
    }
    model
        .validate()
        .map_err(|e| SolveError::InvalidModel(format!("{e:?}")))?;
    if model.thickness <= 0.0 || model.thickness.is_nan() {
        return Err(SolveError::InvalidModel(
            "thickness must be positive".into(),
        ));
    }
    if opts.max_iterations == 0 || !opts.tolerance.is_finite() || opts.tolerance <= 0.0 {
        return Err(SolveError::InvalidModel(
            "NL options need positive max_iterations and tolerance".into(),
        ));
    }
    let ndof = crate::ndof_per_node(model)?;
    let (_, f_full, _) = assembly::assemble(model, ndof, &loads)?;
    let constrained = constrained_dofs(&model.mesh, loads.bcs, ndof)?;
    let prescribed_rows: Vec<(usize, f64)> = prescribed
        .iter()
        .map(|p| {
            let node = model.mesh.node_index(p.node_id).ok_or_else(|| {
                SolveError::InvalidModel(format!("prescribed unknown node {}", p.node_id))
            })?;
            let dof = p.dof;
            if dof == 0 || dof > ndof {
                return Err(SolveError::InvalidModel(format!(
                    "prescribed node {} has unsupported dof {dof}",
                    p.node_id
                )));
            }
            Ok((node * ndof + dof - 1, p.value))
        })
        .collect::<Result<_, _>>()?;
    let is_constrained = |row: usize| constrained.iter().any(|c| c.row == row);
    for (row, _) in &prescribed_rows {
        if !is_constrained(*row) {
            return Err(SolveError::InvalidModel(
                "displacement-controlled J2 dof is not also constrained".into(),
            ));
        }
    }

    let mut elements = build_elements(model, ndof)?;
    let mut u = vec![0.0; f_full.len()];
    let mut frames = Vec::with_capacity(targets.len());

    for &target in targets {
        let f_ext: Vec<f64> = f_full.iter().map(|value| value * target).collect();
        let free_norm = |values: &[f64]| {
            values
                .iter()
                .enumerate()
                .filter(|(row, _)| !is_constrained(*row))
                .map(|(_, value)| value * value)
                .sum::<f64>()
                .sqrt()
        };

        let mut streak = 0usize;
        let mut previous = f64::INFINITY;
        let mut reference_norm = 0.0f64;
        let mut iterations = 0;
        let mut residual = f64::NAN;
        let mut cond = 1.0;
        let mut converged: Option<Evaluation> = None;
        // The candidate FULL increment for this target, measured from
        // the committed state `u`; the element history stays frozen at
        // the increment start until convergence.  Every Newton and
        // line-search iterate is re-evaluated from that start, so the
        // path-dependent return sees one clean strain path per
        // increment — committing intermediate iterates would let an
        // overshoot-and-correct cycle of the free dofs leave permanent
        // plastic strain behind (visible as drift in quadratic
        // elements).
        let mut step_total = vec![0.0; u.len()];
        for iteration in 1..=opts.max_iterations {
            iterations = iteration;
            for (row, value) in &prescribed_rows {
                step_total[*row] = target * value - u[*row];
            }
            let ev = evaluate(model, ndof, &elements, &step_total)?;
            if reference_norm == 0.0 {
                reference_norm = free_norm(&f_ext).max(free_norm(&ev.f_int));
            }
            let norm = reference_norm;
            let mut sum = 0.0;
            for (row, external) in f_ext.iter().enumerate().take(u.len()) {
                if !is_constrained(row) {
                    let delta = external - ev.f_int[row];
                    sum += delta * delta;
                }
            }
            residual = if norm == 0.0 { 0.0 } else { sum.sqrt() / norm };
            if residual.is_finite() && residual <= opts.tolerance {
                streak += 1;
                if streak >= 2 {
                    converged = Some(ev);
                    break;
                }
            } else {
                streak = 0;
            }
            if !residual.is_finite() || residual > 10.0 * previous {
                break;
            }
            previous = residual;

            let mut k = ev.kt.clone();
            let mut rhs: Vec<f64> = f_ext.iter().zip(&ev.f_int).map(|(a, b)| a - b).collect();
            for c in &constrained {
                rhs[c.row] = 0.0;
            }
            apply_dirichlet(&mut k, &mut rhs, &constrained);
            let factorization = match Ldlt::factorize(&k) {
                Ok(factorization) => factorization,
                Err(_) => {
                    // The consistent tangent is used whenever it
                    // factors.  At a very large plastic increment its
                    // strongly softened directional curvature can be
                    // indefinite; an elastic fallback direction is still a
                    // descent direction and preserves incrementation.
                    let elastic = evaluate_with_tangent(model, ndof, &elements, &step_total, true)?;
                    let mut k = elastic.kt;
                    let mut rhs: Vec<f64> = f_ext
                        .iter()
                        .zip(&elastic.f_int)
                        .map(|(a, b)| a - b)
                        .collect();
                    for c in &constrained {
                        rhs[c.row] = 0.0;
                    }
                    apply_dirichlet(&mut k, &mut rhs, &constrained);
                    Ldlt::factorize(&k)?
                }
            };
            cond = factorization.cond_estimate();
            let mut newton_step = factorization.solve(&rhs);
            for (row, _) in &prescribed_rows {
                newton_step[*row] = 0.0;
            }

            let mut accepted: Option<(Vec<f64>, Evaluation)> = None;
            let mut alpha = 1.0;
            for _ in 0..13 {
                let trial: Vec<f64> = step_total
                    .iter()
                    .zip(&newton_step)
                    .map(|(total, value)| total + alpha * value)
                    .collect();
                let trial_ev = evaluate(model, ndof, &elements, &trial)?;
                let norm = reference_norm;
                let mut sum = 0.0;
                for (row, external) in f_ext.iter().enumerate().take(u.len()) {
                    if !is_constrained(row) {
                        let difference = external - trial_ev.f_int[row];
                        sum += difference * difference;
                    }
                }
                let trial_residual = if norm == 0.0 { 0.0 } else { sum.sqrt() / norm };
                if trial_residual <= residual {
                    accepted = Some((trial, trial_ev));
                    break;
                }
                alpha *= 0.5;
            }
            let Some((trial, _)) = accepted else {
                break;
            };
            step_total = trial;
        }

        let Some(ev) = converged else {
            return Err(SolveError::NotConverged {
                iterations,
                residual,
            });
        };
        for (u_i, step_i) in u.iter_mut().zip(&step_total) {
            *u_i += step_i;
        }
        for (element, history) in elements.iter_mut().zip(&ev.history) {
            element.history = history.clone();
        }
        frames.push(NlFrame {
            load_fraction: target,
            solution: stress_solution(
                model,
                ndof,
                &elements,
                &u,
                &ev,
                &f_ext,
                &constrained,
                iterations,
                residual,
                cond,
            ),
        });
    }
    Ok(frames)
}

pub(crate) fn solve_nl(model: &Model, opts: &NlOptions) -> Result<Vec<NlFrame>, SolveError> {
    if opts.max_increments == 0 {
        return Err(SolveError::InvalidModel(
            "NL options need positive max_increments".into(),
        ));
    }
    let loads = model.load_state();
    let ndof = crate::ndof_per_node(model)?;
    let (_, f_full, _) = assembly::assemble(model, ndof, &loads)?;
    if f_full.iter().all(|value| *value == 0.0) {
        let mut elements = build_elements(model, ndof)?;
        for element in &mut elements {
            element.history.fill(HistoryPoint {
                sigma: [0.0; 6],
                equivalent_plastic_strain: 0.0,
            });
        }
        let constrained = constrained_dofs(&model.mesh, loads.bcs, ndof)?;
        let evaluation = evaluate(model, ndof, &elements, &vec![0.0; f_full.len()])?;
        return Ok(vec![NlFrame {
            load_fraction: 1.0,
            solution: stress_solution(
                model,
                ndof,
                &elements,
                &vec![0.0; f_full.len()],
                &evaluation,
                &f_full,
                &constrained,
                0,
                0.0,
                1.0,
            ),
        }]);
    }
    let targets = (1..=opts.max_increments)
        .map(|increment| increment as f64 / opts.max_increments as f64)
        .collect::<Vec<_>>();
    solve_targets(model, opts, &targets, &[])
}

#[cfg(test)]
fn solve_displacement_path(
    model: &Model,
    opts: &NlOptions,
    factors: &[f64],
    prescribed: &[PrescribedDof],
) -> Result<Vec<NlFrame>, SolveError> {
    solve_targets(model, opts, factors, prescribed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BoundaryCondition, ElementSet, Material};

    const E: f64 = 210.0e9;
    const NU: f64 = 0.3;
    const YIELD: f64 = 100.0e6;
    const HARDENING: f64 = 21.0e9;

    fn cube(yield_stress: f64) -> Model {
        let coordinates = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        Model {
            mesh: crate::Mesh {
                nodes: coordinates
                    .iter()
                    .enumerate()
                    .map(|(index, xyz)| crate::Node {
                        id: index as u32 + 1,
                        coords: *xyz,
                    })
                    .collect(),
                elements: vec![crate::Element {
                    id: 1,
                    kind: ElemKind::C3d8,
                    connectivity: (1..=8).collect(),
                }],
            },
            thickness: 1.0,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: E,
                poisson_ratio: NU,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: vec![1],
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            boundary_conditions: vec![
                BoundaryCondition {
                    node_ids: vec![1, 4, 5, 8],
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: vec![2, 3, 6, 7],
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: vec![1, 2, 5, 6],
                    dofs: vec![2],
                },
                BoundaryCondition {
                    node_ids: vec![1, 2, 3, 4],
                    dofs: vec![3],
                },
            ],
            plasticity: Some(PlasticitySpec {
                yield_stress,
                hardening_modulus: HARDENING,
            }),
            ..Model::default()
        }
    }

    fn prescribed() -> Vec<PrescribedDof> {
        [2, 3, 6, 7]
            .into_iter()
            .map(|node_id| PrescribedDof {
                node_id,
                dof: 1,
                value: 1.0,
            })
            .collect()
    }

    fn solve_to_strain(model: &Model, strain: f64) -> f64 {
        let factors = (1..=20)
            .map(|increment| strain * increment as f64 / 20.0)
            .collect::<Vec<_>>();
        let frames =
            solve_displacement_path(model, &NlOptions::default(), &factors, &prescribed()).unwrap();
        frames.last().unwrap().solution.von_mises[0]
    }

    #[test]
    fn radial_return_unloads_to_zero_stress_and_keeps_plastic_strain() {
        let spec = PlasticitySpec {
            yield_stress: YIELD,
            hardening_modulus: HARDENING,
        };
        let elastic = d3_matrix(E, NU);
        let old = HistoryPoint {
            sigma: [0.0; 6],
            equivalent_plastic_strain: 0.0,
        };
        let strain = [0.002, -NU * 0.002, -NU * 0.002, 0.0, 0.0, 0.0];
        let loaded = radial_return(spec, &elastic, E, NU, old, &strain);
        assert!(loaded.equivalent_plastic_strain > 0.0);

        let shear = E / (2.0 * (1.0 + NU));
        let elastic_unload = [
            (loaded.sigma[0] - NU * (loaded.sigma[1] + loaded.sigma[2])) / E,
            (loaded.sigma[1] - NU * (loaded.sigma[0] + loaded.sigma[2])) / E,
            (loaded.sigma[2] - NU * (loaded.sigma[0] + loaded.sigma[1])) / E,
            loaded.sigma[3] / shear,
            loaded.sigma[4] / shear,
            loaded.sigma[5] / shear,
        ];
        let unload = [
            -elastic_unload[0],
            -elastic_unload[1],
            -elastic_unload[2],
            -elastic_unload[3],
            -elastic_unload[4],
            -elastic_unload[5],
        ];
        let unloaded = radial_return(
            spec,
            &elastic,
            E,
            NU,
            HistoryPoint {
                sigma: loaded.sigma,
                equivalent_plastic_strain: loaded.equivalent_plastic_strain,
            },
            &unload,
        );
        assert!(
            unloaded.sigma.iter().all(|sigma| sigma.abs() < 1.0e-7),
            "unloaded sigma={:?}, plastic={} (loaded {:?},{})",
            unloaded.sigma,
            unloaded.equivalent_plastic_strain,
            loaded.sigma,
            loaded.equivalent_plastic_strain
        );
        assert_eq!(
            unloaded.equivalent_plastic_strain,
            loaded.equivalent_plastic_strain
        );
    }

    #[test]
    fn consistent_tangent_matches_finite_difference_of_return_map() {
        let spec = PlasticitySpec {
            yield_stress: YIELD,
            hardening_modulus: HARDENING,
        };
        let elastic = d3_matrix(E, NU);
        let old = HistoryPoint {
            sigma: [0.0; 6],
            equivalent_plastic_strain: 0.0,
        };
        let strain = [0.002, -0.0003, 0.0002, 0.001, -0.0004, 0.0003];
        let returned = radial_return(spec, &elastic, E, NU, old, &strain);
        let h = 1.0e-9;
        for column in 0..6 {
            let mut forward = strain;
            let mut backward = strain;
            forward[column] += h;
            backward[column] -= h;
            let plus = radial_return(spec, &elastic, E, NU, old, &forward);
            let minus = radial_return(spec, &elastic, E, NU, old, &backward);
            for row in 0..6 {
                let finite_difference = (plus.sigma[row] - minus.sigma[row]) / (2.0 * h);
                assert!(
                    (finite_difference - returned.tangent[row][column]).abs()
                        < 1.0e-3 * returned.tangent[row][column].abs().max(1.0),
                    "C[{row}][{column}]={} vs fd={finite_difference}",
                    returned.tangent[row][column]
                );
            }
        }
    }

    #[test]
    fn c3d8_uniaxial_displacement_control_follows_elastic_plastic_curve() {
        let model = cube(YIELD);
        let tangent = E * HARDENING / (E + HARDENING);
        let cases = [
            (0.0003, E * 0.0003),
            (0.001, YIELD + tangent * (0.001 - YIELD / E)),
            (0.002, YIELD + tangent * (0.002 - YIELD / E)),
        ];
        for (strain, expected) in cases {
            let actual = solve_to_strain(&model, strain);
            eprintln!(
                "J2 C3D8: eps={strain:.6}, sigma={actual:.6e} Pa, expected={expected:.6e} Pa"
            );
            assert!(
                (actual - expected).abs() / expected < 0.01,
                "eps={strain}, sigma={actual:.6e}, expected={expected:.6e}"
            );
        }
    }

    #[test]
    fn c3d8_unload_path_reaches_equilibrium() {
        let model = cube(YIELD);
        let loaded_sigma = solve_to_strain(&model, 0.002);
        let permanent_set = (loaded_sigma - YIELD) / HARDENING;
        let mut factors = (1..=20)
            .map(|increment| 0.002 * increment as f64 / 20.0)
            .collect::<Vec<_>>();
        factors.extend(
            (1..=20).map(|increment| 0.002 + (permanent_set - 0.002) * increment as f64 / 20.0),
        );
        let frames =
            solve_displacement_path(&model, &NlOptions::default(), &factors, &prescribed())
                .unwrap();
        assert!(
            (frames[19].solution.von_mises[0] - 1.2909090909e8).abs() < 1.0e6,
            "loaded sigma={:.9e}",
            frames[19].solution.von_mises[0]
        );
        assert!(frames.last().unwrap().solution.final_residual < 1.0e-6);
    }

    #[test]
    fn high_yield_matches_linear_elastic_solution() {
        let mut loaded = cube(1.0e12);
        loaded.boundary_conditions.remove(1);
        loaded.cloads = [2, 3, 6, 7]
            .into_iter()
            .map(|node_id| crate::ConcentratedLoad {
                node_id,
                dof: 1,
                value: 1.2909090909e8 / 4.0,
            })
            .collect();
        let mut linear = loaded.clone();
        linear.plasticity = None;
        let actual = crate::solve(&loaded).unwrap();
        let expected = crate::solve(&linear).unwrap();
        assert!(
            ((actual.displacements[1][0] - expected.displacements[1][0])
                / expected.displacements[1][0])
                .abs()
                < 1.0e-10,
            "high-yield tip={:.9e}, linear tip={:.9e}",
            actual.displacements[1][0],
            expected.displacements[1][0]
        );
        assert_eq!(actual.solver_method, "newton-j2");
    }

    #[test]
    fn b31_plasticity_is_rejected() {
        let mut model = cube(YIELD);
        model.mesh.elements[0].kind = ElemKind::B31;
        model.mesh.elements[0].connectivity.truncate(2);
        model.beam_sections = vec![("all".into(), 1.0, 1.0)];
        match solve_displacement_path(
            &model,
            &NlOptions::default(),
            &[0.001],
            &[PrescribedDof {
                node_id: 2,
                dof: 1,
                value: 1.0,
            }],
        ) {
            Err(SolveError::InvalidModel(message)) => {
                assert!(message.contains("B31"), "{message}")
            }
            other => panic!("expected InvalidModel, got {other:?}"),
        }
    }

    /// Single-element (or two-triangle) plane models under uniaxial
    /// displacement control: ux fixed on x=0, ux prescribed on x=1,
    /// uy fixed on y=0 only, so sigma_yy = 0 settles on the free edge.
    fn plane_model(
        yield_stress: f64,
        coords: &[[f64; 2]],
        elements: Vec<(ElemKind, Vec<u32>)>,
        fixed_x: &[u32],
        prescribed_x: &[u32],
        fixed_y: &[u32],
    ) -> Model {
        let count = elements.len() as u32;
        Model {
            mesh: crate::Mesh {
                nodes: coords
                    .iter()
                    .enumerate()
                    .map(|(index, xy)| crate::Node {
                        id: index as u32 + 1,
                        coords: [xy[0], xy[1], 0.0],
                    })
                    .collect(),
                elements: elements
                    .into_iter()
                    .enumerate()
                    .map(|(index, (kind, connectivity))| crate::Element {
                        id: index as u32 + 1,
                        kind,
                        connectivity,
                    })
                    .collect(),
            },
            thickness: 1.0,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: E,
                poisson_ratio: NU,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=count).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            boundary_conditions: vec![
                BoundaryCondition {
                    node_ids: fixed_x.to_vec(),
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: prescribed_x.to_vec(),
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: fixed_y.to_vec(),
                    dofs: vec![2],
                },
            ],
            plasticity: Some(PlasticitySpec {
                yield_stress,
                hardening_modulus: HARDENING,
            }),
            ..Model::default()
        }
    }

    fn uniaxial_prescribed(nodes: &[u32]) -> Vec<PrescribedDof> {
        nodes
            .iter()
            .map(|&node_id| PrescribedDof {
                node_id,
                dof: 1,
                value: 1.0,
            })
            .collect()
    }

    fn solve_last_von_mises(model: &Model, target: f64, prescribed: &[PrescribedDof]) -> f64 {
        let factors = (1..=20)
            .map(|increment| target * increment as f64 / 20.0)
            .collect::<Vec<_>>();
        let frames =
            solve_displacement_path(model, &NlOptions::default(), &factors, prescribed).unwrap();
        frames.last().unwrap().solution.von_mises[0]
    }

    fn plane_stress_return_for(
        spec: PlasticitySpec,
        old: HistoryPoint,
        strain: [f64; 6],
    ) -> ReturnResult {
        plane_stress_return(spec, E, NU, old, &strain)
    }

    #[test]
    fn plane_stress_return_stays_on_surface_with_zero_sigma_zz() {
        let spec = PlasticitySpec {
            yield_stress: YIELD,
            hardening_modulus: HARDENING,
        };
        let mut old = HistoryPoint {
            sigma: [0.0; 6],
            equivalent_plastic_strain: 0.0,
        };
        // Load into yield, keep loading with shear, unload (elastic),
        // then reload past yield again: every plastic step must sit on
        // sigma_zz = 0 and on the hardened surface.
        let steps = [
            [0.002, -0.0003, 0.0, 0.001, 0.0, 0.0],
            [0.0012, -0.0002, 0.0, -0.0007, 0.0, 0.0],
            [-0.0015, 0.0003, 0.0, 0.0002, 0.0, 0.0],
            [0.0035, -0.0005, 0.0, -0.0018, 0.0, 0.0],
        ];
        for (step, strain) in steps.iter().enumerate() {
            let returned = plane_stress_return_for(spec, old, *strain);
            let plastic = returned.equivalent_plastic_strain > old.equivalent_plastic_strain;
            old = HistoryPoint {
                sigma: returned.sigma,
                equivalent_plastic_strain: returned.equivalent_plastic_strain,
            };
            assert!(
                returned.sigma[2].abs() <= 1.0e-8 * YIELD,
                "step {step}: sigma_zz={:.3e}",
                returned.sigma[2]
            );
            assert_eq!(returned.sigma[4], 0.0);
            assert_eq!(returned.sigma[5], 0.0);
            let surface = YIELD + HARDENING * returned.equivalent_plastic_strain;
            let actual = von_mises(&returned.sigma);
            if plastic {
                assert!(
                    (actual - surface).abs() <= 1.0e-8 * surface,
                    "step {step}: von Mises={actual:.9e}, surface={surface:.9e}"
                );
            } else {
                assert!(
                    actual <= surface * (1.0 + 1.0e-8),
                    "step {step}: {actual:.9e}"
                );
            }
        }
        assert!(old.equivalent_plastic_strain > 0.0);
    }

    #[test]
    fn plane_stress_return_tangent_matches_finite_difference() {
        let spec = PlasticitySpec {
            yield_stress: YIELD,
            hardening_modulus: HARDENING,
        };
        // Elastic: the tangent is the plane-stress elastic matrix.
        let elastic = plane_stress_return_for(
            spec,
            HistoryPoint {
                sigma: [0.0; 6],
                equivalent_plastic_strain: 0.0,
            },
            [0.0002, -0.00005, 0.0, 0.00003, 0.0, 0.0],
        );
        let d = d_matrix(E, NU);
        for (i, &si) in IN_PLANE_SLOTS.iter().enumerate() {
            for (j, &sj) in IN_PLANE_SLOTS.iter().enumerate() {
                assert!(
                    (elastic.tangent[si][sj] - d[i][j]).abs() <= 1.0e-9 * d[i][j].abs().max(1.0),
                    "C[{si}][{sj}]={}",
                    elastic.tangent[si][sj]
                );
            }
        }

        // Plastic: central differences of the return map (the shipped
        // tangent itself uses forward differences, h = 1e-9).
        let old = HistoryPoint {
            sigma: [1.0e8, 2.0e7, 0.0, 3.0e7, 0.0, 0.0],
            equivalent_plastic_strain: 0.0004,
        };
        let strain = [0.0008, -0.0002, 0.0, 0.0005, 0.0, 0.0];
        let returned = plane_stress_return_for(spec, old, strain);
        let h = 1.0e-9;
        for &slot in IN_PLANE_SLOTS.iter() {
            let mut forward = strain;
            let mut backward = strain;
            forward[slot] += h;
            backward[slot] -= h;
            let plus = plane_stress_return_for(spec, old, forward);
            let minus = plane_stress_return_for(spec, old, backward);
            for &row in IN_PLANE_SLOTS.iter() {
                let finite_difference = (plus.sigma[row] - minus.sigma[row]) / (2.0 * h);
                assert!(
                    (finite_difference - returned.tangent[row][slot]).abs()
                        < 1.0e-3 * returned.tangent[row][slot].abs().max(1.0),
                    "C[{row}][{slot}]={} vs fd={finite_difference}",
                    returned.tangent[row][slot]
                );
            }
        }
    }

    #[test]
    fn cps4_uniaxial_displacement_control_follows_elastic_plastic_curve() {
        let model = plane_model(
            YIELD,
            &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            vec![(ElemKind::Cps4, vec![1, 2, 3, 4])],
            &[1, 4],
            &[2, 3],
            &[1, 2],
        );
        let prescribed = uniaxial_prescribed(&[2, 3]);
        let tangent = E * HARDENING / (E + HARDENING);
        let cases = [
            (0.0003, E * 0.0003),
            (0.001, YIELD + tangent * (0.001 - YIELD / E)),
            (0.002, YIELD + tangent * (0.002 - YIELD / E)),
        ];
        for (strain, expected) in cases {
            let actual = solve_last_von_mises(&model, strain, &prescribed);
            eprintln!(
                "J2 CPS4: eps={strain:.6}, sigma={actual:.6e} Pa, expected={expected:.6e} Pa"
            );
            assert!(
                (actual - expected).abs() / expected < 0.015,
                "eps={strain}, sigma={actual:.6e}, expected={expected:.6e}"
            );
        }
    }

    #[test]
    fn cps3_and_cps8_match_the_plane_stress_uniaxial_curve() {
        let tangent = E * HARDENING / (E + HARDENING);
        let expected = YIELD + tangent * (0.002 - YIELD / E);
        let cases = [
            (
                plane_model(
                    YIELD,
                    &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
                    vec![
                        (ElemKind::Cps3, vec![1, 2, 3]),
                        (ElemKind::Cps3, vec![1, 3, 4]),
                    ],
                    &[1, 4],
                    &[2, 3],
                    &[1, 2],
                ),
                "CPS3",
                uniaxial_prescribed(&[2, 3]),
            ),
            (
                plane_model(
                    YIELD,
                    &[
                        [0.0, 0.0],
                        [1.0, 0.0],
                        [1.0, 1.0],
                        [0.0, 1.0],
                        [0.5, 0.0],
                        [1.0, 0.5],
                        [0.5, 1.0],
                        [0.0, 0.5],
                    ],
                    vec![(ElemKind::Cps8, (1..=8).collect())],
                    &[1, 4, 8],
                    &[2, 3, 6],
                    &[1, 2, 5],
                ),
                "CPS8",
                uniaxial_prescribed(&[2, 3, 6]),
            ),
        ];
        for (model, name, loaded) in cases {
            let factors: Vec<f64> = (1..=20)
                .map(|increment| 0.002 * increment as f64 / 20.0)
                .collect();
            let frames =
                solve_displacement_path(&model, &NlOptions::default(), &factors, &loaded).unwrap();
            if name == "CPS8" {
                eprintln!(
                    "{name} final u: {:?}",
                    frames.last().unwrap().solution.displacements
                );
                eprintln!(
                    "{name} residual: {}",
                    frames.last().unwrap().solution.final_residual
                );
            }
            let actual = frames.last().unwrap().solution.von_mises[0];
            assert!(
                (actual - expected).abs() / expected < 0.015,
                "{name}: sigma={actual:.6e}, expected={expected:.6e}"
            );
        }
    }

    #[test]
    fn cps4_unload_path_reaches_zero_stress() {
        let model = plane_model(
            YIELD,
            &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            vec![(ElemKind::Cps4, vec![1, 2, 3, 4])],
            &[1, 4],
            &[2, 3],
            &[1, 2],
        );
        let prescribed = uniaxial_prescribed(&[2, 3]);
        let loaded_sigma = solve_last_von_mises(&model, 0.002, &prescribed);
        let permanent_set = (loaded_sigma - YIELD) / HARDENING;
        let mut factors = (1..=20)
            .map(|increment| 0.002 * increment as f64 / 20.0)
            .collect::<Vec<_>>();
        factors.extend(
            (1..=20).map(|increment| 0.002 + (permanent_set - 0.002) * increment as f64 / 20.0),
        );
        let frames =
            solve_displacement_path(&model, &NlOptions::default(), &factors, &prescribed).unwrap();
        assert!(
            (frames[19].solution.von_mises[0] - 1.2909090909e8).abs() < 1.0e6,
            "loaded sigma={:.9e}",
            frames[19].solution.von_mises[0]
        );
        let last = frames.last().unwrap();
        assert!(
            last.solution.von_mises[0] < 1.0e-3 * YIELD,
            "residual sigma={:.6e}",
            last.solution.von_mises[0]
        );
        assert!(last.solution.final_residual < 1.0e-6);
    }

    #[test]
    fn plane_shear_plasticity_matches_gh_over_3g_plus_h() {
        // Homogeneous pure shear imposed by full Dirichlet: ux = g/2*y,
        // uy = g/2*x (the symmetric split of the engineering shear g).
        // The state sigma = (0, 0, 0, tau) is then exact for both plane
        // stress and plane strain, with tau(g) = G*g while elastic and
        // slope G*H/(3G + H) after yield — q = sqrt(3)*tau = sigma_y +
        // H*p and g = tau/G + sqrt(3)*p.
        let shear = E / (2.0 * (1.0 + NU));
        let sqrt3 = (3.0f64).sqrt();
        let plastic =
            |gamma: f64| (gamma + sqrt3 * YIELD / HARDENING) / (1.0 / shear + 3.0 / HARDENING);
        let cases = [
            (2.0e-4, shear * 2.0e-4),
            (1.5e-3, plastic(1.5e-3)),
            (3.0e-3, plastic(3.0e-3)),
        ];
        for kind in [ElemKind::Cps4, ElemKind::Cpe4] {
            let mut model = plane_model(
                YIELD,
                &[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
                vec![(kind, vec![1, 2, 3, 4])],
                &[1, 2],
                &[3, 4],
                &[1, 4],
            );
            // uy = gamma/2 on the right edge too (full Dirichlet).
            model.boundary_conditions.push(BoundaryCondition {
                node_ids: vec![2, 3],
                dofs: vec![2],
            });
            let prescribed = [
                PrescribedDof {
                    node_id: 3,
                    dof: 1,
                    value: 0.5,
                },
                PrescribedDof {
                    node_id: 4,
                    dof: 1,
                    value: 0.5,
                },
                PrescribedDof {
                    node_id: 2,
                    dof: 2,
                    value: 0.5,
                },
                PrescribedDof {
                    node_id: 3,
                    dof: 2,
                    value: 0.5,
                },
            ];
            for (gamma, tau) in cases {
                let actual = solve_last_von_mises(&model, gamma, &prescribed) / sqrt3;
                eprintln!(
                    "J2 shear {kind:?}: gamma={gamma:.6}, tau={actual:.6e}, expected={tau:.6e}"
                );
                assert!(
                    (actual - tau).abs() / tau < 0.01,
                    "{kind:?}: gamma={gamma}, tau={actual:.6e}, expected={tau:.6e}"
                );
            }
        }
    }
}
