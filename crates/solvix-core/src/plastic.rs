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
//! Element support is deliberately explicit.  C3D8/C3D4 and plane-strain
//! CPE4/CPE8 can use the direct six-component return.  Plane-stress
//! CPS3/CPS4/CPS8 is rejected with `InvalidModel`: the return surface is
//! not closed under a plain six-component return because `sigma_zz = 0`
//! must be enforced while `eps_zz` evolves.  A dedicated three-component
//! fixed-point plane-stress return is the upgrade path.  Plasticity also
//! takes precedence over NLGEOM and uses the small-strain branch documented
//! here.

use crate::assembly::{self, Csr};
use crate::element::{d3_matrix, hex_b, quad8_b, quad_b, tet_gradients, von_mises, G3, W3};
use crate::ldl::Ldlt;
use crate::solver::{apply_dirichlet, constrained_dofs};
use crate::{ElemKind, Model, NlFrame, NlOptions, PlasticitySpec, Reaction, Solution, SolveError};

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
    d: [[f64; 6]; 6],
    youngs_modulus: f64,
    poisson_ratio: f64,
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

fn element_points(
    model: &Model,
    element: &crate::Element,
    node_idx: &[usize],
) -> Result<Vec<IntegrationPoint>, SolveError> {
    let invalid =
        |message: &str| SolveError::InvalidModel(format!("element {}: {message}", element.id));
    match element.kind {
        ElemKind::Cpe4 => {
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
        ElemKind::Cpe8 => {
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
        ElemKind::Cps3 | ElemKind::Cps4 | ElemKind::Cps8 | ElemKind::B31 => {
            Err(SolveError::InvalidModel(format!(
                "J2 plasticity is not yet supported for plane-stress kind {:?}",
                element.kind
            )))
        }
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
            d: d3_matrix(material.youngs_modulus, material.poisson_ratio),
            youngs_modulus: material.youngs_modulus,
            poisson_ratio: material.poisson_ratio,
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
            let mut strain = [0.0; 6];
            for (b_row, strain_component) in point.b.iter().zip(strain.iter_mut()) {
                *strain_component = b_row
                    .iter()
                    .enumerate()
                    .map(|(column, b)| {
                        let node = element.node_idx[column / ndof];
                        b * step[node * ndof + column % ndof]
                    })
                    .sum();
            }
            let returned = radial_return(
                spec,
                &element.d,
                element.youngs_modulus,
                element.poisson_ratio,
                *old,
                &strain,
            );
            new_history.push(HistoryPoint {
                sigma: returned.sigma,
                equivalent_plastic_strain: returned.equivalent_plastic_strain,
            });
            let tangent = if elastic_tangent {
                element.d
            } else {
                returned.tangent
            };

            for (a, force) in element_force.iter_mut().enumerate() {
                for (b_row, sigma) in point.b.iter().zip(returned.sigma) {
                    *force += point.weight * b_row[a] * sigma;
                }
            }
            for a in 0..element_force.len() {
                for c in a..element_force.len() {
                    let mut value = 0.0;
                    for (r, b_row) in point.b.iter().enumerate() {
                        let mut db = 0.0;
                        for (k, b_k) in point.b.iter().enumerate() {
                            db += tangent[r][k] * b_k[c];
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
        let prescribed_delta = |u: &[f64], out: &mut Vec<f64>| {
            for (row, value) in &prescribed_rows {
                out[*row] = target * value - u[*row];
            }
        };
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
        for iteration in 1..=opts.max_iterations {
            iterations = iteration;
            let mut trial_step = vec![0.0; u.len()];
            prescribed_delta(&u, &mut trial_step);
            let ev = evaluate(model, ndof, &elements, &trial_step)?;
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
                    // The closed-form consistent tangent is used whenever
                    // it factors.  At a very large plastic increment its
                    // strongly softened directional curvature can be
                    // indefinite; an elastic fallback direction is still a
                    // descent direction and preserves incrementation.
                    let elastic = evaluate_with_tangent(model, ndof, &elements, &trial_step, true)?;
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
            let mut delta = vec![0.0; u.len()];
            prescribed_delta(&u, &mut delta);

            let mut accepted: Option<(Vec<f64>, Evaluation)> = None;
            let mut alpha = 1.0;
            for _ in 0..13 {
                let trial: Vec<f64> = newton_step
                    .iter()
                    .enumerate()
                    .map(|(row, value)| {
                        if prescribed_rows
                            .iter()
                            .any(|(prescribed, _)| *prescribed == row)
                        {
                            delta[row]
                        } else {
                            alpha * value
                        }
                    })
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
            let Some((trial, trial_ev)) = accepted else {
                break;
            };
            for (u_i, trial_i) in u.iter_mut().zip(trial) {
                *u_i += trial_i;
            }
            for (element, history) in elements.iter_mut().zip(&trial_ev.history) {
                element.history = history.clone();
            }
        }

        let Some(ev) = converged else {
            return Err(SolveError::NotConverged {
                iterations,
                residual,
            });
        };
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
    fn plane_stress_is_explicitly_rejected() {
        let mut model = cube(YIELD);
        model.mesh.elements[0].kind = ElemKind::Cps4;
        model.mesh.elements[0].connectivity.truncate(4);
        model.boundary_conditions[2].dofs = vec![2];
        model.boundary_conditions[3].dofs = vec![2];
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
                assert!(message.contains("plane-stress"), "{message}")
            }
            other => panic!("expected InvalidModel, got {other:?}"),
        }
    }
}
