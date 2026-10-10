//! # Geometric nonlinearity: Updated Lagrangian + full Newton
//!
//! [`solve_nl`] applies the step's total load in equal increments of
//! `load_fraction` and iterates full Newton on
//! `f_ext(lambda) = f_int(u)`.
//!
//! Formulation (v1, documented):
//! - Updated Lagrangian with SEQUENTIAL CONFIGURATION REFRESH: the
//!   reference configuration is the state at the start of each Newton
//!   iteration. Stress history is kept per Gauss point (CPS3:
//!   constant, [sxx, syy, szz, sxy]) and committed along the accepted
//!   Newton / line-search path with the TRUESDELL objective rate
//!   `dsigma = D deps + L sigma + sigma L^T - sigma tr(L)` (L = the
//!   step's displacement gradient) — the rate whose rigid-rotation
//!   behavior is exact and whose consistent tangent is exactly
//!   `K_mat + K_sigma`. The trial internal force is
//!   `f_int = Σ_gp w detJ(x) B_gp(x)^T sigma_gp(du)` and the tangent
//!   is `K_mat(x) + K_sigma(sigma)` with the material stiffness at
//!   the trial geometry and the per-Gauss-point geometric stiffness
//!   reusing the buckling integrals' structure ([`crate::buckle`]) on
//!   current coordinates. Verified on the bench cantilever: Newton
//!   converges in 4-5 iterations per increment up to 25 deg tip
//!   rotations, and the load-deflection curve tracks the exact
//!   elastica within the mesh's linear shear-locking band (~1.5%).
//! - Physical direction (documented deviation from the issue text): a
//!   dead-tip-loaded cantilever STIFFENS, `delta_NL/delta_lin < 1` —
//!   the vertical tip deflection saturates as the beam rotates
//!   (elastica ratio 0.9198 at Pbar = PL^2/EI = 0.914; the issue's
//!   "NL gives larger displacement" has the direction backwards for
//!   this benchmark). Axial stretching softens (carried stress decays
//!   as the section grows), covered by a dedicated test.
//! - Why a moving B never meets a large displacement multiplier: those
//!   products are numerically poisoned on thin bending elements —
//!   `ke(x)·u` cancels large terms (shear-locking conditioning,
//!   ~1e2-1e3 on the bench mesh) and the Jacobian perturbation rides
//!   that amplification into the residual as a scale-free bias
//!   (measured ~7% element-force error at a 5e-4 rad rotation).
//!   A Total-Lagrangian Green-Lagrange variant was also tried and
//!   rejected: the rotation-induced `E22 += ½ u_x,y^2` term produces
//!   a spurious transverse stress `~ E theta^2 / 2` worth ~40%
//!   residual error at `delta/L = 0.7%`. Strain increments are the
//!   small-strain tensor at the reference configuration
//!   ("Green-Lagrange small-term approximation").
//! - Loads are DEAD: assembled once on the original geometry
//!   (follower loads are the upgrade path).
//! - Frame von Mises uses the Gauss-point-averaged carried stress;
//!   nodal extrapolated recovery stays on the original-geometry
//!   linear machinery (v1 ceiling).
//! - Plane elements only (CPS3/CPS4/CPE4); other kinds return
//!   `InvalidModel("NL not yet supported for kind X")`.
//! - Equal load increments v1. Newton carries a 13-halving
//!   backtracking line search; on divergence (residual growing 10x
//!   after the line search, or the iteration budget exhausted) the
//!   increment is cut in half and retried while the increment size
//!   stays above 1/512 of the step, then the solve fails with
//!   [`SolveError::NotConverged`] (the issue asked for a single retry
//!   level; repeated halving is the robustness margin — automatic
//!   incrementation remains the documented upgrade).
//! - ponytail: an LDLT pivot failure mid-Newton (tangent losing
//!   positive-definiteness at a limit point) surfaces as
//!   `InvalidModel` rather than triggering cutback — fine until
//!   snap-through problems matter.

use crate::assembly::{self, Csr};
use crate::buckle::geo_scalar_tri;
use crate::element::{cps3_stiffness, cps4_stiffness, invalid, quad_b, von_mises, Constitutive};
use crate::ldl::Ldlt;
use crate::solver::{apply_dirichlet, constrained_dofs};
use crate::{ElemKind, Model, Reaction, Solution, SolveError};

/// Newton/incrementation controls for [`solve_nl`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NlOptions {
    /// Number of equal load increments over the step (default 20).
    pub max_increments: usize,
    /// Newton iteration budget PER increment (default 16).
    pub max_iterations: usize,
    /// Convergence tolerance on the relative free-dof force residual
    /// `||f_ext - f_int|| / ||f_ext||`, required on two consecutive
    /// iterations (default 1e-6).
    pub tolerance: f64,
}

impl Default for NlOptions {
    fn default() -> Self {
        Self {
            max_increments: 20,
            max_iterations: 16,
            tolerance: 1e-6,
        }
    }
}

/// One converged load increment: the fraction of the step's total load
/// and the solution at that fraction. `solution.displacements` holds
/// TOTAL displacements relative to the original geometry.
#[derive(Debug, Clone, Default)]
pub struct NlFrame {
    pub load_fraction: f64,
    pub solution: Solution,
}

/// 2x2 Gauss abscissae (unit weights), stiffness-loop order.
const GP2: f64 = 0.577_350_269_189_625_8;
const QUAD_GPS: [(f64, f64); 4] = [(-GP2, -GP2), (-GP2, GP2), (GP2, -GP2), (GP2, GP2)];

/// Per-element data for the UL loop (geometry-independent except for
/// the carried strain history).
struct UlElem {
    node_idx: Vec<usize>,
    /// Plane-stress / plane-strain 3x3 constitutive matrix.
    d: [[f64; 3]; 3],
    plane_strain: bool,
    /// Accumulated stress per Gauss point [sxx, syy, szz, sxy],
    /// committed along the accepted path with the Truesdell objective
    /// rate (the rate whose consistent tangent is exactly
    /// `K_mat + K_sigma`).
    sigma_gp: Vec<[f64; 4]>,
}

/// State evaluation at one trial point (reference `u_ref`, trial step
/// `step`).
struct Eval {
    kt: Csr,
    f_int: Vec<f64>,
    /// Stress history including the trial step (to commit on
    /// acceptance).
    sigma_new: Vec<Vec<[f64; 4]>>,
}

/// Closed-form CPS3 B rows at the given coordinates.
fn tri_b(c: &[[f64; 2]; 3]) -> [[f64; 6]; 3] {
    let b = [c[1][1] - c[2][1], c[2][1] - c[0][1], c[0][1] - c[1][1]];
    let cc = [c[2][0] - c[1][0], c[0][0] - c[2][0], c[1][0] - c[0][0]];
    let two_a = b[0] * cc[1] - b[1] * cc[0];
    let inv = 1.0 / two_a;
    let mut bm = [[0.0; 6]; 3];
    for i in 0..3 {
        bm[0][2 * i] = b[i] * inv;
        bm[1][2 * i + 1] = cc[i] * inv;
        bm[2][2 * i] = cc[i] * inv;
        bm[2][2 * i + 1] = b[i] * inv;
    }
    bm
}

/// Evaluate tangent and internal force at the trial point
/// `u_ref + step` (see the [module docs](self) for the formulation
/// and why the reference rows carry the strain increment while the
/// trial geometry carries the force and the tangent).
fn evaluate(
    model: &Model,
    elems: &[UlElem],
    u_ref: &[f64],
    step: &[f64],
) -> Result<Eval, SolveError> {
    let ndof = 2usize; // the kind whitelist in solve_nl guarantees plane-only
    let t = model.thickness;
    let n = model.mesh.nodes.len() * ndof;
    let mut trips = Vec::new();
    for i in 0..n {
        trips.push((i, i, 0.0)); // diagonal entries for the Dirichlet pass
    }
    let mut f_int = vec![0.0; n];
    let mut sigma_new = Vec::with_capacity(elems.len());

    for (ei, e) in model.mesh.elements.iter().enumerate() {
        let el = &elems[ei];
        let node = |k: usize| el.node_idx[k];
        let xy = |k: usize, with_step: bool| {
            let i = node(k);
            let (a, b) = if with_step {
                (step[i * ndof], step[i * ndof + 1])
            } else {
                (0.0, 0.0)
            };
            [
                model.mesh.nodes[i].coords[0] + u_ref[i * ndof] + a,
                model.mesh.nodes[i].coords[1] + u_ref[i * ndof + 1] + b,
            ]
        };
        let due: Vec<f64> = el
            .node_idx
            .iter()
            .flat_map(|&i| (0..ndof).map(move |dd| step[i * ndof + dd]))
            .collect();
        let ref3 = [xy(0, false), xy(1, false), xy(2, false)];
        let ref4 = [xy(0, false), xy(1, false), xy(2, false), xy(3, false)];
        let cur3 = [xy(0, true), xy(1, true), xy(2, true)];
        let cur4 = [xy(0, true), xy(1, true), xy(2, true), xy(3, true)];
        // Strain increments and displacement gradients from the FIXED
        // reference rows (tame products); the stress is updated with
        // the Truesdell objective rate
        // `dsigma = D deps + L sigma + sigma L^T - sigma tr(L)`
        // (L = the step's displacement gradient), whose rigid-rotation
        // behavior is exact and whose consistent tangent is exactly
        // K_mat + K_sigma.
        let mut sig_e: Vec<[f64; 4]> = el.sigma_gp.clone();
        let truesdell = |sig: &mut [f64; 4],
                         de: [f64; 3],
                         l: [[f64; 2]; 2],
                         d: &[[f64; 3]; 3],
                         plane_strain: bool| {
            let trl = l[0][0] + l[1][1];
            sig[0] +=
                d[0][0] * de[0] + d[0][1] * de[1] + 2.0 * (l[0][0] * sig[0] + l[0][1] * sig[3])
                    - sig[0] * trl;
            sig[1] +=
                d[1][0] * de[0] + d[1][1] * de[1] + 2.0 * (l[1][0] * sig[3] + l[1][1] * sig[1])
                    - sig[1] * trl;
            sig[2] += if plane_strain {
                d[0][1] * (de[0] + de[1])
            } else {
                0.0
            } - sig[2] * trl;
            sig[3] += d[2][2] * de[2]
                + (l[0][0] * sig[3] + l[0][1] * sig[1])
                + (sig[0] * l[1][0] + sig[3] * l[1][1])
                - sig[3] * trl;
        };
        let grads = |bm: &[[f64; 8]; 3]| -> ([f64; 3], [[f64; 2]; 2]) {
            let (mut ux_x, mut ux_y, mut uy_x, mut uy_y) = (0.0, 0.0, 0.0, 0.0);
            for i in 0..4 {
                let bx = bm[2][2 * i + 1];
                let by = bm[2][2 * i];
                ux_x += bx * due[2 * i];
                ux_y += by * due[2 * i];
                uy_x += bx * due[2 * i + 1];
                uy_y += by * due[2 * i + 1];
            }
            (
                [
                    bm[0].iter().zip(&due).map(|(a, c)| a * c).sum(),
                    bm[1].iter().zip(&due).map(|(a, c)| a * c).sum(),
                    bm[2].iter().zip(&due).map(|(a, c)| a * c).sum(),
                ],
                [[ux_x, ux_y], [uy_x, uy_y]],
            )
        };
        match e.kind {
            ElemKind::Cps3 => {
                let bmt = tri_b(&ref3);
                let (mut ux_x, mut ux_y, mut uy_x, mut uy_y) = (0.0, 0.0, 0.0, 0.0);
                for i in 0..3 {
                    let bx = bmt[2][2 * i + 1];
                    let by = bmt[2][2 * i];
                    ux_x += bx * due[2 * i];
                    ux_y += by * due[2 * i];
                    uy_x += bx * due[2 * i + 1];
                    uy_y += by * due[2 * i + 1];
                }
                let de = [
                    bmt[0].iter().zip(&due).map(|(a, c)| a * c).sum(),
                    bmt[1].iter().zip(&due).map(|(a, c)| a * c).sum(),
                    bmt[2].iter().zip(&due).map(|(a, c)| a * c).sum(),
                ];
                truesdell(
                    &mut sig_e[0],
                    de,
                    [[ux_x, ux_y], [uy_x, uy_y]],
                    &el.d,
                    el.plane_strain,
                );
            }
            _ => {
                for (k, &(xi, eta)) in QUAD_GPS.iter().enumerate() {
                    let (bm, _) = quad_b(&ref4, xi, eta).ok_or_else(|| {
                        invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                    })?;
                    let (de, l) = grads(&bm);
                    truesdell(&mut sig_e[k], de, l, &el.d, el.plane_strain);
                }
            }
        }
        // Internal force AND per-Gauss-point geometric stiffness in one
        // pass over the trial Gauss points. The GP-resolved K_sigma
        // matters: the exact Jacobian of the trial map carries the
        // GP-resolved geometric part, and an element-mean K_sigma
        // breaks Newton descent on bending-dominated states.
        let mut fe = vec![0.0; due.len()];
        let ks: Vec<Vec<f64>> = match e.kind {
            ElemKind::Cps3 => {
                let bm = tri_b(&cur3);
                let cross = (cur3[1][0] - cur3[0][0]) * (cur3[2][1] - cur3[0][1])
                    - (cur3[1][1] - cur3[0][1]) * (cur3[2][0] - cur3[0][0]);
                let s = [sig_e[0][0], sig_e[0][1], sig_e[0][3]];
                let area = cross.abs() / 2.0;
                for i in 0..3 {
                    fe[2 * i] += t * area * (bm[0][2 * i] * s[0] + bm[2][2 * i] * s[2]);
                    fe[2 * i + 1] += t * area * (bm[1][2 * i + 1] * s[1] + bm[2][2 * i + 1] * s[2]);
                }
                geo_scalar_tri(&cur3, &s, t)
                    .map(|k| k.iter().map(|r| r.to_vec()).collect())
                    .ok_or_else(|| invalid(e.id, "degenerate element in geometric stiffness"))?
            }
            _ => {
                let mut ks4 = [[0.0f64; 4]; 4];
                for (k, &(xi, eta)) in QUAD_GPS.iter().enumerate() {
                    let (bm, detj) = quad_b(&cur4, xi, eta).ok_or_else(|| {
                        invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                    })?;
                    let s = [sig_e[k][0], sig_e[k][1], sig_e[k][3]];
                    for i in 0..4 {
                        fe[2 * i] += t * detj * (bm[0][2 * i] * s[0] + bm[2][2 * i] * s[2]);
                        fe[2 * i + 1] +=
                            t * detj * (bm[1][2 * i + 1] * s[1] + bm[2][2 * i + 1] * s[2]);
                    }
                    let g = |i: usize| (bm[2][2 * i + 1], bm[2][2 * i]);
                    for (i, row) in ks4.iter_mut().enumerate() {
                        let (gx, gy) = g(i);
                        for (j, cell) in row.iter_mut().enumerate() {
                            let (hx, hy) = g(j);
                            let v = gx * (s[0] * hx + s[2] * hy) + gy * (s[2] * hx + s[1] * hy);
                            *cell += t * detj * v;
                        }
                    }
                }
                ks4.iter().map(|r| r.to_vec()).collect()
            }
        };
        for (a, &fv) in fe.iter().enumerate() {
            f_int[node(a / ndof) * ndof + a % ndof] += fv;
        }
        // Tangent: K_mat(trial) + K_sigma at the trial coordinates.
        let ke: Vec<Vec<f64>> = match e.kind {
            ElemKind::Cps3 => {
                let (ke, _) = cps3_stiffness(&cur3, &el.d, t)
                    .ok_or_else(|| invalid(e.id, "degenerate (zero-area) triangle"))?;
                ke.iter().map(|r| r.to_vec()).collect()
            }
            _ => {
                let (ke, _) = cps4_stiffness(&cur4, &el.d, t).ok_or_else(|| {
                    invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                })?;
                ke.iter().map(|r| r.to_vec()).collect()
            }
        };
        for (a, kea) in ke.iter().enumerate() {
            let grow = node(a / ndof) * ndof + a % ndof;
            for (c, &kv) in kea.iter().enumerate() {
                trips.push((grow, node(c / ndof) * ndof + c % ndof, kv));
            }
        }
        for (i, ksi) in ks.iter().enumerate() {
            for (j, &ksv) in ksi.iter().enumerate() {
                if ksv == 0.0 {
                    continue;
                }
                for dd in 0..ndof {
                    trips.push((node(i) * ndof + dd, node(j) * ndof + dd, ksv));
                }
            }
        }
        sigma_new.push(sig_e);
    }

    Ok(Eval {
        kt: assembly::csr_from_triplets(n, trips),
        f_int,
        sigma_new,
    })
}

/// The zero-load shortcut: homogeneous BCs and no load on any free dof
/// means u = 0 is the equilibrium, reported as one full-load frame.
fn zero_solution(model: &Model) -> Solution {
    Solution {
        displacements: vec![[0.0, 0.0, 0.0]; model.mesh.nodes.len()],
        von_mises: vec![0.0; model.mesh.elements.len()],
        nodal_von_mises: vec![0.0; model.mesh.nodes.len()],
        reactions: Vec::new(),
        iterations: 0,
        final_residual: 0.0,
        solver_method: "newton".to_string(),
        cond_estimate: 1.0,
    }
}

/// Solve the model's load state (last step, or the legacy top-level
/// fields) with geometric nonlinearity: Updated Lagrangian kinematics
/// and full Newton; see the [module docs](self) for the v1 formulation
/// and ceilings. Returns one [`NlFrame`] per converged load increment,
/// ending at `load_fraction = 1.0`.
///
/// Errors: the usual [`Model::validate`] failures and `NoConstraints`;
/// `InvalidModel` for non-plane element kinds, bad options, or a
/// non-positive-definite tangent; `NotConverged` when Newton diverges
/// (or exhausts its budget) even after the allowed increment
/// cutbacks.
pub fn solve_nl(model: &Model, opts: &NlOptions) -> Result<Vec<NlFrame>, SolveError> {
    let loads = model.load_state();
    if loads.bcs.is_empty() {
        return Err(SolveError::NoConstraints);
    }
    if let Err(e) = model.validate() {
        return Err(SolveError::InvalidModel(format!("{e:?}")));
    }
    if model.thickness <= 0.0 || model.thickness.is_nan() {
        return Err(SolveError::InvalidModel(
            "thickness must be positive".into(),
        ));
    }
    if opts.max_increments == 0 || opts.max_iterations == 0 {
        return Err(SolveError::InvalidModel(
            "NL options need positive max_increments and max_iterations".into(),
        ));
    }
    if !opts.tolerance.is_finite() || opts.tolerance <= 0.0 {
        return Err(SolveError::InvalidModel(
            "NL tolerance must be positive".into(),
        ));
    }
    for e in &model.mesh.elements {
        match e.kind {
            ElemKind::Cps3 | ElemKind::Cps4 | ElemKind::Cpe4 => {}
            other => {
                return Err(SolveError::InvalidModel(format!(
                    "NL not yet supported for kind {other:?}"
                )))
            }
        }
    }

    let ndof = 2usize;
    // Loads are DEAD: assembled once, on the original geometry.
    let (_, f_full, elems) = assembly::assemble(model, ndof, &loads)?;
    let ngp = |kind: ElemKind| if kind == ElemKind::Cps3 { 1 } else { 4 };
    let mut ul: Vec<UlElem> = elems
        .iter()
        .zip(&model.mesh.elements)
        .map(|(info, e)| UlElem {
            node_idx: info.node_idx.clone(),
            d: match &info.d {
                Constitutive::Plane(d) | Constitutive::PlaneStrain(d) => *d,
                _ => unreachable!("NL kind whitelist keeps only plane elements"),
            },
            plane_strain: matches!(info.d, Constitutive::PlaneStrain(_)),
            sigma_gp: vec![[0.0; 4]; ngp(e.kind)],
        })
        .collect();
    let constrained = constrained_dofs(&model.mesh, loads.bcs, ndof)?;
    let mut is_constrained = vec![false; f_full.len()];
    for c in &constrained {
        is_constrained[c.row] = true;
    }
    // Residuals only make sense on the free dofs: at constrained rows
    // f_int - f_ext IS the reaction, never small.
    let free_norm = |v: &[f64]| {
        v.iter()
            .enumerate()
            .filter(|(i, _)| !is_constrained[*i])
            .map(|(_, x)| x * x)
            .sum::<f64>()
            .sqrt()
    };
    if free_norm(&f_full) == 0.0 {
        return Ok(vec![NlFrame {
            load_fraction: 1.0,
            solution: zero_solution(model),
        }]);
    }

    let residual_of = |ev: &Eval, f_ext: &[f64], f_norm: f64| -> f64 {
        let mut s = 0.0;
        for (i, (a, b)) in f_ext.iter().zip(&ev.f_int).enumerate() {
            if !is_constrained[i] {
                let d = a - b;
                s += d * d;
            }
        }
        s.sqrt() / f_norm
    };

    let mut u = vec![0.0; f_full.len()];
    let zero_step = vec![0.0; f_full.len()];
    let mut frames = Vec::with_capacity(opts.max_increments);
    let mut lam = 0.0f64;
    let mut dlam = 1.0 / opts.max_increments as f64;
    while lam < 1.0 - 1e-12 {
        let target = (lam + dlam).min(1.0);
        let f_ext: Vec<f64> = f_full.iter().map(|v| v * target).collect();
        let f_norm = free_norm(&f_ext);
        let sig_start: Vec<Vec<[f64; 4]>> = ul.iter().map(|e| e.sigma_gp.clone()).collect();
        let u_start = u.clone();
        let mut streak = 0usize;
        let mut prev_res = f64::INFINITY;
        let mut iterations = 0usize;
        let mut residual = f64::NAN;
        let mut cond = 1.0;
        let mut frame_ev: Option<Eval> = None;
        for it in 1..=opts.max_iterations {
            iterations = it;
            // Reference = the committed state; with a zero trial step
            // this residual has no mixed-geometry part at all.
            let ev = evaluate(model, &ul, &u, &zero_step)?;
            residual = residual_of(&ev, &f_ext, f_norm);
            if residual.is_finite() && residual <= opts.tolerance {
                streak += 1;
                if streak >= 2 {
                    frame_ev = Some(ev);
                    break;
                }
            } else {
                streak = 0;
            }
            if !residual.is_finite() || residual > 10.0 * prev_res {
                break; // diverging
            }
            prev_res = residual;
            let mut k = ev.kt.clone();
            let mut rhs: Vec<f64> = f_ext.iter().zip(&ev.f_int).map(|(a, b)| a - b).collect();
            for c in &constrained {
                rhs[c.row] = 0.0;
            }
            apply_dirichlet(&mut k, &mut rhs, &constrained);
            let fact = Ldlt::factorize(&k)?;
            cond = fact.cond_estimate();
            let step = fact.solve(&rhs);
            // Backtracking line search on the trial step; acceptance
            // commits u and the strain history (the reference refresh).
            let mut accepted: Option<(Vec<f64>, Eval)> = None;
            let mut alpha = 1.0f64;
            for _ in 0..13 {
                let alpha_step: Vec<f64> = step.iter().map(|s| alpha * s).collect();
                let ev_t = evaluate(model, &ul, &u, &alpha_step)?;
                let res_t = residual_of(&ev_t, &f_ext, f_norm);
                if res_t <= residual {
                    accepted = Some((alpha_step, ev_t));
                    break;
                }
                alpha *= 0.5;
            }
            let Some((alpha_step, ev_t)) = accepted else {
                break; // no descent along the Newton direction
            };
            for (ui, &d) in u.iter_mut().zip(&alpha_step) {
                *ui += d;
            }
            for (el, e_new) in ul.iter_mut().zip(&ev_t.sigma_new) {
                el.sigma_gp = e_new.clone();
            }
        }
        match frame_ev {
            Some(ev) => {
                lam = target;
                frames.push(NlFrame {
                    load_fraction: target,
                    solution: Solution {
                        displacements: model
                            .mesh
                            .nodes
                            .iter()
                            .enumerate()
                            .map(|(i, _)| [u[i * ndof], u[i * ndof + 1], 0.0])
                            .collect(),
                        von_mises: ul
                            .iter()
                            .map(|el| {
                                let m = el.sigma_gp.len() as f64;
                                let mut s = [0.0f64; 4];
                                for g in &el.sigma_gp {
                                    for c in 0..4 {
                                        s[c] += g[c];
                                    }
                                }
                                von_mises(&[s[0] / m, s[1] / m, s[2] / m, s[3] / m, 0.0, 0.0])
                            })
                            .collect(),
                        nodal_von_mises: crate::stress::nodal_von_mises(
                            &u,
                            ndof,
                            &elems,
                            model.mesh.nodes.len(),
                            None,
                        ),
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
                        solver_method: "newton".to_string(),
                        cond_estimate: cond,
                    },
                });
            }
            None => {
                if dlam <= 1.0 / 512.0 {
                    return Err(SolveError::NotConverged {
                        iterations,
                        residual,
                    });
                }
                // Halve the increment and retry from the last converged
                // state: restore BOTH the strain history and the
                // committed displacements of this increment's start.
                dlam *= 0.5;
                u = u_start.clone();
                for (el, e0) in ul.iter_mut().zip(&sig_start) {
                    el.sigma_gp = e0.clone();
                }
            }
        }
    }
    Ok(frames)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BoundaryCondition, ConcentratedLoad, ElementSet, Material, Method};

    const E: f64 = 210e9;

    /// The bench cantilever: 40x4 CPS4, L = 2, h = t = 0.1, clamped
    /// root, tip load `p` (N) on the right-edge midnode. Returns the
    /// model and the tip node id. P = 1e4 gives the 15.238 mm
    /// Euler-Bernoulli tip deflection.
    fn cantilever(p: f64) -> (Model, u32) {
        let (nx, ny) = (40usize, 4);
        let (l, h) = (2.0, 0.1);
        let mut model = Model {
            mesh: crate::assembly::tests::quad_mesh(nx, ny, l, h),
            thickness: 0.1,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: E,
                poisson_ratio: 0.3,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=(nx * ny) as u32).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            ..Default::default()
        };
        let root: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2],
        }];
        let tip = (2 * (nx + 1) + nx + 1) as u32;
        model.cloads = vec![ConcentratedLoad {
            node_id: tip,
            dof: 2,
            value: p,
        }];
        (model, tip)
    }

    fn tip_uy(model: &Model, tip: u32, sol: &Solution) -> f64 {
        sol.displacements[model.mesh.node_index(tip).unwrap()][1]
    }

    #[test]
    fn cantilever_moderate_load_matches_elastica() {
        // P = 4e5 N (Pbar = PL^2/EI = 0.914): linear FEM tip ~0.553 m,
        // tip rotation ~25 deg. The exact elastica answer for a
        // dead-loaded cantilever is STIFFENING — the vertical tip
        // deflection saturates as the beam rotates (dy/dy_lin =
        // 0.9198 by numerical integration of the elastica; the issue
        // text's "NL gives larger displacement" has the direction
        // backwards for this benchmark). The FEM ratio sits slightly
        // above the elastica ratio because both the linear reference
        // and the NL solution ride the same shear-locked mesh.
        // (The issue's suggested P = 200 N gives delta/L ~ 1e-3 and an
        // NL effect ~1e-5 — below any testable tolerance.)
        let (model, tip) = cantilever(-4.0e5);
        let lin = crate::solve_with_method(&model, Method::Ldlt).unwrap();
        let lin_tip = tip_uy(&model, tip, &lin);
        let frames = solve_nl(&model, &NlOptions::default()).unwrap();
        let nl_tip = tip_uy(&model, tip, &frames.last().unwrap().solution);
        let ratio = nl_tip / lin_tip;
        eprintln!(
            "cantilever P=4e5N (Pbar=0.914): linear tip {lin_tip:.6} m, NL tip {nl_tip:.6} m, ratio {ratio:.4} (elastica 0.9198), {} frames",
            frames.len()
        );
        assert!(
            nl_tip.abs() < lin_tip.abs(),
            "expected large-rotation stiffening (|nl| < |linear|), ratio {ratio:.4}"
        );
        assert!(
            (ratio - 0.9198).abs() < 0.03,
            "NL/elastica mismatch: ratio {ratio:.4} vs 0.9198"
        );
        // Equilibrium: the clamped reactions balance the full tip load.
        let ry: f64 = frames
            .last()
            .unwrap()
            .solution
            .reactions
            .iter()
            .filter(|r| r.dof == 2)
            .map(|r| r.value)
            .sum();
        assert!((ry - 4.0e5).abs() < 1e-4 * 4.0e5, "sum R_y = {ry:.3}");
    }

    #[test]
    fn newton_converges_in_few_iterations() {
        let (model, _) = cantilever(-4.0e5);
        let frames = solve_nl(&model, &NlOptions::default()).unwrap();
        assert_eq!(frames.len(), 20);
        for f in &frames {
            assert!(
                f.solution.iterations <= 6,
                "increment at lambda={:.3} took {} iterations",
                f.load_fraction,
                f.solution.iterations
            );
        }
    }

    #[test]
    fn small_load_nl_matches_linear() {
        // P = 1e4 (the 15.238 mm bench): geometric effects are O((d/L)^2)
        // ~ 1e-4, far inside the 0.5% consistency band.
        let (model, tip) = cantilever(-1.0e4);
        let lin = crate::solve_with_method(&model, Method::Ldlt).unwrap();
        let frames = solve_nl(&model, &NlOptions::default()).unwrap();
        let a = tip_uy(&model, tip, &lin);
        let b = tip_uy(&model, tip, &frames.last().unwrap().solution);
        assert!(
            (a - b).abs() / a.abs() < 5e-3,
            "linear {a:.6e} vs nl {b:.6e}"
        );
    }

    #[test]
    fn load_deflection_curve_is_monotone() {
        let (model, tip) = cantilever(-4.0e5);
        let frames = solve_nl(&model, &NlOptions::default()).unwrap();
        let i = model.mesh.node_index(tip).unwrap();
        let mut prev_frac = 0.0;
        let mut prev_disp = 0.0;
        for f in &frames {
            assert!(f.load_fraction > prev_frac);
            let d = f.solution.displacements[i][1].abs();
            assert!(
                d > prev_disp,
                "tip |uy| not increasing at lambda={:.3}",
                f.load_fraction
            );
            prev_frac = f.load_fraction;
            prev_disp = d;
        }
        assert!((frames.last().unwrap().load_fraction - 1.0).abs() < 1e-12);
    }

    #[test]
    fn unsupported_element_kinds_are_rejected() {
        let mut model =
            crate::tests::solid_model(crate::assembly::tests::hex_mesh(2, 2, 2, 1.0, 1.0, 1.0));
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![1, 2, 3],
        }];
        match solve_nl(&model, &NlOptions::default()) {
            Err(SolveError::InvalidModel(msg)) => {
                assert!(msg.contains("NL not yet supported for kind C3d8"), "{msg}")
            }
            other => panic!("expected InvalidModel, got {other:?}"),
        }
    }

    #[test]
    fn zero_load_step_returns_single_zero_frame() {
        let (model, tip) = cantilever(0.0);
        let frames = solve_nl(&model, &NlOptions::default()).unwrap();
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].load_fraction, 1.0);
        assert_eq!(tip_uy(&model, tip, &frames[0].solution), 0.0);
        assert!(frames[0]
            .solution
            .displacements
            .iter()
            .all(|d| *d == [0.0, 0.0, 0.0]));
    }

    #[test]
    fn axial_bar_large_load_softens() {
        // 20x1 CPS4 bar, fixed left, uniform axial tip load of 2e6 N
        // (nominal sigma = P/A = 8 GPa, linear stretch ~3.8%): the UL
        // machinery is exact for axial stretching and converges at
        // large strain; the sequential refresh softens the response
        // (the carried stress decays as the section grows), so the NL
        // tip exceeds the linear PL/EA value by O(epsilon).
        let (nx, ny) = (20usize, 1);
        let mut model = Model {
            mesh: crate::assembly::tests::quad_mesh(nx, ny, 2.0, 0.05),
            thickness: 0.005,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: E,
                poisson_ratio: 0.3,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=(nx * ny) as u32).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            ..Default::default()
        };
        let root: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2],
        }];
        let p = 2.0e6;
        model.cloads = (0..=ny)
            .map(|j| ConcentratedLoad {
                node_id: (j * (nx + 1) + nx + 1) as u32,
                dof: 1,
                value: p / 2.0,
            })
            .collect();
        let frames = solve_nl(&model, &NlOptions::default()).unwrap();
        let sol = &frames.last().unwrap().solution;
        let tip = sol.displacements[nx][0];
        let lin = p * 2.0 / (E * 0.005 * 0.05);
        let maxit = frames.iter().map(|f| f.solution.iterations).max().unwrap();
        eprintln!("axial: nl tip ux={tip:.6} PL/EA={lin:.6} maxit={maxit}");
        assert!(tip > lin, "expected axial softening, got {tip} vs {lin}");
        assert!(
            tip < lin * 1.06,
            "axial softening too large: {tip} vs {lin}"
        );
        // Reactions balance the full load.
        let rx: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 1)
            .map(|r| r.value)
            .sum();
        assert!((rx + p).abs() < 1e-4 * p, "sum R_x = {rx:.3}");
    }

    #[test]
    fn no_bcs_is_no_constraints() {
        let (mut model, _) = cantilever(-1.0);
        model.boundary_conditions.clear();
        assert!(matches!(
            solve_nl(&model, &NlOptions::default()),
            Err(SolveError::NoConstraints)
        ));
    }
}
