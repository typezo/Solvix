//! # Transient dynamics: implicit Newmark (average acceleration)
//!
//! `M u'' + C u' + K u = f(t)` integrated with the unconditionally
//! stable, second-order-accurate Newmark average-acceleration scheme
//! (beta = 1/4, gamma = 1/2: no numerical damping, conserves energy in
//! undamped free vibration), driven by a step's [`DynamicsSpec`]
//! (`*DYNAMIC` in the deck).
//!
//! v1 choices, deliberate and documented:
//! - Constant time step `dt` for the whole `duration`. The step's loads
//!   (cloads, dloads, gravity of the accumulated state) are a step
//!   function applied at t = 0 and held constant — `*AMPLITUDE` time
//!   curves are the upgrade path.
//! - Rayleigh damping `C = alpha M + beta K` (`*DAMPING`).
//! - Homogeneous Dirichlet only: constrained dofs carry zero
//!   displacement, velocity and acceleration for all time.
//! - Initial conditions: rest (`u = v = 0`) with the consistent initial
//!   acceleration `a_0 = M^-1 f`; `*INITIAL CONDITIONS` is the upgrade
//!   path.
//! - Constant `dt` and constant loads make the effective stiffness
//!   `K_eff = K + a0 M + a1 C` constant too: it is factorized ONCE with
//!   the shared sparse LDLT and reused for every step — cheaper than a
//!   PCG solve per step and exact for this linear problem.

use crate::assembly::{self, Csr};
use crate::ldl::Ldlt;
use crate::modal;
use crate::solver;
use crate::{LoadState, Model, SolveError};

/// Scheme constants for [`solve_dynamic`]: Newmark average acceleration
/// (beta = 1/4, gamma = 1/2). v1 pins the defaults; the struct exists so
/// a later version can expose per-call overrides without an API break.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DynOptions {
    /// Newmark beta; 1/4 = average acceleration.
    pub beta: f64,
    /// Newmark gamma; 1/2 = second-order accurate, no numerical damping.
    pub gamma: f64,
}

impl Default for DynOptions {
    fn default() -> Self {
        Self {
            beta: 0.25,
            gamma: 0.5,
        }
    }
}

/// State at one output time. All vectors use the model dof layout
/// (`node_index * ndof + (dof - 1)`), exactly zero at constrained dofs.
#[derive(Debug, Clone, Default)]
pub struct TimePoint {
    pub time: f64,
    pub displacements: Vec<f64>,
    pub velocities: Vec<f64>,
    pub accelerations: Vec<f64>,
}

/// Transient solution: the state at t = 0 plus one point per time step.
#[derive(Debug, Clone, Default)]
pub struct DynSolution {
    /// `n_steps + 1` points, evenly spaced by `dt` from 0 to `duration`.
    // ponytail: every point stored (3 vectors x ndof per step) — memory
    // ceiling ~ 24 B * dofs * steps; subsample or stream to file if
    // long runs ever need it.
    pub points: Vec<TimePoint>,
    /// Time steps taken (= direct linear solves; the PCG path would
    /// accumulate its iteration counts here instead).
    pub iterations_total: usize,
}

/// `y += a * x`.
fn axpy(y: &mut [f64], x: &[f64], a: f64) {
    for (y, x) in y.iter_mut().zip(x) {
        *y += a * x;
    }
}

/// `k += c * m` in place. Both come from `csr_from_triplets` over the
/// same element scatter with a guaranteed diagonal, and every element
/// stiffness block is structurally full, so m's sparsity pattern is a
/// subset of k's; m-only positions cannot occur (the walk skips them
/// defensively rather than rebuild the pattern). Shared with the
/// nonlinear solver (nl.rs adds K_sigma into the tangent the same way).
pub(crate) fn add_scaled(k: &mut Csr, m: &Csr, c: f64) {
    let n = k.row_ptr.len() - 1;
    for r in 0..n {
        let mut mi = m.row_ptr[r];
        for idx in k.row_ptr[r]..k.row_ptr[r + 1] {
            let col = k.col_idx[idx];
            while mi < m.row_ptr[r + 1] && m.col_idx[mi] < col {
                mi += 1;
            }
            if mi < m.row_ptr[r + 1] && m.col_idx[mi] == col {
                k.vals[idx] += c * m.vals[mi];
                mi += 1;
            }
        }
    }
}

/// Integrate the LAST step carrying [`DynamicsSpec`] (`*DYNAMIC`) with
/// implicit Newmark; see the [module docs](self) for the scheme and v1
/// ceilings. The step's accumulated loads and BCs define the problem;
/// earlier steps are ignored.
///
/// Errors: `InvalidModel` when no step defines dynamics, for the usual
/// [`Model::validate`] failures, non-positive `density`/`thickness`,
/// non-positive `dt`/`duration`, a duration shorter than one step, or a
/// request of more than 5e7 dof-time-points (output memory ceiling);
/// `NoConstraints` for an unconstrained model.
pub fn solve_dynamic(model: &Model) -> Result<DynSolution, SolveError> {
    let opts = DynOptions::default();
    let step = model
        .steps
        .iter()
        .rev()
        .find(|s| s.dynamics.is_some())
        .ok_or_else(|| {
            SolveError::InvalidModel(
                "no step defines *DYNAMIC (solve_dynamic needs a DynamicsSpec)".into(),
            )
        })?;
    let spec = step.dynamics.as_ref().expect("find checked is_some");
    let loads = LoadState {
        bcs: &step.boundary_conditions,
        cloads: &step.cloads,
        dloads: &step.dloads,
        gravity: step.gravity,
    };
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
    if model.density <= 0.0 || model.density.is_nan() {
        return Err(SolveError::InvalidModel("density must be positive".into()));
    }
    if !spec.dt.is_finite() || spec.dt <= 0.0 {
        return Err(SolveError::InvalidModel(
            "*DYNAMIC time increment must be positive".into(),
        ));
    }
    if !spec.duration.is_finite() || spec.duration <= 0.0 {
        return Err(SolveError::InvalidModel(
            "*DYNAMIC duration must be positive".into(),
        ));
    }
    let n_steps = (spec.duration / spec.dt).round();
    if n_steps < 1.0 {
        return Err(SolveError::InvalidModel(
            "*DYNAMIC duration must span at least one time increment".into(),
        ));
    }
    let n_steps = n_steps as usize;

    let ndof = crate::ndof_per_node(model)?;
    let n = model.mesh.nodes.len() * ndof;
    if n_steps.saturating_mul(n) > 50_000_000 {
        return Err(SolveError::InvalidModel(format!(
            "*DYNAMIC asks for {n_steps} steps x {n} dofs of output — refused (ceiling 5e7 dof-time-points; shorten the run or raise dt)"
        )));
    }

    let (k, f, _) = assembly::assemble(model, ndof, &loads)?;
    let m = modal::assemble_mass(model, ndof)?;
    let constrained = solver::constrained_dofs(&model.mesh, loads.bcs, ndof)?;

    let (ra, rb) = (spec.rayleigh_alpha, spec.rayleigh_beta);
    let dt = spec.dt;
    let (beta, gamma) = (opts.beta, opts.gamma);
    let a0 = 1.0 / (beta * dt * dt);
    let a1 = gamma / (beta * dt);
    let a2 = 1.0 / (beta * dt);
    let a3 = 1.0 / (2.0 * beta) - 1.0;
    let a4 = gamma / beta - 1.0;
    let a5 = dt / 2.0 * (gamma / beta - 2.0);

    // K_eff = K + a0 M + a1 C = (1 + a1 rb) K + (a0 + a1 ra) M,
    // then homogeneous Dirichlet and ONE factorization for all steps.
    let mut k_eff = Csr {
        row_ptr: k.row_ptr.clone(),
        col_idx: k.col_idx.clone(),
        vals: k.vals.iter().map(|v| v * (1.0 + a1 * rb)).collect(),
    };
    add_scaled(&mut k_eff, &m, a0 + a1 * ra);
    let mut rhs = f.clone();
    solver::apply_dirichlet(&mut k_eff, &mut rhs, &constrained);
    let k_fact = Ldlt::factorize(&k_eff)?;

    // Consistent initial acceleration a_0 = M^-1 f (rest: u = v = 0).
    let mut m_dir = Csr {
        row_ptr: m.row_ptr.clone(),
        col_idx: m.col_idx.clone(),
        vals: m.vals.clone(),
    };
    solver::apply_dirichlet(&mut m_dir, &mut rhs, &constrained);
    let m_fact = Ldlt::factorize(&m_dir)?;
    let mut u = vec![0.0; n];
    let mut vel = vec![0.0; n];
    let mut acc = m_fact.solve(&rhs);

    // C = ra M + rb K folded into per-vector coefficients.
    let (c_mu, c_ku) = (a0 + a1 * ra, a1 * rb);
    let (c_mv, c_kv) = (a2 + a4 * ra, a4 * rb);
    let (c_ma, c_ka) = (a3 + a5 * ra, a5 * rb);

    let mut points = Vec::with_capacity(n_steps + 1);
    points.push(TimePoint {
        time: 0.0,
        displacements: u.clone(),
        velocities: vel.clone(),
        accelerations: acc.clone(),
    });

    for s in 1..=n_steps {
        let mu = m.mat_vec(&u);
        let ku = k.mat_vec(&u);
        let mv = m.mat_vec(&vel);
        let kv = k.mat_vec(&vel);
        let ma = m.mat_vec(&acc);
        let ka = k.mat_vec(&acc);
        let mut fe = f.clone();
        axpy(&mut fe, &mu, c_mu);
        axpy(&mut fe, &ku, c_ku);
        axpy(&mut fe, &mv, c_mv);
        axpy(&mut fe, &kv, c_kv);
        axpy(&mut fe, &ma, c_ma);
        axpy(&mut fe, &ka, c_ka);
        for c in &constrained {
            fe[c.row] = 0.0;
        }
        let u_next = k_fact.solve(&fe);
        let acc_next: Vec<f64> = u_next
            .iter()
            .zip(u.iter().zip(vel.iter().zip(&acc)))
            .map(|(un, (ui, (vi, ai)))| a0 * (un - ui) - a2 * vi - a3 * ai)
            .collect();
        let vel_next: Vec<f64> = vel
            .iter()
            .zip(acc.iter().zip(&acc_next))
            .map(|(vi, (ai, an))| vi + dt * ((1.0 - gamma) * ai + gamma * an))
            .collect();
        u = u_next;
        vel = vel_next;
        acc = acc_next;
        points.push(TimePoint {
            time: s as f64 * dt,
            displacements: u.clone(),
            velocities: vel.clone(),
            accelerations: acc.clone(),
        });
    }

    Ok(DynSolution {
        points,
        iterations_total: n_steps,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modal::{solve_modal, ModalOptions};
    use crate::{
        solve, BoundaryCondition, ConcentratedLoad, DynamicsSpec, ElementSet, Material,
        Step as StepData,
    };

    const E: f64 = 210e9;
    const NU: f64 = 0.3;
    const RHO: f64 = 7850.0;
    const T: f64 = 0.01;

    fn steel_cantilever(nx: usize, ny: usize, l: f64, h: f64, p: f64) -> Model {
        let mesh = crate::assembly::tests::quad_mesh(nx, ny, l, h);
        let n_elem = mesh.elements.len() as u32;
        let root: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        let tip = (2 * (nx + 1) + nx + 1) as u32;
        Model {
            mesh,
            thickness: T,
            density: RHO,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: E,
                poisson_ratio: NU,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=n_elem).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            steps: vec![StepData {
                name: "dyn".into(),
                boundary_conditions: vec![BoundaryCondition {
                    node_ids: root,
                    dofs: vec![1, 2],
                }],
                cloads: vec![ConcentratedLoad {
                    node_id: tip,
                    dof: 2,
                    value: p,
                }],
                dynamics: Some(DynamicsSpec {
                    dt: 1e-3,
                    duration: 1e-2,
                    rayleigh_alpha: 0.0,
                    rayleigh_beta: 0.0,
                }),
                ..Default::default()
            }],
            ..Default::default()
        }
    }

    fn tip_dof(model: &Model, nx: usize) -> usize {
        let tip = (2 * (nx + 1) + nx + 1) as u32;
        model.mesh.node_index(tip).unwrap() * 2 + 1
    }

    #[test]
    fn forced_response_frequency_and_energy_stability() {
        let (nx, ny, l, h) = (40, 4, 2.0, 0.1);
        let mut model = steel_cantilever(nx, ny, l, h, -1000.0);
        model.steps[0].dynamics = Some(DynamicsSpec {
            dt: 2.5e-4,
            duration: 0.5,
            rayleigh_alpha: 0.0,
            rayleigh_beta: 0.0,
        });
        let f1 = solve_modal(&model, &ModalOptions { num_modes: 3 })
            .unwrap()
            .frequencies_hz[0];

        let sol = solve_dynamic(&model).unwrap();
        assert_eq!(sol.points.len(), 2001);
        assert_eq!(sol.iterations_total, 2000);

        // Energy: KE + SE bounded by the step-load ceiling (max dynamic
        // amplification 2 => E <= 4 x SE_static; average acceleration
        // conserves it exactly in theory, allow roundoff slack).
        let loads = model.load_state();
        let (k, f, _) = assembly::assemble(&model, 2, &loads).unwrap();
        let m = modal::assemble_mass(&model, 2).unwrap();
        let stat = solve(&model).unwrap();
        let u_st: Vec<f64> = stat
            .displacements
            .iter()
            .flat_map(|d| [d[0], d[1]])
            .collect();
        let se_st = 0.5 * modal::dot(&f, &u_st);
        let mut e_max = 0.0f64;
        for p in &sol.points {
            let ke = 0.5 * modal::dot(&p.velocities, &m.mat_vec(&p.velocities));
            let se = 0.5 * modal::dot(&p.displacements, &k.mat_vec(&p.displacements));
            assert!((ke + se).is_finite());
            e_max = e_max.max(ke + se);
        }
        assert!(
            e_max < 4.5 * se_st,
            "energy {e_max:.3e} not bounded by 4.5 x static {se_st:.3e} (blow-up)"
        );

        // Tip oscillation frequency from mean-subtracted zero crossings
        // (constant force: response = static deflection + f1 transient),
        // skipping the first 10 steps of startup.
        let ti = tip_dof(&model, nx);
        let win = &sol.points[10..];
        let mean: f64 = win.iter().map(|p| p.displacements[ti]).sum::<f64>() / win.len() as f64;
        let crossings = win
            .windows(2)
            .filter(|w| {
                let x0 = w[0].displacements[ti] - mean;
                let x1 = w[1].displacements[ti] - mean;
                (x0 > 0.0 && x1 < 0.0) || (x0 < 0.0 && x1 > 0.0)
            })
            .count();
        let window_time = win.last().unwrap().time - win.first().unwrap().time;
        let freq = crossings as f64 / 2.0 / window_time;
        let rel = (freq - f1) / f1;
        eprintln!(
            "forced response: zero-crossing f={freq:.3} Hz vs modal f1={f1:.3} Hz ({:+.2}%), crossings={crossings}, E_max/SE_st={:.3}",
            100.0 * rel,
            e_max / se_st
        );
        assert!(rel.abs() <= 0.10, "frequency off by {:.1}%", 100.0 * rel);
    }

    #[test]
    fn static_limit_under_heavy_damping() {
        let (nx, ny, l, h) = (20, 4, 1.0, 0.1);
        let mut model = steel_cantilever(nx, ny, l, h, -500.0);
        model.steps[0].dynamics = Some(DynamicsSpec {
            dt: 1.0e-4,
            duration: 0.5,
            rayleigh_alpha: 100.0,
            rayleigh_beta: 0.0,
        });
        let sol = solve_dynamic(&model).unwrap();
        let stat = solve(&model).unwrap();
        let ti = tip_dof(&model, nx);
        let u_dyn = sol.points.last().unwrap().displacements[ti];
        let u_stat = stat.displacements[model
            .mesh
            .node_index((2 * (nx + 1) + nx + 1) as u32)
            .unwrap()][1];
        let rel = (u_dyn - u_stat) / u_stat;
        eprintln!(
            "static limit: tip {u_dyn:.6e} vs static {u_stat:.6e} ({:+.3}%)",
            100.0 * rel
        );
        assert!(rel.abs() <= 0.01, "static limit off by {:.3}%", 100.0 * rel);
        // Motion has decayed: residual tip velocity far below the peak
        // of the startup transient (mass damping rate alpha/2 = 50 1/s
        // over 0.5 s; the absolute floor is roundoff of the a0-scaled
        // direct solves, so compare against the observed peak).
        let peak_v = sol
            .points
            .iter()
            .map(|p| p.velocities[ti].abs())
            .fold(0.0f64, f64::max);
        let v_end = sol.points.last().unwrap().velocities[ti].abs();
        eprintln!("static limit: peak |v|={peak_v:.3e}, residual |v|={v_end:.3e}");
        assert!(
            v_end < peak_v / 50.0,
            "residual tip velocity {v_end:.3e} vs peak {peak_v:.3e}"
        );
    }

    #[test]
    fn rayleigh_damping_reduces_peak_response() {
        let (nx, ny, l, h) = (20, 4, 1.0, 0.1);
        let run = |alpha: f64| {
            let mut model = steel_cantilever(nx, ny, l, h, -500.0);
            model.steps[0].dynamics = Some(DynamicsSpec {
                dt: 2.0e-4,
                duration: 0.12,
                rayleigh_alpha: alpha,
                rayleigh_beta: 0.0,
            });
            let sol = solve_dynamic(&model).unwrap();
            let ti = tip_dof(&model, nx);
            let peak = sol
                .points
                .iter()
                .map(|p| p.displacements[ti].abs())
                .fold(0.0, f64::max);
            (peak, model)
        };
        let (undamped, model) = run(0.0);
        let (damped, _) = run(200.0);
        let stat = solve(&model).unwrap();
        let u_stat = stat.displacements[model
            .mesh
            .node_index((2 * (nx + 1) + nx + 1) as u32)
            .unwrap()][1]
            .abs();
        eprintln!("peaks: undamped {undamped:.6e}, alpha=200 {damped:.6e}, static {u_stat:.6e}");
        // Sudden load overshoots the static deflection (amplification ~2)...
        assert!(undamped > 1.5 * u_stat, "undamped peak {undamped:.4e}");
        // ...and mass damping eats a visible share of that peak.
        assert!(
            damped < 0.85 * undamped,
            "damped peak {damped:.4e} not below 85% of undamped {undamped:.4e}"
        );
        assert!(
            damped > u_stat,
            "damped peak {damped:.4e} fell under static {u_stat:.4e}"
        );
    }

    #[test]
    fn dynamics_guards() {
        let mut model = steel_cantilever(4, 2, 1.0, 0.1, -1.0);
        model.steps[0].dynamics = None;
        assert!(matches!(
            solve_dynamic(&model),
            Err(SolveError::InvalidModel(m)) if m.contains("*DYNAMIC")
        ));

        let mut model = steel_cantilever(4, 2, 1.0, 0.1, -1.0);
        model.steps[0].dynamics = Some(DynamicsSpec {
            dt: 0.0,
            duration: 1.0,
            rayleigh_alpha: 0.0,
            rayleigh_beta: 0.0,
        });
        assert!(matches!(
            solve_dynamic(&model),
            Err(SolveError::InvalidModel(m)) if m.contains("increment")
        ));

        let mut model = steel_cantilever(4, 2, 1.0, 0.1, -1.0);
        model.steps[0].dynamics = Some(DynamicsSpec {
            dt: 1.0,
            duration: 0.4,
            rayleigh_alpha: 0.0,
            rayleigh_beta: 0.0,
        });
        assert!(matches!(
            solve_dynamic(&model),
            Err(SolveError::InvalidModel(m)) if m.contains("at least one")
        ));

        let mut model = steel_cantilever(4, 2, 1.0, 0.1, -1.0);
        model.steps[0].boundary_conditions.clear();
        assert!(matches!(
            solve_dynamic(&model),
            Err(SolveError::NoConstraints)
        ));
    }
}
