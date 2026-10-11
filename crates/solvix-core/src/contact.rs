//! # 2D node-to-surface penalty contact (small sliding, optional
//! Coulomb friction)
//!
//! v1 formulation (issue #31), plugged into the Newton machinery of
//! [`crate::nl`] as an active-set strategy:
//!
//! - **Surfaces** are NODE lists ([`ContactPair`]). The master polyline
//!   is ordered by sorting its nodes along the surface's dominant axis
//!   (larger of the x/y bounding-box spans) — exact for straight or
//!   monotone master lines (floors, ramps), the documented v1 ceiling;
//!   genuinely curved masters need connectivity-based ordering later.
//!   Element-edge surfaces (`*SURFACE, TYPE=ELEMENT`) are not supported.
//! - **Gap**: each slave node projects onto every master SEGMENT
//!   (consecutive master-node pair); the projection parameter is
//!   clamped to `[0, 1]` (nodes beyond a segment end do not see it —
//!   small sliding). The segment normal is the left normal of the
//!   sorted polyline tangent, `n = (-ty, tx)`; the gap is
//!   `g = (x_s - p) . n` (negative = penetrating). The penetrating
//!   segment with the deepest gap wins.
//! - **Penalty**: nodal contact force `F = penalty * |g| * t * n`
//!   pushing the slave OUT of the master (t = model thickness), with
//!   master nodes receiving `N_i`-weighted opposite forces. Default
//!   stiffness [`Model::contact_penalty`] = 0 selects the documented
//!   heuristic `penalty = 100 * E_max` (2.1e13 Pa at steel scale):
//!   ~100x a typical adjacent element stiffness, giving ~1e-12 m
//!   penetration under gravity-scale loads without wrecking the
//!   tangent's conditioning.
//! - **Friction** ([`Model::friction_coeff`] = `Some(mu)`, issue #56):
//!   tangential penalty stick/slip on penetrating nodes. STICK: the
//!   tangential spring uses `k_t = penalty/1000` — deliberately below
//!   the normal penalty, which is calibrated for ~1e-12 m penetration;
//!   at the normal value the stick spring shear-lags the surface nodes
//!   (host in-plane stiffness ~E*t) and concentrates the tangential
//!   load at the patch ends far past their `mu*F_n` caps, parking the
//!   corners in a permanent slip state that Newton crawls through. At
//!   `penalty/1000` the stick slip stays ~1e-9 m (zero at contact
//!   scale). The slip s is the relative tangential displacement
//!   `s = (u_s - sum N_i u_mi) . t_hat` from the UNSTRESSED initial
//!   contact point — an EXPLICIT v1 law (total-slip anchor). SLIP when
//!   `|k_t s| > mu * F_n`; the force returns to the friction cone,
//!   `F_t = mu * F_n * sign(s)` opposing the slip direction. The
//!   implicit upgrade (documented, not built): track the slip
//!   incrementally from the point of first contact with a
//!   return-mapping algorithm and the consistent (asymmetric) slip
//!   tangent d|F_t|/du through F_n — the exact tangent is asymmetric
//!   and even its symmetrization is indefinite, unusable under the
//!   symmetric LDLT, so the v1 tangent keeps the STICK blocks
//!   `k_t c_i c_j (t outer t)` active for slipping nodes too (the
//!   standard symmetric-solver approximation). Consequences: explicit
//!   slip states floor the Newton residual at the cap-rounding level
//!   (~5e-6 relative on the bench ramp), and a sliding body finds no
//!   static equilibrium at all — the solve fails with `NotConverged`
//!   after its cutbacks, the reported "slides" signal.
//! - **Cold start**: a node is LINEARIZED as in contact when
//!   `g < activation_gap` (= 1e-9 * mesh diagonal) but only carries
//!   force for `g < 0`. The band exists so a block resting exactly on
//!   a floor (`g = 0 +- rounding noise` at the first iteration) has
//!   normal support from iteration one; without it the first Newton
//!   tangent is vertically singular and the increment is cut forever.
//!   Zero force in the band means no spurious adhesion.
//! - **Active set**: contact states are re-evaluated at every Newton
//!   iteration; convergence with a stable set is just the residual
//!   tolerance. Oscillation is bounded by the MAX 20 SET CHANGES per
//!   increment ceiling in [`crate::nl::solve_nl`] (beyond it the
//!   increment is cut — the documented ceiling before which a sticky
//!   set/Uzawa scheme would be needed).
//! - **Tangent**: penalty normal stiffness only, with the consistent
//!   slave/master coupling blocks `k t c_i c_j n n^T` (`c` from the
//!   projection shape functions; the friction stick blocks add the
//!   matching `k t c_i c_j t t^T`). The normal is held FIXED over an
//!   iteration (small sliding): the `g dn/du` rotation term is dropped,
//!   the standard small-slip node-to-segment linearization.
//! - Contact models ride [`crate::nl::solve_nl`]; plasticity + contact
//!   and thermomechanics + contact are rejected/ignored (v1 ceilings).
//!   A zero-load contact model returns the zero solution (initial
//!   overlap without load is user error).

use crate::{ContactForce, Model, SolveError};

/// Resolved contact problem: node indices, ordered master polyline and
/// the penalty stiffness. Built once per solve from [`Model::contact_pairs`].
pub(crate) struct ContactCtx {
    penalty: f64,
    /// Coulomb friction coefficient when [`crate::Model::friction_coeff`]
    /// is set; `None` keeps contact frictionless.
    friction: Option<f64>,
    /// Nodes linearized as in-contact when g < activation_gap; zero
    /// force unless g < 0 (see module docs).
    activation_gap: f64,
    /// (slave node indices, ordered master polyline node indices)
    pairs: Vec<(Vec<usize>, Vec<usize>)>,
    /// Total slave slots across pairs (forces index into this order).
    n_slave_slots: usize,
}

/// Contact contribution at one trial configuration.
pub(crate) struct ContactTerms {
    /// Symmetric tangent triplets (global dof rows/cols).
    pub(crate) trips: Vec<(usize, usize, f64)>,
    /// Internal-force gradient rows/values: `f_int[row] += value`.
    pub(crate) f: Vec<(usize, f64)>,
    /// Physical contact forces on penetrating slave nodes.
    pub forces: Vec<ContactForce>,
    /// Active-set flags per slave slot (g < activation_gap).
    pub active: Vec<bool>,
}

/// Node ids -> deduplicated node indices, erroring on unknown ids.
fn node_indices(model: &Model, ids: &[u32], what: &str) -> Result<Vec<usize>, SolveError> {
    let mut idx = Vec::with_capacity(ids.len());
    for &id in ids {
        let i = model.mesh.node_index(id).ok_or_else(|| {
            SolveError::InvalidModel(format!("contact pair {what} references unknown node {id}"))
        })?;
        if !idx.contains(&i) {
            idx.push(i);
        }
    }
    Ok(idx)
}

/// Order a master node list into a polyline: dedup, then sort along the
/// dominant bounding-box axis (v1 heuristic, see module docs).
fn sorted_polyline(model: &Model, mut idx: Vec<usize>) -> Vec<usize> {
    if idx.len() < 2 {
        return idx;
    }
    let (mut x0, mut x1, mut y0, mut y1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &i in &idx {
        let c = model.mesh.nodes[i].coords;
        x0 = x0.min(c[0]);
        x1 = x1.max(c[0]);
        y0 = y0.min(c[1]);
        y1 = y1.max(c[1]);
    }
    let by_x = x1 - x0 >= y1 - y0;
    idx.sort_unstable_by(|&a, &b| {
        let (ca, cb) = (model.mesh.nodes[a].coords, model.mesh.nodes[b].coords);
        if by_x {
            ca[0].total_cmp(&cb[0])
        } else {
            ca[1].total_cmp(&cb[1])
        }
    });
    idx
}

/// Build the contact context: resolve/validate ids, order the master
/// polylines, resolve the penalty stiffness (auto = 100 * E_max).
pub(crate) fn build_ctx(model: &Model) -> Result<ContactCtx, SolveError> {
    let penalty = if model.contact_penalty.is_finite() && model.contact_penalty > 0.0 {
        model.contact_penalty
    } else {
        100.0
            * model
                .materials
                .iter()
                .map(|m| m.youngs_modulus)
                .fold(f64::NEG_INFINITY, f64::max)
    };
    if !(penalty.is_finite() && penalty > 0.0) {
        return Err(SolveError::InvalidModel(
            "contact penalty must be positive (materials need a Young's modulus)".into(),
        ));
    }
    let friction = model.friction_coeff;
    if let Some(mu) = friction {
        if !mu.is_finite() || mu < 0.0 {
            return Err(SolveError::InvalidModel(
                "friction coefficient must be finite and non-negative".into(),
            ));
        }
    }
    let mut pairs = Vec::with_capacity(model.contact_pairs.len());
    for (pi, p) in model.contact_pairs.iter().enumerate() {
        let slave = node_indices(model, &p.slave, &format!("{} slave", pi + 1))?;
        let master = sorted_polyline(
            model,
            node_indices(model, &p.master, &format!("{} master", pi + 1))?,
        );
        if master.len() < 2 {
            return Err(SolveError::InvalidModel(format!(
                "contact pair {} master surface needs at least 2 nodes",
                pi + 1
            )));
        }
        if slave.is_empty() {
            return Err(SolveError::InvalidModel(format!(
                "contact pair {} slave surface is empty",
                pi + 1
            )));
        }
        pairs.push((slave, master));
    }
    let n_slave_slots = pairs.iter().map(|(s, _)| s.len()).sum();
    let (mut x0, mut x1, mut y0, mut y1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for n in &model.mesh.nodes {
        x0 = x0.min(n.coords[0]);
        x1 = x1.max(n.coords[0]);
        y0 = y0.min(n.coords[1]);
        y1 = y1.max(n.coords[1]);
    }
    let diag = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
    Ok(ContactCtx {
        penalty,
        friction,
        activation_gap: 1e-9 * diag,
        pairs,
        n_slave_slots,
    })
}

/// Evaluate all contact pairs at the total displacement `u` (two
/// entries per node). See the module docs for the gap/penalty/
/// active-set formulation.
///
/// The gap is evaluated in a cancellation-stabilized association:
/// `g = g_geo + g_disp` with the geometry part built from original
/// coordinates and the displacement part from `u` alone. Evaluating
/// `(x_s - p) . n` directly cancels 0.1-scale coordinates down to the
/// ~1e-12 penetration, so 1-ulp noise is ~1e-5 of the gap and the
/// 1e13-scale penalty amplifies it into a ~1e-4 N residual floor that
/// Newton can never resolve. With the split the coordinate noise is a
/// per-pair CONSTANT (absorbed into the converged penetration) and the
/// u-dependence stays smooth to machine precision.
pub(crate) fn contact_terms(ctx: &ContactCtx, model: &Model, u: &[f64]) -> ContactTerms {
    let mut terms = ContactTerms {
        trips: Vec::new(),
        f: Vec::new(),
        forces: Vec::new(),
        active: vec![false; ctx.n_slave_slots],
    };
    let k = ctx.penalty * model.thickness;
    // Tangential stick stiffness: a fixed 1e-3 fraction of the normal
    // penalty. The normal penalty (100*E) is calibrated for penetration
    // (~1e-12 m); the same value as k_t over-constrains the tangential
    // direction — the host element's in-plane stiffness is ~E*t, so a
    // 100*E stick spring shear-lags the surface nodes and concentrates
    // the tangential load at the patch ends far past their mu*F_n caps
    // (measured: a 2.1e13 k_t parks the corner nodes in a permanent
    // slip state that Newton crawls through). At penalty/1000 the
    // stick slip stays ~1e-9 m (still "zero" at contact scale) and the
    // fair-share distribution sticks cleanly.
    let kt = k * 1e-3;
    let xy0 = |i: usize| {
        let c = model.mesh.nodes[i].coords;
        [c[0] + u[2 * i], c[1] + u[2 * i + 1]]
    };
    let mut slot = 0usize;
    for (slaves, master) in &ctx.pairs {
        'slave: for &s in slaves {
            let xs = xy0(s);
            let xs0 = [model.mesh.nodes[s].coords[0], model.mesh.nodes[s].coords[1]];
            // Deepest penetrating segment wins.
            let mut hit: Option<(f64, usize, usize, f64, [f64; 2])> = None;
            for w in master.windows(2) {
                let (m0, m1) = (w[0], w[1]);
                let a = xy0(m0);
                let b = xy0(m1);
                let t = [b[0] - a[0], b[1] - a[1]];
                let l2 = t[0] * t[0] + t[1] * t[1];
                if l2 == 0.0 {
                    continue; // degenerate segment (duplicate coordinate)
                }
                let l = l2.sqrt();
                let d0 = [xs[0] - a[0], xs[1] - a[1]];
                let mut sp = (d0[0] * t[0] + d0[1] * t[1]) / l2;
                // A slave sitting exactly at a polyline END (stacked
                // blocks sharing corner coordinates) must not flicker
                // in/out of the segment under Poisson-scale drift:
                // accept a small overshoot band, clamp for the shape
                // functions. Truly-beyond-end nodes still skip.
                let sp_tol = 1e-6 * l;
                if sp < -sp_tol || sp > 1.0 + sp_tol {
                    continue; // beyond the segment end (small sliding)
                }
                sp = sp.clamp(0.0, 1.0);
                let n = [-t[1] / l, t[0] / l];
                // Reassociated gap (see the doc comment).
                let (a0, b0) = (
                    [
                        model.mesh.nodes[m0].coords[0],
                        model.mesh.nodes[m0].coords[1],
                    ],
                    [
                        model.mesh.nodes[m1].coords[0],
                        model.mesh.nodes[m1].coords[1],
                    ],
                );
                let geo = [
                    xs0[0] - a0[0] - sp * (b0[0] - a0[0]),
                    xs0[1] - a0[1] - sp * (b0[1] - a0[1]),
                ];
                let dis = [
                    u[2 * s] - u[2 * m0] - sp * (u[2 * m1] - u[2 * m0]),
                    u[2 * s + 1] - u[2 * m0 + 1] - sp * (u[2 * m1 + 1] - u[2 * m0 + 1]),
                ];
                let g = (geo[0] + dis[0]) * n[0] + (geo[1] + dis[1]) * n[1];
                if g < ctx.activation_gap && g < hit.map_or(f64::INFINITY, |h| h.0) {
                    hit = Some((g, m0, m1, sp, n));
                }
            }
            let Some((g, m0, m1, sp, n)) = hit else {
                slot += 1;
                continue 'slave;
            };
            terms.active[slot] = true;
            slot += 1;
            // Internal-force gradient f_int += dPi/du (see module docs
            // for the sign walk-through): slave k*g*n, masters -N_i*k*g*n.
            // PENETRATING nodes only — the activation band is tangent-
            // only (zero force, no adhesion).
            let (w0, w1) = (1.0 - sp, sp);
            let tang = [n[1], -n[0]];
            if g < 0.0 {
                for (row, fx, fy) in [
                    (2 * s, k * g * n[0], k * g * n[1]),
                    (2 * m0, -w0 * k * g * n[0], -w0 * k * g * n[1]),
                    (2 * m1, -w1 * k * g * n[0], -w1 * k * g * n[1]),
                ] {
                    terms.f.push((row, fx));
                    terms.f.push((row + 1, fy));
                }
                // Physical force ON the slave (pushes it out of master).
                let (mut cfx, mut cfy) = (-k * g * n[0], -k * g * n[1]);
                // Coulomb friction (explicit v1, see module docs): stick
                // spring anchored at the initial contact point, slip
                // returns |F_t| to mu * F_n. f_int mirrors the normal
                // pattern: slave +F_t*t_hat, masters -N_i*F_t*t_hat.
                if let Some(mu) = ctx.friction {
                    let slip = (u[2 * s] - w0 * u[2 * m0] - w1 * u[2 * m1]) * tang[0]
                        + (u[2 * s + 1] - w0 * u[2 * m0 + 1] - w1 * u[2 * m1 + 1]) * tang[1];
                    let f_n = -k * g;
                    let sticking = (kt * slip).abs() <= mu * f_n;
                    let ft = if sticking {
                        kt * slip // stick
                    } else {
                        mu * f_n * slip.signum() // slip: back on the cone
                    };
                    for (row, fx, fy) in [
                        (2 * s, ft * tang[0], ft * tang[1]),
                        (2 * m0, -w0 * ft * tang[0], -w0 * ft * tang[1]),
                        (2 * m1, -w1 * ft * tang[0], -w1 * ft * tang[1]),
                    ] {
                        terms.f.push((row, fx));
                        terms.f.push((row + 1, fy));
                    }
                    cfx -= ft * tang[0];
                    cfy -= ft * tang[1];
                }
                terms.forces.push(ContactForce {
                    node_id: model.mesh.nodes[s].id,
                    fx: cfx,
                    fy: cfy,
                    penetration: -g,
                });
            }
            // Consistent tangent k*t*c_i*c_j*(n outer n), fixed normal —
            // active band nodes included (that is the band's purpose:
            // normal support from the first iteration). Friction adds
            // the matching (t outer t) blocks with the same penalty
            // stiffness for EVERY friction-active node, slipping ones
            // included: the exact slip tangent is asymmetric
            // (d|F_t|/du flows through F_n and the penetration), whose
            // symmetrization is indefinite — unusable under the
            // symmetric LDLT. Keeping the stick tangent in slip is the
            // standard symmetric-solver approximation (documented v1
            // ceiling alongside the implicit return-mapping upgrade).
            let blocks = [
                (2 * s, 2 * s, 1.0),
                (2 * s, 2 * m0, -w0),
                (2 * s, 2 * m1, -w1),
                (2 * m0, 2 * s, -w0),
                (2 * m1, 2 * s, -w1),
                (2 * m0, 2 * m0, w0 * w0),
                (2 * m0, 2 * m1, w0 * w1),
                (2 * m1, 2 * m0, w1 * w0),
                (2 * m1, 2 * m1, w1 * w1),
            ];
            let dirs: [[f64; 2]; 2] = if ctx.friction.is_some() {
                [n, tang]
            } else {
                [n, [0.0, 0.0]]
            };
            for (r, c, cf) in blocks {
                for (stiff, dir) in [(k, dirs[0]), (kt, dirs[1])] {
                    for (i, ni) in dir.iter().enumerate() {
                        for (j, nj) in dir.iter().enumerate() {
                            terms.trips.push((r + i, c + j, stiff * cf * ni * nj));
                        }
                    }
                }
            }
        }
    }
    terms
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::tests::quad_mesh;
    use crate::{
        solve, solve_nl, BoundaryCondition, ConcentratedLoad, ContactPair, ElementSet, Material,
        Model, NlOptions, Step,
    };

    const E: f64 = 210e9;
    const RHO: f64 = 7850.0;
    const G: f64 = 9.81;

    fn steel_model(mesh: crate::Mesh, thickness: f64) -> Model {
        let n_elem = mesh.elements.len() as u32;
        Model {
            mesh,
            thickness,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: E,
                poisson_ratio: 0.3,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=n_elem).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            ..Default::default()
        }
    }

    /// Extra FREE master line nodes under y = 0 (floor ids 1000+); the
    /// caller fixes them in the step's boundary conditions.
    fn add_floor(model: &mut Model, xs: &[f64], y: f64) -> Vec<u32> {
        let mut ids = Vec::new();
        for (i, &x) in xs.iter().enumerate() {
            let id = 1000 + i as u32;
            model.mesh.nodes.push(crate::Node {
                id,
                coords: [x, y, 0.0],
            });
            ids.push(id);
        }
        ids
    }

    fn fixed(ids: &[u32]) -> BoundaryCondition {
        BoundaryCondition {
            node_ids: ids.to_vec(),
            dofs: vec![1, 2],
        }
    }

    /// Bottom-edge node ids of a quad_mesh block (y = 0 row).
    fn bottom_edge(_model: &Model, nx: usize) -> Vec<u32> {
        (1..=(nx + 1) as u32).collect()
    }

    #[test]
    fn block_on_rigid_floor_pressure_is_rho_g_h() {
        // (a) CPS4 block (tall, aspect 4: a 1:1 block's own quadratic-
        // in-y gravity field costs ~5% discretization error — the tall
        // geometry drops it below 1%) resting on a fixed floor line at
        // y=0, gravity through a step. ux pinned at ONE top node (kills
        // the free frictionless sliding mode; its Poisson disturbance
        // decays well before the bottom edge).
        let (nx, ny) = (6usize, 24);
        let (w, h) = (0.06, 0.24);
        let mut model = steel_model(quad_mesh(nx, ny, w, h), 1.0);
        let floor = add_floor(&mut model, &[-0.01, 0.02, 0.05, 0.08, 0.11], 0.0);
        let top_mid = (ny * (nx + 1) + nx / 2 + 1) as u32;
        model.contact_pairs = vec![ContactPair {
            slave: bottom_edge(&model, nx),
            master: floor.clone(),
        }];
        model.steps = vec![Step {
            name: "gravity".into(),
            boundary_conditions: vec![
                fixed(&floor),
                BoundaryCondition {
                    node_ids: vec![top_mid],
                    dofs: vec![1],
                },
            ],
            gravity: Some((G, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        // solve() routes through the NL machinery because of the pair.
        let sol = solve(&model).unwrap();
        assert!(!sol.contact_forces.is_empty());
        let p0 = RHO * G * h; // 18.5 kPa at h = 0.24
        let t = model.thickness;
        let dx = w / nx as f64;
        // Uniform interior pressure within 5%: tributary dx, corners dx/2.
        let id_of = |i: usize| (i + 1) as u32;
        let pressure = |i: usize| -> f64 {
            let f = sol
                .contact_forces
                .iter()
                .find(|c| c.node_id == id_of(i))
                .unwrap_or_else(|| panic!("bottom node {} not in contact", id_of(i)))
                .fy;
            let trib = if i == 0 || i == nx { dx / 2.0 } else { dx };
            f / (t * trib)
        };
        let mut worst = 0.0f64;
        for i in 0..=nx {
            let dev = (pressure(i) - p0).abs() / p0;
            eprintln!(
                "floor: node {} pressure {:.1} Pa ({:+.2}%)",
                id_of(i),
                pressure(i),
                100.0 * (pressure(i) / p0 - 1.0)
            );
            worst = worst.max(dev);
        }
        eprintln!("floor: worst interior pressure dev {:.3}%", 100.0 * worst);
        assert!(worst < 0.05, "pressure deviation {worst}");
        // Penetration < h * 1e-6 everywhere (expected ~4e-12 m).
        let max_pen = sol
            .contact_forces
            .iter()
            .map(|c| c.penetration)
            .fold(0.0, f64::max);
        eprintln!("floor: max penetration {max_pen:.3e} m");
        assert!(max_pen < h * 1e-6, "penetration {max_pen}");
        // Global equilibrium: contact holds the whole weight, forces +y.
        let fy: f64 = sol.contact_forces.iter().map(|c| c.fy).sum();
        let weight = RHO * w * h * t * G;
        assert!(
            (fy - weight).abs() / weight < 1e-3,
            "sum fy {fy} vs {weight}"
        );
        assert!(sol.contact_forces.iter().all(|c| c.fy > 0.0));
    }

    #[test]
    fn two_block_stack_load_path_through_interface() {
        // (b) two blocks stacked; master = elastic top edge of the bottom
        // block, slave = bottom edge of the top block. Interface must
        // carry the top block's weight + applied load; the base carries
        // everything.
        let (nx, ny1, ny2) = (8usize, 4, 6);
        let (w, h1, h2) = (0.08, 0.04, 0.06);
        let mut mesh = quad_mesh(nx, ny1, w, h1);
        let off_n = mesh.nodes.len() as u32;
        let off_e = mesh.elements.len() as u32;
        for j in 0..=ny2 {
            for i in 0..=nx {
                mesh.nodes.push(crate::Node {
                    id: off_n + (j * (nx + 1) + i + 1) as u32,
                    coords: [
                        w * i as f64 / nx as f64,
                        h1 + h2 * j as f64 / ny2 as f64,
                        0.0,
                    ],
                });
            }
        }
        for j in 0..ny2 {
            for i in 0..nx {
                let nid = |ii: usize, jj: usize| off_n + (jj * (nx + 1) + ii + 1) as u32;
                mesh.elements.push(crate::Element {
                    id: off_e + (j * nx + i + 1) as u32,
                    kind: crate::ElemKind::Cps4,
                    connectivity: vec![nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)],
                });
            }
        }
        let mut model = steel_model(mesh, 1.0);
        let base: Vec<u32> = (1..=(nx + 1) as u32).collect();
        let base_mid = (nx / 2 + 1) as u32;
        let top_of_bottom: Vec<u32> = (0..=nx).map(|i| (ny1 * (nx + 1) + i + 1) as u32).collect();
        let bottom_of_top: Vec<u32> = (0..=nx).map(|i| off_n + i as u32 + 1).collect();
        let top_of_top: Vec<u32> = (0..=nx)
            .map(|i| off_n + (ny2 * (nx + 1) + i + 1) as u32)
            .collect();
        model.contact_pairs = vec![ContactPair {
            slave: bottom_of_top,
            master: top_of_bottom,
        }];
        // Uniform top load: consistent nodal split of a total downward
        // P over the top edge (corner nodes carry half a segment each).
        let p_total = -2000.0; // N (downward)
                               // Top-block lateral pin (frictionless contact carries no shear:
                               // without it the top block has a free tangential rigid mode).
        let top_pin = off_n + (ny2 * (nx + 1) + nx / 2 + 1) as u32;
        model.steps = vec![Step {
            name: "stack".into(),
            boundary_conditions: vec![
                BoundaryCondition {
                    node_ids: base.clone(),
                    dofs: vec![2],
                },
                BoundaryCondition {
                    node_ids: vec![base_mid],
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: vec![top_pin],
                    dofs: vec![1],
                },
            ],
            cloads: top_of_top
                .iter()
                .enumerate()
                .map(|(i, &nid)| ConcentratedLoad {
                    node_id: nid,
                    dof: 2,
                    value: if i == 0 || i == nx {
                        p_total / 2.0 / nx as f64
                    } else {
                        p_total / nx as f64
                    },
                })
                .collect(),
            gravity: Some((G, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        let sol = solve(&model).unwrap();
        // Interface carries top block weight + full applied load.
        let w_top = RHO * w * h2 * model.thickness * G;
        let fy: f64 = sol.contact_forces.iter().map(|c| c.fy).sum();
        let expect_iface = w_top + p_total.abs();
        assert!(
            (fy - expect_iface).abs() / expect_iface < 0.01,
            "interface fy {fy} vs {expect_iface}"
        );
        // Base reactions carry everything (weight of both + load).
        let ry: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 2)
            .map(|r| r.value)
            .sum();
        let w_all = RHO * w * (h1 + h2) * model.thickness * G;
        let expect_base = w_all + p_total.abs();
        assert!(
            (ry - expect_base).abs() / expect_base < 0.01,
            "base ry {ry} vs {expect_base}"
        );
        // Interface pressure uniform within 5% (interior nodes).
        let p_iface = (w_top + p_total.abs()) / (w * model.thickness);
        let dx = w / nx as f64;
        let mut worst = 0.0f64;
        for i in 1..nx {
            let nid = off_n + i as u32 + 1;
            let f = sol
                .contact_forces
                .iter()
                .find(|c| c.node_id == nid)
                .unwrap()
                .fy;
            let dev = (f / (model.thickness * dx) - p_iface).abs() / p_iface;
            worst = worst.max(dev);
        }
        eprintln!(
            "stack: worst interior interface pressure dev {:.3}%",
            100.0 * worst
        );
        assert!(worst < 0.05, "interface pressure deviation {worst}");
    }

    #[test]
    fn gap_means_zero_contact_forces() {
        // (c) self-supported tall column with a floor 0.02 BELOW: the
        // pair is defined but never activates; no contact forces are
        // reported and the column still solves to the gravity
        // compression rho g h^2 / (2E) (tall so the single lateral pin
        // disturbs only its Saint-Venant neighborhood).
        let (nx, ny) = (4usize, 32);
        let (w, h) = (0.08, 0.32);
        let mut model = steel_model(quad_mesh(nx, ny, w, h), 1.0);
        let floor = add_floor(&mut model, &[-0.01, 0.04, 0.09], -0.02);
        let base = bottom_edge(&model, nx);
        model.contact_pairs = vec![ContactPair {
            slave: base.clone(),
            master: floor.clone(),
        }];
        model.steps = vec![Step {
            name: "gap".into(),
            boundary_conditions: vec![
                fixed(&floor),
                BoundaryCondition {
                    node_ids: base,
                    dofs: vec![2],
                },
                BoundaryCondition {
                    node_ids: vec![1],
                    dofs: vec![1],
                },
            ],
            gravity: Some((G, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        let sol = solve(&model).unwrap();
        assert!(sol.contact_forces.is_empty(), "{:?}", sol.contact_forces);
        let top_left = (ny * (nx + 1) + 1) as u32;
        let i = model.mesh.node_index(top_left).unwrap();
        let expect = -RHO * G * h * h / (2.0 * E);
        let got = sol.displacements[i][1];
        eprintln!("gap: top uy {got:.6e} vs {expect:.6e}");
        assert!(
            (got - expect).abs() / expect.abs() < 0.01,
            "top uy {got:.6e} vs {expect:.6e}"
        );
    }

    #[test]
    fn tilted_ramp_45deg_normal_and_frictionless_resultant() {
        // (d) block on a fixed 45-degree ramp, full gravity. The block
        // hangs tangentially off ONE pinned node; frictionless contact
        // must carry exactly the normal component along the tilted
        // normal n = (-s, c): sum F = m g cos45 * n (components +-W/2),
        // zero tangential resultant, small normal penetration.
        let (nx, ny) = (8usize, 8);
        let (w, h) = (0.1, 0.1);
        let mesh = quad_mesh(nx, ny, w, h);
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let c = s;
        let t_hat = [c, s];
        let n_hat = [-s, c];
        let c0 = [0.0, 0.02];
        // Rotate: local (x, y) -> c0 + x*t + y*n.
        let mut nodes = Vec::with_capacity(mesh.nodes.len() + 5);
        for nd in &mesh.nodes {
            let (lx, ly) = (nd.coords[0], nd.coords[1]);
            nodes.push(crate::Node {
                id: nd.id,
                coords: [
                    c0[0] + lx * t_hat[0] + ly * n_hat[0],
                    c0[1] + lx * t_hat[1] + ly * n_hat[1],
                    0.0,
                ],
            });
        }
        // Ramp master line through c0 along t, slightly wider.
        let ramp_pos = [-0.02, 0.03, 0.06, 0.09, 0.12];
        let mut ramp_ids = Vec::new();
        for (i, &u) in ramp_pos.iter().enumerate() {
            let id = 1000 + i as u32;
            nodes.push(crate::Node {
                id,
                coords: [c0[0] + u * t_hat[0], c0[1] + u * t_hat[1], 0.0],
            });
            ramp_ids.push(id);
        }
        let mut model = steel_model(
            crate::Mesh {
                nodes,
                elements: mesh.elements,
            },
            1.0,
        );
        let slave: Vec<u32> = (1..=(nx + 1) as u32).collect(); // local y=0 edge
        let anchor = (ny / 2 * (nx + 1) + nx / 2 + 1) as u32; // top-face center
        model.contact_pairs = vec![ContactPair {
            slave,
            master: ramp_ids.clone(),
        }];
        model.steps = vec![Step {
            name: "ramp".into(),
            boundary_conditions: vec![
                fixed(&ramp_ids),
                BoundaryCondition {
                    node_ids: vec![anchor],
                    dofs: vec![1, 2],
                },
            ],
            gravity: Some((G, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        let sol = solve(&model).unwrap();
        assert!(!sol.contact_forces.is_empty());
        let w_block = RHO * w * h * model.thickness * G;
        let fx: f64 = sol.contact_forces.iter().map(|cf| cf.fx).sum();
        let fy: f64 = sol.contact_forces.iter().map(|cf| cf.fy).sum();
        let n_comp = fx * n_hat[0] + fy * n_hat[1];
        let t_comp = fx * t_hat[0] + fy * t_hat[1];
        eprintln!(
            "ramp: sum fx {fx:.2} fy {fy:.2}, normal {n_comp:.2} (of m g cos45 = {:.1}), tangential {t_comp:.3e}",
            w_block * c
        );
        // Every contact force points along the tilted normal (-45deg
        // ramp: |fx| == fy within rounding), never axis-aligned.
        for cf in &sol.contact_forces {
            assert!(
                (cf.fx + cf.fy).abs() < 1e-9 * cf.fy,
                "force not along n: {cf:?}"
            );
            assert!(cf.fy > 0.0);
        }
        // Frictionless: the contact resultant has ZERO tangential
        // component; the block slides freely (the one pinned anchor
        // carries whatever tangential share the statically-
        // indeterminate split gives it, so only the direction is exact).
        assert!(t_comp.abs() < 1e-6 * n_comp.abs(), "tangential {t_comp}");
        // Contact is the dominant support: at least half the normal
        // weight component rides the interface.
        assert!(
            n_comp > 0.5 * w_block * c,
            "normal resultant {n_comp} vs m g cos45 {}",
            w_block * c
        );
        // Small normal penetration (measured along n).
        let max_pen = sol
            .contact_forces
            .iter()
            .map(|cf| cf.penetration)
            .fold(0.0, f64::max);
        eprintln!("ramp: max normal penetration {max_pen:.3e} m");
        assert!(max_pen < h * 1e-6, "penetration {max_pen}");
    }

    #[test]
    fn active_set_stabilizes_in_few_newton_iterations() {
        // The floor test's Newton history: every increment converges
        // within the normal budget (contact states settle immediately).
        let (nx, ny) = (6usize, 6);
        let (w, h) = (0.06, 0.06);
        let mut model = steel_model(quad_mesh(nx, ny, w, h), 1.0);
        let floor = add_floor(&mut model, &[-0.005, 0.03, 0.065], 0.0);
        let top_mid = (ny * (nx + 1) + nx / 2 + 1) as u32;
        model.contact_pairs = vec![ContactPair {
            slave: bottom_edge(&model, nx),
            master: floor.clone(),
        }];
        model.steps = vec![Step {
            name: "gravity".into(),
            boundary_conditions: vec![
                fixed(&floor),
                BoundaryCondition {
                    node_ids: vec![top_mid],
                    dofs: vec![1],
                },
            ],
            gravity: Some((G, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        let frames = solve_nl(&model, &NlOptions::default()).unwrap();
        assert_eq!(frames.len(), 20);
        for f in &frames {
            assert!(
                f.solution.iterations <= 8,
                "increment at lambda={:.3} took {} iterations",
                f.load_fraction,
                f.solution.iterations
            );
        }
    }

    #[test]
    fn friction_block_on_30deg_ramp_sticks_or_slides() {
        // Block on a fixed 30-degree ramp under gravity. With
        // mu = 0.8 > tan(30deg) = 0.577 the block STICKS: static
        // equilibrium with (penetration-scale) zero displacement and
        // contact resultant = W*(cos30*n + sin30*t). With mu = 0.4 the
        // friction cone cannot hold the tangential drive (0.4*cos30 =
        // 0.346 < sin30 = 0.5) — the block slides, no static
        // equilibrium exists, and the solve reports NotConverged.
        // Flat slab (w/h = 8): the friction/gravity couple tilts the
        // patch and unloads the downhill corner; a wide flat block
        // keeps every node's mu*F_n above its fair tangential share.
        let (nx, ny) = (12usize, 2);
        let (w, h) = (0.24, 0.03);
        let mesh = quad_mesh(nx, ny, w, h);
        let theta = std::f64::consts::PI / 6.0;
        let (c, s) = (theta.cos(), theta.sin());
        let t_hat = [c, s];
        let n_hat = [-s, c];
        let c0 = [0.0, 0.02];
        let mut nodes = Vec::with_capacity(mesh.nodes.len() + 5);
        for nd in &mesh.nodes {
            let (lx, ly) = (nd.coords[0], nd.coords[1]);
            nodes.push(crate::Node {
                id: nd.id,
                coords: [
                    c0[0] + lx * t_hat[0] + ly * n_hat[0],
                    c0[1] + lx * t_hat[1] + ly * n_hat[1],
                    0.0,
                ],
            });
        }
        let ramp_pos = [-0.04, 0.08, 0.14, 0.2, 0.28];
        let mut ramp_ids = Vec::new();
        for (i, &u) in ramp_pos.iter().enumerate() {
            let id = 1000 + i as u32;
            nodes.push(crate::Node {
                id,
                coords: [c0[0] + u * t_hat[0], c0[1] + u * t_hat[1], 0.0],
            });
            ramp_ids.push(id);
        }
        let mut model = steel_model(
            crate::Mesh {
                nodes,
                elements: mesh.elements,
            },
            1.0,
        );
        let slave: Vec<u32> = (1..=(nx + 1) as u32).collect(); // local y=0 edge
        model.contact_pairs = vec![ContactPair {
            slave,
            master: ramp_ids.clone(),
        }];
        model.steps = vec![Step {
            name: "ramp".into(),
            boundary_conditions: vec![fixed(&ramp_ids)],
            gravity: Some((G, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        let w_block = RHO * w * h * model.thickness * G;

        model.friction_coeff = Some(0.8);
        // Explicit slip states floor the residual at the cap-rounding
        // level (~5e-6 here, one corner node riding its cone): solve at
        // 1e-5 — still far tighter than the physics needs.
        let frames = solve_nl(
            &model,
            &NlOptions {
                tolerance: 1e-5,
                ..NlOptions::default()
            },
        )
        .unwrap();
        let sol = frames.last().unwrap().solution.clone();
        let max_disp = sol
            .displacements
            .iter()
            .map(|d| d[0].abs().max(d[1].abs()))
            .fold(0.0f64, f64::max);
        eprintln!("friction mu=0.8: max |u| = {max_disp:.3e} m");
        // Stick-slip against k_t = penalty/1000 sets the scale (~1e-9 m);
        // 1e-8 m is 4 orders below the block's smallest feature.
        assert!(max_disp < 1e-8, "stuck block moved {max_disp:.3e} m");
        let fx: f64 = sol.contact_forces.iter().map(|cf| cf.fx).sum();
        let fy: f64 = sol.contact_forces.iter().map(|cf| cf.fy).sum();
        let n_comp = fx * n_hat[0] + fy * n_hat[1];
        let t_comp = fx * t_hat[0] + fy * t_hat[1];
        eprintln!(
            "friction mu=0.8: normal {n_comp:.2} (want {:.2}), tangential {t_comp:.2} (want {:.2} = W sin30)",
            w_block * c,
            w_block * s
        );
        assert!((n_comp - w_block * c).abs() / (w_block * c) < 1e-3);
        assert!((t_comp - w_block * s).abs() / (w_block * s) < 1e-3);

        model.friction_coeff = Some(0.4);
        match solve(&model) {
            Err(SolveError::NotConverged { iterations, .. }) => {
                eprintln!("friction mu=0.4: slides (not converged after {iterations} iterations)");
            }
            other => panic!("mu=0.4 on a 30deg ramp must slide, got {other:?}"),
        }
    }

    #[test]
    fn unknown_contact_node_and_short_master_are_invalid() {
        let mut model = steel_model(quad_mesh(2, 2, 1.0, 1.0), 1.0);
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![1, 2],
        }];
        let mut m2 = model.clone();
        m2.contact_pairs = vec![ContactPair {
            slave: vec![999],
            master: vec![1, 2],
        }];
        assert!(matches!(build_ctx(&m2), Err(SolveError::InvalidModel(_))));
        let mut m3 = model;
        m3.contact_pairs = vec![ContactPair {
            slave: vec![1, 2],
            master: vec![5], // one-node master
        }];
        assert!(matches!(build_ctx(&m3), Err(SolveError::InvalidModel(_))));
    }

    #[test]
    fn plasticity_plus_contact_is_rejected() {
        let mut model = steel_model(quad_mesh(2, 2, 0.1, 0.1), 1.0);
        model.plasticity = Some(crate::PlasticitySpec {
            yield_stress: 250e6,
            hardening_modulus: 0.0,
        });
        model.contact_pairs = vec![ContactPair {
            slave: vec![1, 2],
            master: vec![1000, 1001],
        }];
        model.mesh.nodes.push(crate::Node {
            id: 1000,
            coords: [0.0, -0.01, 0.0],
        });
        model.mesh.nodes.push(crate::Node {
            id: 1001,
            coords: [0.1, -0.01, 0.0],
        });
        model.boundary_conditions = vec![
            BoundaryCondition {
                node_ids: vec![1000, 1001],
                dofs: vec![1, 2],
            },
            BoundaryCondition {
                node_ids: vec![3],
                dofs: vec![1],
            },
        ];
        model.steps = vec![Step {
            name: "s".into(),
            boundary_conditions: model.boundary_conditions.clone(),
            gravity: Some((G, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        match solve(&model) {
            Err(SolveError::InvalidModel(m)) => {
                assert!(m.contains("plasticity"), "{m}")
            }
            other => panic!("expected InvalidModel, got {other:?}"),
        }
    }
}
