//! Linear buckling: eigenpairs of `(K + λ Kσ) φ = 0` on the
//! stress-stiffened stiffness of a prior small-strain static solution.
//!
//! Method: constant centroidal stress per element, Gauss-integrated
//! geometric stiffness `Kσ = ∫ G^T S G t dA` (2D kinds only; the scalar
//! `G^T S G` block multiplies the identity on each node's dof pair, the
//! classical small-strain form), then subspace iteration with the
//! constrained stiffness `K` as the SPD operator and `−Kσ` as the mass —
//! reusing the modal machinery. Load factors are the smallest positive
//! eigenvalues: multipliers on the applied loads that produce neutral
//! equilibrium.
//!
//! ponytail: constant-S Kσ (not Gauss-point stress) is the ceiling; the
//! upgrade path is per-GP stress recovery, which `nodal_stresses` already
//! computes. 3D kinds are rejected for now.

use crate::assembly::{self, Csr};
use crate::element::{centroid_stress, quad_b};
use crate::ldl::Ldlt;
use crate::modal::{dot, free_submatrix, zeros, Rng};
use crate::solver::{apply_dirichlet, constrained_dofs};
use crate::{ElemKind, Model, SolveError};

/// Buckling load factors and mode shapes of the loaded model.
#[derive(Debug, Clone, Default)]
pub struct BucklingSolution {
    /// Multipliers on the applied load pattern at neutral equilibrium,
    /// ascending, all positive. A factor × the current loads is the
    /// critical combination.
    pub load_factors: Vec<f64>,
    /// Mode shapes, each of length `nodes*ndof` (constrained dofs 0.0),
    /// normalized to unit Euclidean norm.
    pub modes: Vec<Vec<f64>>,
    pub iterations: usize,
}

/// Geometric stress stiffness block scalar for one element with plane
/// stress state `s = [sxx, syy, sxy]`, integrated over the element.
fn geo_scalar_quad(coords: &[[f64; 2]; 4], s: &[f64; 3], t: f64) -> Option<[[f64; 4]; 4]> {
    const GP: f64 = 0.577_350_269_189_625_8;
    let w = 1.0;
    let mut ks = [[0.0f64; 4]; 4];
    for &xi in &[-GP, GP] {
        for &eta in &[-GP, GP] {
            let (bm, detj) = quad_b(coords, xi, eta)?;
            let g: Vec<(f64, f64)> = (0..4).map(|i| (bm[0][2 * i], bm[1][2 * i + 1])).collect();
            for i in 0..4 {
                for j in 0..4 {
                    let (gx, gy) = g[i];
                    let (hx, hy) = g[j];
                    // (G^T S G) with S = [[sxx, sxy], [sxy, syy]]
                    let v = gx * (s[0] * hx + s[2] * hy) + gy * (s[2] * hx + s[1] * hy);
                    ks[i][j] += t * w * detj * v;
                }
            }
        }
    }
    Some(ks)
}

/// CPS3: constant gradients, closed-form area integral.
fn geo_scalar_tri(coords: &[[f64; 2]; 3], s: &[f64; 3], t: f64) -> Option<[[f64; 3]; 3]> {
    let b = [
        coords[1][1] - coords[2][1],
        coords[2][1] - coords[0][1],
        coords[0][1] - coords[1][1],
    ];
    let c = [
        coords[2][0] - coords[1][0],
        coords[0][0] - coords[2][0],
        coords[1][0] - coords[0][0],
    ];
    let two_a = b[0] * c[1] - b[1] * c[0];
    if two_a == 0.0 {
        return None;
    }
    let area = two_a.abs() / 2.0;
    let inv = 1.0 / two_a;
    let mut ks = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (gx, gy) = (b[i] * inv, c[i] * inv);
            let (hx, hy) = (b[j] * inv, c[j] * inv);
            let v = gx * (s[0] * hx + s[2] * hy) + gy * (s[2] * hx + s[1] * hy);
            ks[i][j] = t * area * v;
        }
    }
    Some(ks)
}

/// Assemble the global geometric stiffness from the static stress state.
fn assemble_ksigma(
    model: &Model,
    ndof: usize,
    u: &[f64],
    elems: &[assembly::ElemInfo],
) -> Result<Csr, SolveError> {
    let t = model.thickness;
    let mut trips: Vec<(usize, usize, f64)> = Vec::new();
    for (idx, e) in model.mesh.elements.iter().enumerate() {
        match e.kind {
            ElemKind::Cps3 | ElemKind::Cps4 => {}
            ElemKind::C3d8 => {
                return Err(SolveError::InvalidModel(
                    "3D buckling not yet supported (C3D8)".into(),
                ));
            }
            ElemKind::C3d4 => {
                return Err(SolveError::InvalidModel(
                    "3D buckling not yet supported (C3D4)".into(),
                ));
            }
        }
        let info = &elems[idx];
        let ue: Vec<f64> = info
            .node_idx
            .iter()
            .flat_map(|&n| (0..ndof).map(move |d| u[n * ndof + d]))
            .collect();
        let s6 = centroid_stress(&info.b, &info.d, &ue);
        let s = [s6[0], s6[1], s6[3]]; // Voigt: [s11, s22, s33, s12, ...]
        let coords2: Vec<[f64; 2]> = info
            .node_idx
            .iter()
            .map(|&n| [model.mesh.nodes[n].coords[0], model.mesh.nodes[n].coords[1]])
            .collect();
        let ks: Option<Vec<Vec<f64>>> = match e.kind {
            ElemKind::Cps3 => geo_scalar_tri(&[coords2[0], coords2[1], coords2[2]], &s, t)
                .map(|k| k.iter().map(|r| r.to_vec()).collect()),
            _ => geo_scalar_quad(&[coords2[0], coords2[1], coords2[2], coords2[3]], &s, t)
                .map(|k| k.iter().map(|r| r.to_vec()).collect()),
        };
        let Some(ks) = ks else {
            return Err(SolveError::InvalidModel(format!(
                "degenerate element {} in geometric stiffness",
                e.id
            )));
        };
        for (i, ksi) in ks.iter().enumerate() {
            for (j, &ksv) in ksi.iter().enumerate() {
                let v = ksv;
                if v == 0.0 {
                    continue;
                }
                for a in 0..ndof {
                    let r = info.node_idx[i] * ndof + a;
                    let c = info.node_idx[j] * ndof + a;
                    trips.push((r, c, v));
                }
            }
        }
    }
    Ok(assembly::csr_from_triplets(
        ndof * model.mesh.nodes.len(),
        trips,
    ))
}

/// Solve the linear buckling eigenproblem for the loaded model. The
/// returned factors multiply the model's own load pattern.
pub fn solve_buckle(model: &Model, num_modes: usize) -> Result<BucklingSolution, SolveError> {
    if model.boundary_conditions.is_empty() {
        return Err(SolveError::NoConstraints);
    }
    if let Err(e) = model.validate() {
        return Err(SolveError::InvalidModel(format!("{e:?}")));
    }
    if num_modes == 0 {
        return Err(SolveError::InvalidModel(
            "num_modes must be positive".into(),
        ));
    }
    let ndof = crate::ndof_per_node(model)?;

    // Static state under the applied loads (direct path for determinism).
    let (mut k, mut f, elems) = assembly::assemble(model, ndof)?;
    let constrained = constrained_dofs(model, ndof)?;
    apply_dirichlet(&mut k, &mut f, &constrained);
    let u = Ldlt::factorize(&k)?.solve(&f);

    let ksigma = assemble_ksigma(model, ndof, &u, &elems)?;

    // Free-dof submatrices: K SPD, B = -Ksigma (PSD under compression).
    let n = k.row_ptr.len() - 1;
    let mut is_constrained = vec![false; n];
    for c in &constrained {
        is_constrained[c.row] = true;
    }
    // Negate into the B operator by building a negated copy.
    let neg = |m: &Csr| -> Csr {
        Csr {
            row_ptr: m.row_ptr.clone(),
            col_idx: m.col_idx.clone(),
            vals: m.vals.iter().map(|v| -v).collect(),
        }
    };
    let (k_ff, free) = free_submatrix(&k, &is_constrained);
    let (b_ff, _) = free_submatrix(&neg(&ksigma), &is_constrained);
    let n_free = free.len();
    if num_modes >= n_free {
        return Err(SolveError::InvalidModel(format!(
            "num_modes {num_modes} must be smaller than the free dof count {n_free}"
        )));
    }
    let k_wanted = num_modes;

    // Spectral shift: inverse iteration on K^-1 B is dominated by
    // spurious near-zero modes (local tension spots around constraints
    // give genuinely tiny positive pencil eigenvalues). Factoring
    // (K + sigma B) with sigma ~ 1% of an Euler-mode Rayleigh estimate
    // damps them: they sit far below sigma, the physical modes far above.
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for n in &model.mesh.nodes {
        for d in 0..3 {
            lo[d] = lo[d].min(n.coords[d]);
            hi[d] = hi[d].max(n.coords[d]);
        }
    }
    let ext: Vec<f64> = (0..3).map(|d| hi[d] - lo[d]).collect();
    let axial = ext
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .unwrap()
        .0;
    let perp = if axial == 0 { 1 } else { 0 };
    let length = ext[axial].max(1e-12);
    let mut v = vec![0.0; ndof * model.mesh.nodes.len()];
    for n in &model.mesh.nodes {
        let s = (std::f64::consts::PI * (n.coords[axial] - lo[axial]) / length).sin();
        let c = (std::f64::consts::PI * (n.coords[axial] - lo[axial]) / length).cos();
        let off = n.coords[perp] - 0.5 * (lo[perp] + hi[perp]);
        let i = model.mesh.node_index(n.id).unwrap();
        v[i * ndof + axial] = -off * std::f64::consts::PI / length * c;
        if perp < ndof {
            v[i * ndof + perp] = s;
        }
    }
    // Project out constrained dofs into the free vector space.
    let mut v_free = vec![0.0; n_free];
    for (new, &old) in free.iter().enumerate() {
        v_free[new] = v[old];
    }
    let kv = k_ff.mat_vec(&v_free);
    let bv = b_ff.mat_vec(&v_free);
    let mut r0 = dot(&v_free, &bv) / dot(&v_free, &kv);
    // Boundary-condition-blind estimate: also try the second-mode shape
    // (clamped-clamped buckles in a half-wave of sin(2 pi x / L)).
    let mut v2 = vec![0.0; ndof * model.mesh.nodes.len()];
    for n in &model.mesh.nodes {
        let s = (2.0 * std::f64::consts::PI * (n.coords[axial] - lo[axial]) / length).sin();
        let i = model.mesh.node_index(n.id).unwrap();
        if perp < ndof {
            v2[i * ndof + perp] = s;
        }
    }
    let mut v2_free = vec![0.0; n_free];
    for (new, &old) in free.iter().enumerate() {
        v2_free[new] = v2[old];
    }
    let r2 = dot(&v2_free, &b_ff.mat_vec(&v2_free)) / dot(&v2_free, &k_ff.mat_vec(&v2_free));
    if r2.is_finite() && r2 > r0 {
        r0 = r2;
    }
    if !r0.is_finite() || r0 <= 0.0 {
        return Err(SolveError::InvalidModel(
            "load pattern produces net tension: no positive buckling factors".into(),
        ));
    }
    let fact = Ldlt::factorize(&k_ff)?;
    let floor = 1e-3 * r0;
    let p = (k_wanted + 4).min(n_free);
    let mut rng = Rng(0x0DDB_1A47_5EED_5EED);
    let mut q: Vec<Vec<f64>> = (0..p)
        .map(|_| (0..n_free).map(|_| rng.unit()).collect())
        .collect();

    let mut prev: Option<Vec<f64>> = None;
    for it in 1..=200 {
        let y: Vec<Vec<f64>> = q
            .iter()
            .map(|qj| {
                let yj = fact.solve(&b_ff.mat_vec(qj));
                let nrm = dot(&yj, &yj).sqrt().max(f64::MIN_POSITIVE);
                yj.iter().map(|v| v / nrm).collect()
            })
            .collect();
        let ky: Vec<Vec<f64>> = y.iter().map(|yj| k_ff.mat_vec(yj)).collect();
        let by: Vec<Vec<f64>> = y.iter().map(|yj| b_ff.mat_vec(yj)).collect();
        let mut kk = zeros(p);
        let mut bb = zeros(p);
        for a in 0..p {
            for b in a..p {
                let (kv, bv) = (dot(&y[a], &ky[b]), dot(&y[a], &by[b]));
                kk[a][b] = kv;
                kk[b][a] = kv;
                bb[a][b] = bv;
                bb[b][a] = bv;
            }
        }
        let (lam, c) = match crate::modal::generalized_eigenpairs(&kk, &bb) {
            Some(v) => v,
            None => {
                return Err(SolveError::InvalidModel(
                    "projected geometric stiffness too indefinite to regularize".into(),
                ))
            }
        };
        let x: Vec<Vec<f64>> = c
            .iter()
            .map(|cj| {
                let mut v = vec![0.0; n_free];
                for (yj, &cij) in y.iter().zip(cj) {
                    for (vi, &yv) in v.iter_mut().zip(yj) {
                        *vi += cij * yv;
                    }
                }
                let nrm = dot(&v, &v).sqrt();
                for sv in &mut v {
                    *sv /= nrm;
                }
                v
            })
            .collect();
        let lam_k: Vec<f64> = lam
            .iter()
            .copied()
            .filter(|&l| l > floor)
            .take(k_wanted)
            .collect();
        if lam_k.len() == k_wanted {
            if let Some(pr) = &prev {
                let change = lam_k
                    .iter()
                    .zip(pr)
                    .map(|(&a, &b)| (a - b).abs() / a.abs().max(f64::MIN_POSITIVE))
                    .fold(0.0f64, f64::max);
                if change < 1e-9 {
                    let modes = x
                        .iter()
                        .zip(lam.iter())
                        .filter(|(_, &l)| l > floor)
                        .take(k_wanted)
                        .map(|(xj, _)| {
                            let mut full = vec![0.0; n];
                            for (new, &old) in free.iter().enumerate() {
                                full[old] = xj[new];
                            }
                            full
                        })
                        .collect();
                    return Ok(BucklingSolution {
                        load_factors: lam_k,
                        modes,
                        iterations: it,
                    });
                }
            }
            prev = Some(lam_k);
        }
        q = x;
    }
    Err(SolveError::NotConverged {
        iterations: 200,
        residual: prev
            .map(|p| p.first().copied().unwrap_or(f64::NAN))
            .unwrap_or(f64::NAN),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BoundaryCondition, ConcentratedLoad, Element, ElementSet, Material, Mesh, Node};

    fn column_model(nx: usize, ny: usize, top_load: f64, pin: bool) -> Model {
        let (l, h, t) = (2.0_f64, 0.1_f64, 0.01_f64);
        let id = |i: usize, j: usize| (j * (nx + 1) + i + 1) as u32;
        let mut nodes = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(Node {
                    id: id(i, j),
                    coords: [i as f64 * l / nx as f64, j as f64 * h / ny as f64, 0.0],
                });
            }
        }
        let mut elements = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                elements.push(Element {
                    id: (j * nx + i + 1) as u32,
                    kind: ElemKind::Cps4,
                    connectivity: vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)],
                });
            }
        }
        // Column lies along x: ends are the SHORT edges i=0 and i=nx.
        let bottom: Vec<u32> = (0..=ny).map(|j| id(0, j)).collect();
        let top: Vec<u32> = (0..=ny).map(|j| id(nx, j)).collect();
        // Pinned-pinned: uy on both ends (roller lines), ux only at one
        // corner — rotation is caught by the far end's uy. Clamped-clamped:
        // both dofs on the whole left end.
        let mut bcs = vec![BoundaryCondition {
            node_ids: bottom.clone(),
            dofs: vec![2],
        }];
        if pin {
            bcs.push(BoundaryCondition {
                node_ids: vec![id(0, 0)],
                dofs: vec![1],
            });
            bcs.push(BoundaryCondition {
                node_ids: top.clone(),
                dofs: vec![2],
            });
        } else {
            bcs.push(BoundaryCondition {
                node_ids: bottom,
                dofs: vec![1],
            });
            bcs.push(BoundaryCondition {
                node_ids: top.clone(),
                dofs: vec![1, 2],
            });
        }
        let load_nodes: Vec<u32> = if pin {
            top.clone()
        } else {
            (0..=ny).map(|j| id(nx - 1, j)).collect()
        };
        Model {
            mesh: Mesh { nodes, elements },
            thickness: t,
            density: 7850.0,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: 210e9,
                poisson_ratio: 0.3,
            }],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=(nx * ny) as u32).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            boundary_conditions: bcs,
            cloads: load_nodes
                .iter()
                .map(|&n| ConcentratedLoad {
                    node_id: n,
                    dof: 1,
                    value: top_load,
                })
                .collect(),
        }
    }

    fn pcr_euler(l_eff: f64) -> f64 {
        let e = 210e9_f64;
        let i = 0.01 * 0.1_f64.powi(3) / 12.0;
        std::f64::consts::PI * std::f64::consts::PI * e * i / (l_eff * l_eff)
    }

    #[test]
    fn pinned_pinned_euler_column() {
        let m = column_model(40, 4, -1.0 / 5.0, true);
        let sol = solve_buckle(&m, 2).unwrap();
        let pcr = pcr_euler(2.0);
        let err = (sol.load_factors[0] - pcr) / pcr * 100.0;
        assert!(
            err.abs() < 10.0,
            "Pcr {:.0} vs Euler {:.0} ({err:+.1}%)",
            sol.load_factors[0],
            pcr
        );
        assert!(sol.load_factors[0] < sol.load_factors[1]);
    }

    #[test]
    fn cantilever_column_half_of_pinned() {
        // Vertical cantilever: L_eff = 2L -> Pcr = Pcr_euler(4.0).
        let (l, h) = (2.0_f64, 0.1_f64);
        let (nx, ny) = (4usize, 40usize);
        let id = |i: usize, j: usize| (j * (nx + 1) + i + 1) as u32;
        let mut nodes = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(Node {
                    id: id(i, j),
                    coords: [i as f64 * h / nx as f64, j as f64 * l / ny as f64, 0.0],
                });
            }
        }
        let mut elements = Vec::new();
        for j in 0..ny {
            for i in 0..nx {
                elements.push(Element {
                    id: (j * nx + i + 1) as u32,
                    kind: ElemKind::Cps4,
                    connectivity: vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)],
                });
            }
        }
        // Vertical column: base j=0 clamped, top j=ny axially loaded (dof 2).
        let bottom: Vec<u32> = (0..=nx).map(|i| id(i, 0)).collect();
        let load_row: Vec<u32> = (0..=nx).map(|i| id(i, ny)).collect();
        let m = Model {
            mesh: Mesh { nodes, elements },
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=(nx * ny) as u32).collect(),
            }],
            cloads: load_row
                .iter()
                .map(|&n| ConcentratedLoad {
                    node_id: n,
                    dof: 2,
                    value: -1.0 / 5.0,
                })
                .collect(),
            boundary_conditions: vec![BoundaryCondition {
                node_ids: bottom,
                dofs: vec![1, 2],
            }],
            ..column_model(1, 1, 0.0, true)
        };
        let sol = solve_buckle(&m, 1).unwrap();
        let expect = pcr_euler(4.0);
        let err = (sol.load_factors[0] - expect) / expect * 100.0;
        assert!(
            err.abs() < 15.0,
            "Pcr {:.0} vs Euler(2L) {:.0} ({err:+.1}%)",
            sol.load_factors[0],
            expect
        );
    }

    #[test]
    fn modes_have_unit_norm_and_zero_constrained() {
        let m = column_model(8, 2, -1.0, true);
        let sol = solve_buckle(&m, 1).unwrap();
        let mode = &sol.modes[0];
        let nrm: f64 = mode.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((nrm - 1.0).abs() < 1e-9);
        assert_eq!(mode[0], 0.0);
        assert_eq!(mode[1], 0.0);
    }

    #[test]
    fn no_bc_errors() {
        let mut m = column_model(4, 2, -1.0, true);
        m.boundary_conditions.clear();
        assert!(matches!(
            solve_buckle(&m, 1),
            Err(SolveError::NoConstraints)
        ));
    }
}
