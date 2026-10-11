//! # Modal analysis: natural frequencies and mode shapes
//!
//! Free vibration, `K phi = lambda M phi`, solved for the smallest `k`
//! eigenpairs with **subspace iteration** (v0.3 choice, documented):
//! inverse iteration of a block of `p > k` Ritz vectors through the
//! sparse LDLT of the iteration matrix, with the small projected
//! pencil `Y^T K Y c = lambda Y^T M Y` diagonalized densely (Cholesky
//! whitening + Jacobi rotations). Subspace iteration is deliberately
//! simpler than shift-invert Lanczos, reuses `ldl::Ldlt` unchanged,
//! and converges robustly for lowest modes of constrained structures
//! (linear rate `lambda_i / lambda_{p+1}`, buffered by the extra
//! `min(4, k)` vectors; relative eigenvalue tolerance 1e-9, at most
//! 100 outer iterations).
//!
//! **Free-free models** (issue #56): when the model carries no
//! boundary conditions, K is singular (the rigid-body modes sit at
//! lambda = 0) and plain inverse iteration cannot separate them from
//! the elastic spectrum. The iteration matrix becomes the SHIFT-INVERT
//! operator `(K + sigma M)^-1 M` with
//! `sigma = 1e-6 * trace(K)/trace(M)`: the trace ratio is a scale
//! estimate of the pencil's mean eigenvalue, so sigma sits far BELOW
//! the first elastic eigenvalue (rigid modes still surface first at
//! ~0 Hz, elastic ones essentially unperturbed — the Ritz values come
//! from the UNshifted projected pencil, sigma only conditions the
//! iteration) and far ABOVE the roundoff-zero of the rigid modes, so
//! `K + sigma M` is safely SPD and factorizes. All `k` modes are
//! returned: the rigid ones first (frequencies clamped at >= 0),
//! then the elastic ones.
//!
//! Constraints are handled by **free-dof extraction** rather than the
//! static-solve zeroing trick: the pencil is restricted to the
//! unconstrained dofs of K and M (rows/columns deleted, renumbered).
//! That is mathematically identical to zeroing constrained rows with
//! `M[i,i] = 1` but keeps the spurious constrained-dof eigenpairs out
//! of the spectrum entirely, and returned modes carry exact zeros at
//! constrained dofs.
//!
//! Consistent (not lumped) mass matrices: CPS3 and C3D4 in closed form,
//! CPS4/CPE4 with 2x2, CPS8/CPE8 with 3x3 and C3D8 with 2x2x2 Gauss
//! quadrature of `rho * N^T N`, B31 in closed form (translational
//! Hermite + axial, no rotary inertia).

use crate::assembly::{self, coords2, scatter, Csr};
use crate::element::{
    b31_mass, hex_b, hex_shapes, invalid, quad8_b, quad8_shapes, quad_b, quad_shapes, G3, W3,
};
use crate::ldl::Ldlt;
use crate::solver;
use crate::{Model, SolveError};

/// Options for [`solve_modal`].
#[derive(Debug, Clone)]
pub struct ModalOptions {
    /// Number of smallest modes to extract.
    pub num_modes: usize,
}

impl Default for ModalOptions {
    fn default() -> Self {
        Self { num_modes: 10 }
    }
}

/// Modal analysis output. `modes[j]` is mass-normalized
/// (`m^T M m = 1`), has one entry per dof
/// (`nodes * ndof_per_node`), and is exactly zero at constrained dofs.
#[derive(Debug, Clone, Default)]
pub struct ModalSolution {
    /// Smallest `num_modes` natural frequencies, ascending, in Hz.
    pub frequencies_hz: Vec<f64>,
    /// Mass-normalized mode shapes, parallel to `frequencies_hz`.
    pub modes: Vec<Vec<f64>>,
    /// Outer subspace-iterations used.
    pub iterations: usize,
}

// ---- element consistent mass matrices ----

/// CPS3 consistent mass, closed form: `rho*t*A/12` times the
/// [[2,1,1],[1,2,1],[1,1,2]] pattern per dof block. `None` for a
/// zero-area triangle (same trigger as the stiffness).
fn cps3_mass(coords: &[[f64; 2]; 3], rho: f64, t: f64) -> Option<[[f64; 6]; 6]> {
    let two_a = (coords[1][0] - coords[0][0]) * (coords[2][1] - coords[0][1])
        - (coords[2][0] - coords[0][0]) * (coords[1][1] - coords[0][1]);
    if two_a == 0.0 {
        return None;
    }
    let w = rho * t * two_a.abs() / 24.0; // rho * t * A / 12
    let mut m = [[0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            let v = w * if i == j { 2.0 } else { 1.0 };
            for d in 0..2 {
                m[2 * i + d][2 * j + d] = v;
            }
        }
    }
    Some(m)
}

/// CPS4 consistent mass: `rho*t*detJ*N^T N` under 2x2 Gauss. `None` if
/// det(J) <= 0 at a Gauss point (same trigger as the stiffness).
fn cps4_mass(coords: &[[f64; 2]; 4], rho: f64, t: f64) -> Option<[[f64; 8]; 8]> {
    let g = 1.0f64 / (3.0f64).sqrt();
    let mut m = [[0.0; 8]; 8];
    for &xi in &[-g, g] {
        for &eta in &[-g, g] {
            let (_, detj) = quad_b(coords, xi, eta)?;
            let n = quad_shapes(xi, eta);
            for i in 0..4 {
                for j in 0..4 {
                    let v = rho * t * detj * n[i] * n[j];
                    for d in 0..2 {
                        m[2 * i + d][2 * j + d] += v;
                    }
                }
            }
        }
    }
    Some(m)
}

/// C3D8 consistent mass: `rho*detJ*N^T N` under 2x2x2 Gauss. `None` if
/// det(J) <= 0 at a Gauss point (same trigger as the stiffness).
fn c3d8_mass(coords: &[[f64; 3]; 8], rho: f64) -> Option<[[f64; 24]; 24]> {
    let g = 1.0f64 / (3.0f64).sqrt();
    let mut m = [[0.0; 24]; 24];
    for &xi in &[-g, g] {
        for &eta in &[-g, g] {
            for &zeta in &[-g, g] {
                let (_, detj) = hex_b(coords, xi, eta, zeta)?;
                let n = hex_shapes(xi, eta, zeta);
                for i in 0..8 {
                    for j in 0..8 {
                        let v = rho * detj * n[i] * n[j];
                        for d in 0..3 {
                            m[3 * i + d][3 * j + d] += v;
                        }
                    }
                }
            }
        }
    }
    Some(m)
}

/// CPS8/CPE8 consistent mass: `rho*t*detJ*N^T N` under 3x3 Gauss
/// (weights 5/9, 8/9). `None` if det(J) <= 0 at a Gauss point (same
/// trigger as the stiffness).
fn cps8_mass(coords: &[[f64; 2]; 8], rho: f64, t: f64) -> Option<[[f64; 16]; 16]> {
    let mut m = [[0.0; 16]; 16];
    for (a, &xi) in [-G3, 0.0, G3].iter().enumerate() {
        for (b, &eta) in [-G3, 0.0, G3].iter().enumerate() {
            let (_, detj) = quad8_b(coords, xi, eta)?;
            let n = quad8_shapes(xi, eta);
            for i in 0..8 {
                for j in 0..8 {
                    let v = rho * t * detj * W3[a] * W3[b] * n[i] * n[j];
                    for d in 0..2 {
                        m[2 * i + d][2 * j + d] += v;
                    }
                }
            }
        }
    }
    Some(m)
}

/// C3D4 consistent mass, closed form: `rho*V/20` times the
/// [[2,1,1,1],...] pattern per dof block. `None` for zero volume.
fn c3d4_mass(coords: &[[f64; 3]; 4], rho: f64) -> Option<[[f64; 12]; 12]> {
    let p = |i: usize| {
        [
            coords[i][0] - coords[3][0],
            coords[i][1] - coords[3][1],
            coords[i][2] - coords[3][2],
        ]
    };
    let (a, b, c) = (p(0), p(1), p(2));
    let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]); // signed 6V
    if det == 0.0 {
        return None;
    }
    let w = rho * det.abs() / 120.0; // rho * V / 20
    let mut m = [[0.0; 12]; 12];
    for i in 0..4 {
        for j in 0..4 {
            let v = w * if i == j { 2.0 } else { 1.0 };
            for d in 0..3 {
                m[3 * i + d][3 * j + d] = v;
            }
        }
    }
    Some(m)
}

// ---- global consistent-mass assembly (mirrors assembly::assemble) ----

/// Assemble the consistent mass matrix (CSR, same dof map and sparsity
/// construction as the stiffness). Assumes `Model::validate()` passed
/// and `ndof = ndof_per_node(model)`.
pub(crate) fn assemble_mass(model: &Model, ndof: usize) -> Result<Csr, SolveError> {
    let n = model.mesh.nodes.len() * ndof;
    let mut trips = Vec::new();
    for i in 0..n {
        trips.push((i, i, 0.0)); // ensure diagonal entries exist
    }
    for (ei, e) in model.mesh.elements.iter().enumerate() {
        let node_idx: Vec<usize> = e
            .connectivity
            .iter()
            .map(|&nid| model.mesh.node_index(nid).unwrap())
            .collect();
        match e.kind {
            crate::ElemKind::Cps3 => {
                let coords = [
                    coords2(model, &node_idx, 0),
                    coords2(model, &node_idx, 1),
                    coords2(model, &node_idx, 2),
                ];
                let me = cps3_mass(&coords, model.density, model.thickness)
                    .ok_or_else(|| invalid(e.id, "degenerate (zero-area) triangle"))?;
                scatter(&mut trips, &node_idx, &me, ndof);
            }
            crate::ElemKind::Cps4 | crate::ElemKind::Cpe4 => {
                let coords = [
                    coords2(model, &node_idx, 0),
                    coords2(model, &node_idx, 1),
                    coords2(model, &node_idx, 2),
                    coords2(model, &node_idx, 3),
                ];
                let me = cps4_mass(&coords, model.density, model.thickness).ok_or_else(|| {
                    invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                })?;
                scatter(&mut trips, &node_idx, &me, ndof);
            }
            crate::ElemKind::Cps8 | crate::ElemKind::Cpe8 => {
                let coords: [[f64; 2]; 8] = core::array::from_fn(|k| coords2(model, &node_idx, k));
                let me = cps8_mass(&coords, model.density, model.thickness).ok_or_else(|| {
                    invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                })?;
                scatter(&mut trips, &node_idx, &me, ndof);
            }
            crate::ElemKind::B31 => {
                let coords = [coords2(model, &node_idx, 0), coords2(model, &node_idx, 1)];
                let (a, _) = model
                    .beam_section_for_element(ei)
                    .ok_or_else(|| invalid(e.id, "no beam section assigned (*BEAM SECTION)"))?;
                let me = b31_mass(&coords, model.density, a)
                    .ok_or_else(|| invalid(e.id, "zero-length beam"))?;
                scatter(&mut trips, &node_idx, &me, ndof);
            }
            crate::ElemKind::C3d8 => {
                let coords: [[f64; 3]; 8] =
                    core::array::from_fn(|k| model.mesh.nodes[node_idx[k]].coords);
                let me = c3d8_mass(&coords, model.density).ok_or_else(|| {
                    invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                })?;
                scatter(&mut trips, &node_idx, &me, ndof);
            }
            crate::ElemKind::C3d4 => {
                let coords: [[f64; 3]; 4] =
                    core::array::from_fn(|k| model.mesh.nodes[node_idx[k]].coords);
                let me = c3d4_mass(&coords, model.density)
                    .ok_or_else(|| invalid(e.id, "degenerate (zero-volume) tetrahedron"))?;
                scatter(&mut trips, &node_idx, &me, ndof);
            }
        }
    }
    Ok(assembly::csr_from_triplets(n, trips))
}

// ---- constraint handling: free-dof principal submatrix ----

/// Principal submatrix of `a` on the dofs where `is_constrained` is
/// false, renumbered to 0..n_free. Also returns the kept old indices,
/// ascending (the map for expanding solution vectors back).
pub(crate) fn free_submatrix(a: &Csr, is_constrained: &[bool]) -> (Csr, Vec<usize>) {
    let free: Vec<usize> = (0..is_constrained.len())
        .filter(|&i| !is_constrained[i])
        .collect();
    let mut new_of = vec![usize::MAX; is_constrained.len()];
    for (new, &old) in free.iter().enumerate() {
        new_of[old] = new;
    }
    let mut trips = Vec::with_capacity(a.vals.len());
    for (new, &r) in free.iter().enumerate() {
        for idx in a.row_ptr[r]..a.row_ptr[r + 1] {
            let c = a.col_idx[idx];
            if !is_constrained[c] {
                trips.push((new, new_of[c], a.vals[idx]));
            }
        }
    }
    (assembly::csr_from_triplets(free.len(), trips), free)
}

// ---- small dense symmetric helpers (p x p, p <= num_modes + 4) ----

pub(crate) fn zeros(p: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; p]; p]
}

pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Lower Cholesky factor of an SPD matrix, `None` when a pivot is not
/// positive (NaN pivots included).
fn cholesky(m: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = m.len();
    let mut g = zeros(n);
    for i in 0..n {
        let mut row = vec![0.0; n];
        for j in 0..=i {
            // Row i is still in `row` (the diagonal case needs its own
            // prefix); earlier rows are already final in g.
            let gj: &Vec<f64> = if j == i { &row } else { &g[j] };
            let s = m[i][j] - row.iter().zip(gj).take(j).map(|(x, y)| x * y).sum::<f64>();
            if i == j {
                if s.is_nan() || s <= 0.0 {
                    return None;
                }
                row[i] = s.sqrt();
            } else {
                row[j] = s / g[j][j];
            }
        }
        g[i] = row;
    }
    Some(g)
}

/// Solve `G y = b` for a lower-triangular `G`.
fn solve_lower(g: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut y = b.to_vec();
    for i in 0..n {
        let mut s = y[i];
        s -= g[i]
            .iter()
            .zip(&y)
            .take(i)
            .map(|(gv, yv)| gv * yv)
            .sum::<f64>();
        y[i] = s / g[i][i];
    }
    y
}

/// Solve `G^T y = b` for a lower-triangular `G`.
fn solve_lower_t(g: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut y = b.to_vec();
    for i in (0..n).rev() {
        let mut s = y[i];
        s -= g
            .iter()
            .zip(&y)
            .skip(i + 1)
            .map(|(grow, yv)| grow[i] * yv)
            .sum::<f64>();
        y[i] = s / g[i][i];
    }
    y
}

/// Cyclic Jacobi diagonalization of a small symmetric matrix.
/// Returns eigenvalues and eigenvector **columns**, unsorted.
fn jacobi_eigen(a0: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a0.len();
    let mut a = a0.to_vec();
    let mut v = zeros(n);
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for (i, row) in a.iter().enumerate() {
            for &e in row.iter().skip(i + 1) {
                off += e * e;
            }
        }
        let scale: f64 = a
            .iter()
            .enumerate()
            .map(|(i, row)| row[i] * row[i])
            .sum::<f64>()
            .sqrt();
        if off.sqrt() <= 1e-14 * scale.max(f64::MIN_POSITIVE) {
            break;
        }
        for p in 0..n - 1 {
            for q in p + 1..n {
                if a[p][q].abs() <= 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let sign = if theta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = (t * t + 1.0).sqrt().recip();
                let s = t * c;
                // Column rotation on both a and v: a <- a J, v <- v J.
                for (ar, vr) in a.iter_mut().zip(v.iter_mut()) {
                    let (akp, akq) = (ar[p], ar[q]);
                    ar[p] = c * akp - s * akq;
                    ar[q] = s * akp + c * akq;
                    let (vkp, vkq) = (vr[p], vr[q]);
                    vr[p] = c * vkp - s * vkq;
                    vr[q] = s * vkp + c * vkq;
                }
                // Row rotation on a: a <- J^T a.
                let newp: Vec<f64> = a[p]
                    .iter()
                    .zip(&a[q])
                    .map(|(&x, &y)| c * x - s * y)
                    .collect();
                let newq: Vec<f64> = a[p]
                    .iter()
                    .zip(&a[q])
                    .map(|(&x, &y)| s * x + c * y)
                    .collect();
                a[p] = newp;
                a[q] = newq;
            }
        }
    }
    (a.iter().enumerate().map(|(i, row)| row[i]).collect(), v)
}

/// Eigenpairs of the dense projected pencil `kk c = lam mm c`
/// (kk, mm symmetric, mm SPD), eigenvalues ascending, eigenvectors as
/// columns. Whitening: `mm = G G^T`, Jacobi on `G^-1 kk G^-T`,
/// `c = G^-T z`. `None` when `mm` stays indefinite even after the
/// ridge below.
pub(crate) fn generalized_eigenpairs(
    kk: &[Vec<f64>],
    mm: &[Vec<f64>],
) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = kk.len();
    let mut mm = mm.to_vec();
    let max_diag = (0..n).map(|i| mm[i][i].abs()).fold(0.0f64, f64::max);
    let mean_diag = max_diag; // scale reference: strongest diagonal entry
    let g = match cholesky(&mm) {
        Some(g) => g,
        None => {
            // Near-converged blocks make cond(Y^T M Y) climb; an adaptive
            // relative ridge keeps the Cholesky alive. Regularizing the
            // projected mass perturbs eigenvalues O(ridge) — 1e-4 relative
            // is the ceiling (raised for buckling pencils with strong local tension).
            let cap = 1e-2 * mean_diag.abs().max(f64::MIN_POSITIVE);
            let mut ridge = 1e-10 * mean_diag.abs();
            loop {
                for (i, row) in mm.iter_mut().enumerate() {
                    row[i] += ridge;
                }
                if let Some(g) = cholesky(&mm) {
                    break g;
                }
                if ridge >= cap {
                    return None; // ridge ceiling exhausted
                }
                ridge *= 100.0;
            }
        }
    };
    // H = G^-1 kk G^-T in three column solves: W = G^-1 kk, then
    // H = G^-1 W^T (kk symmetric), symmetrized against roundoff by
    // averaging with its transpose at assembly.
    let w: Vec<Vec<f64>> = (0..n)
        .map(|j| {
            let col: Vec<f64> = (0..n).map(|i| kk[i][j]).collect();
            solve_lower(&g, &col)
        })
        .collect();
    let h_cols: Vec<Vec<f64>> = (0..n)
        .map(|j| {
            let row_of_w: Vec<f64> = (0..n).map(|i| w[i][j]).collect();
            solve_lower(&g, &row_of_w)
        })
        .collect();
    let mut h = zeros(n);
    for (i, row) in h.iter_mut().enumerate() {
        for (j, hv) in row.iter_mut().enumerate() {
            *hv = 0.5 * (h_cols[j][i] + h_cols[i][j]);
        }
    }
    let (lam, z) = jacobi_eigen(&h);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| lam[a].total_cmp(&lam[b]));
    let lams: Vec<f64> = order.iter().map(|&i| lam[i]).collect();
    let cols: Vec<Vec<f64>> = order
        .iter()
        .map(|&i| {
            let z_col: Vec<f64> = (0..n).map(|k| z[k][i]).collect();
            solve_lower_t(&g, &z_col)
        })
        .collect();
    Some((lams, cols))
}

// ---- deterministic random block (no external rng dependency) ----

/// xorshift64* seeded from a fixed constant: the start block — and so
/// the whole eigen-solution — is reproducible across runs and platforms.
pub(crate) struct Rng(pub(crate) u64);

impl Rng {
    /// Uniform in [-1, 1).
    pub(crate) fn unit(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        (x.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }
}

// ---- public entry point ----

/// Extract the `num_modes` lowest natural frequencies (Hz) and
/// mass-normalized mode shapes via subspace iteration on the CSR
/// stiffness/mass pair; see the [module docs](self) for the algorithm
/// and its v1 ceilings.
///
/// Free-free models (no boundary conditions) are supported: the
/// iteration rides the shift-invert operator `(K + sigma M)^-1 M` (see
/// the module docs for the sigma heuristic), so the rigid-body modes
/// come out near 0 Hz ahead of the elastic ones and all `num_modes`
/// modes are returned.
///
/// Errors: `InvalidModel` for the usual [`Model::validate`] failures
/// plus non-positive `density`/`thickness` and `num_modes` not in
/// `1..free_dof_count`, `NotConverged` when 100 outer iterations do
/// not settle the eigenvalues to 1e-9 relative.
pub fn solve_modal(model: &Model, opts: &ModalOptions) -> Result<ModalSolution, SolveError> {
    let free_free = model.load_state().bcs.is_empty();
    if let Err(e) = model.validate() {
        // The free-free path intentionally runs without boundary
        // conditions; every other validate() failure still applies.
        if !free_free || e != crate::ModelError::NoBoundaryConditions {
            return Err(SolveError::InvalidModel(format!("{e:?}")));
        }
    }
    if model.thickness <= 0.0 || model.thickness.is_nan() {
        return Err(SolveError::InvalidModel(
            "thickness must be positive".into(),
        ));
    }
    if model.density <= 0.0 || model.density.is_nan() {
        return Err(SolveError::InvalidModel("density must be positive".into()));
    }
    let k = opts.num_modes;
    if k == 0 {
        return Err(SolveError::InvalidModel(
            "num_modes must be positive".into(),
        ));
    }

    let ndof = crate::ndof_per_node(model)?;
    let loads = model.load_state();
    let (k_full, _, _) = assembly::assemble(model, ndof, &loads)?;
    let m_full = assemble_mass(model, ndof)?;
    let constrained = solver::constrained_dofs(&model.mesh, loads.bcs, ndof)?;

    let n = k_full.row_ptr.len() - 1;
    let mut is_constrained = vec![false; n];
    for c in &constrained {
        is_constrained[c.row] = true;
    }
    let (k_ff, free) = free_submatrix(&k_full, &is_constrained);
    let (m_ff, _) = free_submatrix(&m_full, &is_constrained);
    let n_free = free.len();
    if k >= n_free {
        return Err(SolveError::InvalidModel(format!(
            "num_modes {k} must be smaller than the free dof count {n_free}"
        )));
    }
    // A singular K_ff here means the constraints leave a rigid-body
    // mode: the same pivot failure the static direct path reports.
    // Free-free: factor the shift-invert matrix K + sigma*M instead
    // (SPD by construction; sigma heuristic in the module docs).
    let fact = if free_free {
        // K and M share the CSR pattern (same element blocks, same
        // diagonal seeds), so entrywise addition is pattern-aligned.
        let diag_sum = |a: &Csr| {
            (0..a.row_ptr.len() - 1)
                .filter_map(|i| {
                    (a.row_ptr[i]..a.row_ptr[i + 1])
                        .find(|&idx| a.col_idx[idx] == i)
                        .map(|idx| a.vals[idx])
                })
                .sum::<f64>()
        };
        let sigma = 1e-6 * diag_sum(&k_ff) / diag_sum(&m_ff);
        let mut trips: Vec<(usize, usize, f64)> = Vec::with_capacity(k_ff.vals.len());
        for r in 0..k_ff.row_ptr.len() - 1 {
            for idx in k_ff.row_ptr[r]..k_ff.row_ptr[r + 1] {
                trips.push((
                    r,
                    k_ff.col_idx[idx],
                    k_ff.vals[idx] + sigma * m_ff.vals[idx],
                ));
            }
        }
        Ldlt::factorize(&assembly::csr_from_triplets(k_ff.row_ptr.len() - 1, trips))?
    } else {
        Ldlt::factorize(&k_ff)?
    };

    let p = (k + k.min(4)).min(n_free);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut q: Vec<Vec<f64>> = (0..p)
        .map(|_| (0..n_free).map(|_| rng.unit()).collect())
        .collect();

    let mut prev: Option<Vec<f64>> = None;
    let mut last_change = f64::INFINITY;
    for it in 1..=100 {
        // Inverse iteration: Y = K^-1 M Q through the shared LDLT.
        let y: Vec<Vec<f64>> = q.iter().map(|qj| fact.solve(&m_ff.mat_vec(qj))).collect();
        // Rayleigh-Ritz projections on the iterated subspace.
        let ky: Vec<Vec<f64>> = y.iter().map(|yj| k_ff.mat_vec(yj)).collect();
        let my: Vec<Vec<f64>> = y.iter().map(|yj| m_ff.mat_vec(yj)).collect();
        let mut kk = zeros(p);
        let mut mm = zeros(p);
        for a in 0..p {
            for b in a..p {
                let (kv, mv) = (dot(&y[a], &ky[b]), dot(&y[a], &my[b]));
                kk[a][b] = kv;
                kk[b][a] = kv;
                mm[a][b] = mv;
                mm[b][a] = mv;
            }
        }
        let (lam, c) = generalized_eigenpairs(&kk, &mm).ok_or_else(|| {
            SolveError::InvalidModel("projected mass matrix is not positive definite".into())
        })?;
        // Ritz vectors X = Y C, then mass-normalize each column.
        let x: Vec<Vec<f64>> = c
            .iter()
            .map(|cj| {
                let mut v = vec![0.0; n_free];
                for (yj, &cij) in y.iter().zip(cj) {
                    for (vi, &yv) in v.iter_mut().zip(yj) {
                        *vi += cij * yv;
                    }
                }
                let nrm = dot(&v, &m_ff.mat_vec(&v)).sqrt();
                for s in &mut v {
                    *s /= nrm;
                }
                v
            })
            .collect();

        let lam_k: Vec<f64> = lam[..k].to_vec();
        if let Some(pr) = &prev {
            // Free-free rigid modes carry roundoff-scale Ritz values;
            // scaling the change by the LARGEST tracked eigenvalue
            // keeps that noise from blocking convergence (for
            // constrained models this is the old per-mode measure
            // with a uniform, tighter denominator).
            let scale = lam_k
                .iter()
                .cloned()
                .fold(0.0f64, f64::max)
                .max(f64::MIN_POSITIVE);
            last_change = lam_k
                .iter()
                .zip(pr)
                .map(|(&a, &b)| (a - b).abs())
                .fold(0.0f64, f64::max)
                / scale;
            if last_change < 1e-9 {
                let frequencies_hz: Vec<f64> = lam_k
                    .iter()
                    .map(|&l| l.max(0.0).sqrt() / (2.0 * std::f64::consts::PI))
                    .collect();
                let modes = x[..k]
                    .iter()
                    .map(|xj| {
                        let mut full = vec![0.0; n];
                        for (new, &old) in free.iter().enumerate() {
                            full[old] = xj[new];
                        }
                        full
                    })
                    .collect();
                return Ok(ModalSolution {
                    frequencies_hz,
                    modes,
                    iterations: it,
                });
            }
        }
        prev = Some(lam_k);
        q = x;
    }
    Err(SolveError::NotConverged {
        iterations: 100,
        residual: last_change,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::tests::{beam_mesh, hex_mesh, quad_mesh};
    use crate::tests::solid_model;
    use crate::{BoundaryCondition, ElementSet, Material, Mesh};

    const E: f64 = 210e9;
    const NU: f64 = 0.3;
    const RHO: f64 = 7850.0;
    const T: f64 = 0.01; // plane-stress thickness

    fn steel_model(mesh: Mesh, thickness: f64) -> Model {
        let n_elem = mesh.elements.len() as u32;
        Model {
            mesh,
            thickness,
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
            ..Default::default()
        }
    }

    /// Slender CPS4 cantilever: nx x ny mesh of an l x h rectangle,
    /// fully clamped on the left edge.
    fn cantilever(nx: usize, ny: usize, l: f64, h: f64) -> Model {
        let mut model = steel_model(quad_mesh(nx, ny, l, h), T);
        let root: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2],
        }];
        model
    }

    /// Euler-Bernoulli first bending frequency of a cantilever.
    fn cantilever_f1_exact(l: f64, h: f64) -> f64 {
        let area = T * h;
        let inertia = T * h * h * h / 12.0;
        1.875f64 * 1.875 / (2.0 * std::f64::consts::PI)
            * (E * inertia / (RHO * area * l * l * l * l)).sqrt()
    }

    fn modes(model: &Model, k: usize) -> ModalSolution {
        solve_modal(model, &ModalOptions { num_modes: k }).unwrap()
    }

    #[test]
    fn cantilever_first_bending_frequency() {
        let (nx, ny, l, h) = (40, 4, 2.0, 0.1);
        let model = cantilever(nx, ny, l, h);
        let sol = modes(&model, 6);
        let exact = cantilever_f1_exact(l, h);
        let rel = (sol.frequencies_hz[0] - exact) / exact;
        eprintln!(
            "cps4 cantilever f1: fem={:.4} Hz, analytical={:.4} Hz ({:+.2}%), iters={}",
            sol.frequencies_hz[0],
            exact,
            100.0 * rel,
            sol.iterations
        );
        assert!(
            rel.abs() <= 0.08,
            "f1 off by {:.1}% (fem={}, exact={})",
            100.0 * rel,
            sol.frequencies_hz[0],
            exact
        );
        assert!(sol.iterations <= 100);
    }

    #[test]
    fn b31_cantilever_first_bending_frequency() {
        // 8 Hermite elements, L = 2, 0.1 x 0.1 section: Euler-Bernoulli
        // f1 = (1.875^2/2pi) sqrt(EI/(rho A L^4)); consistent Hermite
        // mass converges from above, well inside 2% at 8 elements.
        let (n, l) = (8, 2.0);
        let (a, inertia) = (0.01, 1.0e-4 / 12.0);
        let mut model = steel_model(beam_mesh(n, l), 1.0);
        model.beam_sections = vec![("all".into(), a, inertia)];
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![1, 2, 3],
        }];
        let sol = modes(&model, 3);
        let exact = 1.875f64 * 1.875 / (2.0 * std::f64::consts::PI)
            * (E * inertia / (RHO * a * l * l * l * l)).sqrt();
        let rel = (sol.frequencies_hz[0] - exact) / exact;
        eprintln!(
            "b31 cantilever f1: fem={:.4} Hz, analytical={:.4} Hz ({:+.2}%), iters={}",
            sol.frequencies_hz[0],
            exact,
            100.0 * rel,
            sol.iterations
        );
        assert!(
            rel.abs() <= 0.02,
            "f1 off by {:.2}% (fem={}, exact={})",
            100.0 * rel,
            sol.frequencies_hz[0],
            exact
        );
    }

    #[test]
    fn mode1_is_transverse_at_tip() {
        let (nx, ny) = (40, 4);
        let model = cantilever(nx, ny, 2.0, 0.1);
        let sol = modes(&model, 3);
        // Middle node of the free (right) edge.
        let tip = model
            .mesh
            .node_index((2 * (nx + 1) + nx + 1) as u32)
            .unwrap();
        let (ux, uy) = (sol.modes[0][2 * tip], sol.modes[0][2 * tip + 1]);
        assert!(
            uy.abs() > 50.0 * ux.abs(),
            "tip |uy|={uy:.3e} not >> |ux|={ux:.3e}"
        );
    }

    #[test]
    fn modes_are_mass_orthonormal() {
        let model = cantilever(20, 4, 2.0, 0.1);
        let sol = modes(&model, 6);
        let m = assemble_mass(&model, 2).unwrap();
        for i in 0..sol.modes.len() {
            for j in 0..sol.modes.len() {
                let p = dot(&sol.modes[j], &m.mat_vec(&sol.modes[i]));
                let want = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (p - want).abs() < 1e-6,
                    "mode ({i},{j}): m_i^T M m_j = {p:.3e}, want {want}"
                );
            }
        }
    }

    #[test]
    fn frequencies_ascending_and_positive() {
        let sol = modes(&cantilever(20, 4, 2.0, 0.1), 6);
        assert!(sol.frequencies_hz[0] > 0.0);
        for w in sol.frequencies_hz.windows(2) {
            assert!(w[1] > w[0], "{:?}", sol.frequencies_hz);
        }
    }

    #[test]
    fn constrained_dof_entries_are_zero() {
        let (nx, ny) = (20, 4);
        let model = cantilever(nx, ny, 2.0, 0.1);
        let sol = modes(&model, 3);
        let root: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        for mode in &sol.modes {
            for &nid in &root {
                let i = model.mesh.node_index(nid).unwrap();
                assert_eq!(mode[2 * i], 0.0);
                assert_eq!(mode[2 * i + 1], 0.0);
            }
        }
    }

    #[test]
    fn clamped_clamped_first_mode_exceeds_cantilever() {
        let (nx, ny, l, h) = (20, 4, 2.0, 0.1);
        let cant = modes(&cantilever(nx, ny, l, h), 3);
        let mut cc = cantilever(nx, ny, l, h);
        let tip: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + nx + 1) as u32).collect();
        cc.boundary_conditions.push(BoundaryCondition {
            node_ids: tip,
            dofs: vec![1, 2],
        });
        let fixed = modes(&cc, 3);
        // Beam theory ratio (4.730/1.875)^2 = 6.4; allow FEM slack down to 4.
        assert!(
            fixed.frequencies_hz[0] > 4.0 * cant.frequencies_hz[0],
            "clamped-clamped {} not above 4x cantilever {}",
            fixed.frequencies_hz[0],
            cant.frequencies_hz[0]
        );
    }

    #[test]
    fn c3d8_cantilever_first_mode() {
        // 80x2x2 hexes, L = 2, h = b = 0.05: cube elements, L/h = 40 and
        // two elements through the section keep full-integration shear
        // locking modest (one element through is far worse).
        let (nx, ny, nz) = (80, 2, 2);
        let (l, h, b) = (2.0, 0.05, 0.05);
        let mut model = solid_model(hex_mesh(nx, ny, nz, l, h, b));
        let root: Vec<u32> = (0..=nz)
            .flat_map(|k| {
                (0..=ny).map(move |j| (k * (ny + 1) * (nx + 1) + j * (nx + 1) + 1) as u32)
            })
            .collect();
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2, 3],
        }];
        let sol = modes(&model, 3);
        let area = h * b;
        let inertia = b * h * h * h / 12.0;
        let exact = 1.875f64 * 1.875 / (2.0 * std::f64::consts::PI)
            * (E * inertia / (RHO * area * l * l * l * l)).sqrt();
        let rel = (sol.frequencies_hz[0] - exact) / exact;
        eprintln!(
            "c3d8 cantilever f1: fem={:.4} Hz, analytical={:.4} Hz ({:+.2}%), iters={}",
            sol.frequencies_hz[0],
            exact,
            100.0 * rel,
            sol.iterations
        );
        assert!(
            rel.abs() <= 0.15,
            "f1 off by {:.1}% (fem={}, exact={})",
            100.0 * rel,
            sol.frequencies_hz[0],
            exact
        );
    }

    #[test]
    fn consistent_mass_totals_match_rho_times_volume() {
        // Per dof, the mass matrix sums to rho * (volume): a rigid unit
        // translation in one dof carries the full element mass.
        let tri = cps3_mass(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], RHO, T).unwrap();
        let quad = cps4_mass(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], RHO, T).unwrap();
        let hex = c3d8_mass(
            &[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
                [0.0, 1.0, 1.0],
            ],
            RHO,
        )
        .unwrap();
        let tet = c3d4_mass(
            &[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            RHO,
        )
        .unwrap();
        let sum2 = |m: &[[f64; 8]; 8]| {
            (0..4)
                .map(|i| (0..4).map(|j| m[2 * i][2 * j]).sum::<f64>())
                .sum::<f64>()
        };
        let sum6 = |m: &[[f64; 6]; 6]| {
            (0..3)
                .map(|i| (0..3).map(|j| m[2 * i][2 * j]).sum::<f64>())
                .sum::<f64>()
        };
        let sum24 = |m: &[[f64; 24]; 24]| {
            (0..8)
                .map(|i| (0..8).map(|j| m[3 * i][3 * j]).sum::<f64>())
                .sum::<f64>()
        };
        let sum12 = |m: &[[f64; 12]; 12]| {
            (0..4)
                .map(|i| (0..4).map(|j| m[3 * i][3 * j]).sum::<f64>())
                .sum::<f64>()
        };
        let eps = 1e-9 * RHO;
        assert!((sum6(&tri) - RHO * T * 0.5).abs() < eps);
        assert!((sum2(&quad) - RHO * T).abs() < eps);
        assert!((sum24(&hex) - RHO).abs() < eps);
        assert!((sum12(&tet) - RHO / 6.0).abs() < eps);
    }

    #[test]
    fn modal_guards() {
        // No BCs is no longer rejected as NoConstraints: free-free
        // models ride the shift-invert path (the dedicated test below
        // checks the physics on a proper mesh; this stubby locked one
        // may legitimately exhaust the iteration budget).
        let mut model = cantilever(4, 2, 1.0, 0.1);
        model.boundary_conditions.clear();
        assert!(!matches!(
            solve_modal(&model, &ModalOptions { num_modes: 3 }),
            Err(SolveError::NoConstraints)
        ));

        let mut model = cantilever(4, 2, 1.0, 0.1);
        model.density = 0.0;
        assert!(matches!(
            solve_modal(&model, &ModalOptions::default()),
            Err(SolveError::InvalidModel(m)) if m.contains("density")
        ));

        let model = cantilever(4, 2, 1.0, 0.1);
        assert!(matches!(
            solve_modal(&model, &ModalOptions { num_modes: 0 }),
            Err(SolveError::InvalidModel(m)) if m.contains("num_modes")
        ));
        // 4x2 mesh: 15 nodes * 2 dofs - 6 constrained = 24 free dofs;
        // asking for all of them leaves no iteration buffer.
        assert!(matches!(
            solve_modal(&model, &ModalOptions { num_modes: 24 }),
            Err(SolveError::InvalidModel(m)) if m.contains("num_modes")
        ));
    }

    #[test]
    fn free_free_beam_rigid_modes_then_first_bending() {
        // Unconstrained CPS4 strip, L = 1, h = t = 0.1: the 3 planar
        // rigid-body modes (ux, uy, rotz) come out near 0 Hz, then the
        // first elastic mode — free-free beam first bending
        // f1 = (4.730^2/2pi) sqrt(EI/(rho A L^4)) — within 10%.
        let (nx, ny, l, h, t) = (40usize, 4, 1.0, 0.1, 0.1);
        let model = steel_model(quad_mesh(nx, ny, l, h), t);
        let sol = modes(&model, 6);
        let (area, inertia) = (t * h, t * h * h * h / 12.0);
        let exact = 4.730f64 * 4.730 / (2.0 * std::f64::consts::PI)
            * (E * inertia / (RHO * area * l * l * l * l)).sqrt();
        for (i, &f) in sol.frequencies_hz.iter().take(3).enumerate() {
            assert!(f < 0.1, "rigid mode {i}: {f:.6} Hz not ~0");
        }
        let rel = (sol.frequencies_hz[3] - exact) / exact;
        eprintln!(
            "free-free strip: rigid = [{:.3e}, {:.3e}, {:.3e}] Hz, f4: fem={:.2} Hz, exact={:.2} Hz ({:+.2}%), iters={}",
            sol.frequencies_hz[0],
            sol.frequencies_hz[1],
            sol.frequencies_hz[2],
            sol.frequencies_hz[3],
            exact,
            100.0 * rel,
            sol.iterations
        );
        assert!(
            rel.abs() <= 0.10,
            "first bending off by {:.1}% (fem={}, exact={})",
            100.0 * rel,
            sol.frequencies_hz[3],
            exact
        );
        assert!(sol.iterations <= 100);
    }

    #[test]
    fn jacobi_matches_known_spectrum() {
        let a = vec![
            vec![2.0, 1.0, 0.0],
            vec![1.0, 2.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let (lam, v) = jacobi_eigen(&a);
        let mut pairs: Vec<(f64, usize)> = lam
            .iter()
            .copied()
            .enumerate()
            .map(|(i, l)| (l, i))
            .collect();
        pairs.sort_by(|x, y| x.0.total_cmp(&y.0));
        for (got, want) in pairs
            .iter()
            .zip([2.0 - 2.0f64.sqrt(), 2.0, 2.0 + 2.0f64.sqrt()])
        {
            assert!((got.0 - want).abs() < 1e-12, "{} vs {want}", got.0);
        }
        // A v = lam v per column, eigenvalue paired with its column.
        for &(l, j) in &pairs {
            for i in 0..3 {
                let av: f64 = (0..3).map(|k| a[i][k] * v[k][j]).sum();
                assert!((av - l * v[i][j]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn generalized_eigenpairs_solve_the_dense_pencil() {
        // kk = B^T D B with B of full rank, mm = A A^T + 3I (SPD).
        let b = [
            [3.0, 1.0, 0.0, 2.0],
            [0.0, 2.0, 1.0, 1.0],
            [1.0, 0.0, 3.0, 0.5],
        ];
        let d = [5.0, 7.0, 9.0, 11.0];
        let mut kk = zeros(3);
        for r in 0..3 {
            for c in 0..3 {
                kk[r][c] = (0..4).map(|m| b[r][m] * d[m] * b[c][m]).sum();
            }
        }
        let a = [[1.0, 0.3, 0.1], [0.0, 1.0, 0.2], [0.5, 0.0, 1.0]];
        let mut mm = zeros(3);
        for r in 0..3 {
            for c in 0..3 {
                mm[r][c] =
                    (0..3).map(|m| a[m][r] * a[m][c]).sum::<f64>() + if r == c { 3.0 } else { 0.0 };
            }
        }
        let (lam, cols) = generalized_eigenpairs(&kk, &mm).unwrap();
        // kk c = lam mm c, and M-orthonormal columns.
        for (j, c) in cols.iter().enumerate() {
            for i in 0..3 {
                let kc: f64 = (0..3).map(|m| kk[i][m] * c[m]).sum();
                let mc: f64 = (0..3).map(|m| mm[i][m] * c[m]).sum();
                assert!(
                    (kc - lam[j] * mc).abs() < 1e-8 * kc.abs().max(1.0),
                    "row {i} mode {j}: Kc={kc}, lam*Mc={}",
                    lam[j]
                );
            }
        }
        for (i, ci) in cols.iter().enumerate() {
            for (j, cj) in cols.iter().enumerate() {
                let mi: f64 = (0..3)
                    .map(|r| (0..3).map(|s| mm[r][s] * ci[s]).sum::<f64>() * cj[r])
                    .sum();
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((mi - want).abs() < 1e-10, "({i},{j}) {mi}");
            }
        }
    }
}
