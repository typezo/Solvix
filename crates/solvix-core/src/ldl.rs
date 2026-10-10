//! Sparse LDLT factorization with fill-reducing RCM ordering — the
//! direct fallback for systems where Jacobi-PCG exceeds its budget.
//!
//! Pipeline: reverse Cuthill-McKee permutation (BFS from the
//! lowest-degree node of each connected component, neighbors in
//! ascending index order, BFS order reversed) -> elimination tree from
//! the permuted structure -> up-looking LDLT with L stored in CSR rows
//! (the row orientation makes both triangular solves scatter/gather
//! passes over rows, so no CSC conversion is ever needed).

use crate::assembly::Csr;
use crate::SolveError;

/// Factored matrix: `P A P^T = L D L^T` with unit-lower L.
pub(crate) struct Ldlt {
    /// Unit lower triangular L in CSR rows, permuted ordering.
    l_row_ptr: Vec<usize>,
    l_col_idx: Vec<usize>,
    l_vals: Vec<f64>,
    d: Vec<f64>,
    /// new index -> old index.
    invp: Vec<usize>,
}

/// Reverse Cuthill-McKee ordering. Returns perm (old index -> new
/// index).
fn rcm(k: &Csr) -> Vec<usize> {
    let n = k.row_ptr.len() - 1;
    let mut perm = vec![usize::MAX; n];
    let mut order = Vec::with_capacity(n);
    let mut visited = vec![false; n];
    while order.len() < n {
        // Start each component at its minimum-degree node.
        let start = (0..n)
            .filter(|&i| !visited[i])
            .min_by_key(|&i| k.row_ptr[i + 1] - k.row_ptr[i])
            .unwrap();
        visited[start] = true;
        // BFS with a growing Vec as the queue; CSR columns are sorted,
        // so neighbors are visited in ascending index order.
        let mut queue = vec![start];
        let mut head = 0;
        while head < queue.len() {
            let u = queue[head];
            head += 1;
            order.push(u);
            for idx in k.row_ptr[u]..k.row_ptr[u + 1] {
                let v = k.col_idx[idx];
                if !visited[v] {
                    visited[v] = true;
                    queue.push(v);
                }
            }
        }
    }
    for (new, &old) in order.iter().rev().enumerate() {
        perm[old] = new;
    }
    perm
}

/// `out[r][c] = a[invp[r]][invp[c]]`, each row re-sorted ascending.
fn permuted(k: &Csr, perm: &[usize], invp: &[usize]) -> Csr {
    let n = perm.len();
    let mut row_ptr = Vec::with_capacity(n + 1);
    row_ptr.push(0);
    let mut col_idx = Vec::with_capacity(k.col_idx.len());
    let mut vals = Vec::with_capacity(k.vals.len());
    for &old in invp {
        let mut row: Vec<(usize, f64)> = (k.row_ptr[old]..k.row_ptr[old + 1])
            .map(|p| (perm[k.col_idx[p]], k.vals[p]))
            .collect();
        row.sort_unstable_by_key(|&(c, _)| c);
        for (c, v) in row {
            col_idx.push(c);
            vals.push(v);
        }
        row_ptr.push(col_idx.len());
    }
    Csr {
        row_ptr,
        col_idx,
        vals,
    }
}

/// Elimination tree (parent array, roots = `usize::MAX`) from the
/// strictly lower triangle of a symmetric CSR matrix.
fn elimination_tree(a: &Csr) -> Vec<usize> {
    let n = a.row_ptr.len() - 1;
    let mut parent = vec![usize::MAX; n];
    let mut ancestor = vec![usize::MAX; n];
    for i in 0..n {
        for p in a.row_ptr[i]..a.row_ptr[i + 1] {
            let mut node = a.col_idx[p];
            // Path from a column j < i up the tree toward i, with
            // path compression via `ancestor`.
            while node < i {
                let next = ancestor[node];
                ancestor[node] = i;
                if next == usize::MAX {
                    parent[node] = i;
                    break;
                }
                node = next;
            }
        }
    }
    parent
}

/// Row pattern of L for pivot `i`: the elimination-tree reach of row
/// `i` of A (cs_ereach). Uses `seen` as a generation-marked scratch;
/// returns unsorted indices k < i.
fn ereach(
    a: &Csr,
    parent: &[usize],
    seen: &mut [usize],
    gen: usize,
    i: usize,
    pattern: &mut Vec<usize>,
) {
    pattern.clear();
    seen[i] = gen;
    for p in a.row_ptr[i]..a.row_ptr[i + 1] {
        let mut node = a.col_idx[p];
        while node < i && seen[node] != gen {
            pattern.push(node);
            seen[node] = gen;
            match parent[node] {
                usize::MAX => break,
                par => node = par,
            }
        }
    }
}

impl Ldlt {
    /// Factorize a symmetric matrix (after Dirichlet elimination this is
    /// the SPD stiffness matrix). Fails with
    /// `SolveError::InvalidModel` when a pivot collapses, which is how
    /// structural singularities (underconstrained models) surface on
    /// the direct path.
    pub(crate) fn factorize(k: &Csr) -> Result<Ldlt, SolveError> {
        let n = k.row_ptr.len() - 1;
        let perm = rcm(k);
        let mut invp = vec![0; n];
        for (old, &new) in perm.iter().enumerate() {
            invp[new] = old;
        }
        let a = permuted(k, &perm, &invp);
        let parent = elimination_tree(&a);

        let mut l_row_ptr = vec![0usize];
        let mut l_col_idx = Vec::new();
        let mut l_vals = Vec::new();
        let mut d = vec![0.0; n];
        let mut x = vec![0.0; n];
        let mut seen = vec![usize::MAX; n];
        let mut pattern = Vec::new();

        for i in 0..n {
            let gen = i;
            ereach(&a, &parent, &mut seen, gen, i, &mut pattern);
            pattern.sort_unstable();
            // Dense scatter of row i of the permuted matrix.
            let mut aii = 0.0;
            let mut row_max = 0.0f64;
            for p in a.row_ptr[i]..a.row_ptr[i + 1] {
                let c = a.col_idx[p];
                x[c] = a.vals[p];
                row_max = row_max.max(a.vals[p].abs());
                if c == i {
                    aii = a.vals[p];
                }
            }
            // Sparse forward solve L w = a restricted to the row
            // pattern: x[k] holds w_k once every pattern entry below k
            // is final, because rows of L only reach smaller indices.
            for &kk in &pattern {
                let mut w = x[kk];
                for p in l_row_ptr[kk]..l_row_ptr[kk + 1] {
                    w -= l_vals[p] * x[l_col_idx[p]];
                }
                x[kk] = w; // store the final w_kk for later gathers
                let lik = w / d[kk];
                aii -= lik * w;
                l_col_idx.push(kk);
                l_vals.push(lik);
            }
            l_row_ptr.push(l_col_idx.len());
            // Reset the whole workspace for this row: the pattern
            // entries hold intermediate w values (fill included) and
            // the scatter leaves the structural entries dirty; either
            // going stale would poison a later row.
            for &kk in &pattern {
                x[kk] = 0.0;
            }
            for p in a.row_ptr[i]..a.row_ptr[i + 1] {
                x[a.col_idx[p]] = 0.0;
            }
            // ponytail: un-pivoted LDLT — valid for SPD stiffness
            // matrices; indefinite systems would need 2x2 Bunch-Kaufman
            // pivoting (extend the pivot step and the solve, keep this
            // threshold as the detector).
            if !aii.is_finite() || aii <= 1e-13 * row_max.max(f64::MIN_POSITIVE) {
                return Err(SolveError::InvalidModel(format!(
                    "stiffness matrix is not positive definite (pivot {} at equation {i})",
                    if aii < 0.0 { "negative" } else { "near-zero" }
                )));
            }
            d[i] = aii;
        }
        Ok(Ldlt {
            l_row_ptr,
            l_col_idx,
            l_vals,
            d,
            invp,
        })
    }

    /// Solve `A x = b` (both sides in the original ordering).
    pub(crate) fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.d.len();
        // L y = P b, unit diagonal.
        let mut y: Vec<f64> = (0..n).map(|i| b[self.invp[i]]).collect();
        for i in 0..n {
            let mut s = y[i];
            for p in self.l_row_ptr[i]..self.l_row_ptr[i + 1] {
                s -= self.l_vals[p] * y[self.l_col_idx[p]];
            }
            y[i] = s;
        }
        // D z = y.
        for (z, &dv) in y.iter_mut().zip(&self.d) {
            *z /= dv;
        }
        // L^T x = z: finalize x[i] top-down, then scatter its row.
        for i in (0..n).rev() {
            let xi = y[i];
            for p in self.l_row_ptr[i]..self.l_row_ptr[i + 1] {
                y[self.l_col_idx[p]] -= self.l_vals[p] * xi;
            }
        }
        let mut x = vec![0.0; n];
        for (&old, &v) in self.invp.iter().zip(&y) {
            x[old] = v;
        }
        x
    }

    /// Diagonal-based condition surrogate max|D_ii| / min|D_ii|. NOT a
    /// true condition number (that needs the full spectrum), but a
    /// cheap indicator of how hard the direct solve found the system.
    pub(crate) fn cond_estimate(&self) -> f64 {
        let mut lo = f64::INFINITY;
        let mut hi = 0.0f64;
        for &dv in &self.d {
            let m = dv.abs();
            lo = lo.min(m);
            hi = hi.max(m);
        }
        hi / lo
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly;
    use crate::assembly::tests::{hex_mesh, model_with_mesh, quad_mesh};
    use crate::solver;
    use crate::tests::solid_model;

    fn bandwidth(k: &Csr) -> usize {
        (0..k.row_ptr.len() - 1)
            .flat_map(|r| {
                (k.row_ptr[r]..k.row_ptr[r + 1]).map(move |p| r as isize - k.col_idx[p] as isize)
            })
            .map(|d| d.unsigned_abs())
            .max()
            .unwrap()
    }

    fn grid_csr(nx: usize, ny: usize) -> Csr {
        let model = model_with_mesh(quad_mesh(nx, ny, 1.0, 1.0));
        assembly::assemble(&model, 2).unwrap().0
    }

    #[test]
    fn rcm_is_a_valid_permutation_that_cuts_bandwidth() {
        let k = grid_csr(20, 5);
        let n = k.row_ptr.len() - 1;
        let perm = rcm(&k);
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..n).collect::<Vec<_>>(), "not a permutation");
        let mut invp = vec![0; n];
        for (old, &new) in perm.iter().enumerate() {
            invp[new] = old;
        }
        let reordered = permuted(&k, &perm, &invp);
        assert!(
            bandwidth(&reordered) < bandwidth(&k),
            "rcm bandwidth {} not smaller than natural {}",
            bandwidth(&reordered),
            bandwidth(&k)
        );
    }

    fn csr_from_rows(rows: &[&[(usize, f64)]]) -> Csr {
        let mut row_ptr = vec![0usize];
        let mut col_idx = Vec::new();
        let mut vals = Vec::new();
        for row in rows {
            for &(c, v) in *row {
                col_idx.push(c);
                vals.push(v);
            }
            row_ptr.push(col_idx.len());
        }
        Csr {
            row_ptr,
            col_idx,
            vals,
        }
    }

    #[test]
    fn factorize_and_solve_small_spd_system() {
        // [4 1 0; 1 3 1; 0 1 2]
        let k = csr_from_rows(&[
            &[(0, 4.0), (1, 1.0)],
            &[(0, 1.0), (1, 3.0), (2, 1.0)],
            &[(1, 1.0), (2, 2.0)],
        ]);
        let b = [1.0, 2.0, 3.0];
        let f = Ldlt::factorize(&k).unwrap();
        let x = f.solve(&b);
        let ku = k.mat_vec(&x);
        for (i, e) in b.iter().enumerate() {
            assert!((ku[i] - e).abs() < 1e-12, "row {i}: {} vs {e}", ku[i]);
        }
        assert!(f.cond_estimate() > 1.0);
    }

    #[test]
    fn factorize_rejects_singular_matrix() {
        // [1 1; 1 1] — second pivot is exactly zero.
        let k = csr_from_rows(&[&[(0, 1.0), (1, 1.0)], &[(0, 1.0), (1, 1.0)]]);
        match Ldlt::factorize(&k) {
            Err(SolveError::InvalidModel(msg)) => assert!(msg.contains("pivot"), "{msg}"),
            Err(other) => panic!("expected pivot failure, got {other:?}"),
            Ok(_) => panic!("expected pivot failure, got a factor"),
        }
    }

    #[test]
    fn fill_and_timing_3d_cantilever() {
        use crate::solve_with_method;
        use crate::Method;
        use std::time::Instant;

        let (nx, ny, nz) = (40usize, 4, 4);
        let (l, h, b) = (2.0, 0.1, 0.1);
        let mut model = solid_model(hex_mesh(nx, ny, nz, l, h, b));
        let root: Vec<u32> = (0..=nz)
            .flat_map(|k| {
                (0..=ny).map(move |j| (k * (ny + 1) * (nx + 1) + j * (nx + 1) + 1) as u32)
            })
            .collect();
        model.boundary_conditions.push(crate::BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2, 3],
        });
        let tip = ((nz / 2) * (ny + 1) * (nx + 1) + (ny / 2) * (nx + 1) + nx + 1) as u32;
        model.cloads = vec![crate::ConcentratedLoad {
            node_id: tip,
            dof: 2,
            value: -10000.0,
        }];

        let ndof = crate::ndof_per_node(&model).unwrap();
        let (mut k, mut f, _) = assembly::assemble(&model, ndof).unwrap();
        let constrained = solver::constrained_dofs(&model, ndof).unwrap();
        solver::apply_dirichlet(&mut k, &mut f, &constrained);

        let t0 = Instant::now();
        let fact = Ldlt::factorize(&k).unwrap();
        let t_factor = t0.elapsed();
        let t0 = Instant::now();
        let u = fact.solve(&f);
        let t_solve = t0.elapsed();
        let r = k.mat_vec(&u);
        let res = f
            .iter()
            .zip(&r)
            .map(|(a, c)| a - c)
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            / f.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(res < 1e-8, "ldlt residual {res}");
        eprintln!(
            "3d cantilever 40x4x4: n={} nnz(K)={} nnz(L)={} factor={:?} solve={:?} cond_est={:.3e}",
            f.len(),
            k.vals.len(),
            fact.l_vals.len(),
            t_factor,
            t_solve,
            fact.cond_estimate()
        );

        // End-to-end through the public API, both paths, for timing.
        for method in [Method::Auto, Method::Ldlt] {
            let t0 = Instant::now();
            let sol = solve_with_method(&model, method).unwrap();
            eprintln!(
                "{method:?} solve_with_method: {:?} residual {:.3e} cond {:.3e} iters {}",
                t0.elapsed(),
                sol.final_residual,
                sol.cond_estimate,
                sol.iterations
            );
            // The direct residual sits at the eps*cond floor (~1e-10
            // here); allow one extra decade over the PCG bar.
            let tol = if method == Method::Ldlt { 1e-8 } else { 1e-10 };
            assert!(sol.final_residual <= tol);
        }
    }
}
