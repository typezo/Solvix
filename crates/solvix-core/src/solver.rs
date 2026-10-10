//! Dirichlet constraints and Jacobi-preconditioned conjugate gradient.

use crate::assembly::Csr;
use crate::{Model, SolveError};

/// One constrained dof, kept for reaction output.
pub(crate) struct ConstrainedDof {
    pub(crate) row: usize,
    pub(crate) node_id: u32,
    pub(crate) dof: usize,
}

/// Expand boundary conditions to constrained dof rows. `ndof` is the
/// model's dof count per node (2 or 3).
pub(crate) fn constrained_dofs(
    model: &Model,
    ndof: usize,
) -> Result<Vec<ConstrainedDof>, SolveError> {
    let mut out = Vec::new();
    for bc in &model.boundary_conditions {
        for &nid in &bc.node_ids {
            let idx = model.mesh.node_index(nid).ok_or_else(|| {
                SolveError::InvalidModel(format!("boundary condition on unknown node {nid}"))
            })?;
            for &dof in &bc.dofs {
                if dof == 0 || dof > ndof {
                    return Err(SolveError::InvalidModel(format!(
                        "boundary condition on node {nid} has unsupported dof {dof}"
                    )));
                }
                out.push(ConstrainedDof {
                    row: idx * ndof + dof - 1,
                    node_id: nid,
                    dof,
                });
            }
        }
    }
    Ok(out)
}

/// Homogeneous Dirichlet: zero row and column, K[i,i] = 1, f[i] = 0.
/// Single O(nnz) pass: an entry dies if its row or column is constrained.
pub(crate) fn apply_dirichlet(k: &mut Csr, f: &mut [f64], constrained: &[ConstrainedDof]) {
    let n = f.len();
    let mut is_constrained = vec![false; n];
    for c in constrained {
        is_constrained[c.row] = true;
    }
    for r in 0..n {
        for idx in k.row_ptr[r]..k.row_ptr[r + 1] {
            let c = k.col_idx[idx];
            if is_constrained[r] || is_constrained[c] {
                k.vals[idx] = if r == c { 1.0 } else { 0.0 };
            }
        }
    }
    for c in constrained {
        f[c.row] = 0.0;
    }
}

/// Jacobi-preconditioned CG on an SPD system. Converges when the relative
/// residual ||r||/||f|| <= `tol`. Returns (u, iterations, final residual).
pub(crate) fn pcg(
    k: &Csr,
    f: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<(Vec<f64>, usize, f64), SolveError> {
    let n = f.len();
    let mut u = vec![0.0; n];
    let mut r = f.to_vec();
    let b_norm = f.iter().map(|v| v * v).sum::<f64>().sqrt();
    if b_norm == 0.0 {
        return Ok((u, 0, 0.0));
    }
    let inv_diag: Vec<f64> = (0..n)
        .map(|i| {
            let d = (k.row_ptr[i]..k.row_ptr[i + 1])
                .find(|&idx| k.col_idx[idx] == i)
                .map_or(0.0, |idx| k.vals[idx]);
            if d == 0.0 {
                0.0
            } else {
                1.0 / d
            }
        })
        .collect();
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let rel_res = |r: &[f64]| dot(r, r).sqrt() / b_norm;

    let mut z: Vec<f64> = r.iter().zip(&inv_diag).map(|(&x, &m)| x * m).collect();
    let mut p = z.clone();
    let mut rz = dot(&r, &z);

    for it in 1..=max_iter {
        let kp = k.mat_vec(&p);
        let pap = dot(&p, &kp);
        if pap <= 0.0 || pap.is_nan() {
            return Err(SolveError::NotConverged {
                iterations: it,
                residual: rel_res(&r),
            });
        }
        let alpha = rz / pap;
        for (ui, pi) in u.iter_mut().zip(&p) {
            *ui += alpha * pi;
        }
        for (ri, ki) in r.iter_mut().zip(&kp) {
            *ri -= alpha * ki;
        }
        let res = rel_res(&r);
        if res <= tol {
            return Ok((u, it, res));
        }
        for (zi, (ri, &m)) in z.iter_mut().zip(r.iter().zip(&inv_diag)) {
            *zi = ri * m;
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for (pi, zi) in p.iter_mut().zip(&z) {
            *pi = zi + beta * *pi;
        }
    }
    Err(SolveError::NotConverged {
        iterations: max_iter,
        residual: rel_res(&r),
    })
}
