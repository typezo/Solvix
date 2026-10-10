//! Post-processing: von Mises stress at element centroids and reactions.

use crate::assembly::{Csr, ElemInfo};
use crate::element::centroid_von_mises;
use crate::solver::ConstrainedDof;
use crate::Reaction;

/// Von Mises per element, parallel to `model.mesh.elements`.
pub(crate) fn von_mises_at_centroids(u: &[f64], elems: &[ElemInfo]) -> Vec<f64> {
    elems
        .iter()
        .map(|info| {
            let ue: Vec<f64> = info
                .node_idx
                .iter()
                .flat_map(|&n| [u[2 * n], u[2 * n + 1]])
                .collect();
            centroid_von_mises(&info.b, &info.d, &ue)
        })
        .collect()
}

/// Reactions R = K_orig * u - f_orig at the constrained dofs.
/// Must use the unmodified stiffness and load vector: after Dirichlet
/// modification those rows read 1*u_i - 0 = 0.
pub(crate) fn reactions(
    k_orig: &Csr,
    u: &[f64],
    f_orig: &[f64],
    constrained: &[ConstrainedDof],
) -> Vec<Reaction> {
    let ku = k_orig.mat_vec(u);
    constrained
        .iter()
        .map(|c| Reaction {
            node_id: c.node_id,
            dof: c.dof,
            value: ku[c.row] - f_orig[c.row],
        })
        .collect()
}
