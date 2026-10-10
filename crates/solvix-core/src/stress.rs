//! Post-processing: von Mises stress at element centroids, nodal stress
//! recovery (Gauss extrapolation + nodal averaging), and reactions.

use crate::assembly::{Csr, ElemInfo};
use crate::element::{
    centroid_von_mises, extrapolate_to_corners, stress_components, von_mises, CentroidB,
};
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

/// Von Mises per node, parallel to `model.mesh.nodes`. Per element, the
/// stress components at the 2x2 Gauss points are extrapolated to the
/// corner nodes with bilinear shape functions (CPS4) or taken constant
/// (CPS3), then averaged over all contributing elements.
/// v1 choice: simple unweighted average per node; area/detJ weighting of
/// contributions is a possible later refinement.
pub(crate) fn nodal_von_mises(u: &[f64], elems: &[ElemInfo], n_nodes: usize) -> Vec<f64> {
    let mut sum = vec![[0.0f64; 3]; n_nodes];
    let mut count = vec![0u32; n_nodes];
    for info in elems {
        let ue: Vec<f64> = info
            .node_idx
            .iter()
            .flat_map(|&n| [u[2 * n], u[2 * n + 1]])
            .collect();
        let mut scatter = |node: usize, s: [f64; 3]| {
            for k in 0..3 {
                sum[node][k] += s[k];
            }
            count[node] += 1;
        };
        match &info.b {
            CentroidB::Tri(bm) => {
                let s = stress_components(bm, &info.d, &ue);
                for &n in &info.node_idx {
                    scatter(n, s);
                }
            }
            CentroidB::Quad { gauss, .. } => {
                let mut gp = [[0.0f64; 3]; 4];
                for (j, bm) in gauss.iter().enumerate() {
                    gp[j] = stress_components(bm, &info.d, &ue);
                }
                for (i, s) in extrapolate_to_corners(&gp).into_iter().enumerate() {
                    scatter(info.node_idx[i], s);
                }
            }
        }
    }
    (0..n_nodes)
        .map(|i| {
            if count[i] == 0 {
                0.0 // node referenced by no element
            } else {
                let c = count[i] as f64;
                von_mises(&[sum[i][0] / c, sum[i][1] / c, sum[i][2] / c])
            }
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
