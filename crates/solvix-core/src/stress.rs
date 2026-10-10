//! Post-processing: von Mises stress at element centroids, nodal stress
//! recovery (Gauss extrapolation + nodal averaging), and reactions.

use crate::assembly::{Csr, ElemInfo};
use crate::element::{centroid_stress, nodal_stresses, von_mises};
use crate::solver::ConstrainedDof;
use crate::Reaction;

/// Element nodal displacements gathered from the global vector with the
/// model's dof layout (`ndof` dofs per node).
fn gather_ue(u: &[f64], node_idx: &[usize], ndof: usize) -> Vec<f64> {
    node_idx
        .iter()
        .flat_map(|&n| (0..ndof).map(move |d| u[n * ndof + d]))
        .collect()
}

/// Von Mises per element, parallel to `model.mesh.elements`. With
/// `thermal`, each element's recovered stress has its (constant)
/// thermal stress vector subtracted first: sigma = D(B u) - sigma_th,
/// so the value includes thermal strain.
pub(crate) fn von_mises_at_centroids(
    u: &[f64],
    ndof: usize,
    elems: &[ElemInfo],
    thermal: Option<&[[f64; 6]]>,
) -> Vec<f64> {
    elems
        .iter()
        .enumerate()
        .map(|(ei, info)| {
            let ue = gather_ue(u, &info.node_idx, ndof);
            let mut s = centroid_stress(&info.b, &info.d, &ue);
            if let Some(th) = thermal {
                for (a, b) in s.iter_mut().zip(&th[ei]) {
                    *a -= b;
                }
            }
            von_mises(&s)
        })
        .collect()
}

/// Von Mises per node, parallel to `model.mesh.nodes`. Per element, the
/// stress components at the Gauss points are extrapolated to the corner
/// nodes (bilinear for quads, trilinear for hexes; Abaqus-style) or taken
/// constant (triangles, tets), then averaged over contributing elements.
/// With `thermal`, each element's thermal stress vector is subtracted
/// from its own contributions BEFORE the nodal average (it is constant
/// per element but differs between neighbors).
/// v1 choice: simple unweighted average per node; volume weighting of
/// contributions is a possible later refinement.
pub(crate) fn nodal_von_mises(
    u: &[f64],
    ndof: usize,
    elems: &[ElemInfo],
    n_nodes: usize,
    thermal: Option<&[[f64; 6]]>,
) -> Vec<f64> {
    let mut sum = vec![[0.0f64; 6]; n_nodes];
    let mut count = vec![0u32; n_nodes];
    for (ei, info) in elems.iter().enumerate() {
        let ue = gather_ue(u, &info.node_idx, ndof);
        for (i, mut s) in nodal_stresses(&info.b, &info.d, &ue)
            .into_iter()
            .enumerate()
        {
            if let Some(th) = thermal {
                for (a, b) in s.iter_mut().zip(&th[ei]) {
                    *a -= b;
                }
            }
            let node = info.node_idx[i];
            for k in 0..6 {
                sum[node][k] += s[k];
            }
            count[node] += 1;
        }
    }
    (0..n_nodes)
        .map(|i| {
            if count[i] == 0 {
                0.0 // node referenced by no element
            } else {
                let c = count[i] as f64;
                let avg = [
                    sum[i][0] / c,
                    sum[i][1] / c,
                    sum[i][2] / c,
                    sum[i][3] / c,
                    sum[i][4] / c,
                    sum[i][5] / c,
                ];
                von_mises(&avg)
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
