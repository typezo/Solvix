//! Global stiffness assembly: element triplets to CSR, plus the load vector.
//! Dof map: `node_index * 2 + (dof - 1)`.

use crate::element::{cps3_stiffness, cps4_stiffness, d_matrix, invalid, CentroidB};
use crate::{Model, SolveError};

/// Compressed sparse row matrix (symmetric pattern, sorted column indices).
pub(crate) struct Csr {
    pub(crate) row_ptr: Vec<usize>,
    pub(crate) col_idx: Vec<usize>,
    pub(crate) vals: Vec<f64>,
}

impl Csr {
    pub(crate) fn mat_vec(&self, u: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.row_ptr.len() - 1];
        for (r, o) in out.iter_mut().enumerate() {
            for idx in self.row_ptr[r]..self.row_ptr[r + 1] {
                *o += self.vals[idx] * u[self.col_idx[idx]];
            }
        }
        out
    }
}

/// Per-element data kept for centroid stress recovery.
pub(crate) struct ElemInfo {
    pub(crate) node_idx: Vec<usize>,
    pub(crate) b: CentroidB,
    pub(crate) d: [[f64; 3]; 3],
}

fn csr_from_triplets(n: usize, trips: Vec<(usize, usize, f64)>) -> Csr {
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (r, c, v) in trips {
        rows[r].push((c, v));
    }
    let mut row_ptr = Vec::with_capacity(n + 1);
    row_ptr.push(0);
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    for row in &mut rows {
        row.sort_unstable_by_key(|&(c, _)| c);
        let mut i = 0;
        while i < row.len() {
            let c = row[i].0;
            let mut v = 0.0;
            while i < row.len() && row[i].0 == c {
                v += row[i].1;
                i += 1;
            }
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

/// Scatter an element stiffness into global triplets. Local dof a maps to
/// global `node_idx[a/2] * 2 + a%2`.
fn scatter<const N: usize>(
    trips: &mut Vec<(usize, usize, f64)>,
    node_idx: &[usize],
    ke: &[[f64; N]; N],
) {
    for a in 0..N {
        let gr = node_idx[a / 2] * 2 + a % 2;
        for c in 0..N {
            let gc = node_idx[c / 2] * 2 + c % 2;
            trips.push((gr, gc, ke[a][c]));
        }
    }
}

/// Assemble K (CSR) and f. Assumes `Model::validate()` already passed.
pub(crate) fn assemble(model: &Model) -> Result<(Csr, Vec<f64>, Vec<ElemInfo>), SolveError> {
    let ndof = model.mesh.nodes.len() * 2;
    let mut trips = Vec::new();
    for i in 0..ndof {
        trips.push((i, i, 0.0)); // ensure diagonal entries exist
    }
    let mut f = vec![0.0; ndof];
    let mut elems = Vec::with_capacity(model.mesh.elements.len());

    for (ei, e) in model.mesh.elements.iter().enumerate() {
        let mat = model
            .material_for_element(ei)
            .ok_or_else(|| invalid(e.id, "no material assigned"))?;
        if mat.youngs_modulus <= 0.0 || mat.poisson_ratio.abs() >= 1.0 {
            return Err(invalid(e.id, "invalid material constants"));
        }
        let d = d_matrix(mat.youngs_modulus, mat.poisson_ratio);
        let node_idx: Vec<usize> = e
            .connectivity
            .iter()
            .map(|&nid| model.mesh.node_index(nid).unwrap())
            .collect();
        let coords_of = |k: usize| {
            let c = model.mesh.nodes[node_idx[k]].coords;
            [c[0], c[1]]
        };
        let b = match e.kind {
            crate::ElemKind::Cps3 => {
                let coords = [coords_of(0), coords_of(1), coords_of(2)];
                let (ke, b) = cps3_stiffness(&coords, &d, model.thickness)
                    .ok_or_else(|| invalid(e.id, "degenerate (zero-area) triangle"))?;
                scatter(&mut trips, &node_idx, &ke);
                b
            }
            crate::ElemKind::Cps4 => {
                let coords = [coords_of(0), coords_of(1), coords_of(2), coords_of(3)];
                let (ke, b) = cps4_stiffness(&coords, &d, model.thickness).ok_or_else(|| {
                    invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                })?;
                scatter(&mut trips, &node_idx, &ke);
                b
            }
        };
        elems.push(ElemInfo { node_idx, b, d });
    }

    for cl in &model.cloads {
        if cl.dof != 1 && cl.dof != 2 {
            return Err(SolveError::InvalidModel(format!(
                "load on node {} has unsupported dof {}",
                cl.node_id, cl.dof
            )));
        }
        let row = model
            .mesh
            .node_index(cl.node_id)
            .map(|i| i * 2 + cl.dof - 1)
            .ok_or_else(|| {
                SolveError::InvalidModel(format!("load on unknown node {}", cl.node_id))
            })?;
        f[row] += cl.value;
    }

    Ok((csr_from_triplets(ndof, trips), f, elems))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::{BoundaryCondition, ElemKind, Element, Material, Mesh, Node};

    // Structured nx*ny quad mesh of a width x height rectangle, row-major
    // node ids starting at 1.
    pub(crate) fn quad_mesh(nx: usize, ny: usize, width: f64, height: f64) -> Mesh {
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push(Node {
                    id: (j * (nx + 1) + i + 1) as u32,
                    coords: [
                        width * i as f64 / nx as f64,
                        height * j as f64 / ny as f64,
                        0.0,
                    ],
                });
            }
        }
        let mut elements = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let nid = |ii: usize, jj: usize| (jj * (nx + 1) + ii + 1) as u32;
                elements.push(Element {
                    id: (j * nx + i + 1) as u32,
                    kind: ElemKind::Cps4,
                    connectivity: vec![nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)],
                });
            }
        }
        Mesh { nodes, elements }
    }

    fn model_with_mesh(mesh: Mesh) -> Model {
        Model {
            mesh,
            thickness: 1.0,
            materials: vec![Material {
                name: "mat".into(),
                youngs_modulus: 1.0,
                poisson_ratio: 0.25,
            }],
            elsets: vec![crate::ElementSet {
                name: "all".into(),
                element_ids: (1..=4).collect(),
            }],
            mat_assignments: vec![("all".into(), "mat".into())],
            boundary_conditions: vec![BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1, 2],
            }],
            cloads: vec![],
        }
    }

    #[test]
    fn free_free_stiffness_is_symmetric_and_kills_translation() {
        let model = model_with_mesh(quad_mesh(2, 2, 1.0, 1.0));
        let (k, _, _) = assemble(&model).unwrap();

        let entry = |i: usize, j: usize| -> f64 {
            (k.row_ptr[i]..k.row_ptr[i + 1])
                .find(|&idx| k.col_idx[idx] == j)
                .map_or(0.0, |idx| k.vals[idx])
        };
        let n = k.row_ptr.len() - 1;
        for i in 0..n {
            for idx in k.row_ptr[i]..k.row_ptr[i + 1] {
                let j = k.col_idx[idx];
                assert!(
                    (k.vals[idx] - entry(j, i)).abs() < 1e-12,
                    "K[{i},{j}] != K[{j},{i}]"
                );
            }
        }
        let u: Vec<f64> = model.mesh.nodes.iter().flat_map(|_| [1.0, 0.5]).collect();
        let ku = k.mat_vec(&u);
        assert!(ku.iter().all(|&v| v.abs() < 1e-10));
    }
}
