//! Global stiffness assembly: element triplets to CSR, plus the load vector.
//! Dof map: `node_index * ndof + (dof - 1)` with ndof = 2 for plane-only
//! models, 3 once any solid element is present.

use crate::element::{
    c3d4_stiffness, c3d8_stiffness, cps3_stiffness, cps4_stiffness, d3_matrix, d_matrix, invalid,
    CentroidB, Constitutive,
};
use crate::{ElemKind, Model, SolveError};

/// Compressed sparse row matrix (symmetric pattern, sorted column indices).
#[derive(Clone)]
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

/// Per-element data kept for stress recovery.
pub(crate) struct ElemInfo {
    pub(crate) node_idx: Vec<usize>,
    pub(crate) b: CentroidB,
    pub(crate) d: Constitutive,
}

/// Triplets to CSR: entries summed per (row, col), columns sorted.
/// Shared by the stiffness and the consistent-mass assembly.
pub(crate) fn csr_from_triplets(n: usize, trips: Vec<(usize, usize, f64)>) -> Csr {
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

/// Scatter an element matrix (stiffness or mass) into global triplets.
/// Local dof a maps to global `node_idx[a/ndof] * ndof + a%ndof`.
pub(crate) fn scatter<const N: usize>(
    trips: &mut Vec<(usize, usize, f64)>,
    node_idx: &[usize],
    ke: &[[f64; N]; N],
    ndof: usize,
) {
    for a in 0..N {
        let gr = node_idx[a / ndof] * ndof + a % ndof;
        for c in 0..N {
            let gc = node_idx[c / ndof] * ndof + c % ndof;
            trips.push((gr, gc, ke[a][c]));
        }
    }
}

/// Assemble K (CSR) and f. Assumes `Model::validate()` already passed and
/// `ndof = ndof_per_node(model)` (2 or 3).
pub(crate) fn assemble(
    model: &Model,
    ndof: usize,
) -> Result<(Csr, Vec<f64>, Vec<ElemInfo>), SolveError> {
    let ndof_total = model.mesh.nodes.len() * ndof;
    let mut trips = Vec::new();
    for i in 0..ndof_total {
        trips.push((i, i, 0.0)); // ensure diagonal entries exist
    }
    let mut f = vec![0.0; ndof_total];
    let mut elems = Vec::with_capacity(model.mesh.elements.len());

    for (ei, e) in model.mesh.elements.iter().enumerate() {
        let mat = model
            .material_for_element(ei)
            .ok_or_else(|| invalid(e.id, "no material assigned"))?;
        if mat.youngs_modulus <= 0.0 || mat.poisson_ratio.abs() >= 1.0 {
            return Err(invalid(e.id, "invalid material constants"));
        }
        let solid = matches!(e.kind, ElemKind::C3d8 | ElemKind::C3d4);
        if solid && mat.poisson_ratio >= 0.5 {
            // 1-2*nu hits the incompressible limit; the plane check above
            // already rejects |nu| >= 1.
            return Err(invalid(e.id, "poisson ratio must be < 0.5 for solids"));
        }
        let node_idx: Vec<usize> = e
            .connectivity
            .iter()
            .map(|&nid| model.mesh.node_index(nid).unwrap())
            .collect();
        let (b, d) = match e.kind {
            ElemKind::Cps3 => {
                let coords = [
                    coords2(model, &node_idx, 0),
                    coords2(model, &node_idx, 1),
                    coords2(model, &node_idx, 2),
                ];
                let d = d_matrix(mat.youngs_modulus, mat.poisson_ratio);
                let (ke, b) = cps3_stiffness(&coords, &d, model.thickness)
                    .ok_or_else(|| invalid(e.id, "degenerate (zero-area) triangle"))?;
                scatter(&mut trips, &node_idx, &ke, ndof);
                (b, Constitutive::Plane(d))
            }
            ElemKind::Cps4 => {
                let coords = [
                    coords2(model, &node_idx, 0),
                    coords2(model, &node_idx, 1),
                    coords2(model, &node_idx, 2),
                    coords2(model, &node_idx, 3),
                ];
                let d = d_matrix(mat.youngs_modulus, mat.poisson_ratio);
                let (ke, b) = cps4_stiffness(&coords, &d, model.thickness).ok_or_else(|| {
                    invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                })?;
                scatter(&mut trips, &node_idx, &ke, ndof);
                (b, Constitutive::Plane(d))
            }
            ElemKind::C3d8 => {
                let coords: [[f64; 3]; 8] =
                    core::array::from_fn(|k| model.mesh.nodes[node_idx[k]].coords);
                let d = d3_matrix(mat.youngs_modulus, mat.poisson_ratio);
                let (ke, b) = c3d8_stiffness(&coords, &d).ok_or_else(|| {
                    invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                })?;
                scatter(&mut trips, &node_idx, &ke, ndof);
                (b, Constitutive::Solid(Box::new(d)))
            }
            ElemKind::C3d4 => {
                let coords: [[f64; 3]; 4] =
                    core::array::from_fn(|k| model.mesh.nodes[node_idx[k]].coords);
                let d = d3_matrix(mat.youngs_modulus, mat.poisson_ratio);
                let (ke, b) = c3d4_stiffness(&coords, &d)
                    .ok_or_else(|| invalid(e.id, "degenerate (zero-volume) tetrahedron"))?;
                scatter(&mut trips, &node_idx, &ke, ndof);
                (b, Constitutive::Solid(Box::new(d)))
            }
        };
        elems.push(ElemInfo { node_idx, b, d });
    }

    for cl in &model.cloads {
        if cl.dof == 0 || cl.dof > ndof {
            return Err(SolveError::InvalidModel(format!(
                "load on node {} has unsupported dof {}",
                cl.node_id, cl.dof
            )));
        }
        let row = model
            .mesh
            .node_index(cl.node_id)
            .map(|i| i * ndof + cl.dof - 1)
            .ok_or_else(|| {
                SolveError::InvalidModel(format!("load on unknown node {}", cl.node_id))
            })?;
        f[row] += cl.value;
    }

    Ok((csr_from_triplets(ndof_total, trips), f, elems))
}

/// In-plane coordinates of the k-th element node (plane elements).
pub(crate) fn coords2(model: &Model, node_idx: &[usize], k: usize) -> [f64; 2] {
    let c = model.mesh.nodes[node_idx[k]].coords;
    [c[0], c[1]]
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

    // Structured nx*ny*nz hex mesh of an lx*ly*lz box, k-major node ids
    // starting at 1. C3D8 ordering: bottom face CCW, then top face CCW.
    pub(crate) fn hex_mesh(nx: usize, ny: usize, nz: usize, lx: f64, ly: f64, lz: f64) -> Mesh {
        let nid =
            |i: usize, j: usize, k: usize| (k * (ny + 1) * (nx + 1) + j * (nx + 1) + i + 1) as u32;
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1));
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    nodes.push(Node {
                        id: nid(i, j, k),
                        coords: [
                            lx * i as f64 / nx as f64,
                            ly * j as f64 / ny as f64,
                            lz * k as f64 / nz as f64,
                        ],
                    });
                }
            }
        }
        let mut elements = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    elements.push(Element {
                        id: (k * nx * ny + j * nx + i + 1) as u32,
                        kind: ElemKind::C3d8,
                        connectivity: vec![
                            nid(i, j, k),
                            nid(i + 1, j, k),
                            nid(i + 1, j + 1, k),
                            nid(i, j + 1, k),
                            nid(i, j, k + 1),
                            nid(i + 1, j, k + 1),
                            nid(i + 1, j + 1, k + 1),
                            nid(i, j + 1, k + 1),
                        ],
                    });
                }
            }
        }
        Mesh { nodes, elements }
    }

    pub(crate) fn model_with_mesh(mesh: Mesh) -> Model {
        let n_elem = mesh.elements.len();
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
                element_ids: (1..=n_elem as u32).collect(),
            }],
            mat_assignments: vec![("all".into(), "mat".into())],
            density: 1.0,
            boundary_conditions: vec![BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1, 2],
            }],
            cloads: vec![],
        }
    }

    fn assert_symmetric(k: &Csr) {
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
    }

    #[test]
    fn free_free_stiffness_is_symmetric_and_kills_translation() {
        let model = model_with_mesh(quad_mesh(2, 2, 1.0, 1.0));
        let (k, _, _) = assemble(&model, 2).unwrap();

        assert_symmetric(&k);
        let u: Vec<f64> = model.mesh.nodes.iter().flat_map(|_| [1.0, 0.5]).collect();
        let ku = k.mat_vec(&u);
        assert!(ku.iter().all(|&v| v.abs() < 1e-10));
    }

    #[test]
    fn free_free_solid_stiffness_is_symmetric_and_kills_translations() {
        let model = model_with_mesh(hex_mesh(2, 2, 2, 1.0, 1.0, 1.0));
        let (k, _, _) = assemble(&model, 3).unwrap();

        assert_symmetric(&k);
        let u: Vec<f64> = model
            .mesh
            .nodes
            .iter()
            .flat_map(|_| [1.0, 0.5, -0.25])
            .collect();
        let ku = k.mat_vec(&u);
        assert!(ku.iter().all(|&v| v.abs() < 1e-10));
    }
}
