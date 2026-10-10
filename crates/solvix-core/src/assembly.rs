//! Global stiffness assembly: element triplets to CSR, plus the load vector.
//! Dof map: `node_index * ndof + (dof - 1)` with ndof = 2 for plane-only
//! models, 3 once any solid element is present.

use std::thread;

use crate::element::{
    b31_stiffness, c3d4_stiffness, c3d8_stiffness, cps3_stiffness, cps4_stiffness, cps8_stiffness,
    d3_matrix, d_matrix, dpe_matrix, invalid, quad8_b, quad8_shapes, CentroidB, Constitutive, G3,
    W3,
};
use crate::{DistributedLoad, ElemKind, LoadState, Model, SolveError};

/// Worker count for the threaded hot paths: `available_parallelism`
/// capped at 8 (more workers only add contention on bandwidth-bound
/// SpMV loops).
pub(crate) fn worker_threads() -> usize {
    thread::available_parallelism()
        .map_or(1, |n| n.get())
        .min(8)
}

/// Below this element count the threaded assembly costs more in thread
/// startup than it saves.
const PARALLEL_ASSEMBLY_MIN_ELEMS: usize = 256;
/// Below this row count an SpMV is too small to amortize thread startup.
const PARALLEL_SPMV_MIN_ROWS: usize = 20_000;

/// Compressed sparse row matrix (symmetric pattern, sorted column indices).
#[derive(Clone)]
pub(crate) struct Csr {
    pub(crate) row_ptr: Vec<usize>,
    pub(crate) col_idx: Vec<usize>,
    pub(crate) vals: Vec<f64>,
}

impl Csr {
    pub(crate) fn mat_vec(&self, u: &[f64]) -> Vec<f64> {
        self.mat_vec_with(u, worker_threads())
    }

    /// Row-chunked SpMV: rows are independent and each output row keeps
    /// its exact serial summation order, so the result is bit-identical
    /// to the serial product for any chunking. Parallel only when the
    /// thread count and row count both clear their gates.
    pub(crate) fn mat_vec_with(&self, u: &[f64], threads: usize) -> Vec<f64> {
        let n = self.row_ptr.len() - 1;
        let mut out = vec![0.0; n];
        if threads > 1 && n > PARALLEL_SPMV_MIN_ROWS {
            let chunk = n.div_ceil(threads);
            thread::scope(|s| {
                for (ci, rows) in out.chunks_mut(chunk).enumerate() {
                    let start = ci * chunk;
                    s.spawn(move || {
                        for (i, o) in rows.iter_mut().enumerate() {
                            let r = start + i;
                            for idx in self.row_ptr[r]..self.row_ptr[r + 1] {
                                *o += self.vals[idx] * u[self.col_idx[idx]];
                            }
                        }
                    });
                }
            });
        } else {
            for (r, o) in out.iter_mut().enumerate() {
                for idx in self.row_ptr[r]..self.row_ptr[r + 1] {
                    *o += self.vals[idx] * u[self.col_idx[idx]];
                }
            }
        }
        out
    }
}

/// One (row, col, value) contribution destined for the CSR builder.
type Triplet = (usize, usize, f64);

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

/// One element's stiffness triplets (scattered into `trips`) plus the
/// data stress recovery needs. Extracted per element so the serial and
/// the threaded assembly run byte-identical element code.
fn element_entry(
    trips: &mut Vec<(usize, usize, f64)>,
    model: &Model,
    ei: usize,
    ndof: usize,
) -> Result<ElemInfo, SolveError> {
    let e = &model.mesh.elements[ei];
    let mat = model
        .material_for_element(ei)
        .ok_or_else(|| invalid(e.id, "no material assigned"))?;
    if mat.youngs_modulus <= 0.0 || mat.poisson_ratio.abs() >= 1.0 {
        return Err(invalid(e.id, "invalid material constants"));
    }
    let solid = matches!(e.kind, ElemKind::C3d8 | ElemKind::C3d4);
    let plane_strain = matches!(e.kind, ElemKind::Cpe4 | ElemKind::Cpe8);
    if (solid || plane_strain) && mat.poisson_ratio >= 0.5 {
        // 1-2*nu hits the incompressible limit; the plane check above
        // already rejects |nu| >= 1.
        return Err(invalid(
            e.id,
            "poisson ratio must be < 0.5 for solids and plane strain",
        ));
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
            scatter(trips, &node_idx, &ke, ndof);
            (b, Constitutive::Plane(d))
        }
        ElemKind::Cps4 | ElemKind::Cpe4 => {
            let coords = [
                coords2(model, &node_idx, 0),
                coords2(model, &node_idx, 1),
                coords2(model, &node_idx, 2),
                coords2(model, &node_idx, 3),
            ];
            // Same geometry/stiffness code; only the constitutive
            // matrix switches between the plane families.
            let strain = e.kind == ElemKind::Cpe4;
            let d = if strain {
                dpe_matrix(mat.youngs_modulus, mat.poisson_ratio)
            } else {
                d_matrix(mat.youngs_modulus, mat.poisson_ratio)
            };
            let (ke, b) = cps4_stiffness(&coords, &d, model.thickness).ok_or_else(|| {
                invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
            })?;
            scatter(trips, &node_idx, &ke, ndof);
            (
                b,
                if strain {
                    Constitutive::PlaneStrain(d)
                } else {
                    Constitutive::Plane(d)
                },
            )
        }
        ElemKind::Cps8 | ElemKind::Cpe8 => {
            let coords: [[f64; 2]; 8] = core::array::from_fn(|k| coords2(model, &node_idx, k));
            let strain = e.kind == ElemKind::Cpe8;
            let d = if strain {
                dpe_matrix(mat.youngs_modulus, mat.poisson_ratio)
            } else {
                d_matrix(mat.youngs_modulus, mat.poisson_ratio)
            };
            let (ke, b) = cps8_stiffness(&coords, &d, model.thickness).ok_or_else(|| {
                invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
            })?;
            scatter(trips, &node_idx, &ke, ndof);
            (
                b,
                if strain {
                    Constitutive::PlaneStrain(d)
                } else {
                    Constitutive::Plane(d)
                },
            )
        }
        ElemKind::B31 => {
            let coords = [coords2(model, &node_idx, 0), coords2(model, &node_idx, 1)];
            let (a, inertia) = model
                .beam_section_for_element(ei)
                .ok_or_else(|| invalid(e.id, "no beam section assigned (*BEAM SECTION)"))?;
            if a <= 0.0 || inertia <= 0.0 {
                return Err(invalid(e.id, "beam section needs positive A and I"));
            }
            let (ke, b) = b31_stiffness(&coords, mat.youngs_modulus, a, inertia)
                .ok_or_else(|| invalid(e.id, "zero-length beam"))?;
            scatter(trips, &node_idx, &ke, ndof);
            (
                b,
                Constitutive::Beam {
                    e: mat.youngs_modulus,
                    a,
                    i: inertia,
                },
            )
        }
        ElemKind::C3d8 => {
            let coords: [[f64; 3]; 8] =
                core::array::from_fn(|k| model.mesh.nodes[node_idx[k]].coords);
            let d = d3_matrix(mat.youngs_modulus, mat.poisson_ratio);
            let (ke, b) = c3d8_stiffness(&coords, &d).ok_or_else(|| {
                invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
            })?;
            scatter(trips, &node_idx, &ke, ndof);
            (b, Constitutive::Solid(Box::new(d)))
        }
        ElemKind::C3d4 => {
            let coords: [[f64; 3]; 4] =
                core::array::from_fn(|k| model.mesh.nodes[node_idx[k]].coords);
            let d = d3_matrix(mat.youngs_modulus, mat.poisson_ratio);
            let (ke, b) = c3d4_stiffness(&coords, &d)
                .ok_or_else(|| invalid(e.id, "degenerate (zero-volume) tetrahedron"))?;
            scatter(trips, &node_idx, &ke, ndof);
            (b, Constitutive::Solid(Box::new(d)))
        }
    };
    Ok(ElemInfo { node_idx, b, d })
}

/// Assemble K (CSR) and f. Assumes `Model::validate()` already passed and
/// `ndof = ndof_per_node(model)` (2 or 3). The load vector collects the
/// state's concentrated loads, edge pressures and gravity.
pub(crate) fn assemble(
    model: &Model,
    ndof: usize,
    loads: &LoadState,
) -> Result<(Csr, Vec<f64>, Vec<ElemInfo>), SolveError> {
    assemble_with(model, ndof, loads, worker_threads())
}

/// [`assemble`] with the worker-thread count pinned: `1` forces the
/// serial path (the determinism and speedup tests use that baseline).
pub(crate) fn assemble_with(
    model: &Model,
    ndof: usize,
    loads: &LoadState,
    threads: usize,
) -> Result<(Csr, Vec<f64>, Vec<ElemInfo>), SolveError> {
    let ndof_total = model.mesh.nodes.len() * ndof;
    let n_elem = model.mesh.elements.len();
    let mut trips = Vec::with_capacity(ndof_total + n_elem * 24);
    for i in 0..ndof_total {
        trips.push((i, i, 0.0)); // ensure diagonal entries exist
    }
    let mut f = vec![0.0; ndof_total];
    let mut elems = Vec::with_capacity(n_elem);

    // Threaded element pass: 2*threads CONTIGUOUS chunks interleaved
    // over the workers (coarser steal without a work-stealing runtime).
    // Chunks are contiguous element ranges merged in chunk order, so
    // the concatenated triplet stream equals the serial one exactly;
    // csr_from_triplets then sums per (r,c) in the serial order and the
    // assembled matrix is bit-identical to the serial assembly.
    let n_chunks = if threads > 1 && n_elem >= PARALLEL_ASSEMBLY_MIN_ELEMS {
        2 * threads
    } else {
        1
    };
    let chunk_len = n_elem.div_ceil(n_chunks);
    type ChunkOut = Result<(Vec<Triplet>, Vec<ElemInfo>), SolveError>;
    let run_chunk = |c: usize| -> ChunkOut {
        let a = c * chunk_len;
        let b = (a + chunk_len).min(n_elem);
        let mut ct = Vec::with_capacity((b - a) * 24);
        let mut ce = Vec::with_capacity(b - a);
        for ei in a..b {
            ce.push(element_entry(&mut ct, model, ei, ndof)?);
        }
        Ok((ct, ce))
    };
    if n_chunks == 1 {
        let (ct, ce) = run_chunk(0)?;
        trips.extend(ct);
        elems.extend(ce);
    } else {
        let per_thread = thread::scope(|s| {
            let run = &run_chunk;
            let handles: Vec<_> = (0..threads)
                .map(|i| {
                    s.spawn(move || {
                        let mut outs = Vec::new();
                        let mut a = i;
                        while a < n_chunks {
                            outs.push((a, run(a)));
                            a += threads;
                        }
                        outs
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        // Merge in chunk order: the first failure is the first bad
        // element in enumeration order, exactly like the serial loop.
        let mut chunks: Vec<_> = per_thread.into_iter().flatten().collect();
        chunks.sort_unstable_by_key(|&(c, _)| c);
        for (_, r) in chunks {
            let (ct, ce) = r?;
            trips.extend(ct);
            elems.extend(ce);
        }
    }

    for cl in loads.cloads {
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

    edge_pressure_forces(model, loads.dloads, ndof, &mut f)?;
    gravity_forces(model, loads.gravity, ndof, &mut f)?;

    Ok((csr_from_triplets(ndof_total, trips), f, elems))
}

/// Consistent nodal forces for edge pressure (`*DLOAD` Pn). Linear edge:
/// the total `magnitude * thickness * edge_length` splits 1/2 - 1/2 over
/// the two end nodes (quadratic edges are v0.4). With CCW connectivity
/// the outward normal is the edge tangent crossed with +z; positive
/// pressure acts INTO the element.
fn edge_pressure_forces(
    model: &Model,
    dloads: &[DistributedLoad],
    ndof: usize,
    f: &mut [f64],
) -> Result<(), SolveError> {
    for dl in dloads {
        let e = model
            .mesh
            .elements
            .iter()
            .find(|e| e.id == dl.element_id)
            .ok_or_else(|| {
                SolveError::InvalidModel(format!(
                    "distributed load on unknown element {}",
                    dl.element_id
                ))
            })?;
        let max_face = match e.kind {
            ElemKind::Cps3 => 3,
            ElemKind::Cps4 | ElemKind::Cpe4 => 4,
            _ => {
                return Err(SolveError::InvalidModel(format!(
                    "edge pressure on element {} ({:?}) unsupported (linear quad/tri edges only)",
                    e.id, e.kind
                )))
            }
        };
        if dl.face == 0 || dl.face > max_face {
            return Err(SolveError::InvalidModel(format!(
                "element {} ({:?}) has no face {}",
                e.id, e.kind, dl.face
            )));
        }
        let n = e.connectivity.len();
        let na = e.connectivity[dl.face - 1];
        let nb = e.connectivity[dl.face % n];
        let (ia, ib) = (
            model.mesh.node_index(na).unwrap(),
            model.mesh.node_index(nb).unwrap(),
        );
        let xy = |i: usize| [model.mesh.nodes[i].coords[0], model.mesh.nodes[i].coords[1]];
        let (pa, pb) = (xy(ia), xy(ib));
        let d = [pb[0] - pa[0], pb[1] - pa[1]];
        let len = (d[0] * d[0] + d[1] * d[1]).sqrt();
        if len == 0.0 {
            return Err(SolveError::InvalidModel(format!(
                "zero-length pressure edge on element {}",
                e.id
            )));
        }
        let n_in = [-d[1] / len, d[0] / len];
        let fmag = dl.magnitude * model.thickness * len / 2.0;
        for &i in &[ia, ib] {
            f[i * ndof] += fmag * n_in[0];
            f[i * ndof + 1] += fmag * n_in[1];
        }
    }
    Ok(())
}

/// Gravity body force: per element `f = rho * t * A * g * dir`, lumped by
/// the shape-function integral (equal quarters for the linear quads,
/// thirds for CPS3, 3x3 Gauss weights for the quadratic quads, equal
/// halves for beams). Plane elements and B31 beams only.
fn gravity_forces(
    model: &Model,
    gravity: Option<(f64, [f64; 3])>,
    ndof: usize,
    f: &mut [f64],
) -> Result<(), SolveError> {
    let Some((g, dir)) = gravity else {
        return Ok(());
    };
    let dir2 = [dir[0], dir[1]];
    if dir2 == [0.0, 0.0] {
        return Err(SolveError::InvalidModel(
            "gravity direction has no in-plane component".into(),
        ));
    }
    for (ei, e) in model.mesh.elements.iter().enumerate() {
        let idx: Vec<usize> = e
            .connectivity
            .iter()
            .map(|&nid| model.mesh.node_index(nid).unwrap())
            .collect();
        // Nodal force magnitudes (before the direction split).
        let per: Vec<f64> = match e.kind {
            ElemKind::Cps3 => {
                let p = [
                    coords2(model, &idx, 0),
                    coords2(model, &idx, 1),
                    coords2(model, &idx, 2),
                ];
                let cross = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                    - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]);
                vec![model.density * model.thickness * cross.abs() / 2.0 * g / 3.0; 3]
            }
            ElemKind::Cps4 | ElemKind::Cpe4 => {
                let p: Vec<[f64; 2]> = (0..4).map(|k| coords2(model, &idx, k)).collect();
                let mut shoelace = 0.0;
                for k in 0..4 {
                    let a = p[k];
                    let b = p[(k + 1) % 4];
                    shoelace += a[0] * b[1] - b[0] * a[1];
                }
                vec![model.density * model.thickness * shoelace.abs() / 2.0 * g / 4.0; 4]
            }
            ElemKind::Cps8 | ElemKind::Cpe8 => {
                // Consistent lumping: rho*t*g*sum(w*detJ*N_i) over 3x3.
                let coords: [[f64; 2]; 8] = core::array::from_fn(|k| coords2(model, &idx, k));
                let mut w_i = [0.0; 8];
                for (a, &xi) in [-G3, 0.0, G3].iter().enumerate() {
                    for (b, &eta) in [-G3, 0.0, G3].iter().enumerate() {
                        let (_, detj) = quad8_b(&coords, xi, eta).ok_or_else(|| {
                            invalid(e.id, "non-positive Jacobian determinant at a Gauss point")
                        })?;
                        let wgt = W3[a] * W3[b] * model.density * model.thickness * g;
                        for (wi, &ni) in w_i.iter_mut().zip(&quad8_shapes(xi, eta)) {
                            *wi += wgt * detj * ni;
                        }
                    }
                }
                w_i.to_vec()
            }
            ElemKind::B31 => {
                let (a, _) = model
                    .beam_section_for_element(ei)
                    .ok_or_else(|| invalid(e.id, "no beam section assigned (*BEAM SECTION)"))?;
                let p0 = model.mesh.nodes[idx[0]].coords;
                let p1 = model.mesh.nodes[idx[1]].coords;
                let l = ((p1[0] - p0[0]).powi(2) + (p1[1] - p0[1]).powi(2)).sqrt();
                vec![model.density * a * l * g / 2.0; 2]
            }
            _ => {
                return Err(SolveError::InvalidModel(
                    "gravity on solid elements unsupported (plane and beam only)".into(),
                ))
            }
        };
        for (&i, &p) in idx.iter().zip(&per) {
            f[i * ndof] += p * dir2[0];
            f[i * ndof + 1] += p * dir2[1];
        }
    }
    Ok(())
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

    // Structured nx*ny quad8 mesh of a width x height rectangle over a
    // (2nx+1)x(2ny+1) grid, SKIPPING the (odd, odd) points — those are
    // element centers, not serendipity nodes, and would sit unreferenced
    // (zero stiffness rows). Ids are sequential in row-major order of the
    // kept nodes; connectivity is corners CCW, then midsides 1-2, 2-3,
    // 3-4, 4-1.
    pub(crate) fn quad8_mesh(nx: usize, ny: usize, width: f64, height: f64) -> Mesh {
        let (gnx, gny) = (2 * nx + 1, 2 * ny + 1);
        let mut id_of = vec![0u32; gnx * gny];
        let mut nodes = Vec::new();
        let mut next = 0u32;
        for j in 0..gny {
            for i in 0..gnx {
                if i % 2 == 1 && j % 2 == 1 {
                    continue;
                }
                next += 1;
                id_of[j * gnx + i] = next;
                nodes.push(Node {
                    id: next,
                    coords: [
                        width * i as f64 / (2 * nx) as f64,
                        height * j as f64 / (2 * ny) as f64,
                        0.0,
                    ],
                });
            }
        }
        let nid = |i: usize, j: usize| id_of[j * gnx + i];
        let mut elements = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let c = (2 * i, 2 * j);
                elements.push(Element {
                    id: (j * nx + i + 1) as u32,
                    kind: ElemKind::Cps8,
                    connectivity: vec![
                        nid(c.0, c.1),
                        nid(c.0 + 2, c.1),
                        nid(c.0 + 2, c.1 + 2),
                        nid(c.0, c.1 + 2),
                        nid(c.0 + 1, c.1),
                        nid(c.0 + 2, c.1 + 1),
                        nid(c.0 + 1, c.1 + 2),
                        nid(c.0, c.1 + 1),
                    ],
                });
            }
        }
        Mesh { nodes, elements }
    }

    // n-element B31 chain along x over [0, l], node/element ids from 1.
    pub(crate) fn beam_mesh(n: usize, l: f64) -> Mesh {
        Mesh {
            nodes: (0..=n)
                .map(|i| Node {
                    id: (i + 1) as u32,
                    coords: [l * i as f64 / n as f64, 0.0, 0.0],
                })
                .collect(),
            elements: (1..=n)
                .map(|i| Element {
                    id: i as u32,
                    kind: ElemKind::B31,
                    connectivity: vec![i as u32, i as u32 + 1],
                })
                .collect(),
        }
    }

    pub(crate) fn model_with_mesh(mesh: Mesh) -> Model {
        let n_elem = mesh.elements.len();
        Model {
            steps: Vec::new(),
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
            beam_sections: Vec::new(),
            density: 1.0,
            boundary_conditions: vec![BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1, 2],
            }],
            cloads: vec![],
            ..Model::default()
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
        let (k, _, _) = assemble(&model, 2, &model.load_state()).unwrap();

        assert_symmetric(&k);
        let u: Vec<f64> = model.mesh.nodes.iter().flat_map(|_| [1.0, 0.5]).collect();
        let ku = k.mat_vec(&u);
        assert!(ku.iter().all(|&v| v.abs() < 1e-10));
    }

    #[test]
    fn free_free_solid_stiffness_is_symmetric_and_kills_translations() {
        let model = model_with_mesh(hex_mesh(2, 2, 2, 1.0, 1.0, 1.0));
        let (k, _, _) = assemble(&model, 3, &model.load_state()).unwrap();

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

    #[test]
    fn parallel_assembly_is_bit_identical_to_serial() {
        // 12^3 hexes (1728 elements) and 20x20 quad8s (400 elements):
        // both clear the 256-element parallel gate.
        let hex = model_with_mesh(hex_mesh(12, 12, 12, 1.0, 1.0, 1.0));
        let quad8 = model_with_mesh(quad8_mesh(20, 20, 1.0, 1.0));
        for (model, ndof) in [(&hex, 3usize), (&quad8, 2usize)] {
            let loads = model.load_state();
            let (ks, fs, es) = assemble_with(model, ndof, &loads, 1).unwrap();
            let (kp, fp, ep) = assemble_with(model, ndof, &loads, worker_threads()).unwrap();
            assert_eq!(ks.row_ptr, kp.row_ptr);
            assert_eq!(ks.col_idx, kp.col_idx);
            assert!(ks
                .vals
                .iter()
                .zip(&kp.vals)
                .all(|(a, b)| a.to_bits() == b.to_bits()));
            assert!(fs.iter().zip(&fp).all(|(a, b)| a.to_bits() == b.to_bits()));
            let serial_idx: Vec<_> = es.iter().map(|e| e.node_idx.clone()).collect();
            assert!(ep.iter().zip(serial_idx).all(|(e, i)| e.node_idx == i));
        }
    }

    #[test]
    fn parallel_spmv_is_bit_identical_to_serial() {
        // 40x12x12 hexes -> 20_787 rows, past the 20_000-row gate.
        let model = model_with_mesh(hex_mesh(40, 12, 12, 2.0, 0.1, 0.1));
        let loads = model.load_state();
        let (k, _, _) = assemble_with(&model, 3, &loads, 1).unwrap();
        let n = k.row_ptr.len() - 1;
        assert!(n > PARALLEL_SPMV_MIN_ROWS);
        // Deterministic LCG vector in [-0.5, 0.5).
        let mut u = vec![0.0; n];
        let mut s = 0x853c_49e6_748f_ea9bu64;
        for x in &mut u {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *x = ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5;
        }
        let a = k.mat_vec_with(&u, 1);
        let b = k.mat_vec_with(&u, worker_threads());
        assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
    }
}
