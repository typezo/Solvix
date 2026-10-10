//! # solvix-core
//!
//! FEM kernel for Solvix: mesh data model, linear-elastic elements,
//! global assembly and the linear solver.
//!
//! This file defines the **public type contract** for the whole workspace.
//! `solvix-io` parses Abaqus decks into `Model`; `solvix-cli` orchestrates.
//!
//! Implementation modules: [`element`] (CPS3/CPS4/C3D8/C3D4 stiffness),
//! [`assembly`] (CSR assembly), [`solver`] (Dirichlet + preconditioned
//! CG), [`stress`] (von Mises + reactions).

mod assembly;
mod element;
mod solver;
mod stress;

/// A mesh node. Plane elements use x, y; solid elements use all three.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node {
    pub id: u32,
    pub coords: [f64; 3],
}

/// Supported element kinds: 2D plane stress (2 dofs/node) and 3D linear
/// solids (3 dofs/node). Plane and solid kinds cannot be mixed in one
/// model (rejected with `InvalidModel("mixed dimensionality")`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElemKind {
    /// 3-node linear triangle, plane stress.
    Cps3,
    /// 4-node bilinear quadrilateral, plane stress.
    Cps4,
    /// 8-node trilinear hexahedron, 2x2x2 Gauss.
    C3d8,
    /// 4-node constant-strain tetrahedron.
    C3d4,
}

impl ElemKind {
    /// Number of nodes per element.
    pub fn node_count(&self) -> usize {
        match self {
            ElemKind::Cps3 => 3,
            ElemKind::Cps4 => 4,
            ElemKind::C3d8 => 8,
            ElemKind::C3d4 => 4,
        }
    }
}

/// A finite element. `connectivity` holds **node ids** (not indices),
/// ordered counter-clockwise in the element's own numbering.
#[derive(Debug, Clone, PartialEq)]
pub struct Element {
    pub id: u32,
    pub kind: ElemKind,
    pub connectivity: Vec<u32>,
}

/// The mesh: nodes and elements, each carrying its own id.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mesh {
    pub nodes: Vec<Node>,
    pub elements: Vec<Element>,
}

impl Mesh {
    /// Position of `id` in `nodes`, if present.
    pub fn node_index(&self, id: u32) -> Option<usize> {
        self.nodes.iter().position(|n| n.id == id)
    }

    /// Position of `id` in `elements`, if present.
    pub fn element_index(&self, id: u32) -> Option<usize> {
        self.elements.iter().position(|e| e.id == id)
    }
}

/// Isotropic linear-elastic material (SI: Pa).
#[derive(Debug, Clone, PartialEq)]
pub struct Material {
    pub name: String,
    pub youngs_modulus: f64,
    pub poisson_ratio: f64,
}

/// A named set of element ids (`*ELSET` / `*ELSET, ELSET=name`).
#[derive(Debug, Clone, PartialEq)]
pub struct ElementSet {
    pub name: String,
    pub element_ids: Vec<u32>,
}

/// Homogeneous boundary condition on a list of nodes.
/// `dofs` uses 1-based local dof numbers: 1 = x, 2 = y, 3 = z
/// (3D models only).
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryCondition {
    pub node_ids: Vec<u32>,
    pub dofs: Vec<usize>,
}

/// Concentrated load (`*CLOAD`), already expanded to single nodes.
/// `dof` is 1-based: 1 = x, 2 = y, 3 = z (3D models only). Value in Newtons.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConcentratedLoad {
    pub node_id: u32,
    pub dof: usize,
    pub value: f64,
}

/// A complete, solvable v0.1 model: one static step.
#[derive(Debug, Clone, Default)]
pub struct Model {
    pub mesh: Mesh,
    /// Out-of-plane thickness for plane-stress elements (default 1.0).
    pub thickness: f64,
    pub materials: Vec<Material>,
    pub elsets: Vec<ElementSet>,
    /// (elset name, material name) pairs from `*SOLID SECTION`.
    pub mat_assignments: Vec<(String, String)>,
    pub boundary_conditions: Vec<BoundaryCondition>,
    pub cloads: Vec<ConcentratedLoad>,
}

impl Model {
    /// Material assigned to the element at `elem_index`, if any.
    /// Resolution: elset containing the element id -> assigned material.
    pub fn material_for_element(&self, elem_index: usize) -> Option<&Material> {
        let elem = self.mesh.elements.get(elem_index)?;
        for (elset, mat) in &self.mat_assignments {
            if let Some(set) = self.elsets.iter().find(|s| &s.name == elset) {
                if set.element_ids.contains(&elem.id) {
                    return self.materials.iter().find(|m| &m.name == mat);
                }
            }
        }
        None
    }

    /// Structural checks shared by parser tests and the solver.
    pub fn validate(&self) -> Result<(), ModelError> {
        if self.mesh.nodes.is_empty() || self.mesh.elements.is_empty() {
            return Err(ModelError::EmptyMesh);
        }
        if self.materials.is_empty() {
            return Err(ModelError::NoMaterial);
        }
        if self.boundary_conditions.is_empty() {
            return Err(ModelError::NoBoundaryConditions);
        }
        for e in &self.mesh.elements {
            if e.connectivity.len() != e.kind.node_count() {
                return Err(ModelError::MalformedElement(e.id));
            }
            for &n in &e.connectivity {
                if self.mesh.node_index(n).is_none() {
                    return Err(ModelError::UnknownNode(n));
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ModelError {
    EmptyMesh,
    NoMaterial,
    NoBoundaryConditions,
    MalformedElement(u32),
    UnknownNode(u32),
}

/// Reaction force at a constrained dof (Newtons).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Reaction {
    pub node_id: u32,
    pub dof: usize,
    pub value: f64,
}

/// Solver output. Arrays are parallel to `mesh.nodes` / `mesh.elements`.
#[derive(Debug, Clone, Default)]
pub struct Solution {
    /// Displacement per node; z component is 0.0 for plane models.
    pub displacements: Vec<[f64; 3]>,
    /// Von Mises stress at each element centroid (Pa).
    pub von_mises: Vec<f64>,
    /// Von Mises stress per node (Pa): element Gauss-point stresses
    /// extrapolated to nodes and averaged over contributing elements.
    /// Parallel to `mesh.nodes`.
    pub nodal_von_mises: Vec<f64>,
    pub reactions: Vec<Reaction>,
    pub iterations: usize,
    pub final_residual: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SolveError {
    EmptyMesh,
    NoConstraints,
    InvalidModel(String),
    NotConverged { iterations: usize, residual: f64 },
}

/// Solve the linear-static problem `K u = f` with homogeneous Dirichlet
/// constraints. Conjugate gradient on CSR with Jacobi preconditioning.
/// Plane-only models solve with 2 dofs/node, models containing solid
/// elements with 3 (see [`ndof_per_node`]).
///
/// Relative tolerance 1e-10, `max_iter = max(1000, 10 * ndof)`.
pub fn solve(model: &Model) -> Result<Solution, SolveError> {
    // Free-body check first: validate() reports NoBoundaryConditions as a
    // ModelError, but the contract wants SolveError::NoConstraints.
    if model.boundary_conditions.is_empty() {
        return Err(SolveError::NoConstraints);
    }
    if let Err(e) = model.validate() {
        return Err(SolveError::InvalidModel(format!("{e:?}")));
    }
    if model.thickness <= 0.0 || model.thickness.is_nan() {
        return Err(SolveError::InvalidModel(
            "thickness must be positive".into(),
        ));
    }

    let ndof = ndof_per_node(model)?;
    let (mut k, mut f, elems) = assembly::assemble(model, ndof)?;
    let constrained = solver::constrained_dofs(model, ndof)?;

    // Keep the unmodified system for reaction recovery.
    let k_orig = assembly::Csr {
        row_ptr: k.row_ptr.clone(),
        col_idx: k.col_idx.clone(),
        vals: k.vals.clone(),
    };
    let f_orig = f.clone();

    solver::apply_dirichlet(&mut k, &mut f, &constrained);
    let ndof_total = f.len();
    let max_iter = 1000.max(10 * ndof_total);
    let (u, iterations, final_residual) = solver::pcg(&k, &f, 1e-10, max_iter)?;

    let displacements = model
        .mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let b = i * ndof;
            [u[b], u[b + 1], if ndof == 3 { u[b + 2] } else { 0.0 }]
        })
        .collect();
    let von_mises = stress::von_mises_at_centroids(&u, ndof, &elems);
    let nodal_von_mises = stress::nodal_von_mises(&u, ndof, &elems, model.mesh.nodes.len());
    let reactions = stress::reactions(&k_orig, &u, &f_orig, &constrained);

    Ok(Solution {
        displacements,
        von_mises,
        nodal_von_mises,
        reactions,
        iterations,
        final_residual,
    })
}

/// Dofs per node for the model: 2 while every element is a plane kind,
/// 3 once any C3D8/C3D4 solid is present. Mixing plane and solid kinds
/// in one model is rejected: there is no formulation that shares a
/// stiffness matrix between them.
pub(crate) fn ndof_per_node(model: &Model) -> Result<usize, SolveError> {
    let mut plane = false;
    let mut solid = false;
    for e in &model.mesh.elements {
        match e.kind {
            ElemKind::Cps3 | ElemKind::Cps4 => plane = true,
            ElemKind::C3d8 | ElemKind::C3d4 => solid = true,
        }
    }
    match (plane, solid) {
        (true, true) => Err(SolveError::InvalidModel("mixed dimensionality".into())),
        (false, true) => Ok(3),
        (true, false) => Ok(2),
        (false, false) => Ok(2), // unreachable after validate() (EmptyMesh)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contract_compiles_and_mesh_lookup_works() {
        let mesh = Mesh {
            nodes: vec![Node {
                id: 7,
                coords: [0.0, 0.0, 0.0],
            }],
            elements: vec![Element {
                id: 1,
                kind: ElemKind::Cps4,
                connectivity: vec![7, 8, 9, 10],
            }],
        };
        let model = Model {
            mesh,
            thickness: 1.0,
            materials: vec![Material {
                name: "steel".into(),
                youngs_modulus: 210e9,
                poisson_ratio: 0.3,
            }],
            ..Default::default()
        };
        let model = Model {
            boundary_conditions: vec![BoundaryCondition {
                node_ids: vec![7],
                dofs: vec![1, 2],
            }],
            ..model
        };
        assert_eq!(model.mesh.node_index(7), Some(0));
        assert_eq!(model.mesh.element_index(1), Some(0));
        assert!(matches!(model.validate(), Err(ModelError::UnknownNode(8))));
    }

    const E: f64 = 210e9;
    const NU: f64 = 0.3;

    fn steel() -> Material {
        Material {
            name: "steel".into(),
            youngs_modulus: E,
            poisson_ratio: NU,
        }
    }

    fn quad_model(nx: usize, ny: usize, width: f64, height: f64) -> Model {
        let n_elem = (nx * ny) as u32;
        Model {
            mesh: assembly::tests::quad_mesh(nx, ny, width, height),
            thickness: 1.0,
            materials: vec![steel()],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=n_elem).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            ..Default::default()
        }
    }

    fn solid_model(mesh: Mesh) -> Model {
        let n_elem = mesh.elements.len() as u32;
        Model {
            mesh,
            thickness: 1.0, // unused by solids
            materials: vec![steel()],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=n_elem).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            ..Default::default()
        }
    }

    /// Uniform sigma_z tension on a unit cube: w fixed on the bottom face,
    /// ux at node 1 and uy at nodes 1-2 (y=0) to pin rigid sliding without
    /// fighting the Poisson contraction.
    fn cube_tension_bcs(model: &mut Model) {
        model.boundary_conditions = vec![
            BoundaryCondition {
                node_ids: vec![1, 2, 3, 4],
                dofs: vec![3],
            },
            BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1],
            },
            BoundaryCondition {
                node_ids: vec![1, 2],
                dofs: vec![2],
            },
        ];
        let sigma = 1e6;
        model.cloads = [5, 6, 7, 8]
            .iter()
            .map(|&nid| ConcentratedLoad {
                node_id: nid,
                dof: 3,
                value: sigma / 4.0,
            })
            .collect();
    }

    fn assert_uniform_tension(model: &Model, sol: &Solution) {
        let sigma = 1e6;
        for (i, &vm) in sol.von_mises.iter().enumerate() {
            assert!((vm - sigma).abs() / sigma < 1e-6, "element {i}: vM={vm}");
        }
        assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
        for (i, &nvm) in sol.nodal_von_mises.iter().enumerate() {
            assert!((nvm - sigma).abs() / sigma < 1e-6, "node {i}: vM={nvm}");
        }
        let w = sigma / E;
        for &nid in &[5, 6, 7, 8] {
            let i = model.mesh.node_index(nid).unwrap();
            assert!(
                (sol.displacements[i][2] - w).abs() / w < 1e-6,
                "node {nid}: w={}",
                sol.displacements[i][2]
            );
        }
        let rz: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 3)
            .map(|r| r.value)
            .sum();
        assert!((rz + sigma).abs() < 1e-6 * sigma);
    }

    fn fix_left_edge_x(model: &mut Model, nx: usize, ny: usize) {
        // ux on the left edge, uy pinned at one node only: clamping uy
        // everywhere would forbid the Poisson contraction and spoil the
        // uniform-stress patch solution.
        let ids: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        model.boundary_conditions.push(BoundaryCondition {
            node_ids: ids,
            dofs: vec![1],
        });
        model.boundary_conditions.push(BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![2],
        });
    }

    #[test]
    fn cps4_patch_single_quad() {
        let mut model = quad_model(1, 1, 1.0, 1.0);
        fix_left_edge_x(&mut model, 1, 1);
        let sigma = 1e6;
        model.cloads = vec![
            ConcentratedLoad {
                node_id: 2,
                dof: 1,
                value: sigma / 2.0,
            },
            ConcentratedLoad {
                node_id: 4,
                dof: 1,
                value: sigma / 2.0,
            },
        ];
        let sol = solve(&model).unwrap();
        // Bilinear quad reproduces constant strain exactly.
        assert!((sol.von_mises[0] - sigma).abs() / sigma < 1e-6);
        // Uniform field: nodal recovery is exact at every node and agrees
        // with the centroid value (extrapolation patch test).
        assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
        for (i, &nvm) in sol.nodal_von_mises.iter().enumerate() {
            assert!((nvm - sigma).abs() / sigma < 1e-6, "node {i}: vM={nvm}");
            assert!((nvm - sol.von_mises[0]).abs() / sigma < 1e-6);
        }
        let ux_right = sigma / E;
        assert!((sol.displacements[1][0] - ux_right).abs() / ux_right < 1e-6);
        assert!((sol.displacements[3][0] - ux_right).abs() / ux_right < 1e-6);
        // Total reaction balances the applied load.
        let rx: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 1)
            .map(|r| r.value)
            .sum();
        assert!((rx + sigma).abs() < 1e-6 * sigma);
        assert!(sol.final_residual <= 1e-10);
    }

    #[test]
    fn cps4_patch_2x2() {
        let mut model = quad_model(2, 2, 1.0, 1.0);
        fix_left_edge_x(&mut model, 2, 2);
        let sigma = 1e6;
        // Consistent nodal loads for uniform traction on the right edge.
        let s = sigma / 4.0;
        model.cloads = vec![
            ConcentratedLoad {
                node_id: 3,
                dof: 1,
                value: s,
            },
            ConcentratedLoad {
                node_id: 6,
                dof: 1,
                value: 2.0 * s,
            },
            ConcentratedLoad {
                node_id: 9,
                dof: 1,
                value: s,
            },
        ];
        let sol = solve(&model).unwrap();
        for (i, &vm) in sol.von_mises.iter().enumerate() {
            assert!((vm - sigma).abs() / sigma < 2e-2, "element {i}: vM={vm}");
        }
        // Shared (interior) nodes average several contributions and stay
        // uniform: nodal == centroid == exact within the same tolerance.
        assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
        for (i, &nvm) in sol.nodal_von_mises.iter().enumerate() {
            assert!((nvm - sigma).abs() / sigma < 2e-2, "node {i}: vM={nvm}");
        }
        let rx: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 1)
            .map(|r| r.value)
            .sum();
        assert!((rx + sigma).abs() < 2e-2 * sigma);
    }

    #[test]
    fn cps3_patch_two_triangles() {
        // Unit square split diagonally; uniform tension is exact for CST.
        let model = Model {
            mesh: Mesh {
                nodes: vec![
                    Node {
                        id: 1,
                        coords: [0.0, 0.0, 0.0],
                    },
                    Node {
                        id: 2,
                        coords: [1.0, 0.0, 0.0],
                    },
                    Node {
                        id: 3,
                        coords: [1.0, 1.0, 0.0],
                    },
                    Node {
                        id: 4,
                        coords: [0.0, 1.0, 0.0],
                    },
                ],
                elements: vec![
                    Element {
                        id: 1,
                        kind: ElemKind::Cps3,
                        connectivity: vec![1, 2, 3],
                    },
                    Element {
                        id: 2,
                        kind: ElemKind::Cps3,
                        connectivity: vec![1, 3, 4],
                    },
                ],
            },
            thickness: 1.0,
            materials: vec![steel()],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: vec![1, 2],
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            boundary_conditions: vec![
                BoundaryCondition {
                    node_ids: vec![1, 4],
                    dofs: vec![1],
                },
                BoundaryCondition {
                    node_ids: vec![1],
                    dofs: vec![2],
                },
            ],
            cloads: vec![
                ConcentratedLoad {
                    node_id: 2,
                    dof: 1,
                    value: 5e5,
                },
                ConcentratedLoad {
                    node_id: 3,
                    dof: 1,
                    value: 5e5,
                },
            ],
        };
        let sigma = 1e6;
        let sol = solve(&model).unwrap();
        for (i, &vm) in sol.von_mises.iter().enumerate() {
            assert!((vm - sigma).abs() / sigma < 1e-6, "element {i}: vM={vm}");
        }
        // Constant-stress triangle: nodal value equals the exact stress.
        for (i, &nvm) in sol.nodal_von_mises.iter().enumerate() {
            assert!((nvm - sigma).abs() / sigma < 1e-6, "node {i}: vM={nvm}");
        }
        let ux_right = sigma / E;
        assert!((sol.displacements[1][0] - ux_right).abs() / ux_right < 1e-6);
        let uy_top = -NU * sigma / E;
        assert!((sol.displacements[2][1] - uy_top).abs() / uy_top.abs() < 1e-6);
        assert!(sol.iterations > 0);
        assert!(sol.final_residual <= 1e-10);
    }

    #[test]
    fn cantilever_tip_deflection() {
        // 20x4, not the suggested 10x4: fully-integrated bilinear quads
        // lock ~28% stiff at 10x4 for L/h=10; 20x4 lands inside 15%.
        let nx = 20;
        let ny = 4;
        let (l, h) = (1.0, 0.1);
        let mut model = quad_model(nx, ny, l, h);
        // Clamped root: ux and uy on the whole left edge.
        let root: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        model.boundary_conditions.push(BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2],
        });
        let p = -1000.0; // downward tip load
        let tip = (2 * (nx + 1) + nx + 1) as u32; // middle of the right edge
        model.cloads = vec![ConcentratedLoad {
            node_id: tip,
            dof: 2,
            value: p,
        }];
        let sol = solve(&model).unwrap();
        assert!(sol.final_residual <= 1e-10);
        assert!(sol.iterations < 1000.max(10 * sol.displacements.len() * 2));
        let i = model.mesh.node_index(tip).unwrap();
        let delta = sol.displacements[i][1];
        let inertia = h * h * h / 12.0;
        let beam = p * l * l * l / (3.0 * E * inertia);
        let rel = (delta - beam) / beam;
        assert!(
            rel.abs() <= 0.15,
            "tip deflection off by {:.1}% (fem={delta:.3e}, beam={beam:.3e})",
            100.0 * rel
        );
        // Vertical reactions balance the applied load (sum R_y = -P).
        let ry: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 2)
            .map(|r| r.value)
            .sum();
        assert!((ry + p).abs() < 1e-6 * p.abs());
    }

    #[test]
    fn cantilever_root_nodal_von_mises_peak() {
        // 40x4; beam theory peak 6PL/h^2 = 120 MPa at the clamped edge.
        // Extrapolated nodal values capture the boundary peak that the
        // element-centroid average smears (centroid reads ~half of it).
        let nx = 40;
        let ny = 4;
        let (l, h) = (1.0, 0.1);
        let mut model = quad_model(nx, ny, l, h);
        let root: Vec<u32> = (0..=ny).map(|j| (j * (nx + 1) + 1) as u32).collect();
        model.boundary_conditions.push(BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2],
        });
        let p = 200_000.0; // 6PL/h^2 = 120 MPa
        let tip = (2 * (nx + 1) + nx + 1) as u32;
        model.cloads = vec![ConcentratedLoad {
            node_id: tip,
            dof: 2,
            value: -p,
        }];
        let sol = solve(&model).unwrap();
        assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
        let peak = sol.nodal_von_mises.iter().cloned().fold(0.0, f64::max);
        assert!(
            peak > 100e6,
            "max nodal von Mises {peak:.3e} Pa, want > 100 MPa (beam peak 120 MPa)"
        );
    }

    #[test]
    fn kirsch_quarter_plate_nodal_hole_edge() {
        // Radial mesh matching examples/plate_hole.inp: 8x8 CPS4 over the
        // quarter annulus, r = 0.2..1.0, theta = 0..pi/2, t = 0.01.
        // Nominal far-field 1 MPa -> Kirsch hole-edge peak 3 MPa.
        let (nr, nt) = (8, 8);
        let nid = |i: usize, j: usize| (j * (nt + 1) + i + 1) as u32;
        let mut nodes = Vec::new();
        for j in 0..=nr {
            for i in 0..=nt {
                let theta = i as f64 * std::f64::consts::FRAC_PI_2 / nt as f64;
                let r = 0.2 + 0.8 * j as f64 / nr as f64;
                nodes.push(Node {
                    id: nid(i, j),
                    coords: [r * theta.cos(), r * theta.sin(), 0.0],
                });
            }
        }
        let elements: Vec<Element> = (0..nr)
            .flat_map(|j| {
                (0..nt).map(move |i| Element {
                    id: (j * nt + i + 1) as u32,
                    kind: ElemKind::Cps4,
                    connectivity: vec![nid(i, j), nid(i, j + 1), nid(i + 1, j + 1), nid(i + 1, j)],
                })
            })
            .collect();
        let mut model = Model {
            mesh: Mesh { nodes, elements },
            thickness: 0.01,
            materials: vec![steel()],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=(nr * nt) as u32).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            ..Default::default()
        };
        // y-symmetry on theta=0, x-symmetry on theta=pi/2.
        model.boundary_conditions = vec![
            BoundaryCondition {
                node_ids: (0..=nr).map(|j| nid(0, j)).collect(),
                dofs: vec![2],
            },
            BoundaryCondition {
                node_ids: (0..=nr).map(|j| nid(nt, j)).collect(),
                dofs: vec![1],
            },
        ];
        // Far-edge x-traction, nominal 1 MPa: consistent nodal loads
        // sigma*t*ds*cos(theta), half weight at the two end nodes.
        let ds = std::f64::consts::FRAC_PI_2 / nt as f64;
        model.cloads = (0..=nt)
            .filter_map(|i| {
                let theta = i as f64 * std::f64::consts::FRAC_PI_2 / nt as f64;
                let w = if i == 0 || i == nt { 0.5 } else { 1.0 };
                let f = 1e6 * 0.01 * ds * theta.cos() * w;
                (f > 0.0).then(|| ConcentratedLoad {
                    node_id: nid(i, nr),
                    dof: 1,
                    value: f,
                })
            })
            .collect();
        let sol = solve(&model).unwrap();
        assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
        let peak = sol.nodal_von_mises.iter().cloned().fold(0.0, f64::max);
        assert!(
            peak > 2.4e6,
            "max nodal von Mises {peak:.4e} Pa, want > 2.4 MPa (Kirsch 3.0 MPa)"
        );
    }

    #[test]
    fn missing_constraints_is_no_constraints() {
        let mut model = quad_model(1, 1, 1.0, 1.0);
        model.cloads = vec![ConcentratedLoad {
            node_id: 2,
            dof: 1,
            value: 1.0,
        }];
        assert!(matches!(solve(&model), Err(SolveError::NoConstraints)));
    }

    #[test]
    fn c3d8_patch_single_cube() {
        let mut model = solid_model(assembly::tests::hex_mesh(1, 1, 1, 1.0, 1.0, 1.0));
        cube_tension_bcs(&mut model);
        let sol = solve(&model).unwrap();
        assert!(sol.final_residual <= 1e-10);
        assert_uniform_tension(&model, &sol);
    }

    #[test]
    fn c3d4_patch_cube() {
        // Unit cube split into 6 Kuhn tets around the 1-7 diagonal (the
        // minimal conforming tet mesh of a cube; two tets cannot tile it).
        // Constant-strain tets reproduce the uniform-tension solution
        // exactly and nodal averaging stays uniform across the shared
        // interior faces.
        let corners = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let nodes: Vec<Node> = corners
            .iter()
            .enumerate()
            .map(|(i, &c)| Node {
                id: i as u32 + 1,
                coords: c,
            })
            .collect();
        let tets = [
            [1, 2, 3, 7],
            [1, 3, 4, 7],
            [1, 4, 8, 7],
            [1, 8, 5, 7],
            [1, 5, 6, 7],
            [1, 6, 2, 7],
        ];
        let elements = tets
            .iter()
            .enumerate()
            .map(|(i, t)| Element {
                id: i as u32 + 1,
                kind: ElemKind::C3d4,
                connectivity: t.to_vec(),
            })
            .collect();
        let mut model = solid_model(Mesh { nodes, elements });
        // Consistent nodal loads for uniform sigma_z: constant traction on
        // a triangular face splits in equal thirds, and only two of the
        // surface triangles lie on z=1 ((5,6,7) and (5,7,8)) -> corner
        // nodes 5 and 7 carry sigma/3, nodes 6 and 8 sigma/6.
        model.boundary_conditions = vec![
            BoundaryCondition {
                node_ids: vec![1, 2, 3, 4],
                dofs: vec![3],
            },
            BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1],
            },
            BoundaryCondition {
                node_ids: vec![1, 2],
                dofs: vec![2],
            },
        ];
        let sigma = 1e6;
        model.cloads = [(5u32, 2.0), (6, 1.0), (7, 2.0), (8, 1.0)]
            .iter()
            .map(|&(nid, w)| ConcentratedLoad {
                node_id: nid,
                dof: 3,
                value: sigma * w / 6.0,
            })
            .collect();
        let sol = solve(&model).unwrap();
        assert!(sol.final_residual <= 1e-10);
        assert_uniform_tension(&model, &sol);
    }

    #[test]
    fn c3d8_cantilever_tip_deflection() {
        // 40x4x4 hexes, L=2, square section h=b=0.1: Euler-Bernoulli
        // tip deflection PL^3/(3EI) = 15.238 mm. Fully integrated linear
        // hexes shear-lock in bending; 40 elements along the length land
        // inside the same 15% band as the 2D quad test.
        let (nx, ny, nz) = (40, 4, 4);
        let (l, h, b) = (2.0, 0.1, 0.1);
        let mut model = solid_model(assembly::tests::hex_mesh(nx, ny, nz, l, h, b));
        // Clamped root: all dofs on the x = 0 face.
        let root: Vec<u32> = (0..=nz)
            .flat_map(|k| {
                (0..=ny).map(move |j| (k * (ny + 1) * (nx + 1) + j * (nx + 1) + 1) as u32)
            })
            .collect();
        model.boundary_conditions.push(BoundaryCondition {
            node_ids: root,
            dofs: vec![1, 2, 3],
        });
        let p = -10000.0; // downward tip load -> 15.238 mm
        let tip = ((nz / 2) * (ny + 1) * (nx + 1) + (ny / 2) * (nx + 1) + nx + 1) as u32;
        model.cloads = vec![ConcentratedLoad {
            node_id: tip,
            dof: 2,
            value: p,
        }];
        let sol = solve(&model).unwrap();
        assert!(sol.final_residual <= 1e-10);
        let i = model.mesh.node_index(tip).unwrap();
        let delta = sol.displacements[i][1];
        let inertia = b * h * h * h / 12.0;
        let beam = p * l * l * l / (3.0 * E * inertia);
        let rel = (delta - beam) / beam;
        assert!(
            rel.abs() <= 0.15,
            "tip deflection off by {:.1}% (fem={delta:.3e}, beam={beam:.3e})",
            100.0 * rel
        );
        let ry: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 2)
            .map(|r| r.value)
            .sum();
        assert!((ry + p).abs() < 1e-6 * p.abs());
    }

    #[test]
    fn mixed_plane_and_solid_is_invalid() {
        let mut model = solid_model(assembly::tests::hex_mesh(1, 1, 1, 1.0, 1.0, 1.0));
        model.mesh.elements.push(Element {
            id: 2,
            kind: ElemKind::Cps4,
            connectivity: vec![1, 2, 3, 4],
        });
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![1, 2, 3],
        }];
        match solve(&model) {
            Err(SolveError::InvalidModel(m)) => assert!(m.contains("mixed dimensionality")),
            other => panic!("expected InvalidModel(mixed dimensionality), got {other:?}"),
        }
    }

    #[test]
    fn inverted_hex_is_invalid_model() {
        let mut model = solid_model(assembly::tests::hex_mesh(1, 1, 1, 1.0, 1.0, 1.0));
        model.mesh.elements[0].connectivity.swap(1, 2); // clockwise bottom face
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![1, 2, 3],
        }];
        match solve(&model) {
            Err(SolveError::InvalidModel(m)) => assert!(m.contains("Jacobian")),
            other => panic!("expected InvalidModel(Jacobian), got {other:?}"),
        }
    }
}
