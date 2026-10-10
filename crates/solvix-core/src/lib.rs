//! # solvix-core
//!
//! FEM kernel for Solvix: mesh data model, linear-elastic elements,
//! global assembly and the linear solver.
//!
//! This file defines the **public type contract** for the whole workspace.
//! `solvix-io` parses Abaqus decks into `Model`; `solvix-cli` orchestrates.
//!
//! Implementation modules: [`element`] (CPS3/CPS4 stiffness), [`assembly`]
//! (CSR assembly), [`solver`] (Dirichlet + preconditioned CG), [`stress`]
//! (von Mises + reactions).

mod assembly;
mod element;
mod solver;
mod stress;

/// A mesh node. `coords[2]` is kept for the 3D roadmap; plane elements use x, y.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node {
    pub id: u32,
    pub coords: [f64; 3],
}

/// Supported element kinds (v0.1: 2D plane stress).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElemKind {
    /// 3-node linear triangle, plane stress.
    Cps3,
    /// 4-node bilinear quadrilateral, plane stress.
    Cps4,
}

impl ElemKind {
    /// Number of nodes per element.
    pub fn node_count(&self) -> usize {
        match self {
            ElemKind::Cps3 => 3,
            ElemKind::Cps4 => 4,
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
/// `dofs` uses 1-based local dof numbers: 1 = x, 2 = y.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryCondition {
    pub node_ids: Vec<u32>,
    pub dofs: Vec<usize>,
}

/// Concentrated load (`*CLOAD`), already expanded to single nodes.
/// `dof` is 1-based: 1 = x, 2 = y. Value in Newtons.
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
    /// Displacement per node; z component is 0.0 in v0.1.
    pub displacements: Vec<[f64; 3]>,
    /// Von Mises stress at each element centroid (Pa).
    pub von_mises: Vec<f64>,
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

    let (mut k, mut f, elems) = assembly::assemble(model)?;
    let constrained = solver::constrained_dofs(model)?;

    // Keep the unmodified system for reaction recovery.
    let k_orig = assembly::Csr {
        row_ptr: k.row_ptr.clone(),
        col_idx: k.col_idx.clone(),
        vals: k.vals.clone(),
    };
    let f_orig = f.clone();

    solver::apply_dirichlet(&mut k, &mut f, &constrained);
    let ndof = f.len();
    let max_iter = 1000.max(10 * ndof);
    let (u, iterations, final_residual) = solver::pcg(&k, &f, 1e-10, max_iter)?;

    let displacements = model
        .mesh
        .nodes
        .iter()
        .enumerate()
        .map(|(i, _)| [u[2 * i], u[2 * i + 1], 0.0])
        .collect();
    let von_mises = stress::von_mises_at_centroids(&u, &elems);
    let reactions = stress::reactions(&k_orig, &u, &f_orig, &constrained);

    Ok(Solution {
        displacements,
        von_mises,
        reactions,
        iterations,
        final_residual,
    })
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
    fn missing_constraints_is_no_constraints() {
        let mut model = quad_model(1, 1, 1.0, 1.0);
        model.cloads = vec![ConcentratedLoad {
            node_id: 2,
            dof: 1,
            value: 1.0,
        }];
        assert!(matches!(solve(&model), Err(SolveError::NoConstraints)));
    }
}
