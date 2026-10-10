//! # solvix-core
//!
//! FEM kernel for Solvix: mesh data model, linear-elastic elements,
//! global assembly and the linear solver.
//!
//! This file defines the **public type contract** for the whole workspace.
//! `solvix-io` parses Abaqus decks into `Model`; `solvix-cli` orchestrates.
//! Implementation of `solve` is issue #1.

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
/// Implementation owned by issue #1.
pub fn solve(_model: &Model) -> Result<Solution, SolveError> {
    todo!("issue #1: FEM kernel")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contract_compiles_and_mesh_lookup_works() {
        let mesh = Mesh {
            nodes: vec![Node { id: 7, coords: [0.0, 0.0, 0.0] }],
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
        assert!(matches!(
            model.validate(),
            Err(ModelError::UnknownNode(8))
        ));
    }
}
