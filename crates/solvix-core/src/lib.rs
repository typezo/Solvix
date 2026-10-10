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
//! CG), [`stress`] (von Mises + reactions). [`solve`] tries PCG first
//! and falls back to the RCM-ordered sparse LDLT in `ldl` when the
//! iteration budget runs out; [`solve_with_method`] pins the method.
//! [`solve_modal`] (backed by `modal`) extracts the smallest natural
//! frequencies and mass-normalized mode shapes. [`solve_thermal`] and
//! [`solve_thermomech`] (backed by `heat`) run steady-state conduction
//! and the sequential thermal-stress pass. [`solve_dynamic`] (backed by
//! the `dyn.rs` module) runs implicit Newmark transient dynamics with
//! Rayleigh damping. [`solve_nl`] (backed by `nl`) runs geometric
//! nonlinearity for plane models: Updated Lagrangian kinematics with
//! full Newton (NLGEOM steps).

mod assembly;
pub mod buckle;
#[path = "dyn.rs"]
mod dynamics;
mod element;
mod heat;
mod ldl;
mod modal;
mod nl;
mod plastic;
mod solver;
mod stress;

/// A mesh node. Plane elements use x, y; solid elements use all three.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Node {
    pub id: u32,
    pub coords: [f64; 3],
}

/// Supported element kinds: plane stress/strain (2 dofs/node), 3D
/// linear solids and 2D beams (3 dofs/node). Plane and solid kinds
/// cannot be mixed in one model (rejected with
/// `InvalidModel("mixed dimensionality")`); beams cannot be mixed with
/// any continuum kind (v1: B31-only models — their third dof is rotz,
/// which no continuum kind shares).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElemKind {
    /// 3-node linear triangle, plane stress.
    Cps3,
    /// 4-node bilinear quadrilateral, plane stress.
    Cps4,
    /// 8-node serendipity quadrilateral, plane stress, 3x3 Gauss.
    Cps8,
    /// 4-node bilinear quadrilateral, plane strain.
    Cpe4,
    /// 8-node serendipity quadrilateral, plane strain, 3x3 Gauss.
    Cpe8,
    /// 8-node trilinear hexahedron, 2x2x2 Gauss.
    C3d8,
    /// 4-node constant-strain tetrahedron.
    C3d4,
    /// 2-node Euler-Bernoulli beam in the xy-plane (ux, uy, rotz);
    /// section A and I from [`Model::beam_sections`].
    B31,
}

impl ElemKind {
    /// Number of nodes per element.
    pub fn node_count(&self) -> usize {
        match self {
            ElemKind::Cps3 => 3,
            ElemKind::Cps4 | ElemKind::Cpe4 => 4,
            ElemKind::Cps8 | ElemKind::Cpe8 | ElemKind::C3d8 => 8,
            ElemKind::C3d4 => 4,
            ElemKind::B31 => 2,
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

/// Distributed load (`*DLOAD`): pressure on one element edge, in Pa.
/// `face` is 1-based Abaqus edge numbering (P1..P4 on CPS4, P1..P3 on
/// CPS3): edge `n` connects local nodes `n` and `n+1` (cyclic). Positive
/// magnitude acts INTO the element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistributedLoad {
    pub element_id: u32,
    pub face: usize,
    pub magnitude: f64,
}

/// Prescribed temperature (`*TEMPERATURE`) on a list of nodes. Units are
/// caller's choice (only differences enter the stress pass, against
/// T_ref = 0).
#[derive(Debug, Clone, PartialEq)]
pub struct TempBc {
    pub node_ids: Vec<u32>,
    pub temperature: f64,
}

/// Surface heat flux (`*DFLUX` S-face) on one element edge, in W/m^2.
/// `face` uses the same 1-based Abaqus edge numbering as
/// [`DistributedLoad`] (S1..S4 on quads, S1..S3 on CPS3). Positive
/// magnitude is heat flux INTO the element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatFlux {
    pub element_id: u32,
    pub face: usize,
    pub magnitude: f64,
}

/// Transient dynamics definition of a step (`*DYNAMIC`): implicit
/// Newmark average-acceleration integration with a constant time step.
#[derive(Debug, Clone, PartialEq)]
pub struct DynamicsSpec {
    /// Constant time increment in seconds.
    pub dt: f64,
    /// Total step time in seconds (`duration/dt` steps are taken).
    pub duration: f64,
    /// Rayleigh mass damping alpha: `C = alpha M + beta K` (`*DAMPING`).
    pub rayleigh_alpha: f64,
    /// Rayleigh stiffness damping beta: `C = alpha M + beta K`.
    pub rayleigh_beta: f64,
}

/// Linear-hardening J2 plasticity parameters reduced from an Abaqus
/// `*PLASTIC` table. `yield_stress` is the first-row von Mises yield
/// stress and `hardening_modulus` is the slope between the first and last
/// `(sigma_y, equivalent_plastic_strain)` rows.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlasticitySpec {
    /// Initial von Mises yield stress in Pa.
    pub yield_stress: f64,
    /// Isotropic hardening modulus `d sigma_y / d eps_pl` in Pa.
    pub hardening_modulus: f64,
}

/// One analysis step: its boundary conditions and loads. Steps built by
/// the `solvix-io` parser carry the ACCUMULATED state at the end of that
/// step (linear static = superposition): OP=NEW/MOD semantics are
/// resolved when the steps are produced, so each step's loads are its
/// own plus everything still active from earlier steps.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Step {
    pub name: String,
    pub boundary_conditions: Vec<BoundaryCondition>,
    pub cloads: Vec<ConcentratedLoad>,
    pub dloads: Vec<DistributedLoad>,
    /// Gravity: (magnitude in m/s^2, direction vector). Body force per
    /// unit mass is `magnitude * dir`; `dir` need not be normalized.
    /// Plane models use only the x/y components.
    pub gravity: Option<(f64, [f64; 3])>,
    /// Prescribed temperatures (`*TEMPERATURE`), already expanded to
    /// node lists. Thermal BCs live on steps only (no legacy top-level
    /// counterpart).
    pub temp_bcs: Vec<TempBc>,
    /// Surface heat fluxes (`*DFLUX` S1..S4).
    pub heat_fluxes: Vec<HeatFlux>,
    /// Volumetric heat generation in W/m^3 (`*DFLUX` BN), one value for
    /// the whole model.
    pub heat_source: Option<f64>,
    /// Transient dynamics (`*DYNAMIC`): when set, the step's loads are
    /// applied as a step function at t = 0 and held constant, and
    /// [`solve_dynamic`] integrates the step. Steps only (no legacy
    /// top-level counterpart).
    pub dynamics: Option<DynamicsSpec>,
    /// Geometric nonlinearity (`*STEP, NLGEOM=YES`): when the LAST step
    /// carries this flag, [`solve`] routes through [`solve_nl`]
    /// (Updated Lagrangian + full Newton) instead of the linear static
    /// path; [`solve_with_method`] always stays linear. Default false.
    pub nlgeom: bool,
}

/// A complete, solvable model.
///
/// Step duality: loads/BCs live either in [`Model::steps`] (multi-step
/// decks, v0.3) or in the legacy top-level [`Model::boundary_conditions`]
/// and [`Model::cloads`] fields (the implicit "step 0", single-step decks
/// and programmatic models). [`solve`] normalizes through
/// [`Model::load_state`]: `steps` wins when non-empty (the last step's
/// accumulated state), otherwise the legacy fields.
#[derive(Debug, Clone)]
pub struct Model {
    pub mesh: Mesh,
    /// Out-of-plane thickness for plane-stress elements (default 1.0).
    pub thickness: f64,
    /// Mass density in kg/m^3 (feeds `*DENSITY`, gravity loads and
    /// [`solve_modal`]). Default 7850.0 (steel); one value covers the
    /// whole model (per-material density is v0.4).
    pub density: f64,
    /// Isotropic thermal conductivity in W/(m·K) (`*CONDUCTIVITY`).
    /// Default 45.0 (steel); one value covers the whole model
    /// (per-material thermal constants are a later refinement).
    pub conductivity: f64,
    /// Coefficient of thermal expansion in 1/K (`*EXPANSION`). Default
    /// 1.2e-5 (steel); one value covers the whole model.
    pub thermal_expansion: f64,
    /// Model-wide J2 plasticity (`*PLASTIC`). v0.5 supports one active
    /// plastic law; per-material plasticity is a later refinement.
    pub plasticity: Option<PlasticitySpec>,
    pub materials: Vec<Material>,
    pub elsets: Vec<ElementSet>,
    /// (elset name, material name) pairs from `*SOLID SECTION` and
    /// `*BEAM SECTION`.
    pub mat_assignments: Vec<(String, String)>,
    /// (elset name, area A, second moment I) pairs from `*BEAM SECTION`
    /// (B31 elements).
    pub beam_sections: Vec<(String, f64, f64)>,
    pub boundary_conditions: Vec<BoundaryCondition>,
    pub cloads: Vec<ConcentratedLoad>,
    pub steps: Vec<Step>,
}

impl Default for Model {
    /// Non-zero defaults: `thickness` 1.0, `density` 7850.0,
    /// `conductivity` 45.0 and `thermal_expansion` 1.2e-5 (steel); the
    /// remaining fields default as their own types do.
    fn default() -> Self {
        Self {
            mesh: Mesh::default(),
            thickness: 1.0,
            density: 7850.0,
            conductivity: 45.0,
            thermal_expansion: 1.2e-5,
            plasticity: None,
            materials: Vec::new(),
            elsets: Vec::new(),
            mat_assignments: Vec::new(),
            beam_sections: Vec::new(),
            boundary_conditions: Vec::new(),
            cloads: Vec::new(),
            steps: Vec::new(),
        }
    }
}

/// Borrowed view of a model's effective BC/load state after step
/// normalization (see [`Model::load_state`]).
pub(crate) struct LoadState<'a> {
    pub(crate) bcs: &'a [BoundaryCondition],
    pub(crate) cloads: &'a [ConcentratedLoad],
    pub(crate) dloads: &'a [DistributedLoad],
    pub(crate) gravity: Option<(f64, [f64; 3])>,
}

/// Equivalent thermal loads for [`solve_thermomech`]: an extra RHS
/// vector (already in the model's dof layout, added to the assembled
/// mechanical loads before reactions are recovered) and per-element
/// thermal stress vectors subtracted during stress recovery.
pub(crate) struct ThermalLoads {
    pub(crate) f: Vec<f64>,
    pub(crate) sigma_th: Vec<[f64; 6]>,
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

    /// Beam section (area, second moment) assigned to the element at
    /// `elem_index`, if any. Same elset resolution as
    /// [`Model::material_for_element`].
    pub fn beam_section_for_element(&self, elem_index: usize) -> Option<(f64, f64)> {
        let elem = self.mesh.elements.get(elem_index)?;
        for (elset, a, i) in &self.beam_sections {
            if let Some(set) = self.elsets.iter().find(|s| &s.name == elset) {
                if set.element_ids.contains(&elem.id) {
                    return Some((*a, *i));
                }
            }
        }
        None
    }

    /// The model's effective BC/load state after step normalization: the
    /// LAST step's accumulated state when `steps` is non-empty (linear
    /// superposition; OP=NEW/MOD were resolved when the steps were
    /// built), otherwise the legacy top-level fields.
    pub(crate) fn load_state(&self) -> LoadState<'_> {
        match self.steps.last() {
            Some(s) => LoadState {
                bcs: &s.boundary_conditions,
                cloads: &s.cloads,
                dloads: &s.dloads,
                gravity: s.gravity,
            },
            None => LoadState {
                bcs: &self.boundary_conditions,
                cloads: &self.cloads,
                dloads: &[],
                gravity: None,
            },
        }
    }

    /// Structural checks shared by parser tests and the solver.
    pub fn validate(&self) -> Result<(), ModelError> {
        if self.mesh.nodes.is_empty() || self.mesh.elements.is_empty() {
            return Err(ModelError::EmptyMesh);
        }
        if self.materials.is_empty() {
            return Err(ModelError::NoMaterial);
        }
        if self.load_state().bcs.is_empty() {
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
        if let Some(p) = &self.plasticity {
            if !p.yield_stress.is_finite()
                || p.yield_stress <= 0.0
                || !p.hardening_modulus.is_finite()
                || p.hardening_modulus < 0.0
            {
                return Err(ModelError::InvalidPlasticity);
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
    InvalidPlasticity,
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
    /// Which linear solver produced the displacements: "pcg" or "ldlt".
    /// `Auto` (the [`solve`] default) reports the method that actually
    /// delivered the answer — "ldlt" means PCG missed its budget and
    /// the direct fallback took over.
    pub solver_method: String,
    /// Diagonal-based conditioning indicator, NOT a true condition
    /// number: on the LDLT path the ratio of largest to smallest
    /// |D_ii| pivot; on the PCG path the ratio of extreme Jacobi
    /// diagonal entries of K. Treat it as a coarse ill-conditioning
    /// warning scale, not as ||K|| * ||K^-1||.
    pub cond_estimate: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SolveError {
    EmptyMesh,
    NoConstraints,
    InvalidModel(String),
    NotConverged { iterations: usize, residual: f64 },
}

/// Linear solver selection for [`solve_with_method`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// PCG first; on `NotConverged` (or a PCG breakdown) fall back to
    /// the sparse LDLT. The default behavior of [`solve`].
    Auto,
    /// Jacobi-PCG only — surface the iteration failure, no fallback.
    Pcg,
    /// RCM-ordered sparse LDLT only.
    Ldlt,
}

pub use dynamics::{solve_dynamic, DynOptions, DynSolution, TimePoint};
pub use heat::{solve_thermal, solve_thermomech, ThermalSolution};
pub use modal::{solve_modal, ModalOptions, ModalSolution};
pub use nl::{solve_nl, NlFrame, NlOptions};

/// Solve the linear-static problem `K u = f` with homogeneous Dirichlet
/// constraints. [`Method::Auto`]: conjugate gradient on CSR with Jacobi
/// preconditioning, falling back to a direct sparse LDLT (RCM ordering)
/// when PCG exceeds its iteration budget. Plane-only models solve with
/// 2 dofs/node, models containing solid elements with 3 (see
/// [`ndof_per_node`]).
///
/// When the LAST step sets `nlgeom` (`*STEP, NLGEOM=YES`), the solve
/// routes through [`solve_nl`] instead and returns the final frame's
/// solution (total displacements under the full load). Models without
/// steps, or whose last step is linear, take the linear path below.
/// A model with [`Model::plasticity`] takes the small-strain J2 branch of
/// that NL path first; NLGEOM is intentionally ignored for that material
/// solve.
///
/// Relative tolerance 1e-10, `max_iter = max(1000, 10 * ndof)`.
pub fn solve(model: &Model) -> Result<Solution, SolveError> {
    if model.plasticity.is_some() {
        let frames = nl::solve_nl(model, &NlOptions::default())?;
        return Ok(frames
            .last()
            .map(|f| f.solution.clone())
            .unwrap_or_default());
    }
    if model.steps.last().is_some_and(|s| s.nlgeom) {
        let frames = nl::solve_nl(model, &NlOptions::default())?;
        return Ok(frames
            .last()
            .map(|f| f.solution.clone())
            .unwrap_or_default());
    }
    solve_with_method(model, Method::Auto)
}

/// [`solve`] with the linear solver pinned; see [`Method`]. Multi-step
/// models solve the accumulated state at the end of the LAST step
/// (linear static superposition).
pub fn solve_with_method(model: &Model, method: Method) -> Result<Solution, SolveError> {
    solve_with_thermal(model, method, None)
}

/// The shared static-solve pipeline: [`solve_with_method`] passes no
/// thermal loads; `heat::solve_thermomech` adds the equivalent thermal
/// load vector and per-element thermal stress for recovery.
pub(crate) fn solve_with_thermal(
    model: &Model,
    method: Method,
    thermal: Option<ThermalLoads>,
) -> Result<Solution, SolveError> {
    // Free-body check first: validate() reports NoBoundaryConditions as a
    // ModelError, but the contract wants SolveError::NoConstraints.
    let loads = model.load_state();
    if loads.bcs.is_empty() {
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
    let (mut k, mut f, elems) = assembly::assemble(model, ndof, &loads)?;
    if let Some(th) = &thermal {
        for (a, b) in f.iter_mut().zip(&th.f) {
            *a += b;
        }
    }
    let constrained = solver::constrained_dofs(&model.mesh, loads.bcs, ndof)?;

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

    let jacobi_cond = || {
        (0..ndof_total)
            .filter_map(|i| {
                (k.row_ptr[i]..k.row_ptr[i + 1])
                    .find(|&idx| k.col_idx[idx] == i)
                    .map(|idx| k.vals[idx].abs())
            })
            .fold((0.0f64, f64::INFINITY), |(hi, lo), m| {
                (hi.max(m), lo.min(m))
            })
    };
    let direct = || -> Result<(Vec<f64>, f64, f64), SolveError> {
        let fact = ldl::Ldlt::factorize(&k)?;
        let u = fact.solve(&f);
        let b_norm = f.iter().map(|v| v * v).sum::<f64>().sqrt();
        let res = if b_norm == 0.0 {
            0.0
        } else {
            let r = k.mat_vec(&u);
            f.iter()
                .zip(&r)
                .map(|(a, c)| a - c)
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt()
                / b_norm
        };
        Ok((u, res, fact.cond_estimate()))
    };

    let (u, iterations, final_residual, solver_method, cond_estimate) = match method {
        Method::Pcg => {
            let (u, iterations, final_residual) = solver::pcg(&k, &f, 1e-10, max_iter)?;
            let (hi, lo) = jacobi_cond();
            (u, iterations, final_residual, "pcg".to_string(), hi / lo)
        }
        Method::Ldlt => {
            let (u, final_residual, cond) = direct()?;
            (u, 0, final_residual, "ldlt".to_string(), cond)
        }
        Method::Auto => match solver::pcg(&k, &f, 1e-10, max_iter) {
            Ok((u, iterations, final_residual)) => {
                let (hi, lo) = jacobi_cond();
                (u, iterations, final_residual, "pcg".to_string(), hi / lo)
            }
            Err(_) => {
                // PCG missed its budget (or broke down); try the direct
                // path and let a real factorization failure (singular
                // matrix) surface with the sharper diagnostic.
                let (u, final_residual, cond) = direct()?;
                (u, 0, final_residual, "ldlt".to_string(), cond)
            }
        },
    };

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
    let thermal_sigma = thermal.as_ref().map(|t| t.sigma_th.as_slice());
    let von_mises = stress::von_mises_at_centroids(&u, ndof, &elems, thermal_sigma);
    let nodal_von_mises =
        stress::nodal_von_mises(&u, ndof, &elems, model.mesh.nodes.len(), thermal_sigma);
    let reactions = stress::reactions(&k_orig, &u, &f_orig, &constrained);

    Ok(Solution {
        displacements,
        von_mises,
        nodal_von_mises,
        reactions,
        iterations,
        final_residual,
        solver_method,
        cond_estimate,
    })
}

/// Dofs per node for the model: 2 while every element is a plane kind,
/// 3 once any C3D8/C3D4 solid or B31 beam is present. Mixing plane and
/// solid kinds in one model is rejected: there is no formulation that
/// shares a stiffness matrix between them. Beams must likewise live in
/// B31-only models: their third dof is rotz, not a translation, so
/// sharing the 3-dof layout with solids or the 2-dof one with plane
/// elements is meaningless.
pub(crate) fn ndof_per_node(model: &Model) -> Result<usize, SolveError> {
    let mut plane = false;
    let mut solid = false;
    let mut beam = false;
    for e in &model.mesh.elements {
        match e.kind {
            ElemKind::Cps3 | ElemKind::Cps4 | ElemKind::Cps8 | ElemKind::Cpe4 | ElemKind::Cpe8 => {
                plane = true
            }
            ElemKind::C3d8 | ElemKind::C3d4 => solid = true,
            ElemKind::B31 => beam = true,
        }
    }
    if beam && (plane || solid) {
        return Err(SolveError::InvalidModel(
            "mixed beam and continuum elements (v1: B31-only models)".into(),
        ));
    }
    if beam {
        return Ok(3);
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

    pub(crate) fn solid_model(mesh: Mesh) -> Model {
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

    /// Every solve-path assertion in this suite runs through both the
    /// default (PCG with LDLT fallback) and the direct LDLT path.
    const BOTH_METHODS: [Method; 2] = [Method::Auto, Method::Ldlt];

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
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
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
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
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
    }

    #[test]
    fn cps3_patch_two_triangles() {
        // Unit square split diagonally; uniform tension is exact for CST.
        let model = Model {
            steps: Vec::new(),
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
            density: 7850.0,
            materials: vec![steel()],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: vec![1, 2],
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            beam_sections: Vec::new(),
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
            ..Model::default()
        };
        let sigma = 1e6;
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
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
            assert!(sol.iterations > 0 || method == Method::Ldlt);
            assert!(sol.final_residual <= 1e-10);
        }
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
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            assert!(sol.final_residual <= 1e-10);
            assert!(sol.iterations < 1000.max(10 * sol.displacements.len() * 2));
            let i = model.mesh.node_index(tip).unwrap();
            let delta = sol.displacements[i][1];
            let inertia = h * h * h / 12.0;
            let beam = p * l * l * l / (3.0 * E * inertia);
            let rel = (delta - beam) / beam;
            assert!(
                rel.abs() <= 0.15,
                "{method:?}: tip deflection off by {:.1}% (fem={delta:.3e}, beam={beam:.3e})",
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
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
            let peak = sol.nodal_von_mises.iter().cloned().fold(0.0, f64::max);
            assert!(
                peak > 100e6,
                "{method:?}: max nodal von Mises {peak:.3e} Pa, want > 100 MPa (beam peak 120 MPa)"
            );
        }
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
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
            let peak = sol.nodal_von_mises.iter().cloned().fold(0.0, f64::max);
            assert!(
                peak > 2.4e6,
                "{method:?}: max nodal von Mises {peak:.4e} Pa, want > 2.4 MPa (Kirsch 3.0 MPa)"
            );
        }
    }

    fn quad8_model(nx: usize, ny: usize, width: f64, height: f64) -> Model {
        let n_elem = (nx * ny) as u32;
        Model {
            mesh: assembly::tests::quad8_mesh(nx, ny, width, height),
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

    /// 0.1 x 0.1 square section: A = 0.01, I = b*h^3/12 = 1e-4/12.
    const BEAM_A: f64 = 0.01;
    const BEAM_I: f64 = 1.0e-4 / 12.0;

    fn beam_model(n: usize, l: f64) -> Model {
        Model {
            mesh: assembly::tests::beam_mesh(n, l),
            thickness: 1.0, // unused by beams
            density: 7850.0,
            materials: vec![steel()],
            elsets: vec![ElementSet {
                name: "all".into(),
                element_ids: (1..=n as u32).collect(),
            }],
            mat_assignments: vec![("all".into(), "steel".into())],
            beam_sections: vec![("all".into(), BEAM_A, BEAM_I)],
            ..Default::default()
        }
    }

    /// Node id at (x, y) — patch tests address nodes by coordinate, not
    /// by the quad8 mesh's renumbered ids.
    fn node_at(model: &Model, x: f64, y: f64) -> u32 {
        model
            .mesh
            .nodes
            .iter()
            .find(|n| (n.coords[0] - x).abs() < 1e-12 && (n.coords[1] - y).abs() < 1e-12)
            .map(|n| n.id)
            .unwrap_or_else(|| panic!("no node at ({x}, {y})"))
    }

    #[test]
    fn cps8_patch_2x2() {
        // 2x2 quad8 mesh: every emitted node is referenced (element
        // centers are skipped by the mesh helper). Uniform x-traction
        // via consistent quadratic-edge loads on the right edge.
        let mut model = quad8_model(2, 2, 1.0, 1.0);
        let left: Vec<u32> = (0..=4)
            .map(|j| node_at(&model, 0.0, j as f64 * 0.25))
            .collect();
        model.boundary_conditions = vec![
            BoundaryCondition {
                node_ids: left,
                dofs: vec![1],
            },
            BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![2],
            },
        ];
        // Two quadratic edges of length 1/2: edge corners sigma*t*L/6,
        // midsides 2*sigma*t*L/3; the shared middle corner takes two
        // corner shares.
        let sigma = 1e6;
        let (c, m) = (sigma / 12.0, sigma / 3.0);
        let loads = [
            (1.0f64, 0.0, c),
            (1.0, 0.25, m),
            (1.0, 0.5, 2.0 * c),
            (1.0, 0.75, m),
            (1.0, 1.0, c),
        ];
        model.cloads = loads
            .iter()
            .map(|&(x, y, v)| ConcentratedLoad {
                node_id: node_at(&model, x, y),
                dof: 1,
                value: v,
            })
            .collect();
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            // Quadratic serendipity reproduces constant strain exactly
            // in every element.
            for (i, &vm) in sol.von_mises.iter().enumerate() {
                assert!((vm - sigma).abs() / sigma < 1e-8, "element {i}: vM={vm}");
            }
            // Constant stress field: nodal extrapolation exact at every
            // referenced node.
            assert_eq!(sol.nodal_von_mises.len(), model.mesh.nodes.len());
            for (i, &nvm) in sol.nodal_von_mises.iter().enumerate() {
                assert!((nvm - sigma).abs() / sigma < 1e-8, "node {i}: vM={nvm}");
            }
            let ux_right = sigma / E;
            for j in 0..=4 {
                let nid = node_at(&model, 1.0, j as f64 * 0.25);
                let i = model.mesh.node_index(nid).unwrap();
                assert!(
                    (sol.displacements[i][0] - ux_right).abs() / ux_right < 1e-8,
                    "node {nid}: ux={}",
                    sol.displacements[i][0]
                );
            }
            let rx: f64 = sol
                .reactions
                .iter()
                .filter(|r| r.dof == 1)
                .map(|r| r.value)
                .sum();
            assert!((rx + sigma).abs() < 1e-6 * sigma);
        }
    }

    #[test]
    fn cps8_cantilever_beats_cps4() {
        // Same L=2, h=0.1, t=0.1 cantilever, P = 1e4 at the tip:
        // Euler-Bernoulli 15.238 mm. 20x2 CPS8 (40 nodes-worth of dofs)
        // must beat 40x4 CPS4 (41 nodes) — quadratics do not shear-lock.
        let (l, h) = (2.0, 0.1);
        let p = -1.0e4;
        let inertia = 0.1 * h * h * h / 12.0;
        let beam = p * l * l * l / (3.0 * E * inertia);
        let tip_deflection = |model: &Model, tip: u32| {
            let sol = solve(model).unwrap();
            let i = model.mesh.node_index(tip).unwrap();
            (sol.displacements[i][1] - beam) / beam
        };
        // CPS8 20x2.
        let mut m8 = quad8_model(20, 2, l, h);
        m8.thickness = 0.1;
        let root8: Vec<u32> = (0..=4)
            .map(|j| node_at(&m8, 0.0, j as f64 * h / 4.0))
            .collect();
        m8.boundary_conditions = vec![BoundaryCondition {
            node_ids: root8,
            dofs: vec![1, 2],
        }];
        let tip8 = node_at(&m8, l, h / 2.0);
        m8.cloads = vec![ConcentratedLoad {
            node_id: tip8,
            dof: 2,
            value: p,
        }];
        let rel8 = tip_deflection(&m8, tip8);
        // CPS4 40x4 (the c3d8-test geometry, plane stress).
        let mut m4 = quad_model(40, 4, l, h);
        m4.thickness = 0.1;
        let root4: Vec<u32> = (0..=4)
            .map(|j| node_at(&m4, 0.0, j as f64 * h / 4.0))
            .collect();
        m4.boundary_conditions = vec![BoundaryCondition {
            node_ids: root4,
            dofs: vec![1, 2],
        }];
        let tip4 = node_at(&m4, l, h / 2.0);
        m4.cloads = vec![ConcentratedLoad {
            node_id: tip4,
            dof: 2,
            value: p,
        }];
        let rel4 = tip_deflection(&m4, tip4);
        eprintln!(
            "cantilever L=2: cps8 20x2 rel={rel8:+.3}%  cps4 40x4 rel={rel4:+.3}% (beam {beam:.6e} m)"
        );
        assert!(
            rel8.abs() < rel4.abs(),
            "CPS8 20x2 ({rel8:.3}%) must beat CPS4 40x4 ({rel4:.3}%)"
        );
        assert!(rel8.abs() < 0.03, "CPS8 error {rel8:.3}% not < 3%");
    }

    #[test]
    fn cpe4_patch_matches_plane_strain_theory() {
        // Traction-controlled uniaxial stress: sigma_xx = sigma in both
        // plane stress and plane strain, but the stiffness differs.
        let sigma = 1e6;
        let mut cpe = quad_model(1, 1, 1.0, 1.0);
        for e in &mut cpe.mesh.elements {
            e.kind = ElemKind::Cpe4;
        }
        fix_left_edge_x(&mut cpe, 1, 1);
        cpe.cloads = vec![
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
        // ... and the plane-stress twin for contrast.
        let mut cps = quad_model(1, 1, 1.0, 1.0);
        fix_left_edge_x(&mut cps, 1, 1);
        cps.cloads = cpe.cloads.clone();
        for method in BOTH_METHODS {
            let spe = solve_with_method(&cpe, method).unwrap();
            // sigma_zz = nu*sigma -> vM = sigma*sqrt(1 - nu + nu^2).
            let vm_exact = sigma * (1.0 - NU + NU * NU).sqrt();
            assert!(
                (spe.von_mises[0] - vm_exact).abs() / vm_exact < 1e-6,
                "vM={} vs {}",
                spe.von_mises[0],
                vm_exact
            );
            // Stiffer in-plane: uniaxial STRESS compliance from the
            // plane-strain D inverse: eps_x = sigma(1-nu^2)/E,
            // eps_y = -nu(1+nu) sigma/E.
            let ux = sigma * (1.0 - NU * NU) / E;
            for &nid in &[2u32, 4] {
                let i = cpe.mesh.node_index(nid).unwrap();
                assert!(
                    (spe.displacements[i][0] - ux).abs() / ux < 1e-6,
                    "node {nid}: ux={}",
                    spe.displacements[i][0]
                );
            }
            // eps_y = -nu(1+nu)*sigma/E (vs -nu*sigma/E for plane stress).
            let uy = -NU * (1.0 + NU) * sigma / E;
            let top = cpe.mesh.node_index(4).unwrap();
            assert!((spe.displacements[top][1] - uy).abs() / uy.abs() < 1e-6);
            let sps = solve_with_method(&cps, method).unwrap();
            assert!((sps.von_mises[0] - sigma).abs() / sigma < 1e-6);
            assert!((sps.von_mises[0] - spe.von_mises[0]).abs() / sigma > 1e-3);
        }
    }

    #[test]
    fn b31_cantilever_tip_deflection_exact() {
        // One Hermite element: exact for a tip point load.
        let l = 2.0;
        let p = -1.0e4;
        let mut model = beam_model(1, l);
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![1, 2, 3],
        }];
        model.cloads = vec![ConcentratedLoad {
            node_id: 2,
            dof: 2,
            value: p,
        }];
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            let beam = p * l * l * l / (3.0 * E * BEAM_I);
            let delta = sol.displacements[1][1];
            assert!(
                ((delta - beam) / beam).abs() < 1e-6,
                "{method:?}: tip={delta:.6e}, beam={beam:.6e}"
            );
            // Rotations match v' = PL^2/2EI at the tip; the nodal fiber
            // stress at the clamp is M*c/I = PL*(h/2)/I = 120 MPa (the
            // element centroid value is the mid-span one, PL/2 -> half).
            let rot = p * l * l / (2.0 * E * BEAM_I);
            assert!((sol.displacements[1][2] - rot).abs() / rot.abs() < 1e-6);
            let sigma = p.abs() * l * 0.05 / BEAM_I;
            assert!(
                (sol.nodal_von_mises[0] - sigma).abs() / sigma < 1e-6,
                "clamp vM={} vs {sigma}",
                sol.nodal_von_mises[0]
            );
            assert!((sol.von_mises[0] - sigma / 2.0).abs() / sigma < 1e-6);
            // Reactions balance: R_y = -P, M = PL.
            let ry: f64 = sol
                .reactions
                .iter()
                .filter(|r| r.dof == 2)
                .map(|r| r.value)
                .sum();
            assert!((ry + p).abs() < 1e-6 * p.abs());
            let rm: f64 = sol
                .reactions
                .iter()
                .filter(|r| r.dof == 3)
                .map(|r| r.value)
                .sum();
            assert!((rm - p.abs() * l).abs() < 1e-6 * p.abs() * l);
        }
    }

    #[test]
    fn mixed_beam_and_plane_is_invalid() {
        let mut model = beam_model(1, 1.0);
        model.mesh.nodes.push(Node {
            id: 3,
            coords: [1.0, 1.0, 0.0],
        });
        model.mesh.nodes.push(Node {
            id: 4,
            coords: [0.0, 1.0, 0.0],
        });
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
            Err(SolveError::InvalidModel(m)) => assert!(m.contains("B31-only")),
            other => panic!("expected InvalidModel(B31-only), got {other:?}"),
        }
    }

    #[test]
    fn inverted_quad8_is_invalid_model() {
        let mut model = quad8_model(1, 1, 1.0, 1.0);
        model.mesh.elements[0].connectivity.swap(1, 3); // clockwise corners
        model.boundary_conditions = vec![BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![1, 2],
        }];
        match solve(&model) {
            Err(SolveError::InvalidModel(m)) => assert!(m.contains("Jacobian")),
            other => panic!("expected InvalidModel(Jacobian), got {other:?}"),
        }
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
    fn edge_pressure_matches_equivalent_cload() {
        // Pressure P on edge 2 (right, nodes 2-3) of one CPS4 == the total
        // P*len split 1/2-1/2 as CLOADs: consistent linear-edge nodal forces.
        let mut m = quad_model(1, 1, 1.0, 1.0);
        fix_left_edge_x(&mut m, 1, 1);
        let p = 1.0e6;
        let mut pressurized = m.clone();
        pressurized.steps = vec![Step {
            name: "pressure".into(),
            boundary_conditions: m.boundary_conditions.clone(),
            dloads: vec![DistributedLoad {
                element_id: 1,
                face: 2,
                magnitude: p,
            }],
            ..Step::default()
        }];
        m.cloads = vec![
            ConcentratedLoad {
                node_id: 2,
                dof: 1,
                value: -p / 2.0,
            },
            ConcentratedLoad {
                node_id: 4,
                dof: 1,
                value: -p / 2.0,
            },
        ];
        let a = solve(&pressurized).unwrap();
        let b = solve(&m).unwrap();
        for (i, (ua, ub)) in a.displacements.iter().zip(&b.displacements).enumerate() {
            for d in 0..3 {
                assert!(
                    (ua[d] - ub[d]).abs() < 1e-9,
                    "node {i} dof {d}: dload={ua:?} cload={ub:?}"
                );
            }
        }
        // Uniform uniaxial compression: the patch stress is exactly P.
        for method in BOTH_METHODS {
            let sol = solve_with_method(&pressurized, method).unwrap();
            assert!((sol.von_mises[0] - p).abs() / p < 1e-6);
            let rx: f64 = sol
                .reactions
                .iter()
                .filter(|r| r.dof == 1)
                .map(|r| r.value)
                .sum();
            assert!((rx - p).abs() < 1e-6 * p);
        }
    }

    #[test]
    fn gravity_column_base_reactions_balance_weight() {
        // Column w=0.1 x L=1 x t=0.01 under g=9.81 downward: the summed
        // vertical base reactions equal rho * A * L * g exactly.
        let (w, l, t) = (0.1, 1.0, 0.01);
        let mut m = quad_model(1, 8, w, l);
        m.thickness = t;
        let bcs = vec![
            BoundaryCondition {
                node_ids: vec![1, 2],
                dofs: vec![2],
            },
            BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1],
            },
        ];
        m.steps = vec![Step {
            name: "gravity".into(),
            boundary_conditions: bcs,
            gravity: Some((9.81, [0.0, -1.0, 0.0])),
            ..Step::default()
        }];
        let sol = solve(&m).unwrap();
        let ry: f64 = sol
            .reactions
            .iter()
            .filter(|r| r.dof == 2)
            .map(|r| r.value)
            .sum();
        let expect = 7850.0 * w * t * l * 9.81;
        assert!(
            (ry - expect).abs() < 1e-9 * expect,
            "sum R_y={ry:.6}, weight={expect:.6}"
        );
    }

    #[test]
    fn c3d8_patch_single_cube() {
        let mut model = solid_model(assembly::tests::hex_mesh(1, 1, 1, 1.0, 1.0, 1.0));
        cube_tension_bcs(&mut model);
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            assert!(sol.final_residual <= 1e-10);
            assert_uniform_tension(&model, &sol);
        }
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
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            assert!(sol.final_residual <= 1e-10);
            assert_uniform_tension(&model, &sol);
        }
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
        for method in BOTH_METHODS {
            let sol = solve_with_method(&model, method).unwrap();
            // eps*cond(K) ~ 1e-10 for this mesh: the direct path rides
            // that floor, give it one extra decade.
            let tol = if method == Method::Ldlt { 1e-8 } else { 1e-10 };
            assert!(sol.final_residual <= tol);
            let i = model.mesh.node_index(tip).unwrap();
            let delta = sol.displacements[i][1];
            let inertia = b * h * h * h / 12.0;
            let beam = p * l * l * l / (3.0 * E * inertia);
            let rel = (delta - beam) / beam;
            assert!(
                rel.abs() <= 0.15,
                "{method:?}: tip deflection off by {:.1}% (fem={delta:.3e}, beam={beam:.3e})",
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
    }

    #[test]
    fn slender_column_pcg_budget_exhaustion_falls_back_to_ldlt() {
        // 1x40 quad column of aspect 1000 (w = 1e-3, L = 1):
        // bending-dominated and ill-conditioned enough that Jacobi-PCG
        // burns its whole iteration budget (10 * ndof) without
        // reaching 1e-10. The Auto path must fall back to LDLT and
        // still nail the analytical axial stretch PL/(EA). (A 2-element
        // column needs element aspect ~1e8 to stall PCG, and there the
        // element stiffness itself is numerically rank-deficient in
        // f64; the tall column hits the same failure mode with
        // well-conditioned elements.)
        let (nx, ny) = (1, 40);
        let (w, l) = (1e-3, 1.0);
        let mut model = quad_model(nx, ny, w, l);
        // Pin uy on the base (y=0) and ux on the left base node only,
        // so the Poisson contraction stays free: the uniform-stress
        // solution remains exact for CPS4.
        model.boundary_conditions = vec![
            BoundaryCondition {
                node_ids: vec![1, 2],
                dofs: vec![2],
            },
            BoundaryCondition {
                node_ids: vec![1],
                dofs: vec![1],
            },
        ];
        let top_left = (ny * (nx + 1) + 1) as u32;
        let top_right = ((ny + 1) * (nx + 1)) as u32;
        let p = 1000.0;
        model.cloads = vec![
            ConcentratedLoad {
                node_id: top_left,
                dof: 2,
                value: p / 2.0,
            },
            ConcentratedLoad {
                node_id: top_right,
                dof: 2,
                value: p / 2.0,
            },
        ];
        // PCG alone must fail: that is what triggers the fallback.
        match solve_with_method(&model, Method::Pcg) {
            Err(SolveError::NotConverged {
                iterations,
                residual,
            }) => {
                eprintln!("pcg stalled: it={iterations} res={residual:.3e}");
            }
            Ok(sol) => panic!(
                "pcg converged (it={}, res={:.3e}); column not ill-conditioned enough",
                sol.iterations, sol.final_residual
            ),
            Err(other) => panic!("unexpected {other:?}"),
        }
        for method in [Method::Auto, Method::Ldlt] {
            let sol = solve_with_method(&model, method).unwrap();
            assert_eq!(
                sol.solver_method, "ldlt",
                "{method:?}: expected the direct path"
            );
            let stretch = p * l / (E * w); // A = w * t, t = 1
            for &nid in &[top_left, top_right] {
                let i = model.mesh.node_index(nid).unwrap();
                let uy = sol.displacements[i][1];
                assert!(
                    (uy - stretch).abs() / stretch < 1e-6,
                    "{method:?}: node {nid} uy={uy:.6e}, exact {stretch:.6e}"
                );
            }
            let ry: f64 = sol
                .reactions
                .iter()
                .filter(|r| r.dof == 2)
                .map(|r| r.value)
                .sum();
            assert!((ry + p).abs() < 1e-6 * p);
            assert!(
                sol.cond_estimate > 1e4,
                "{method:?}: cond estimate {}",
                sol.cond_estimate
            );
        }
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
