# Solvix Architecture

Solvix is deliberately split into a small FEM kernel, format adapters, and a
presentation-facing CLI. This document describes the v0.1 architecture, the
data contract shared by the crates, and the conventions to preserve as features
are added.

## Goals

- **Correct and inspectable FEM:** implement a deliberately finite v0.1 scope
  (2D linear-static plane stress) rather than pretending broad Abaqus coverage.
- **Agent-operable:** JSON and exit codes are public API, not decorative.
- **Abaqus-compatible where it matters:** use familiar input keywords and node,
  element, and set semantics while making unsupported behavior explicit.
- **Auditable supply chain:** runtime crates use no third-party Rust packages.
- **Separation of concerns:** numerical code does not parse files or print CLI
  tables; I/O code does not decide solver strategy.

## Crate map

| Crate | Owns | Does not own |
|-------|------|--------------|
| `solvix-core` | `Model`, mesh/material/loads/constraints, element kinematics and constitutive behavior, assembly, Jacobi-PCG solver, `Solution`, numerical errors | Input syntax, files, command-line flags, JSON presentation |
| `solvix-io` | `.inp` lexer/parser, keyword/set expansion, parse report, VTU serialization | Element equations, solver selection, process exit codes |
| `solvix-cli` | Argument parsing, command orchestration, stdout/stderr policy, JSON envelope, exit-code mapping | New physics or parser rules |

The compile-time dependency direction is:

```text
solvix-cli
   ├── solvix-io
   │     └── solvix-core
   └── solvix-core
```

`solvix-core` must not import either adapter. `solvix-io` must not import the
CLI. This prevents file-format details or presentation concerns from becoming
prerequisites for testing the kernel.

## End-to-end flow

```text
Abaqus .inp
   |
   v
solvix-io::parse_inp
   |  Model + ParseReport
   v
solvix-core::Model::validate
   |
   v
solvix-core::solve
   |  Solution
   v
solvix-io::write_vtu
   |
   v
JSON stdout + VTU artifact
```

The CLI owns orchestration, but each transformation has a testable library
entry point:

- `parse_inp(text) -> (Model, ParseReport)` or `ParseError`;
- `Model::validate() -> Result<(), ModelError>`;
- `solve(model) -> Solution` or `SolveError`;
- `write_vtu(path, model, solution) -> io::Result<()>`.

## The `Model` contract

`Model` is the complete, solvable v0.1 problem. It is an immutable input to the
solver.

Important invariants:

- Nodes and elements retain their deck IDs; connectivity stores node IDs, not
  vector indices.
- Plane elements use the X/Y coordinates. A third coordinate is retained in the
  data model for the 3D roadmap and is not consumed by CPS3/CPS4 behavior.
- CPS3 has three nodes; CPS4 has four. Element connectivity is counter-clockwise
  in the element's own numbering.
- Materials are isotropic linear-elastic with Young's modulus and Poisson's
  ratio.
- Element sets associate elements with a material through `*SOLID SECTION`.
- Plane-stress thickness defaults to `1.0` when no section thickness is given.
- Boundary DOFs and concentrated-load DOFs are 1-based: `1` is X and `2` is Y.
- Loads are expanded to individual node/DOF/value records before `solve`.
- `validate()` rejects an empty mesh, missing material, absent boundary
  conditions, malformed element connectivity, and references to unknown nodes.

The model is unit-neutral in representation but not in meaning: all input
quantities must use one self-consistent system. Bundled examples use meters,
pascals, and newtons.

## The `Solution` contract

`solve` accepts a reference to a validated `Model` and returns a `Solution`
without mutating the model.

| Field | Contract |
|-------|----------|
| `displacements` | One `[f64; 3]` per node, in `Model::mesh::nodes` order; Z is zero for v0.1 |
| `von_mises` | One scalar per element, in `Model::mesh::elements` order, evaluated at the element centroid |
| `reactions` | Forces at constrained node/DOF pairs |
| `iterations` | Number of PCG iterations performed |
| `final_residual` | Final normalized residual used to assess convergence |

Result arrays must remain parallel to mesh vectors. Writers and callers must
not assume IDs are contiguous, because Abaqus permits sparse node and element
numbering.

The current solver is conjugate gradient on the assembled sparse system with
Jacobi (diagonal) preconditioning. Homogeneous Dirichlet boundary conditions
are imposed before solution. A non-converged iteration is a `SolveError`, not a
silently approximate success.

## Error and exit-code boundaries

Library errors describe semantic classes; only the CLI maps them to process
exit codes.

| Boundary | Library type | CLI exit |
|----------|--------------|---------:|
| Invalid deck syntax/model structure | `ParseError`, `ModelError` | `2` |
| Bad command/path/output directory | CLI usage or `std::io::Error` before parsing | `1` |
| Failed assembly/solve/convergence/output during solve | `SolveError` or solve-phase I/O error | `3` |
| Successful command | `Result::Ok` | `0` |

This split lets library tests assert precise error values while agents receive
coarse, dependable process states.

## Zero-dependency policy

The runtime crates currently declare no third-party dependencies. They use the
Rust standard library and internal workspace crates only.

The policy exists to:

- keep cold builds and cross-compilation predictable;
- minimize transitive supply-chain risk;
- avoid serialization/framework choices leaking into the FEM kernel;
- make numerical code and byte-oriented I/O directly reviewable.

Rules for future changes:

1. Prefer standard-library types and explicit, tested algorithms.
2. Do not add a dependency to work around a local API design issue.
3. JSON and XML/VTU output are generated explicitly; do not introduce a
   general serialization framework without a separate architectural decision.
4. A proposed dependency must state its security posture, MSRV/build impact,
   feature set, alternatives, and migration path.
5. Development-only tooling must remain separate from the runtime dependency
   graph.

## Adding an element kind

Use the following checklist so a new element changes every required boundary in
the same change.

### 1. Define topology and semantics

- Add the `ElemKind` variant and its node count.
- Document node ordering and required orientation in the variant docs.
- Decide whether it is plane stress, plane strain, axisymmetric, or 3D.
- Identify its DOFs and whether existing DOF helpers and boundary/load
  validation support them.

### 2. Implement element mathematics in `solvix-core`

- Define shape functions and their natural-coordinate derivatives.
- Implement the Jacobian and its sign/orientation check.
- Implement the strain-displacement matrix for the element's strain definition.
- Integrate the isotropic constitutive matrix using the element's integration
  rule and consistent thickness treatment.
- Form the local stiffness matrix and expose only the data needed by assembly.
- Add constant-strain and patch tests with known analytical displacement or
  stress fields.
- Exercise non-contiguous IDs and at least one element with non-rectangular
  geometry where relevant.

### 3. Extend assembly and the solver contract

- Add element DOF mapping without hard-coded vector-index assumptions.
- Include the new variant in global assembly and constraint/load projection.
- Confirm reaction recovery for all constrained DOFs.
- Confirm expected PCG convergence on a well-conditioned reference model.
- Update `Model::validate` and error messages for topology and orientation
  failures.

### 4. Add parser support

- Map the Abaqus `TYPE=` spelling to the new variant.
- Reject misspellings and unsupported element families explicitly.
- Add fixture decks covering minimum and representative connectivity.
- Preserve set membership and material assignment behavior.
- Record useful parse warnings without turning unsupported input into a
  different successful model.

### 5. Add result writing

- Map the element to the correct VTK cell type and number of points.
- Preserve element result ordering.
- Verify displacement component count and von Mises semantics.
- Round-trip the VTU through ParaView or a dependency-free XML check in tests.

### 6. Update public behavior and documentation

- Update the README keyword/element support table.
- Update CLI `info` model summaries and JSON tests.
- Add an example only when it clarifies behavior beyond existing fixtures.
- Record unit assumptions and any limitations that an agent could otherwise
  mistake for solver failure.

Element additions should land with their tests and fixtures rather than behind
long-lived partially implemented variants.

## Source organization guidance

As the single-file v0.1 contracts grow into implementations, split modules by
responsibility, not by syntax keyword:

- `solvix-core`: `model`, `mesh`, `element`, `material`, `assembly`, `solver`;
- `solvix-io`: `inp`, parse report, `vtu`;
- `solvix-cli`: argument model, commands, JSON rendering, process policy.

Public types should remain narrow. Prefer moving complete data through `Model`
and `Solution` rather than adding format-specific fields to the kernel.
