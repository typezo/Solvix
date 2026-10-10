# Solvix

> **The agent-era, open-source Abaqus.**

Solvix is a Rust finite element solver built for both engineers and AI agents.
It reads Abaqus-style `.inp` decks, solves 2D linear-static plane-stress
problems, and emits machine-readable JSON plus ParaView-ready VTU results.

## Why Solvix?

- **Abaqus-style input:** reuse familiar keyword decks and meshing workflows.
- **Linear static, concretely scoped:** v0.1 supports CPS3 and CPS4 plane-stress
  elements with a Jacobi-preconditioned conjugate-gradient solver.
- **Agent-first operation:** every operational command accepts `--json`,
  successful payloads are structured, and exit codes distinguish usage, parse,
  and solve failures.
- **Human-friendly visualization:** nodal displacement and element von Mises
  results are written as VTK UnstructuredGrid `.vtu` files for ParaView.
- **Small and inspectable:** the runtime crates have no third-party Rust
  dependencies, keeping the build and supply chain easy to audit.

## Quickstart

```console
$ cargo build --release
$ ./target/release/solvix run examples/cantilever.inp --json
```

The command solves the model, writes one JSON result object to stdout, and
creates a ParaView-loadable VTU file. Use `-o <directory>` to choose a separate
output directory.

Inspect a deck before spending time on a solve:

```console
$ ./target/release/solvix info examples/cantilever.inp --json
```

## CLI contract for agents

`--json` is the machine interface. JSON payloads go to stdout; progress and
diagnostic text go to stderr. Agents should check both the process exit code
and the payload's `ok` field.

| Exit code | Meaning | Recommended agent response |
|---------:|---------|-----------------------------|
| `0` | Success | Read the JSON payload and downstream artifacts |
| `1` | Usage or I/O error | Correct arguments/paths, then retry |
| `2` | `.inp` parse/validation error | Repair the deck; inspect `line` and `message` |
| `3` | Solve error | Inspect model constraints and solver diagnostics |

See the [Agent Cookbook](docs/AGENT_COOKBOOK.md) for copy-ready recipes.

## Crate architecture

| Crate | Responsibility | Depends on |
|-------|----------------|------------|
| [`solvix-core`](crates/solvix-core) | Mesh and model types, CPS3/CPS4 element behavior, global assembly, Jacobi-PCG linear solver | Rust standard library only |
| [`solvix-io`](crates/solvix-io) | Abaqus `.inp` parser, parse reports/errors, VTU writer | `solvix-core` |
| [`solvix-cli`](crates/solvix-cli) | `solvix` binary, command dispatch, JSON output, exit-code policy | `solvix-core`, `solvix-io` |

The intended dependency direction is strict:

```text
solvix-cli -> solvix-io -> solvix-core
```

The FEM kernel therefore remains independent of file formats and command-line
presentation. More detail is in [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## v0.1 keyword support

Solvix uses a self-consistent unit system; it does not silently convert units.
The bundled examples use SI coordinates in meters, loads in newtons, and
elastic properties in pascals.

| Keyword | v0.1 support | Notes |
|---------|--------------|-------|
| `*NODE` | Supported | Node ID, X, Y; Z may be represented internally but plane elements use X/Y |
| `*ELEMENT` | Supported | `TYPE=CPS3` and `TYPE=CPS4`; an `ELSET=` name may be supplied |
| `CPS3` | Supported | 3-node linear triangle, plane stress |
| `CPS4` | Supported | 4-node bilinear quadrilateral, plane stress |
| `*NSET` | Supported | Named node sets used by boundaries and loads |
| `*ELSET` | Supported | Explicit element-ID lists and named sets |
| `GENERATE` | Supported | Generated IDs for `*NSET` and `*ELSET` |
| `*MATERIAL` | Supported | Named isotropic material |
| `*ELASTIC` | Supported | Young's modulus and Poisson's ratio |
| `*SOLID SECTION` | Supported | Assigns material to an element set and supplies thickness |
| `*BOUNDARY` | Supported | Prescribed DOFs `1` (X) and `2` (Y); set or node references |
| `*CLOAD` | Supported | Concentrated nodal load on DOF `1` or `2` |
| `*STEP` | Supported | One static step per v0.1 model |
| `*STATIC` | Supported | Linear-static analysis declaration |

Unsupported or malformed keywords are reported rather than silently changing
the meaning of a solve.

## Roadmap

v0.2 development is tracked on the
[Solvix issue board](https://github.com/typezo/Solvix/issues):

- 3D linear hexahedral elements (`C3D8`)
- Eigenvalue/modal analysis
- An MCP server so agents can call Solvix through a standard tool interface

## Benchmark (QA #4 will fill)

## License

Solvix is released under the [Apache License 2.0](LICENSE).

**FIGHT FOR OPEN SOURCE.**
