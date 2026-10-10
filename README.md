# Solvix

**The agent-era, open-source Abaqus. A finite element solver in Rust.**

Solvix reads Abaqus-style `.inp` decks, solves linear-static problems, and
answers in JSON — built to be driven by humans *and* AI agents.

```console
$ solvix run examples/cantilever.inp --json
{"ok": true, "job": "cantilever", "nodes": 861, "elements": 800, "iterations": 234,
 "residual": 8.7e-13, "max_disp_mm": 4.73, "max_von_mises_mpa": 38.2, "outputs": ["cantilever.vtu"]}
```

## Status

| Version | Scope |
|---------|-------|
| v0.1 (in flight) | 2D plane stress (CPS3/CPS4), `.inp` parser, VTU output, JSON CLI, benchmarked |
| v0.2 (planned)   | 3D hexahedra (C3D8), eigenvalue/modal, MCP server for agents |

See the [issue board](https://github.com/typezo/Solvix/issues) for the roadmap.

## Layout

- `crates/solvix-core` — FEM kernel: mesh, elements, assembly, sparse solver
- `crates/solvix-io` — Abaqus `.inp` parser, VTU writer
- `crates/solvix-cli` — `solvix` binary, JSON-first output

## Build

```console
$ cargo build --release
$ cargo test --workspace
```

Apache-2.0. FIGHT FOR OPEN SOURCE.
