# Examples gallery

This gallery covers every deck currently in [`examples/`](../examples).
The commands below are the exact debug-workspace invocations used for
verification. Run them from the repository root:

```console
cargo run -q -p solvix-cli -- run examples/cantilever.inp --json
```

Successful runs print one JSON object on stdout, write progress and
diagnostics to stderr, and create a ParaView-loadable `<job>.vtu` file in the
current directory by default. Use `-o <directory>` to keep generated VTU files
elsewhere. Timings (`solve_ms`) are machine-dependent; the other values shown
here are stable QA results.

| Deck | Physics / purpose | Exact run command | Expected QA signal |
|------|-------------------|-------------------|--------------------|
| [`cantilever.inp`](../examples/cantilever.inp) | 2D plane-stress static bending with fully integrated CPS4 quads | `cargo run -q -p solvix-cli -- run examples/cantilever.inp --json` | 13.832 mm maximum tip displacement versus the 15.238 mm Euler-Bernoulli reference |
| [`plate_hole.inp`](../examples/plate_hole.inp) | Kirsch stress concentration around a circular hole in a quarter-symmetric plate | `cargo run -q -p solvix-cli -- run examples/plate_hole.inp --json` | 2.927 MPa maximum extrapolated nodal von Mises stress versus the 3.0 MPa analytical value |
| [`broken.inp`](../examples/broken.inp) | Intentional parser-failure diagnostic | `cargo run -q -p solvix-cli -- run examples/broken.inp --json` | Exit code 2 with `line 7: malformed coordinate 'x.y'` |
| [`cantilever3d.inp`](../examples/cantilever3d.inp) | 3D static bending with fully integrated C3D8 hexahedra | `cargo run -q -p solvix-cli -- run examples/cantilever3d.inp --json` | 13.719 mm maximum tip displacement versus the 15.238 mm Euler-Bernoulli reference |
| [`thermal_bar.inp`](../examples/thermal_bar.inp) | Steady-state heat conduction followed by sequential thermal stress | `cargo run -q -p solvix-cli -- run examples/thermal_bar.inp --json` | `thermal: true`, maximum temperature 100.000, and a linear 100-to-0 gradient |
| [`dyn_cantilever.inp`](../examples/dyn_cantilever.inp) | Undamped transient response to a suddenly applied load using implicit Newmark integration | `cargo run -q -p solvix-cli -- run examples/dyn_cantilever.inp --json` | 27.468 mm peak displacement, approximately 1.8 times static equilibrium |

## `cantilever.inp` - plane-stress static bending

This is the bundled 2 m by 0.1 m steel cantilever: 40x4 CPS4 elements,
0.01 m thickness, a fully fixed left edge, and a total 1,000 N downward tip
load. It exercises linear plane stress, concentrated nodal loads, fixed
displacement boundaries, Jacobi-preconditioned conjugate gradient, and stress
recovery.

```console
cargo run -q -p solvix-cli -- run examples/cantilever.inp --json
```

Verified stable result: 205 nodes, 160 elements, 390 PCG iterations,
residual `6.933e-11`, `max_disp_mm` **13.832**, and maximum von Mises stress
**79.483 MPa**. The Euler-Bernoulli beam reference is 15.238 mm, so this
mesh is 9.23% stiff. That offset is expected parasitic shear locking in
fully integrated bilinear quadrilaterals; span and through-thickness
refinement moves the result toward the reference. See the
[convergence plot](benchmarks/cantilever-convergence.svg) and
[benchmark discussion](benchmarks/cantilever.md) for the refinement history.

## `plate_hole.inp` - Kirsch stress concentration

This quarter-symmetric unit square has a central circular hole of radius
0.2 m. Symmetry constraints replace the omitted quarters, and consistent arc
tractions totaling approximately 10,000 N impose a nominal 1 MPa uniaxial
far-field load. The 8x8 radial CPS4 mesh demonstrates curved quadrilateral
geometry, symmetry boundary conditions, distributed edge loading represented
by consistent nodal forces, and extrapolated nodal stress recovery.

```console
cargo run -q -p solvix-cli -- run examples/plate_hole.inp --json
```

Verified stable result: 81 nodes, 64 elements, 77 PCG iterations, residual
`4.313e-11`, `max_disp_mm` 0.006, and `max_von_mises_mpa` 2.085. The JSON's
stress field is the maximum element-centroid value, which smooths the steep
hole-boundary gradient. The QA comparison metric is the maximum extrapolated
**nodal** von Mises stress: **2.927 MPa**, or 2.45% below the Kirsch
analytical stress-concentration result of **3.0 MPa**. See the
[Kirsch benchmark plot](benchmarks/kirsch-stress.svg) and
[analysis](benchmarks/kirsch.md) for this distinction.

## `broken.inp` - parser-failure smoke test

This deck contains one deliberately malformed node coordinate (`x.y`) in an
otherwise plausible one-element static model. It does not demonstrate a
physics solve; it verifies that parsing failures remain machine-readable and
point to the offending input line.

```console
cargo run -q -p solvix-cli -- run examples/broken.inp --json
```

The expected payload is:

```json
{"ok":false,"stage":"parse","error":"line 7: malformed coordinate 'x.y'"}
```

The command intentionally exits with code **2**, writes the JSON error on
stdout, and produces no VTU file. Under `set -e`, a shell will stop on this
expected nonzero status unless the caller handles it explicitly.

## `cantilever3d.inp` - solid-element static bending

This 3D cantilever is 2 m long with a 0.1 m square cross-section. Its
40x4x4 mesh uses 640 C3D8 linear hexahedral elements, a fully fixed root
face, and a total 10,000 N downward tip load distributed over 25 tip nodes.
It demonstrates three-dimensional geometry, solid-element assembly, and
3-displacement boundary conditions and results.

```console
cargo run -q -p solvix-cli -- run examples/cantilever3d.inp --json
```

Verified stable result: 1,025 nodes, 640 elements, 562 PCG iterations,
residual `9.367e-11`, `max_disp_mm` **13.719**, and maximum von Mises stress
**80.559 MPa**. The result is 9.96% below the 15.238 mm Euler-Bernoulli
reference, consistent with shear locking in the fully integrated linear
bricks at this aspect ratio. A rendered C3D8 viewport is available in
[`docs/screenshots/ui-3d.png`](screenshots/ui-3d.png).

## `thermal_bar.inp` - sequential heat transfer and thermal stress

This 0.2 m by 0.04 m steel strip uses 20x4 CPS4 elements. The left edge is
held at 100 °C and the right edge at 0 °C; the same left edge is mechanically
clamped. The deck demonstrates steady-state scalar heat conduction, Dirichlet
temperature boundaries, thermal expansion, and the sequential thermal-stress
solve.

```console
cargo run -q -p solvix-cli -- run examples/thermal_bar.inp --json
```

Verified stable result: 105 nodes, 80 elements, 52 PCG iterations, residual
`7.159e-11`, `max_temperature_c` **100.000**, `max_disp_mm` 0.123, and
maximum von Mises stress 193.603 MPa. With no internal source, the expected
one-dimensional profile is linear:

```text
T(x) = 100 * (1 - x / 0.2)
```

Thus the maximum occurs at x = 0, the minimum is zero at x = 0.2 m, and the
profile is independent of conductivity.

## `dyn_cantilever.inp` - transient Newmark dynamics

This deck reuses the 2D cantilever geometry and static tip load but adds
steel density. A 1,000 N downward load is applied as a step at time zero and
held. The solve takes 1,000 implicit Newmark steps at `dt = 1e-4 s` over
0.1 s, covering roughly two first-bending-mode periods without damping.

```console
cargo run -q -p solvix-cli -- run examples/dyn_cantilever.inp --json
```

Verified stable result: 205 nodes, 160 elements, `solver_method` `newmark`,
`dynamic: true`, `steps` 1,000, `dt` `1.0e-4`, `duration` 0.1,
`max_disp_mm` **27.468**, and `final_disp_mm` 9.327. The peak is 1.80 times
the 15.238 mm static equilibrium displacement, as expected for an undamped
suddenly applied load approaching the textbook factor-of-two asymptote. The
`0.000 MPa` von Mises value is expected in v1: dynamic stress recovery is not
part of this deck's final-time VTU output.

## Python and `abaqus2solvix` inventory

There are currently no `examples/*.py` files. The public
[roadmap](../ROADMAP.md) mentions the `abaqus2solvix` translator, but its
Python fixtures are internal crate tests rather than public examples, so they
are intentionally omitted from this gallery.
