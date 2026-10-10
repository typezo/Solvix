# Cantilever convergence benchmark

This benchmark exercises Solvix's CPS4 plane-stress elements on the bundled
2 m cantilever. The left edge is fully fixed, the total downward tip load is
1,000 N, and the deck uses SI units. The benchmark's Euler-Bernoulli target is
**15.238 mm**.

The runs below used `cargo run -q -p solvix-cli` (a debug build) with a
Jacobi-preconditioned conjugate-gradient relative residual tolerance of
`1e-10`. `solve_ms` measures assembly, solve, and stress recovery; it is a
single-run timing and will vary by machine. The 10×4, 20×4, 80×8, and 160×8
decks were generated from the same row-major rectangular pattern as
`examples/cantilever.inp`. Five equally spaced tip nodes receive `-200 N` on
four-row meshes; nine receive `-111.111 N` on eight-row meshes, for a nominal
1,000 N resultant. Active DOFs exclude the two constrained directions on every
node of the fixed edge.

| Mesh (CPS4) | Nodes | Active DOFs | Iterations | `solve_ms` | `max_disp_mm` | Error vs. 15.238 mm |
|------------:|------:|------------:|-----------:|-----------:|--------------:|--------------------:|
| 10×4 | 55 | 100 | 201 | 16 | 5.985 | −60.72% |
| 20×4 | 105 | 200 | 315 | 47 | 10.953 | −28.12% |
| 40×4 | 205 | 400 | 390 | 119 | 13.832 | −9.23% |
| 80×8 | 729 | 1,440 | 838 | 911 | 14.877 | −2.37% |
| 160×8 | 1,449 | 2,880 | 943 | 2,040 | 15.147 | −0.60% |

The CLI reports the magnitude of the maximum nodal displacement. For this
load case that maximum occurs at the tip, so `max_disp_mm` is the tip
deflection.

![CPS4 tip displacement versus span refinement](cantilever-convergence.svg)

## Shear-locking trend

Every simulated result lies below the reference and approaches it from below
as the mesh is refined. This is the expected signature of parasitic shear in
fully integrated bilinear quadrilaterals: a rectangular CPS4 element cannot
represent the ideal pure-bending displacement field without generating a
spurious shear strain. That extra shear energy makes the mesh too stiff and
underpredicts deflection.

The effect is strongest on the coarse, high-aspect-ratio elements. Refining
the span from 10 to 40 elements while retaining four rows improves the result
from 60.72% low to 9.23% low. Increasing to eight rows and continuing span
refinement reduces the error further, to 0.60% at 160×8. The monotonic
increase in displacement—and the residual below `1e-10` on every run—supports
convergence of the solver while documenting the known element behavior.
