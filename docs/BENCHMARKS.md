# Solvix parity scoreboard

Analytical regression gates for the solver (issue #34): every case
compares one metric of a solved deck against a closed-form reference
with a per-case tolerance. Cases live in `benchmarks/manifest.toml`
(decks in `examples/` and `benchmarks/decks/`); `solvix bench` runs
the built-in suite plus the manifest.

## Workflow

Run the suite locally (same command CI runs in
`.github/workflows/bench.yml` on every push/PR — any FAIL exits 3
and fails the build):

```console
cargo run -q -p solvix-cli -- bench --manifest benchmarks/manifest.toml
```

JSON for agents: add `--json`. Regenerate this file after a
deliberate solver change (human-run, not CI):

```console
cargo run -q -p solvix-cli -- bench --manifest benchmarks/manifest.toml --markdown > docs/BENCHMARKS.md
```

## Scoreboard (last local run)

| case | metric | value | expected | error | verdict |
|------|--------|-------|----------|-------|---------|
| cantilever40x4 | max_disp_mm | 13.832 mm | 15.238 | -9.23% | PASS |
| kirsch | max_nodal_von_mises_mpa | 2.927 MPa | 3.000 | -2.45% | PASS |
| patch_cps4 | max_centroid_von_mises_mpa | 1.000 MPa | 1.000 | +0.00% | PASS |
| cantilever40x4 | max_disp_mm | 13.832 mm | 15.238 | -9.23% | PASS |
| kirsch8x8 | max_nodal_von_mises_mpa | 2.927 MPa | 3.000 | -2.45% | PASS |
| patch_cps4 | max_centroid_von_mises_mpa | 1.000 MPa | 1.000 | +0.00% | PASS |
| cantilever3d | max_disp_mm | 13.719 mm | 15.238 | -9.97% | PASS |
| thermal_bar | max_temperature_c | 100.000 C | 100.000 | +0.00% | PASS |
| b31_cantilever | max_uy_mm | 158.730 mm | 158.730 | +0.00% | PASS |
| modal_cantilever | freq_1_hz | 21.916 Hz | 20.886 | +4.93% | PASS |
| euler_column | pcr_n | 421312.490 N | 432081.000 | -2.49% | PASS |

11/11 passed.

## References

- **cantilever40x4** — Euler-Bernoulli tip deflection P*L^3/(3*E*I) = 15.238 mm
- **kirsch8x8** — Kirsch: 3*sigma at the hole edge for 1 MPa far field
- **patch_cps4** — distorted CPS4 patch under uniform traction: exact 1 MPa
- **cantilever3d** — same cantilever in C3D8 solids (bending locks, wider band)
- **thermal_bar** — steady 1D gradient T(x) = 100*(1 - x/L): max T = 100 at the clamp
- **b31_cantilever** — B31 Hermite cantilever: P*L^3/(3*E*I) = 158.730 mm (exact)
- **modal_cantilever** — first bending mode 1.875^2/(2*pi)*sqrt(E*I/(rho*A))/L^2 = 20.886 Hz
- **euler_column** — pin-ended Euler column Pcr = pi^2*E*I/L^2 = 432081 N
