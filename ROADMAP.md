# Solvix Roadmap

**North star:** an open-source solver infinitely close to Abaqus in capability,
and better than Abaqus at being driven — by scripts, by AI agents, by anyone.
Every release advances both pillars:

- **Parity pillar** — what Abaqus does: elements, analysis types, materials,
  output database, CAE-style GUI.
- **AI-native pillar** — what Abaqus can't: JSON/REST everything, MCP tools,
  self-describing APIs, self-verifying releases, zero-install single binary.

Policies (all versions): zero external deps in workspace crates · JSON-first UX ·
benchmark-gated releases (QA agent) · one issue owns one directory ·
`cargo test/clippy/fmt` green at every merge.

## Version ladder

### v0.1.0 — Linear static 2D + agent CLI + web UI ✅ (2026-10-10)
CPS3/CPS4 plane stress · `.inp` subset parser · VTU out · `solvix run/info --json` ·
`solvix serve` REST + Abaqus-style canvas UI · cantilever + Kirsch benchmarks
(13.832→14.877 mm converging to 15.238 Euler–Bernoulli; exit codes 0/1/2/3).

### v0.1.x — Harden the flagship
- [x] README/UI screenshots + UI integration smoke in CI (#16)
- [x] `solvix bench` built-in verification suite (#17)
- [x] Cross-platform release binaries: linux/windows/macos (#33)
- [ ] `solvix --version`, job naming polish, VTU binary format option

### v0.2.0 — Into 3D
- [x] C3D8 hex + C3D4 tet, 3 dof/node (#6)
- [x] 3D BCs/loads, VTU volume cells, patch + cube benchmarks
- [x] Direct sparse LDLT solver fallback + condition estimate (#18)
- [x] UI 3D viewport: software renderer, rotate/pan/zoom (#20)

### v0.2.x — Real Abaqus decks
- [x] `*DLOAD` edge pressure, gravity, multi-step OP=NEW/MOD, *DENSITY (#26)
- [x] Gmsh `.msh` import (#19)
- [x] MCP server crate: tools run/info/model/results/bench (#8)
- [ ] Keyword coverage matrix auto-published to docs (#32)

### v0.3.0 — Dynamics of structures
- [x] Modal: subspace iteration, consistent mass, mode shapes in UI + CLI (#23 #40)
- [x] Buckling: linear eigenvalue on prestressed K (#24)
- [ ] Restart: binary state files, resume jobs

### v0.3.x — Element library breadth
- [x] B31 beam, CPS8/CPE8 quadratic plane strain (#25)
- [ ] Plane strain/axisymmetric element families
- [ ] Orthotropic + thermal expansion material models (#22)

### v0.4.0 — Multiphysics I: heat
- [x] Steady-state conduction (scalar quads/hexes/tets), flux BCs (#27)
- [x] Thermal stress coupling (sequential, E*a*dT exact to 5e-11) (#27)

### v0.4.x — Multiphysics II: time
- [x] Transient dynamics: Newmark implicit, Rayleigh damping (#28)
- [ ] Thermal transient (v0.5 candidate)

### v0.5.0 — Nonlinearity
- [x] Geometric: Updated Lagrangian + Newton–Raphson (#29)
- [x] Material: J2 plasticity, radial return (#30)
- [x] Contact 2D: penalty node-to-surface (#31)

### v0.6.0 — Output database + history
- [ ] `.solvix` binary ODB: frames, field/history outputs, random access
- [ ] UI: XY plots, animations, PNG export (#24)
- [ ] Multithreaded assembly + PCG (std threads)

### v0.7.0 — Scale + hardening (shipped ahead of plan)

- [x] Threaded assembly + SpMV, bit-identical (#42)
- [x] NAFEMS-style manifest benchmarks + CI regression gate (#34)
- [x] Binary VTU appended format (#43)
- [x] Abaqus `.frd` reader (#37) · keyword matrix + error taxonomy (#22)
- [x] Gallery site (#10) · examples gallery (#44)

### v0.7.x — Parked-solution upgrades

- [ ] Plane-stress J2 return, Timoshenko, friction, free-free shift (#55 #56)
- [ ] *AMPLITUDE, *INITIAL CONDITIONS, restart files (#57)

### v0.8.0 — CAE Part module parity
- [ ] Parametric geometry: 2D sketch booleans, extrude
- [ ] Meshing: Delaunay 2D + structured 3D, size control
- [ ] Assembly instances, display groups in UI

### v0.8.0 — Migration path for Abaqus users
- [x] `abaqus2solvix`: script translator with --validate (#35, pulled ahead from v0.8)
- [x] `.frd` reader for post-processing existing results (#37)
- [ ] Error taxonomy: agent-readable failure codes (#32)

### v0.9.0 — Scale
- [ ] Direct sparse solver with pivoting, out-of-core option
- [ ] Domain decomposition + threaded solve at scale
- [x] NAFEMS-style benchmark suite as CI regression gate (#34)

### v1.0.0 — Parity declaration
Full docs, tutorial gallery, published benchmark comparisons vs Abaqus
tutorial results, LTS support policy. **FIGHT FOR OPEN SOURCE.**

## The AI-native stack (continuous)

| Layer | Mechanism | Status |
|---|---|---|
| CLI | `--json` everywhere, exit codes 0/1/2/3 | ✅ v0.1 |
| REST | `/api/*` + `/api/agent` self-description | ✅ v0.1 |
| UI console | REPL proxying the API, copy-JSON buttons | ✅ v0.1 |
| MCP | `solvix-mcp` server crate | v0.2.x |
| Self-test | `solvix bench` verification suite | v0.1.x |
| Failure codes | structured error taxonomy in every layer | v0.8 |
