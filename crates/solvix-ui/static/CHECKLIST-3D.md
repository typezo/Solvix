# 3D viewport verification checklist

Issue #20 — Abaqus-style 3D viewport with an orthographic Canvas renderer.

## Scope

- [x] Implementation changes are confined to `crates/solvix-ui/static/`.
- [x] Existing `mock/model.json` and `mock/results.json` are unchanged.
- [x] Added `mock/model3d.json` and `mock/results3d.json`; no other mock files
      were modified.
- [x] No dependency, package manifest, bundler, or Rust-source change was added.
- [x] No version-control mutation was performed.

## Mock data

- [x] `model3d.json` parses and contains a 4 × 1 × 1 m C3D8 beam with 30 nodes
      (5 × 2 × 3 grid), 8 hexahedra, fixed root nodes, and six downward tip
      loads totaling 10,000 N.
- [x] `results3d.json` contains 30 three-component displacement rows and 8
      element von Mises values.
- [x] Root displacement is zero and downward displacement magnitude increases
      monotonically from the root toward the tip.
- [x] Reported maxima agree with the arrays: 0.012190 mm and 0.238 MPa.
- [x] The embedded final fallback is the same 30-node mock data, so `file://`
      startup still works when a browser denies local `fetch` access.

## 3D detection and data path

- [x] A model switches to the 3D renderer when any node has nonzero `z` or any
      element is `c3d8`/`c3d4`.
- [x] A model with all-zero `z` and plane element kinds remains on the existing
      2D renderer and retains 2D pan/cursor/zoom behavior.
- [x] Startup defaults to job `cantilever3d`.
- [x] The normal API pipeline remains in use:
      `GET /api/model?job=cantilever3d`, then `GET /api/results/cantilever3d`.
- [x] On a fresh server with no cached result, startup detects the missing
      result and invokes `POST /api/run {"path":"examples/cantilever3d.inp"}`;
      the browser-less smoke test reached `model: api · results: api`.
- [x] The real deck test loaded 1,025 nodes and 640 C3D8 elements and rendered
      its API result in deformed von Mises mode.
- [x] API transport failure selects the matching 3D mock pair. A compact mock
      result is never attached to the full API-sized model.
- [x] The current API emits legacy `[id,x,y]` nodes for C3D8 decks. A
      connectivity-based compatibility path reconstructs the missing depth levels
      from each hexahedron's bottom/top face pairs; on the real deck it produces
      a centered z range of −0.05…+0.05 m. New payloads containing `[id,x,y,z]`
      use their z coordinates directly.

## Renderer math and viewport

- [x] Matrix multiplication, transpose, vector multiplication, axis rotations,
      unit-vector rotations, and arcball sphere mapping are pure functions.
- [x] Orthographic projection and fit-to-view calculations are pure functions.
- [x] Pure-function checks passed for identity multiplication, a 90° rotation,
      center projection, fit bounds, pan, centered zoom, and unit arcball vector.
- [x] Ordinary pointer drag performs arcball rotation in 3D mode.
- [x] Shift-drag, middle-drag, or right-drag pans in 3D; right-click does not
      open the browser context menu.
- [x] Wheel zoom is anchored at the pointer and clamped to a finite scale.
- [x] Fit resets around the complete 3D bounds and accounts for legend space.
- [x] The 3D axis triad projects X/Y/Z with the active rotation.
- [x] Element edges are painter-sorted by average camera z before drawing.
- [x] C3D8 draws exactly 12 edges and C3D4 exactly 6 edges (both counts tested).
- [x] Deformed geometry uses the existing displacement indexing and 3D bbox
      diagonal; the Results module retains the default 5% deformation setting.
- [x] Von Mises mode colors the element edges and reuses the six-tick jet
      legend with MPa units.
- [x] Displacement-magnitude mode adds depth-sorted colored nodal point
      overlays and an appropriately labeled legend.
- [x] The dashed undeformed 3D wireframe can be toggled independently.
- [x] 3D node/element labels, selected nset rings/elset edges, fixed-DOF
      glyphs, and CLOAD arrows for DOFs 1–3 are supported.
- [x] V1 intentionally uses colored edges rather than shaded element faces.
      Displacement vector arrows (optional in the issue) are not included.

## UI compatibility

- [x] The module bar, model tree, JSON console/REPL, result chips, and status
      bar remain present.
- [x] The viewport label distinguishes `3D · Deformed/Undeformed` from the
      existing 2D labels.
- [x] The status bar reports renderer dimension and element-kind counts, for
      example `3D · C3D8 × 8` and `2D · CPS4 × 160`.
- [x] The job selector includes both `cantilever3d` and the original plane
      `cantilever` jobs; loading the latter automatically returns to 2D.
- [x] The plane API regression loaded 205 nodes / 160 CPS4 elements and its
      205-row result, then rendered the existing 2D deformed von Mises path.

## Static and repository checks

- [x] `node --check crates/solvix-ui/static/app.js` passes with Node v24.21.0.
- [x] All four mock JSON files parse successfully.
- [x] HTML parser check passes: 47 unique required IDs, no unclosed elements.
- [x] CSS structural check passes: 97 balanced declaration blocks.
- [x] Local server returned HTTP 200 for the HTML, JavaScript, CSS, both new
      3D mocks, and both unchanged plane mocks.
- [x] Browser-less VM smoke tests passed for file-style embedded 3D fallback,
      real API 3D data, and real API 2D data.
- [x] No graphical browser executable was available, so interactions were
      verified through pure-function tests plus DOM/canvas stub event/render
      traces rather than a screenshot.
- [x] `git diff --check` passes.
- [x] `cargo fmt --check` passes.
- [x] `cargo clippy --workspace --all-targets -- -D warnings` passes.
- [x] `cargo test --workspace` passes: 97 tests (47 CLI, 24 core, 8 I/O,
      18 UI) plus workspace doc-test suites.
