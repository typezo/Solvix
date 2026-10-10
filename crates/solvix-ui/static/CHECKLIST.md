# Solvix UI verification checklist

Issue #14 — Abaqus-CAE-style frontend, vanilla JavaScript + Canvas.

## Implementation and scope

- [x] Added only `index.html`, `style.css`, `app.js`, and this checklist under
      `crates/solvix-ui/static/`.
- [x] No npm dependency, package file, bundler configuration, or build step.
- [x] Both files under `static/mock/` remain byte-for-byte unchanged.

## Static/browser-less checks

- [x] `node --check crates/solvix-ui/static/app.js` passes with Node v24.21.0.
- [x] HTML parse check: all required elements resolve, IDs are unique, and no
      element is left open.
- [x] CSS structural check: comments/strings excluded, all 96 `{` blocks have a
      matching `}`.
- [x] `index.html`, `app.js`, `style.css`, and both mock JSON files all return
      HTTP 200 from a local static file server.
- [x] Both mock JSON files parse successfully.
- [x] Ran the startup path in a Node VM with lightweight DOM/Canvas stubs and a
      missing API server: it fetched `static/mock/model.json` and
      `static/mock/results.json`, normalized 2 nodes and 1 element, rendered the
      no-API model path, and reported both sources as mock data.
- [x] Served the UI from the workspace `solvix-ui` server and verified the HTML,
      JavaScript, CSS, `/api/model?job=cantilever`, and `/api/agent` routes.
- [x] No graphical Chrome/Firefox executable was available in this environment,
      so pixel output was not screenshotted. Interaction and rendering branches
      were checked by code trace plus the browser-less startup/render smoke test.

## Layout and model database

- [x] Left model-tree panel contains Model DB groups for Nodes, Elements, Sets,
      Materials, BCs, Loads, Steps, and Jobs.
- [x] Top module bar exposes Model, Mesh, Job, and Results states.
- [x] Bottom console/REPL and persistent status bar are present.
- [x] Status bar reports active job, node/element counts, m/Pa units, cursor
      coordinates, and zoom in px/m.
- [x] The viewport uses the requested dark `#1a1a2e` background with Abaqus-like
      tan/yellow mesh strokes and dark control panels.

## Viewport behavior

- [x] Draws element edges from element connectivity; node symbols are optional.
- [x] Supports node-number and element-ID text overlays.
- [x] Fit calculates model bounds with margins and also accounts for the legend.
- [x] Supports pointer-capture dragging for pan.
- [x] Wheel zoom keeps the world coordinate under the cursor fixed and clamps the
      scale to a finite range.
- [x] Draws +X/+Y coordinate axes.
- [x] Draws a scale bar with a round physical length in m/mm/µm/nm.
- [x] Renders fixed-DOF glyphs as orange squares.
- [x] Renders CLOAD glyphs as red arrows in the positive/negative DOF 1 or DOF 2
      direction. Missing referenced nodes are skipped rather than crashing the
      renderer (the tiny mock deck intentionally references absent node 205).
- [x] Nset/elset checkboxes highlight selected node rings and element regions.
- [x] Undeformed wireframe can be toggled independently from the deformed mesh.

## Results behavior

- [x] Result displacements are indexed in mesh node order.
- [x] Deformation factor is `slider_percent / 100 × model_bbox_diagonal /
      max_|u|`; default Results-module value is 5%.
- [x] Deformation slider range is 0–20% and its label shows the active value.
- [x] Element von Mises contours use element-constant values.
- [x] `|U|` contour mode uses nodal displacement magnitudes and interpolates each
      quadrilateral/triangle region.
- [x] Both contour modes use the implemented jet colormap.
- [x] Contour legend has six numeric ticks and displays MPa for von Mises or mm
      for displacement magnitude.
- [x] Result chips expose maximum displacement, maximum von Mises, iteration
      count, active contour, and deformation percentage.
- [x] Length mismatches between model/result arrays produce a JSON data-consistency
      warning instead of an uncaught exception.

## API, job, console, and JSON UX

- [x] Initial model request is `GET /api/model?job=NAME`; initial result request
      is `GET /api/results/NAME`.
- [x] Fetch/HTTP errors fall back to `static/mock/model.json` and
      `static/mock/results.json`.
- [x] `file://` startup has an embedded copy of the exact mock payload as a final
      fallback for browsers that deny local `fetch` access.
- [x] Deck path input and Run button submit `{"path":"..."}` with `POST /api/run`.
- [x] A successful run response adds its job to the dropdown and loads model and
      result payloads for that job.
- [x] Console pretty-prints every request/response as JSON.
- [x] REPL supports `run <path>`, `model <job>`, `results <job>`, `agent`,
      `help`, and `clear`.
- [x] REPL provides command buttons, Tab completion, and arrow-key history.
- [x] Copy-JSON controls cover model tree sections, complete model/result data,
      job/run state, viewport state, console transcript, and individual console
      entries. Clipboard has a non-secure-context fallback.

## Repository verification

- [x] `cargo test --workspace` — all 52 tests pass.
- [x] `cargo fmt --check` — passes.
- [x] `cargo clippy --workspace --all-targets -- -D warnings` — passes.
- [x] `git diff --check` for files under `crates/solvix-ui/static/` — passes.
