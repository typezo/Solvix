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

## Issue #50 — ODB animations, XY plots, probe, PNG export

### API (main.rs)

- [x] `POST /api/run` now mirrors the CLI's `solve_auto` routing: NL decks
      (NLGEOM/plasticity/contact) call `solve_nl` directly so the converged
      increments are kept; dynamic runs already hold their output points.
- [x] Frame summaries (`FrameSummary { t, load_fraction, displacements,
      max_abs_u_mm }`) are built by pure functions `nl_frame_summaries`
      (t = load fraction) and `dyn_frame_summaries` (t = seconds,
      load_fraction = t/duration), both subsampled to a 200-frame cap with
      first/last always included; full per-node triples stay in memory in
      the job cache.
- [x] `GET /api/frames/{job}` →
      `{"ok":true,"frames":N,"points":[{"t":x,"max_abs_u_mm":y}],"has_full":bool}`
      (empty points and `has_full:false` for plain static jobs; 404 for
      unknown jobs).
- [x] `GET /api/frame/{job}/{i}` →
      `{"ok":true,"t":x,"load_fraction":f,"displacements":[[x,y,z] per node]}`;
      out-of-range `i` is a 404 carrying the frame count.
- [x] `/api/agent` documents both new routes; its `body` field is now
      JSON-escaped, fixing pre-existing invalid JSON when a body string
      contained quotes (verified with a strict JSON parse over HTTP).
- [x] Tests: summary value math, 300→200 subsampling with monotone t for
      both builders, an end-to-end NL run (`*STEP, NLGEOM=YES` → 20 frames,
      t = 0.05..1.0) through `/api/run` → `/api/frames` → `/api/frame`
      (including out-of-range 404 and re-run dropping frames), a dynamic
      run capped at 200 frames with the final time preserved, and router
      dispatch for the new paths (`..`/`/`/non-numeric index rejected).

### Frontend (static/)

- [x] After an API results load, `GET /api/frames/{job}` populates the
      player; frame triples are fetched lazily from `/api/frame/{job}/{i}`
      and cached per job. The issue #47 dynamic frames remain the fallback.
- [x] Animation controls: play/pause button + slider at ~20 fps
      (`setInterval` 50 ms, wrapping), reusing the existing deform pipeline
      (`activeDisplacements` prefers the selected frame's triples).
- [x] Hand-rolled SVG XY plot of |u|max vs t below the viewport: frame box,
      4 labeled ticks per axis, curve polyline, cyan current-frame marker
      line that tracks the slider; hidden when fewer than 2 points.
- [x] Click-probe in BOTH renderers: a click (pointer down/up with < 4 px
      movement) runs the shared screen-space `nearestOf` search over the
      displayed nodes (2D `worldToScreen` or 3D projected nodes); the
      tooltip draws node id, |u| in mm, and nodal von Mises in MPa when
      present, follows rotation/deformation, and clears on empty-space
      clicks. Dragging/rotating never triggers a probe.
- [x] PNG export button (`PNG`) downloads the current canvas via
      `toDataURL("image/png")`; the filename carries the job and frame.
- [x] Frame label distinguishes `t = X s` (transient jobs) from `λ = X`
      (NL load fraction); the viewport label, result chips, and the
      viewport copy-JSON payload include frame/playing/probe state.
- [x] All copy-JSON buttons remain; new pure helpers (`xyPlotSVG`,
      `nearestOf`, `frameLabelAt`) are exposed under
      `window.__SOLVIX_EXPOSE_INTERNALS__`.

### Browser-less and live checks

- [x] `node --check crates/solvix-ui/static/app.js` passes (Node v24.21.0).
- [x] HTML id uniqueness and CSS brace balance re-checked (52 ids, all
      blocks closed).
- [x] Node VM startup smoke test with DOM/Canvas stubs: embedded-mock
      startup renders, `xyPlotSVG` emits the curve/marker/8 tick labels
      and returns `""` for short series, `nearestOf` picks the nearest
      entry and handles the empty list.
- [x] Live server over HTTP: static assets 200, `dyn_cantilever` run →
      200 capped frames with correct point JSON, `/api/agent` parses as
      strict JSON.
- [x] `cargo test --workspace` — all green; `cargo fmt --check` — passes;
      `cargo clippy --workspace --all-targets -- -D warnings` — CLEAN.
