# Keyword coverage matrix and error taxonomy

Source of truth: `crates/solvix-io/src/parser.rs`. This file is kept in
sync with the parser by the `keyword_docs_in_sync` test in
`crates/solvix-io/tests/keywords_doc.rs` — every keyword documented here
must have a match arm in the parser and vice versa.

Line model: `**…` comments and blank lines are dropped; every remaining
line is a keyword line (`*NAME, KEY=VALUE, …`, case-insensitive) or a
comma-separated data line owned by the preceding keyword. Set names,
material names, and surface names are matched case-insensitively and
stored uppercased.

Status legend:

- **supported** — parsed into the model; exit is unaffected.
- **supported (partial)** — parsed, but some parameters or forms are
  accepted-with-warning.
- **warned + skipped** — block is dropped and a warning is recorded;
  parsing continues (exit stays 0 unless something else fails).
- **skipped** — data lines consumed silently, nothing recorded.

Any keyword not in the table hits the catch-all arm: warning
`unsupported keyword '*NAME' — block skipped`.

## Keyword matrix

| Keyword | Status | Parameters accepted | Semantics | Notes |
|---------|--------|---------------------|-----------|-------|
| *HEADING | skipped | none | Data lines consumed, no effect. | Safe to keep deck titles. |
| *NODE | supported | none | Data: `id, x, y[, z]`; z defaults to 0. | 3 or 4 fields required; duplicate ids and malformed numbers are fatal. |
| *ELEMENT | supported | `TYPE=`, `ELSET=` | Data: `id, n1, …` (one element per record; a trailing comma continues connectivity on the next line). `ELSET=` upserts the element into a set. | TYPE limited to CPS3, CPS4, CPS8, CPE4, CPE8, C3D8, C3D4, B31. Any other TYPE (or a missing TYPE) warns and skips the whole block. Field-count mismatches and duplicate ids are fatal. |
| *ELSET | supported | `ELSET=<name>`, `GENERATE` | Data: id list, or `start, end[, step]` with GENERATE (step defaults to 1). | Missing ELSET= name warns and skips the block. GENERATE step 0 and bad field counts are fatal. |
| *NSET | supported | `NSET=<name>`, `GENERATE` | Same as *ELSET for node sets; later blocks append. | Same fatal/warn rules as *ELSET. |
| *MATERIAL | supported | `NAME=<name>` | Opens a material block; property keywords (*ELASTIC, *DENSITY, *PLASTIC, *CONDUCTIVITY, *EXPANSION, *DAMPING) apply to it until the next keyword outside the block. | Missing NAME warns and ignores the block. A material that reaches the next keyword without *ELASTIC warns. |
| *ELASTIC | supported | none | Data: `youngs_modulus, poisson_ratio`. | Exactly 2 fields or fatal. Outside *MATERIAL warns and is ignored. Extra rows warn (temperature dependence is future work). |
| *DENSITY | supported | none | Data: `density`. One value covers the whole model (needed for gravity, dynamics, modal). | Outside *MATERIAL warns. Overwriting a different value warns; extra fields/rows warn (temperature dependence is future work). |
| *PLASTIC | supported | none | Data rows: `sigma, eps_pl`. Reduced to J2 linear hardening: yield from the first row, hardening modulus H from first/last rows (single row means H = 0). | Requires sigma > 0, eps_pl >= 0, non-decreasing eps_pl (violations fatal). First/last eps_pl must differ when multiple rows. Outside *MATERIAL warns and is ignored; multiple blocks warn and overwrite one model-wide law. |
| *CONDUCTIVITY | supported | none | Data: `conductivity`. One value covers the whole model (thermal solves). | Outside *MATERIAL warns. Overwrite/extra-field/extra-row warnings like *DENSITY. |
| *EXPANSION | supported | none | Data: `thermal expansion coefficient`. One value covers the whole model (thermomechanical strain). | Same warnings as *CONDUCTIVITY. |
| *DAMPING | supported | `ALPHA=`, `BETA=` | One Rayleigh pair for the whole model. Either both parameters, or a data line `alpha, beta`. Step-level placement overrides the material-level default; applies to *DYNAMIC steps. | Only ALPHA=/BETA= are honored; other parameters warn. ALPHA= without BETA= (or neither form present) is fatal. In a step without *DYNAMIC it warns and is dropped. |
| *SOLID SECTION | supported | `ELSET=`, `MATERIAL=` | Assigns the material to an element set. Optional data line: `thickness` (defaults to 1.0 when absent). | Missing ELSET/MATERIAL warns and skips. Unknown ELSET/MATERIAL names are fatal. Thickness data must be exactly 1 field. |
| *BEAM SECTION | supported | `ELSET=`, `MATERIAL=` | Assigns the material to a B31 element set. Required data: `area A, second moment I`. | Missing ELSET/MATERIAL warns and skips; unknown names fatal; missing or malformed data line fatal. |
| *SURFACE | supported | `NAME=`, `TYPE=NODE` | Data: node-set names and/or bare node ids; the union becomes the surface's node list (v0.5: NODE surfaces only). | Missing NAME, or any TYPE other than NODE (e.g. ELEMENT), warns and skips the block. Empty data fields are fatal. |
| *CONTACT PAIR | supported | `INTERACTION=` (warned, ignored) | Data: `slave_surface, master_surface`; each resolves to a *SURFACE name or an NSET name. Pairs live at model level, active in all steps, frictionless. | Field count must be 2. Unknown surface names are fatal. |
| *BOUNDARY | supported | `OP=NEW`/`OP=MOD` | Data: `node/set, dof1[, dof2][, value]` — homogeneous constraints; dof ranges expand. Inside *STEP it accumulates into the step's total BC state; OP=NEW clears prior BCs in the step. | Only homogeneous values: a nonzero prescribed value warns and is ignored. OP outside NEW/MOD warns (treated as MOD). Outside steps after the first, or OP=NEW outside *STEP, warns and is ignored. |
| *CLOAD | supported | `OP=NEW`/`OP=MOD` | Data: `node/set, dof, value`; set references expand to one load per node. Step accumulation and OP semantics as *BOUNDARY. | Exactly 3 fields. Unknown node sets fatal. |
| *DLOAD | supported | `OP=NEW`/`OP=MOD` | Inside *STEP only. Data: `element/elset, P1..P4, magnitude` (face pressure on plane elements) or `GRAV, magnitude, dir1[, dir2][, dir3]`. OP=NEW clears prior dloads and gravity in the step. | Outside *STEP warns and is ignored. P faces beyond the element's face count, unknown elements, unknown targets, and non-P/GRAV types warn per line. |
| *TEMPERATURE | supported | `OP=NEW`/`OP=MOD` | Inside *STEP only. Data: `node/set, temperature` — nodal temperature BCs for thermal/thermomechanical solves. | Outside *STEP warns and is ignored. Exactly 2 fields; unknown node sets fatal. |
| *DFLUX | supported | `OP=NEW`/`OP=MOD` | Inside *STEP only. Data: `element/elset, S1..S4, magnitude` (surface heat flux on plane element edges) or `element, BN, magnitude` (one body flux for the whole model). | Outside *STEP warns and is ignored. Non-S/BN labels, flux on non-plane elements, unknown elements/targets warn per line; BN overwrite warns. |
| *STEP | supported | `NAME=`, `OP=NEW`/`OP=MOD`, `NLGEOM=YES`/`NO` | Opens a step. MOD (default) inherits the previous step's total BC/load state; NEW starts clean. NAME defaults to `Step-N`. NLGEOM=YES routes to the nonlinear solver for that step only — never inherited. | A missing *END STEP is tolerated (EOF closes the step). Bad OP or NLGEOM values warn (treated as MOD / NO). |
| *STATIC | supported | none (any warn) | Marks the step as linear-static. Data lines consumed. | All parameters warn (v0.1: linear static only). |
| *DYNAMIC | supported | none (all warn) | Inside *STEP only. Data line 1: `dt, duration` (implicit Newmark, constant dt). Optional data line 2: `alpha, beta` — the Rayleigh pair inline. | Every keyword parameter warns (e.g. EXPLICIT). Field-count violations fatal. Extra rows warn. Outside *STEP warns and is ignored. |
| *END STEP | supported | none | Closes the open step. | Data lines consumed. |
| *HEAT TRANSFER | supported (partial) | `STEADY STATE` | Marks the step as a steady-state thermal solve. | Any other parameter (e.g. TRANSIENT) warns and is ignored. Data lines consumed. |
| *SURFACE INTERACTION | warned + skipped | none | Not supported (v0.5: frictionless contact only). | Block dropped with a warning. |
| *FRICTION | warned + skipped | none | Not supported (v0.5: frictionless contact only). | Block dropped with a warning. |
| *CFLUX | warned + skipped | none | Not supported (v0.4: *TEMPERATURE and *DFLUX S/BN only). | Block dropped with a warning. |

Coverage summary: 29 keywords in the parser match — 24 supported
(including *DYNAMIC and *DAMPING with reduced-parameter forms),
1 supported with partial parameters (*HEAT TRANSFER), 3 warned + skipped
(*SURFACE INTERACTION, *FRICTION, *CFLUX), 1 skipped silently
(*HEADING). Every unrecognized keyword warns and is skipped via the
catch-all arm. 8 element types parse: CPS3, CPS4, CPS8, CPE4, CPE8,
C3D8, C3D4, B31.

## Error taxonomy

### Parse stage (exit code 2)

`ParseError { line, message }` — `line` is the 1-based deck line.
Agents should read `line`, jump to it in the deck, and repair.

| Class | Typical message | Agent recipe |
|-------|-----------------|--------------|
| Malformed number | `malformed <what> '<token>'` | The field is not a number in the expected type. Check for typos, units pasted as text, or a comma mis-split; fix the field on `line`. |
| Field-count violation | `<keyword> expects …; got N field(s)` | The data line has too few/many comma fields for the keyword. Compare against the "Semantics" column above and re-split the row. |
| Duplicate id | `duplicate node id N` / `duplicate element id N` | Two mesh entities share an id. Renumber one of them (and any references in sets/connectivity). |
| Unknown reference | `unknown node set 'X'` / `unknown surface 'X'` / `*SOLID SECTION references unknown ELSET 'X'` / `… unknown MATERIAL 'X'` | The name is not defined before use (sets must precede use; surfaces must be TYPE=NODE). Define the set/surface/material earlier or fix the spelling — names are case-insensitive. |
| Missing required data | `*DENSITY expects a density value`, `*PLASTIC expects data rows …`, `*BEAM SECTION expects a data line 'A, I'`, `*DAMPING needs ALPHA=/BETA= or a data line …` | The keyword block has no (or incomplete) data. Add the data line in the documented shape. |
| Value-domain violation | `*PLASTIC needs positive finite sigma …`, `*PLASTIC eps_pl rows must be non-decreasing`, `GENERATE step must be >= 1`, `*DAMPING needs both ALPHA= and BETA=` | Values violate documented invariants. Correct the numbers (e.g. sort *PLASTIC rows by eps_pl). |

Non-fatal notes (`ParseReport.warnings`, exit still 0): skipped
unsupported keywords/blocks, ignored parameters, ignored nonzero
prescribed displacements, overwrite notices. Agents should surface
these to the user because they change what the solve means even though
it runs.

### Solve stage (exit code 3)

`SolveError` (Debug-formatted in diagnostics):

| Class | Raised when | Agent recipe |
|-------|-------------|--------------|
| `EmptyMesh` | The model has no elements. | The deck defined no supported elements — check for *ELEMENT blocks that were skipped as unsupported TYPE (see parse warnings). |
| `NoConstraints` | No constrained DOFs (singular system). | Add or repair *BOUNDARY entries; verify the referenced NSET is defined and non-empty. |
| `InvalidModel(String)` | Model validation failed (bad material data, mixed plane/solid dimensionality, missing section assignment, NL-unsupported element kind, non-positive density, bad modal/buckling options, …). | The message names the exact invariant; fix the named deck field. These are deck problems, not solver bugs. |
| `NotConverged { iterations, residual }` | Iterative process exhausted its budget: PCG (after the LDLT fallback path), Newton increments in NL solves, or eigen iterations. | Inspect load magnitudes vs. stiffness and constraint quality; for NL steps reduce load step severity or refine the mesh. |

Exit codes: 0 success, 1 usage/IO, 2 parse (`ParseError`), 3 solve
(`SolveError`). JSON payloads go to stdout; diagnostics to stderr. See
the [Agent Cookbook](AGENT_COOKBOOK.md) for copy-ready recipes.
