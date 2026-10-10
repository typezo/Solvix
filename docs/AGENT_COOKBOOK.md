# Agent Cookbook

Recipes for AI agents and scripts that drive the `solvix` binary safely. The
CLI's machine contract is simple: ask for `--json`, capture the exit code, keep
stdout reserved for the JSON payload, and treat VTU files as downstream
artifacts named by that payload.

## Core rules

1. Pass `--json` explicitly for every operational command.
2. Never parse stdout as human text when `--json` is set.
3. Check the exit code first, then validate that stdout contains one JSON
   object with the expected `ok` value.
4. Do not use `set -e` around Solvix calls when a non-zero status is part of
   your control flow. Capture `$?` and branch on it.
5. Run `info` before `run` for an unfamiliar or externally generated deck.
6. Treat solver output files as opaque until the run exits; do not read a VTU
   file while the solver may still be writing it.

## Exit-code handling

| Code | Class | Agent action |
|-----:|-------|--------------|
| `0` | Success | Parse JSON and inspect `outputs` |
| `1` | Usage or I/O | Check subcommand syntax, input path, permissions, and output directory |
| `2` | Parse/validation | Use the reported line/message to patch the `.inp`; do not retry unchanged |
| `3` | Solve | Check constraints, loads, material assignment, and convergence diagnostics |

### Bash recipe

```bash
#!/usr/bin/env bash
set -uo pipefail

deck="$1"
stdout_file="$(mktemp)"
stderr_file="$(mktemp)"

./target/release/solvix run "$deck" --json \
  >"$stdout_file" 2>"$stderr_file"
status=$?

case "$status" in
  0) python3 -m json.tool "$stdout_file" ;;
  1) echo "Fix invocation or file access; see $stderr_file" >&2 ;;
  2) echo "Fix the input deck; see $stderr_file" >&2 ;;
  3) echo "Inspect model and solver diagnostics; see $stderr_file" >&2 ;;
  *) echo "Unexpected Solvix exit code: $status" >&2 ;;
esac

rm -f "$stdout_file" "$stderr_file"
exit "$status"
```

The exit status can be made actionable in one line when a shell pipeline is not
required:

```bash
solvix run model.inp --json >result.json 2>run.log
status=$?
```

### JSON error contract

The exact fields are command-specific, but failures retain the `ok`/error
envelope so agents need one parsing strategy:

```json
{
  "ok": false,
  "error": {
    "kind": "parse",
    "line": 42,
    "message": "CPS4 element 7 references unknown node 108"
  }
}
```

Usage and solve failures use the same shape with `"kind": "usage"` or
`"kind": "solve"` respectively. Error diagnostics can also be mirrored to
stderr; stdout remains valid JSON when `--json` was requested.

## Info before run

Use `solvix info` as a cheap preflight. It parses and reports the model without
performing the linear solve.

```bash
solvix info examples/cantilever.inp --json >info.json 2>info.log
```

An illustrative success payload looks like:

```json
{
  "ok": true,
  "job": "cantilever",
  "nodes": 861,
  "elements": 800,
  "element_types": {
    "CPS4": 800
  },
  "materials": ["STEEL"],
  "boundary_conditions": 1,
  "concentrated_loads": 2,
  "warnings": []
}
```

The values above show payload shape; they are not a benchmark claim.

### Preflight decision table

| Preflight result | Next action |
|------------------|-------------|
| Exit `0`, no warnings | Run the solve |
| Exit `0` with skipped-keyword warnings | Decide whether the skipped keyword affects physics before running |
| Exit `2` | Repair the deck; a solve cannot safely proceed |
| Missing node/element/material counts | Treat as a failed contract and inspect stderr |
| Empty loads or constraints | Ask the user or inspect the deck rather than launching a doomed solve |

### Full agent workflow

```bash
deck="examples/cantilever.inp"

solvix info "$deck" --json >model-info.json 2>model-info.log
info_status=$?
if [ "$info_status" -ne 0 ]; then
  echo "Model preflight failed with exit $info_status" >&2
  exit "$info_status"
fi

mkdir -p artifacts
solvix run "$deck" --json -o artifacts >run-result.json 2>run-result.log
run_status=$?
if [ "$run_status" -ne 0 ]; then
  echo "Solve failed with exit $run_status" >&2
  exit "$run_status"
fi
```

Only after the final command exits `0` should an agent inspect files listed in
`run-result.json`.

## Run result payload

```bash
solvix run examples/cantilever.inp --json -o artifacts
```

An illustrative success payload is:

```json
{
  "ok": true,
  "job": "cantilever",
  "nodes": 861,
  "elements": 800,
  "iterations": 234,
  "residual": 8.7e-13,
  "max_disp": 0.00473,
  "max_von_mises": 38200000.0,
  "outputs": ["artifacts/cantilever.vtu"]
}
```

Field names may gain additional details, but agents should rely on:

- `ok` for command status;
- `outputs` for artifact paths;
- `iterations` and `residual` for solver progress/convergence;
- model-count fields to correlate the result with the preflight model.

Do not infer benchmark performance from one run. Systematic QA belongs to the
benchmark suite and must be reported only through its defined methodology.

## VTU post-processing

`outputs` names VTU files containing an inline-XML VTK UnstructuredGrid.
Point data includes displacement; cell data includes element von Mises stress.

### ParaView

1. Open the path from `outputs`.
2. Apply the reader.
3. Select the displacement array for glyph or deformation visualization.
4. Select von Mises stress for cell-centered coloring.
5. Use warp/glyph scale factors appropriate for the model's displacement and
   coordinate units.

### Dependency-free XML inspection

This Python snippet uses only the standard library to list data arrays and
their value ranges:

```python
#!/usr/bin/env python3
import sys
import xml.etree.ElementTree as ET


def as_floats(text):
    return [float(value) for value in text.split()]


def inspect(path):
    root = ET.parse(path).getroot()
    piece = root.find(".//UnstructuredGrid/Piece")
    print(f"file={path}")
    print(f"points={piece.get('NumberOfPoints')}")
    print(f"cells={piece.get('NumberOfCells')}")

    for section in ("PointData", "CellData"):
        parent = piece.find(section)
        if parent is None:
            continue
        for array in parent.findall("DataArray"):
            values = as_floats(array.text or "")
            print(
                f"{section[:-4]} name={array.get('Name')!r} "
                f"components={array.get('NumberOfComponents')} "
                f"min={min(values):.12g} max={max(values):.12g}"
            )


if __name__ == "__main__":
    inspect(sys.argv[1])
```

Run it after the solver exits:

```bash
python3 inspect_vtu.py artifacts/cantilever.vtu
```

### Connecting JSON and VTU

Agents should use the JSON payload as the manifest:

1. Confirm exit code `0` and `ok: true`.
2. Read each path in `outputs`.
3. Report scalar summaries from JSON when answering status questions.
4. Use VTU when the user asks for fields, contours, images, or detailed
   post-processing.
5. Preserve unit context from the `.inp` deck; VTU arrays do not rewrite the
   model into a different unit system.

## Practical troubleshooting

| Symptom | Likely class | First check |
|---------|--------------|-------------|
| Unknown flag or missing path | Exit `1` | Command usage and path spelling |
| Cannot open `.inp` | Exit `1` | Current working directory and permissions |
| Error cites an input line | Exit `2` | That keyword and its comma-separated fields |
| Unknown node/element/set | Exit `2` | Set expansion, `GENERATE` ranges, and ID spelling |
| Singularity or constraint failure | Exit `3` | Whether rigid-body motion is constrained |
| No convergence | Exit `3` | Mesh quality, material positivity, and load/boundary data |
| JSON parse fails on stdout | Contract failure | Check whether another program wrote to stdout or `--json` was omitted |

## Driving Solvix via MCP (v0.2)

`solvix mcp` speaks Model Context Protocol over stdio (newline-delimited
JSON-RPC 2.0). Register it with any MCP client:

```json
{
  "mcpServers": {
    "solvix": { "command": "/path/to/solvix", "args": ["mcp"] }
  }
}
```

Tools exposed: `solvix_info` (deck stats), `solvix_run` (solve + summary),
`solvix_results` (displacements mm + von Mises MPa arrays), `solvix_bench`
(verification suite — exit-equivalent green means the binary self-checks).

Example raw session:

```console
$ echo '{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"solvix_run","arguments":{"path":"examples/cantilever3d.inp"}}}' | solvix mcp
{"jsonrpc":"2.0","id":1,"result":{"content":[{"type":"text","text":"{\"ok\":true,...\"max_disp_mm\":13.719...}"}],"isError":false}}
```

Agent pattern: `solvix_info` first (cheap), then `solvix_run`, then
`solvix_results` only if you need arrays. `solvix_bench` after install or
upgrade to verify the binary on your machine.
