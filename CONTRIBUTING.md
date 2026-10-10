# Contributing to Solvix

Thanks for helping build Solvix. This project favors small, reviewable, and
well-scoped changes that preserve a dependency-free Rust workspace.

## Start on the issue board

Work begins on the [issue board](https://github.com/typezo/Solvix/issues).
Before making changes, find or open an issue that describes the intended
outcome, alternatives considered, and acceptance criteria. Link the issue in
your pull request so reviewers have its context.

## One issue, one crate directory

An issue owns exactly one crate directory:

- numerical behavior and solver internals belong in `crates/solvix-core`;
- input parsing and result serialization belong in `crates/solvix-io`;
- command behavior and JSON output belong in `crates/solvix-cli`.

Keep code changes within the directory named by the issue. If a feature needs
coordinated changes in another crate, split the work into a separate issue and
pull request. Do not take simultaneous ownership of multiple crate directories
under one issue.

Tests and narrowly scoped workspace configuration required by an issue may live
where they are most useful, but review will reject unrelated drive-by changes.

## Zero external dependencies

Workspace crates use the Rust standard library only. Do not add dependencies,
dev-dependencies, or build-dependencies from crates.io or any other external
registry. New algorithms, parsing, serialization, and test tooling should be
implemented in Rust within the appropriate crate.

If an issue cannot be completed without an external dependency, raise that
constraint on the issue before implementation. Such a change requires an
explicit decision to revise this policy.

## Verify your work

All three checks must pass before requesting review:

```console
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

Add or update tests that exercise the behavior in the issue. For parser,
solver, and CLI changes, include a focused regression case wherever practical.

## JSON-first UX

Solvix is designed to be driven by humans and AI agents. Preserve JSON-first
command behavior:

- successful machine-readable results are emitted as valid JSON on stdout;
- diagnostics intended for humans are emitted on stderr;
- JSON schemas, fields, and exit semantics are changed deliberately and
  documented in the linked issue;
- output remains deterministic unless the operation inherently requires
  additional context.

## Commit style

Use [Conventional Commits](https://www.conventionalcommits.org/) with a scope
matching the owning crate:

```text
feat(core): assemble sparse systems with element matrices
fix(io): reject duplicate node labels
feat(cli): emit machine-readable diagnostics
docs: explain the dependency policy
test(core): cover triangular element stiffness
chore: update release workflow
```

Keep commit subjects in the imperative mood. Do not amend or rewrite commits
that another maintainer has already integrated.

AI coding agents must not run `git commit`, push commits, or otherwise perform
version-control operations. Leave changes in the working tree; the maintainer
agent owns version control.
