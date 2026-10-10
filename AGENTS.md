# Agent Guide

This guide applies to AI coding agents working in the Solvix repository.

## Operating rules

1. Start at the [issue board](https://github.com/typezo/Solvix/issues). Work on
   one accepted issue at a time and preserve the scope stated in that issue.
2. One issue owns one directory. Keep implementation changes inside the crate
   directory identified by the issue (`crates/solvix-core`,
   `crates/solvix-io`, or `crates/solvix-cli`).
3. Workspace crates have a zero external-dependency policy. Do not add runtime,
   development, or build dependencies from external registries.
4. Preserve Solvix's JSON-first UX: valid JSON for machine-readable results on
   stdout, human diagnostics on stderr, and deliberate changes to field names
   or exit semantics.
5. Do not run `git commit` or perform any other version-control mutation. Leave
   the completed changes in the working tree; the maintainer agent handles
   version control.

## Verification

Before reporting an implementation as complete, run:

```console
cargo test --workspace
```

Also follow any issue-specific test or benchmark instructions. For work that
changes Rust source, formatting and lint checks are expected to pass:

```console
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
```
