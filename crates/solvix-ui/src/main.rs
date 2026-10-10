//! # solvix-ui — the Abaqus-CAE-style web app for Solvix
//!
//! `solvix serve` starts an HTTP server: static UI under `/` and `/static/`,
//! JSON REST API under `/api/*` (issue #13). The same API serves the browser
//! UI (issue #14) and AI agents (issue #15).
//!
//! API v1 contract:
//! - `GET  /api/model?job=NAME`      → mesh+model JSON (see issue #13)
//! - `POST /api/run`  body {"path": "examples/cantilever.inp"}
//!                                   → run summary JSON (CLI schema)
//! - `GET  /api/results/{job}`       → displacements + von Mises arrays
//! - `GET  /api/agent`               → machine-readable API description
//!
//! Implementation owned by issue #13.

fn main() {
    todo!("issue #13: solvix serve")
}
