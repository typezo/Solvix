//! # solvix-io
//!
//! Abaqus `.inp` reader and results writers.
//!
//! Implementation owned by issue #2 (parser) — this file fixes the interface.

use solvix_core::{Model, Solution};

/// Non-fatal notes collected while parsing (skipped keywords, defaults used).
#[derive(Debug, Clone, Default)]
pub struct ParseReport {
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ParseError {
    /// 1-based line number in the input deck.
    pub line: usize,
    pub message: String,
}

/// Parse an Abaqus input deck (v0.1 keyword subset, see issue #2).
///
/// Implementation owned by issue #2.
pub fn parse_inp(_text: &str) -> Result<(Model, ParseReport), ParseError> {
    todo!("issue #2: .inp parser")
}

/// Write nodal displacements and element von Mises stress as an
/// ASCII inline-XML VTK UnstructuredGrid `.vtu` file (ParaView-ready).
///
/// Implementation owned by issue #2.
pub fn write_vtu(
    _path: &std::path::Path,
    _model: &Model,
    _solution: &Solution,
) -> std::io::Result<()> {
    todo!("issue #2: VTU writer")
}
