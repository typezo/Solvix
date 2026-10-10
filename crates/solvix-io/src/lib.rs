//! # solvix-io
//!
//! Abaqus `.inp` reader and results writers.
//!
//! Implementation owned by issue #2 (parser) — this file fixes the interface.

mod msh;
mod parser;
mod vtu;

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
pub fn parse_inp(text: &str) -> Result<(Model, ParseReport), ParseError> {
    parser::parse_inp(text)
}

/// Parse a Gmsh ASCII `.msh` mesh (MSH 2.2 or 4.1).
///
/// Mesh import creates a default material and element sets from physical
/// groups. It intentionally creates no boundary conditions or loads; add
/// those to the returned model before solving.
pub fn parse_msh(text: &str) -> Result<(Model, ParseReport), ParseError> {
    msh::parse_msh(text)
}

/// Write nodal displacements and element von Mises stress as an
/// ASCII inline-XML VTK UnstructuredGrid `.vtu` file (ParaView-ready).
pub fn write_vtu(
    path: &std::path::Path,
    model: &Model,
    solution: &Solution,
) -> std::io::Result<()> {
    vtu::write_vtu(path, model, solution)
}
