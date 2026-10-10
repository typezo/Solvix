//! # solvix-tools
//!
//! Migration tools for Solvix users arriving from Abaqus (issue #35).
//!
//! `abaqus2solvix` (see `src/main.rs`) translates a bounded subset of
//! Abaqus Python scripts into Solvix `.inp` decks:
//!
//! * [`python_parse`] — hand-rolled tokenizer + recursive-descent
//!   parser for the scripting subset (attribute chains, keyword-arg
//!   calls, numbers/strings/lists/tuples, assignments; control flow
//!   blocks are skipped whole with a warning);
//! * [`translate`] — translation rules mapping the AST onto
//!   `solvix_core::Model` (structured mesh with poor-man's circular
//!   holes, materials, sections, steps/BCs/loads);
//! * [`meshgen`] — the structured CPS4 grid generator;
//! * [`deck`] — the `.inp` writer (inverse of `solvix_io`'s parser).

pub mod deck;
pub mod meshgen;
pub mod python_parse;
pub mod translate;
