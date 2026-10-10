//! Abaqus ASCII `.frd` results reader (issue #37).
//!
//! Solvix reads the post-processing subset it can represent without an
//! element library: `NODE` coordinate blocks, `DISP` nodal displacement
//! blocks and `STRESS` nodal stress blocks (six components in Voigt order
//! xx, yy, zz, xy, xz, yz). Element records, other result quantities and
//! step metadata are skipped; each unsupported *block* is reported once.
//!
//! Record layout: Abaqus writes fixed-width FORTRAN records, while
//! hand-written and converted files commonly use commas. Both are
//! accepted. Fixed-width records follow `A3,I10,E12.5…` (the three
//! character marker may be omitted); comma records are `id,value,…`.
//! Whitespace-separated records also work because they are unambiguous.
//!
//! Markers: `-1`, `-2`, `-3` and `-4` are record markers when a line
//! contains nothing else (otherwise `-1`/`-2` prefix a data record and
//! `-2` continues a STRESS record). `1C`, `1P` and `100N` structural
//! records are tolerated (`100N` opens the node block in Abaqus-written
//! files), `9999` ends a block, and `-5` ends the file.

use std::str::FromStr;

use super::{ParseError, ParseReport};

/// Results extracted from an Abaqus ASCII `.frd` file.
///
/// Ids are model-independent: callers map them onto their own mesh (an
/// `.frd` file can contain foreign or renumbered nodes).
#[derive(Debug, Clone, Default)]
pub struct FrdResult {
    pub nodes: Vec<(u32, [f64; 3])>,
    pub displacements: Vec<(u32, [f64; 3])>,
    pub stresses6: Vec<(u32, [f64; 6])>,
    /// Non-fatal notes collected while parsing (unsupported records and
    /// stray lines that were skipped).
    pub report: ParseReport,
}

#[derive(Clone, Copy)]
enum Mode {
    Nodes,
    Displacements,
    Stresses,
    /// A recognized but unsupported block: its records are skipped.
    Skipping,
}

#[derive(Clone, Copy)]
struct Line<'a> {
    number: usize,
    text: &'a str,
}

impl Line<'_> {
    /// A three-character FORTRAN marker such as `-1 ` or `-2 `.
    fn marker(&self) -> &str {
        let end = self
            .text
            .char_indices()
            .nth(2)
            .map(|(index, character)| index + character.len_utf8())
            .unwrap_or(self.text.len());
        &self.text[..end]
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum RecordMarker {
    Start,
    Continuation,
}

/// The `A3` FORTRAN marker at the start of a data record. Comma-written
/// files sometimes glue the marker to its first comma (`-2,`), which is
/// accepted too.
fn record_marker(line: &Line<'_>) -> Option<RecordMarker> {
    let marker = line.marker().trim_end();
    let marker = marker.strip_suffix(',').unwrap_or(marker);
    match marker {
        "-1" => Some(RecordMarker::Start),
        "-2" => Some(RecordMarker::Continuation),
        _ => None,
    }
}

struct PendingStress {
    line: usize,
    id: u32,
    values: Vec<f64>,
}

struct Parsed {
    id: Option<u32>,
    values: Vec<f64>,
}

#[derive(Clone, Copy)]
enum ErrorClass {
    Empty,
    Id,
    Value,
}

struct FieldFailure {
    class: ErrorClass,
    error: ParseError,
}

impl FieldFailure {
    fn empty(line: usize, message: String) -> Self {
        Self::new(ErrorClass::Empty, line, message)
    }

    fn id(line: usize, message: String) -> Self {
        Self::new(ErrorClass::Id, line, message)
    }

    fn value(line: usize, message: String) -> Self {
        Self::new(ErrorClass::Value, line, message)
    }

    fn new(class: ErrorClass, line: usize, message: String) -> Self {
        Self {
            class,
            error: ParseError { line, message },
        }
    }
}

#[derive(Clone, Copy)]
enum Field {
    Node,
    Displacement,
    Stress,
}

impl Field {
    fn record(self) -> &'static str {
        match self {
            Field::Node => "NODE",
            Field::Displacement => "DISP",
            Field::Stress => "STRESS",
        }
    }

    fn value(self) -> &'static str {
        match self {
            Field::Node => "node coordinate",
            Field::Displacement => "displacement component",
            Field::Stress => "stress component",
        }
    }
}

enum Control {
    StepHeader,
    Block(Mode),
    Unsupported(String),
}

/// Parse the Solvix-supported subset of an Abaqus ASCII `.frd` results
/// file: node coordinates, nodal displacements and six-component nodal
/// stresses. Ids are returned unchanged so callers can map them onto
/// their own model; unsupported blocks are skipped and reported.
pub fn parse_frd(text: &str) -> Result<FrdResult, ParseError> {
    let mut result = FrdResult::default();
    let mut mode = None;
    let mut pending_stress = None;

    for line in lines(text) {
        let trimmed = line.text.trim();
        if trimmed.is_empty() {
            continue;
        }

        // `-5` terminates the file; anything after it is ignored.
        if trimmed == "-5" {
            flush_pending(&mut pending_stress, &mut result.stresses6)?;
            break;
        }

        if let Some(control) = control_record(trimmed) {
            flush_pending(&mut pending_stress, &mut result.stresses6)?;
            mode = match control {
                Control::StepHeader => None,
                Control::Block(block) => Some(block),
                Control::Unsupported(key) => {
                    result.report.warnings.push(format!(
                        "line {}: unsupported '{}' block skipped",
                        line.number, key
                    ));
                    Some(Mode::Skipping)
                }
            };
            continue;
        }

        // `9999` closes a block; the `-1`/`-2`/`-3`/`-4` record markers
        // close a record but leave the block active for its siblings.
        if trimmed == "9999" {
            flush_pending(&mut pending_stress, &mut result.stresses6)?;
            mode = None;
            continue;
        }
        if matches!(trimmed, "-1" | "-2" | "-3" | "-4") {
            flush_pending(&mut pending_stress, &mut result.stresses6)?;
            continue;
        }

        let Some(active) = mode else {
            result.report.warnings.push(format!(
                "line {}: skipped record outside NODE/DISP/STRESS block",
                line.number
            ));
            continue;
        };

        match active {
            Mode::Nodes => {
                let parsed = parse_fields(&line, true, Field::Node)?;
                if parsed.values.len() != 3 {
                    return err(
                        line.number,
                        format!(
                            "NODE record expects id, x, y, z — got {} coordinate(s)",
                            parsed.values.len()
                        ),
                    );
                }
                let id = parsed.id.expect("id requested");
                let values: [f64; 3] = parsed.values.try_into().expect("three coordinates");
                result.nodes.push((id, values));
            }
            Mode::Displacements => {
                let parsed = parse_fields(&line, true, Field::Displacement)?;
                if parsed.values.len() != 3 {
                    return err(
                        line.number,
                        format!(
                            "DISP record expects id, ux, uy, uz — got {} component(s)",
                            parsed.values.len()
                        ),
                    );
                }
                let id = parsed.id.expect("id requested");
                let values: [f64; 3] = parsed.values.try_into().expect("three components");
                result.displacements.push((id, values));
            }
            Mode::Stresses => {
                if pending_stress.is_some() {
                    let parsed = parse_fields(&line, false, Field::Stress)?;
                    append_stress_continuation(&mut pending_stress, parsed, line.number)?;
                    flush_pending(&mut pending_stress, &mut result.stresses6)?;
                } else if record_marker(&line) == Some(RecordMarker::Continuation) {
                    return err(
                        line.number,
                        "STRESS continuation has no matching record".to_string(),
                    );
                } else {
                    let parsed = parse_fields(&line, true, Field::Stress)?;
                    let id = parsed.id.expect("id requested");
                    match parsed.values.len() {
                        0 => {
                            return err(line.number, "STRESS record has no components".to_string())
                        }
                        1..=5 => {
                            pending_stress = Some(PendingStress {
                                line: line.number,
                                id,
                                values: parsed.values,
                            })
                        }
                        6 => result
                            .stresses6
                            .push((id, parsed.values.try_into().expect("six components"))),
                        _ => {
                            return err(
                                line.number,
                                format!(
                                    "STRESS record expects at most six components — got {}",
                                    parsed.values.len()
                                ),
                            )
                        }
                    }
                }
            }
            Mode::Skipping => {}
        }
    }

    flush_pending(&mut pending_stress, &mut result.stresses6)?;
    Ok(result)
}

fn control_record(text: &str) -> Option<Control> {
    let trimmed = text.trim_start();
    // `100N` is Abaqus's own node-block header; its count and format
    // payload is ignored, but the `-1` records that follow are nodes.
    // Requiring the key letter keeps numeric data records (for example a
    // whitespace node id `100 1.0 …`) out of the control path.
    if let Some(rest) = trimmed
        .strip_prefix("100")
        .filter(|rest| rest.starts_with(|key: char| key.is_ascii_uppercase()))
    {
        return Some(if rest.starts_with('N') {
            Control::Block(Mode::Nodes)
        } else {
            Control::StepHeader
        });
    }
    for marker in ["1C", "1U"] {
        let Some(rest) = trimmed.strip_prefix(marker) else {
            continue;
        };
        // A bare structural `1C` record opens the file in real Abaqus
        // output; it carries no block to read (or skip).
        let Some(key) = rest.split_whitespace().next() else {
            return Some(Control::StepHeader);
        };
        let key = key.to_ascii_uppercase();
        return Some(match key.as_str() {
            "NODE" => Control::Block(Mode::Nodes),
            "DISP" => Control::Block(Mode::Displacements),
            "STRESS" => Control::Block(Mode::Stresses),
            other => Control::Unsupported(other.to_string()),
        });
    }
    trimmed.starts_with("1P").then_some(Control::StepHeader)
}

/// Parse one data record, accepting comma, fixed-column and whitespace
/// layouts. `with_id` distinguishes a primary record (`id, values…`) from
/// a continuation (`values…`).
fn parse_fields(line: &Line<'_>, with_id: bool, field: Field) -> Result<Parsed, ParseError> {
    let marker = record_marker(line);
    let body = if marker.is_some() {
        &line.text[line.marker().len()..]
    } else {
        line.text
    };

    let mut candidates = Vec::new();
    if body.contains(',') {
        candidates.push(comma_fields(body));
    } else {
        candidates.push(column_fields(body, with_id.then_some(10)));
        if marker.is_some() {
            // Tolerate a hand-written record that omits the A3 marker.
            candidates.push(column_fields(line.text, with_id.then_some(10)));
        }
        candidates.push(line.text.split_whitespace().map(str::to_string).collect());
    }

    // If no layout parses, report the most diagnostic error: a bad value
    // beats a bad id, which beats an empty fixed-column candidate.
    let mut value_error = None;
    let mut other_error = None;
    for candidate in &candidates {
        match try_fields(candidate, with_id, line.number, field) {
            Ok(parsed) => return Ok(parsed),
            Err(failure) => {
                match failure.class {
                    ErrorClass::Value => value_error.get_or_insert(failure.error),
                    ErrorClass::Id | ErrorClass::Empty => other_error.get_or_insert(failure.error),
                };
            }
        }
    }

    Err(value_error.or(other_error).unwrap_or_else(|| ParseError {
        line: line.number,
        message: format!("{} record is empty", field.record()),
    }))
}

fn comma_fields(text: &str) -> Vec<String> {
    let mut fields: Vec<String> = text.split(',').map(str::trim).map(str::to_string).collect();
    while fields.last().is_some_and(String::is_empty) {
        fields.pop();
    }
    fields
}

/// Split a FORTRAN record into `I10` followed by `E12.5` fields. Field
/// widths matter because Abaqus omits the separator between a negative
/// value and the value before it.
fn column_fields(text: &str, id_width: Option<usize>) -> Vec<String> {
    if !text.is_ascii() {
        // Numeric records are ASCII; let the whitespace candidate report
        // the offending field instead of slicing through a character.
        return Vec::new();
    }
    let mut fields = Vec::new();
    let mut start = 0usize;
    let mut width = id_width.unwrap_or(12);
    while start < text.len() {
        let end = (start + width).min(text.len());
        fields.push(text[start..end].trim().to_string());
        start = end;
        width = 12;
    }
    while fields.last().is_some_and(String::is_empty) {
        fields.pop();
    }
    fields
}

fn try_fields(
    fields: &[String],
    with_id: bool,
    line: usize,
    field: Field,
) -> Result<Parsed, FieldFailure> {
    let (id, values) = if with_id {
        let Some(first) = fields.first() else {
            return Err(FieldFailure::empty(
                line,
                format!("{} record is empty", field.record()),
            ));
        };
        let id = first
            .parse::<u32>()
            .map_err(|_| FieldFailure::id(line, format!("malformed node id '{first}'")))?;
        (Some(id), &fields[1..])
    } else {
        (None, fields)
    };

    let mut numbers = Vec::with_capacity(values.len());
    for value in values {
        let number = f64::from_str(value).map_err(|_| {
            FieldFailure::value(line, format!("malformed {} '{value}'", field.value()))
        })?;
        if !number.is_finite() {
            return Err(FieldFailure::value(
                line,
                format!("{} '{value}' is not finite", field.value()),
            ));
        }
        numbers.push(number);
    }
    Ok(Parsed {
        id,
        values: numbers,
    })
}

fn append_stress_continuation(
    pending: &mut Option<PendingStress>,
    parsed: Parsed,
    line: usize,
) -> Result<(), ParseError> {
    let Some(record) = pending.as_mut() else {
        return err(
            line,
            "STRESS continuation has no matching record".to_string(),
        );
    };
    if parsed.values.is_empty() {
        return err(line, "STRESS continuation has no components".to_string());
    }
    let total = record.values.len() + parsed.values.len();
    if total > 6 {
        return err(
            line,
            format!("STRESS continuation overflows six components ({total})"),
        );
    }
    record.values.extend_from_slice(&parsed.values);
    Ok(())
}

fn flush_pending(
    pending: &mut Option<PendingStress>,
    stresses: &mut Vec<(u32, [f64; 6])>,
) -> Result<(), ParseError> {
    let Some(record) = pending.take() else {
        return Ok(());
    };
    if record.values.len() != 6 {
        return err(
            record.line,
            format!(
                "incomplete STRESS record for node {} — got {} of six components",
                record.id,
                record.values.len()
            ),
        );
    }
    stresses.push((record.id, record.values.try_into().expect("six components")));
    Ok(())
}

fn lines(text: &str) -> Vec<Line<'_>> {
    text.lines()
        .enumerate()
        .map(|(index, raw)| Line {
            number: index + 1,
            text: raw.trim_end_matches(['\r', ' ', '\t']),
        })
        .collect()
}

fn err<T>(line: usize, message: String) -> Result<T, ParseError> {
    Err(ParseError { line, message })
}
