//! Gmsh ASCII `.msh` reader for MSH 2.2 (and the compatible MSH 2.0
//! legacy layout) and MSH 4.1.
//!
//! MSH files describe a mesh, not an analysis. Consequently the imported
//! [`Model`] has no boundary conditions or loads: callers add those
//! programmatically before invoking the solver. A bare imported model fails
//! `Model::validate` with `NoBoundaryConditions` (and `solve` reports that
//! condition as `SolveError::NoConstraints`).
//!
//! MSH also has no material blocks. To make the resulting model usable, a
//! default linear-elastic material is created and assigned to every element.
//! Physical groups on the same dimension as supported elements become element
//! sets. Lower-dimensional groups would be node sets; the current shared
//! `Model` contract has no persistent NSET field, so they are not mis-encoded
//! as element sets.

use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use solvix_core::{ElemKind, Element, ElementSet, Material, Model, Node};

use super::{ParseError, ParseReport};

#[derive(Clone, Copy)]
struct Line<'a> {
    number: usize,
    text: &'a str,
}

impl Line<'_> {
    fn fields(&self) -> Vec<&str> {
        self.text.split_whitespace().collect()
    }
}

#[derive(Debug, Clone)]
struct PhysicalName {
    dimension: usize,
    tag: u32,
    name: String,
}

#[derive(Debug, Clone)]
struct RawElement {
    id: u32,
    kind: ElemKind,
    connectivity: Vec<u32>,
    entity_dimension: usize,
    entity_tag: u32,
    physical_tags: Vec<u32>,
    line: usize,
}

#[derive(Clone, Copy, PartialEq)]
enum MshVersion {
    Legacy,
    V41,
}

pub(super) fn parse_msh(text: &str) -> Result<(Model, ParseReport), ParseError> {
    let lines = nonempty_lines(text);
    let mut report = ParseReport::default();
    let mut version = None;
    let mut nodes = None;
    let mut elements = None;
    let mut physical_names = Vec::new();
    let mut entity_groups = HashMap::new();

    let mut i = 0usize;
    while i < lines.len() {
        let opening = lines[i];
        if !opening.text.starts_with('$') {
            return err(
                opening.number,
                format!("expected MSH section header, got '{}'", opening.text),
            );
        }

        let section = opening.text[1..].to_ascii_lowercase();
        if section.starts_with("end") {
            return err(
                opening.number,
                format!("unexpected section end '{}'", opening.text),
            );
        }

        let end = format!("$end{section}");
        let Some(offset) = lines[i + 1..]
            .iter()
            .position(|line| line.text.eq_ignore_ascii_case(&end))
        else {
            return err(
                opening.number,
                format!("missing '{end}' for section '${section}'"),
            );
        };
        let end_index = i + 1 + offset;
        let contents = &lines[i + 1..end_index];

        match section.as_str() {
            "meshformat" => {
                if version.is_some() {
                    return err(opening.number, "duplicate $MeshFormat section".to_string());
                }
                version = Some(parse_mesh_format(contents, opening.number)?);
            }
            "physicalnames" => {
                physical_names = parse_physical_names(contents, opening.number)?;
            }
            "entities" => {
                entity_groups = parse_entities(contents, opening.number)?;
            }
            "nodes" => {
                if nodes.is_some() {
                    return err(opening.number, "duplicate $Nodes section".to_string());
                }
                let selected_version = version.ok_or_else(|| ParseError {
                    line: opening.number,
                    message: "$Nodes appears before $MeshFormat".to_string(),
                })?;
                nodes = Some(if selected_version == MshVersion::Legacy {
                    parse_legacy_nodes(contents, opening.number)?
                } else {
                    parse_v41_nodes(contents, opening.number)?
                });
            }
            "elements" => {
                if elements.is_some() {
                    return err(opening.number, "duplicate $Elements section".to_string());
                }
                let selected_version = version.ok_or_else(|| ParseError {
                    line: opening.number,
                    message: "$Elements appears before $MeshFormat".to_string(),
                })?;
                if nodes.is_none() {
                    return err(
                        opening.number,
                        "$Elements appears before $Nodes".to_string(),
                    );
                }
                elements = Some(if selected_version == MshVersion::Legacy {
                    parse_legacy_elements(contents, opening.number, &mut report)?
                } else {
                    parse_v41_elements(contents, opening.number, &mut report)?
                });
            }
            "comments" => {}
            other => report.warnings.push(format!(
                "line {}: unsupported section '${other}' skipped",
                opening.number
            )),
        }

        i = end_index + 1;
    }

    let last_line = lines.last().map_or(1, |line| line.number);
    if version.is_none() {
        return err(last_line, "missing $MeshFormat section".to_string());
    }
    let Some(mesh_nodes) = nodes else {
        return err(last_line, "missing $Nodes section".to_string());
    };
    let Some(mesh_elements) = elements else {
        return err(last_line, "missing $Elements section".to_string());
    };
    if mesh_nodes.is_empty() || mesh_elements.is_empty() {
        return err(
            last_line,
            "MSH must contain at least one node and one supported element".to_string(),
        );
    }

    let mut node_indices = HashMap::with_capacity(mesh_nodes.len());
    for (index, node) in mesh_nodes.iter().enumerate() {
        if node_indices.insert(node.id, index).is_some() {
            return err(last_line, format!("duplicate node tag {}", node.id));
        }
    }
    let mut element_ids = HashSet::with_capacity(mesh_elements.len());
    for element in &mesh_elements {
        if !element_ids.insert(element.id) {
            return err(
                element.line,
                format!("duplicate element tag {}", element.id),
            );
        }
        for node in &element.connectivity {
            if !node_indices.contains_key(node) {
                return err(
                    element.line,
                    format!("element {} references unknown node {}", element.id, node),
                );
            }
        }
    }

    let mut model = Model {
        thickness: 1.0,
        ..Model::default()
    };
    let elsets =
        physical_element_sets(&mesh_elements, &physical_names, &entity_groups, &mut report);
    model.mesh.nodes = mesh_nodes;
    model.mesh.elements = mesh_elements
        .into_iter()
        .map(|raw| Element {
            id: raw.id,
            kind: raw.kind,
            connectivity: raw.connectivity,
        })
        .collect();
    model.elsets = elsets;

    // MSH has no material definition. This is the parser-side equivalent of
    // MATERIAL name=DEFAULT / ELASTIC 210e9,0.3 plus SOLID SECTION, ELSET=ALL.
    model.materials.push(Material {
        name: "DEFAULT".to_string(),
        youngs_modulus: 210.0e9,
        poisson_ratio: 0.3,
    });
    let all_name = if model.elsets.iter().any(|set| set.name == "ALL") {
        "SOLVIX_ALL"
    } else {
        "ALL"
    };
    model.elsets.push(ElementSet {
        name: all_name.to_string(),
        element_ids: model.mesh.elements.iter().map(|e| e.id).collect(),
    });
    model
        .mat_assignments
        .push((all_name.to_string(), "DEFAULT".to_string()));
    report.warnings.push(format!(
        "no materials in MSH; created DEFAULT (E=210e9, nu=0.3) and assigned all elements to {all_name}"
    ));

    Ok((model, report))
}

fn parse_mesh_format(contents: &[Line<'_>], opening: usize) -> Result<MshVersion, ParseError> {
    let Some(header) = contents.first() else {
        return err(opening, "$MeshFormat is empty".to_string());
    };
    let fields = header.fields();
    if fields.len() != 3 {
        return err(
            header.number,
            format!(
                "$MeshFormat expects version file-type data-size; got {} field(s)",
                fields.len()
            ),
        );
    }

    let version = number::<f64>(fields[0], header.number, "mesh format version")?;
    let file_type = number::<u8>(fields[1], header.number, "mesh file type")?;
    let data_size = number::<u8>(fields[2], header.number, "mesh data size")?;
    if file_type != 0 {
        return err(
            header.number,
            "binary MSH files are not supported".to_string(),
        );
    }
    if data_size != 8 {
        return err(header.number, "expected MSH data size 8".to_string());
    }

    if version == 2.0 || version == 2.2 {
        Ok(MshVersion::Legacy)
    } else if version == 4.1 {
        Ok(MshVersion::V41)
    } else {
        err(
            header.number,
            format!("unsupported MSH version {version}; expected 2.2 or 4.1"),
        )
    }
}

fn parse_legacy_nodes(contents: &[Line<'_>], opening: usize) -> Result<Vec<Node>, ParseError> {
    let Some(header) = contents.first() else {
        return err(opening, "$Nodes is empty".to_string());
    };
    let header_fields = header.fields();
    if header_fields.len() != 1 {
        return err(
            header.number,
            format!(
                "legacy $Nodes expects node count; got {} field(s)",
                header_fields.len()
            ),
        );
    }
    let count = number::<usize>(header_fields[0], header.number, "legacy node count")?;
    if contents.len() != count + 1 {
        return err(
            header.number,
            format!(
                "legacy $Nodes declares {count} nodes but contains {}",
                contents.len() - 1
            ),
        );
    }

    let mut nodes = Vec::with_capacity(count);
    for line in &contents[1..] {
        let fields = line.fields();
        if fields.len() != 4 {
            return err(
                line.number,
                format!(
                    "legacy node expects tag x y z; got {} field(s)",
                    fields.len()
                ),
            );
        }
        let tag = number::<u32>(fields[0], line.number, "node tag")?;
        let x = number::<f64>(fields[1], line.number, "node coordinate")?;
        let y = number::<f64>(fields[2], line.number, "node coordinate")?;
        let z = number::<f64>(fields[3], line.number, "node coordinate")?;
        push_node(&mut nodes, tag, [x, y, z], line.number)?;
    }
    Ok(nodes)
}

fn parse_v41_nodes(contents: &[Line<'_>], opening: usize) -> Result<Vec<Node>, ParseError> {
    let Some(header) = contents.first() else {
        return err(opening, "$Nodes is empty".to_string());
    };
    let fields = header.fields();
    if fields.len() != 4 {
        return err(
            header.number,
            format!(
                "MSH 4.1 $Nodes expects 4 header fields; got {}",
                fields.len()
            ),
        );
    }
    let block_count = number::<usize>(fields[0], header.number, "node block count")?;
    let total = number::<usize>(fields[1], header.number, "node count")?;
    let min_tag = number::<u32>(fields[2], header.number, "minimum node tag")?;
    let max_tag = number::<u32>(fields[3], header.number, "maximum node tag")?;
    if total == 0 || min_tag == 0 || max_tag < min_tag {
        return err(header.number, "invalid MSH 4.1 node header".to_string());
    }

    let mut nodes = Vec::with_capacity(total);
    let mut cursor = 1usize;
    for _ in 0..block_count {
        let Some(block_header) = contents.get(cursor) else {
            return err(header.number, "truncated MSH 4.1 node block".to_string());
        };
        cursor += 1;
        let block_fields = block_header.fields();
        if block_fields.len() != 4 {
            return err(
                block_header.number,
                format!(
                    "MSH 4.1 node block header expects dimension tag parametric count; got {} fields",
                    block_fields.len()
                ),
            );
        }
        let dimension = number::<usize>(block_fields[0], block_header.number, "entity dimension")?;
        let _entity_tag = number::<u32>(block_fields[1], block_header.number, "entity tag")?;
        let parametric = number::<u8>(block_fields[2], block_header.number, "parametric flag")?;
        let count = number::<usize>(block_fields[3], block_header.number, "nodes in block")?;
        if parametric != 0 || dimension > 3 || count == 0 {
            return err(
                block_header.number,
                "invalid MSH 4.1 node block header".to_string(),
            );
        }

        let mut tags = Vec::with_capacity(count);
        for _ in 0..count {
            let Some(tag_line) = contents.get(cursor) else {
                return err(block_header.number, "truncated node tag list".to_string());
            };
            let tag_fields = tag_line.fields();
            if tag_fields.len() != 1 {
                return err(
                    tag_line.number,
                    "expected one node tag per line".to_string(),
                );
            }
            tags.push(number::<u32>(tag_fields[0], tag_line.number, "node tag")?);
            cursor += 1;
        }
        for tag in tags {
            let Some(coordinate_line) = contents.get(cursor) else {
                return err(
                    block_header.number,
                    "truncated node coordinates".to_string(),
                );
            };
            let coordinate_fields = coordinate_line.fields();
            if coordinate_fields.len() != 3 {
                return err(
                    coordinate_line.number,
                    format!(
                        "node coordinates expect x y z; got {} fields",
                        coordinate_fields.len()
                    ),
                );
            }
            let x = number::<f64>(
                coordinate_fields[0],
                coordinate_line.number,
                "node coordinate",
            )?;
            let y = number::<f64>(
                coordinate_fields[1],
                coordinate_line.number,
                "node coordinate",
            )?;
            let z = number::<f64>(
                coordinate_fields[2],
                coordinate_line.number,
                "node coordinate",
            )?;
            push_node(&mut nodes, tag, [x, y, z], coordinate_line.number)?;
            cursor += 1;
        }
    }

    if cursor != contents.len() {
        return err(
            contents[cursor].number,
            "unexpected trailing data in MSH 4.1 $Nodes".to_string(),
        );
    }
    if nodes.len() != total {
        return err(
            header.number,
            format!(
                "MSH 4.1 $Nodes declares {total} nodes but parsed {}",
                nodes.len()
            ),
        );
    }
    Ok(nodes)
}

fn parse_legacy_elements(
    contents: &[Line<'_>],
    opening: usize,
    report: &mut ParseReport,
) -> Result<Vec<RawElement>, ParseError> {
    let Some(header) = contents.first() else {
        return err(opening, "$Elements is empty".to_string());
    };
    let header_fields = header.fields();
    if header_fields.len() != 1 {
        return err(
            header.number,
            format!(
                "legacy $Elements expects element count; got {} fields",
                header_fields.len()
            ),
        );
    }
    let count = number::<usize>(header_fields[0], header.number, "legacy element count")?;
    if contents.len() != count + 1 {
        return err(
            header.number,
            format!(
                "legacy $Elements declares {count} elements but contains {}",
                contents.len() - 1
            ),
        );
    }

    let mut elements = Vec::with_capacity(count);
    for line in &contents[1..] {
        let fields = line.fields();
        if fields.len() < 3 {
            return err(
                line.number,
                format!(
                    "legacy element record expects at least 3 fields; got {}",
                    fields.len()
                ),
            );
        }
        let id = number::<u32>(fields[0], line.number, "element tag")?;
        let element_type = number::<u32>(fields[1], line.number, "element type")?;
        let tag_count = number::<usize>(fields[2], line.number, "element tag count")?;
        let tag_end = 3 + tag_count;
        if tag_count > fields.len().saturating_sub(3) || fields.len() < tag_end {
            return err(
                line.number,
                format!("element {id} has an invalid tag count"),
            );
        }
        let entity_tag = if tag_count == 0 {
            0
        } else {
            number::<u32>(fields[3], line.number, "entity tag")?
        };
        let physical_tags = if tag_count > 1 {
            fields[4..tag_end]
                .iter()
                .map(|field| number::<u32>(field, line.number, "physical tag"))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            Vec::new()
        };

        let Some(kind) = supported_kind(element_type) else {
            report.warnings.push(format!(
                "line {}: unsupported MSH element type {element_type} skipped (supported: triangle3, quad4, tetra4, hexa8)",
                line.number
            ));
            continue;
        };
        let node_count = kind.node_count();
        if fields.len() != tag_end + node_count {
            return err(
                line.number,
                format!(
                    "MSH element type {element_type} expects {node_count} nodes; record has {}",
                    fields.len().saturating_sub(tag_end)
                ),
            );
        }
        let connectivity = fields[tag_end..]
            .iter()
            .map(|field| number::<u32>(field, line.number, "element node tag"))
            .collect::<Result<Vec<_>, _>>()?;
        elements.push(RawElement {
            id,
            kind,
            connectivity,
            entity_dimension: element_dimension(kind),
            entity_tag,
            physical_tags,
            line: line.number,
        });
    }
    Ok(elements)
}

fn parse_v41_elements(
    contents: &[Line<'_>],
    opening: usize,
    report: &mut ParseReport,
) -> Result<Vec<RawElement>, ParseError> {
    let Some(header) = contents.first() else {
        return err(opening, "$Elements is empty".to_string());
    };
    let fields = header.fields();
    if fields.len() != 4 {
        return err(
            header.number,
            format!(
                "MSH 4.1 $Elements expects 4 header fields; got {}",
                fields.len()
            ),
        );
    }
    let block_count = number::<usize>(fields[0], header.number, "element block count")?;
    let total = number::<usize>(fields[1], header.number, "element count")?;
    let min_tag = number::<u32>(fields[2], header.number, "minimum element tag")?;
    let max_tag = number::<u32>(fields[3], header.number, "maximum element tag")?;
    if total == 0 || min_tag == 0 || max_tag < min_tag {
        return err(header.number, "invalid MSH 4.1 element header".to_string());
    }

    let mut elements = Vec::with_capacity(total);
    let mut cursor = 1usize;
    let mut skipped = 0usize;
    for _ in 0..block_count {
        let Some(block_header) = contents.get(cursor) else {
            return err(header.number, "truncated MSH 4.1 element block".to_string());
        };
        cursor += 1;
        let block_fields = block_header.fields();
        if block_fields.len() != 4 {
            return err(
                block_header.number,
                format!(
                    "MSH 4.1 element block header expects dimension tag type count; got {} fields",
                    block_fields.len()
                ),
            );
        }
        let dimension = number::<usize>(block_fields[0], block_header.number, "entity dimension")?;
        let entity_tag = number::<u32>(block_fields[1], block_header.number, "entity tag")?;
        let element_type = number::<u32>(block_fields[2], block_header.number, "element type")?;
        let count = number::<usize>(block_fields[3], block_header.number, "elements in block")?;
        if dimension > 3 || count == 0 {
            return err(
                block_header.number,
                "invalid MSH 4.1 element block header".to_string(),
            );
        }

        let Some(kind) = supported_kind(element_type) else {
            report.warnings.push(format!(
                "line {}: unsupported MSH element type {element_type} block skipped (supported: triangle3, quad4, tetra4, hexa8)",
                block_header.number
            ));
            if cursor + count > contents.len() {
                return err(
                    block_header.number,
                    "truncated unsupported element block".to_string(),
                );
            }
            cursor += count;
            skipped += count;
            continue;
        };

        let node_count = kind.node_count();
        for _ in 0..count {
            let Some(element_line) = contents.get(cursor) else {
                return err(
                    block_header.number,
                    "truncated MSH 4.1 element records".to_string(),
                );
            };
            let element_fields = element_line.fields();
            if element_fields.len() != 1 + node_count {
                return err(
                    element_line.number,
                    format!(
                        "MSH element type {element_type} expects tag plus {node_count} nodes; got {} fields",
                        element_fields.len()
                    ),
                );
            }
            let id = number::<u32>(element_fields[0], element_line.number, "element tag")?;
            let connectivity = element_fields[1..]
                .iter()
                .map(|field| number::<u32>(field, element_line.number, "element node tag"))
                .collect::<Result<Vec<_>, _>>()?;
            elements.push(RawElement {
                id,
                kind,
                connectivity,
                entity_dimension: dimension,
                entity_tag,
                physical_tags: Vec::new(),
                line: element_line.number,
            });
            cursor += 1;
        }
    }

    if cursor != contents.len() {
        return err(
            contents[cursor].number,
            "unexpected trailing data in MSH 4.1 $Elements".to_string(),
        );
    }
    if elements.len() + skipped != total {
        return err(
            header.number,
            format!(
                "MSH 4.1 $Elements declares {total} elements but parsed {} supported and {skipped} unsupported",
                elements.len()
            ),
        );
    }
    Ok(elements)
}

fn parse_physical_names(
    contents: &[Line<'_>],
    opening: usize,
) -> Result<Vec<PhysicalName>, ParseError> {
    let Some(header) = contents.first() else {
        return err(opening, "$PhysicalNames is empty".to_string());
    };
    let header_fields = header.fields();
    if header_fields.len() != 1 {
        return err(
            header.number,
            format!(
                "$PhysicalNames expects count; got {} fields",
                header_fields.len()
            ),
        );
    }
    let count = number::<usize>(header_fields[0], header.number, "physical name count")?;
    if contents.len() != count + 1 {
        return err(
            header.number,
            format!(
                "$PhysicalNames declares {count} names but contains {}",
                contents.len() - 1
            ),
        );
    }

    let mut names = Vec::with_capacity(count);
    let mut seen = HashSet::new();
    let mut seen_names = HashSet::new();
    for line in &contents[1..] {
        let Some((dimension, remainder)) = line.text.split_once(char::is_whitespace) else {
            return err(
                line.number,
                "physical name expects dimension tag \"name\"".to_string(),
            );
        };
        let Some((tag, raw_name)) = remainder.split_once(char::is_whitespace) else {
            return err(
                line.number,
                "physical name expects dimension tag \"name\"".to_string(),
            );
        };
        let dimension = number::<usize>(dimension.trim(), line.number, "physical group dimension")?;
        let tag = number::<u32>(tag.trim(), line.number, "physical group tag")?;
        let raw_name = raw_name.trim();
        let Some(name) = raw_name
            .strip_prefix('"')
            .and_then(|value| value.strip_suffix('"'))
        else {
            return err(
                line.number,
                "physical group name must be double-quoted".to_string(),
            );
        };
        if dimension > 3 || name.trim().is_empty() {
            return err(line.number, "invalid physical group definition".to_string());
        }
        let name = name.trim().to_ascii_uppercase();
        if !seen.insert((dimension, tag)) || !seen_names.insert(name.clone()) {
            return err(
                line.number,
                "duplicate physical group tag or name".to_string(),
            );
        }
        names.push(PhysicalName {
            dimension,
            tag,
            name,
        });
    }
    Ok(names)
}

fn parse_entities(
    contents: &[Line<'_>],
    opening: usize,
) -> Result<HashMap<(usize, u32), Vec<u32>>, ParseError> {
    let Some(header) = contents.first() else {
        return err(opening, "$Entities is empty".to_string());
    };
    let header_fields = header.fields();
    if header_fields.len() != 4 {
        return err(
            header.number,
            format!(
                "$Entities expects point curve surface volume counts; got {} fields",
                header_fields.len()
            ),
        );
    }
    let point_count = number::<usize>(header_fields[0], header.number, "point entity count")?;
    let curve_count = number::<usize>(header_fields[1], header.number, "curve entity count")?;
    let surface_count = number::<usize>(header_fields[2], header.number, "surface entity count")?;
    let volume_count = number::<usize>(header_fields[3], header.number, "volume entity count")?;
    let expected = 1 + point_count + curve_count + surface_count + volume_count;
    if contents.len() != expected {
        return err(
            header.number,
            format!(
                "$Entities declares {expected} records but contains {}",
                contents.len() - 1
            ),
        );
    }

    let mut groups = HashMap::new();
    let mut cursor = 1usize;
    let counts = [point_count, curve_count, surface_count, volume_count];
    for (dimension, &count) in counts.iter().enumerate() {
        for _ in 0..count {
            let line = contents[cursor];
            cursor += 1;
            let fields = line.fields();
            let Some(tag_text) = fields.first() else {
                return err(line.number, "empty entity record".to_string());
            };
            let tag = number::<u32>(tag_text, line.number, "entity tag")?;
            let coordinate_count = if dimension == 0 { 3 } else { 6 };
            let physical_offset = 1 + coordinate_count;
            if fields.len() <= physical_offset {
                return err(
                    line.number,
                    format!("truncated {dimension}-D entity record"),
                );
            }
            for field in &fields[1..physical_offset] {
                number::<f64>(field, line.number, "entity bounding coordinate")?;
            }
            let physical_count = number::<usize>(
                fields[physical_offset],
                line.number,
                "entity physical tag count",
            )?;
            let physical_end = physical_offset + 1 + physical_count;
            if fields.len() <= physical_end {
                return err(
                    line.number,
                    "entity has an invalid physical tag count".to_string(),
                );
            }
            let physical_tags = fields[physical_offset + 1..physical_end]
                .iter()
                .map(|field| number::<u32>(field, line.number, "entity physical tag"))
                .collect::<Result<Vec<_>, _>>()?;

            if dimension == 0 {
                if fields.len() != physical_end {
                    return err(
                        line.number,
                        "point entity record has trailing fields".to_string(),
                    );
                }
            } else {
                let boundary_count =
                    number::<usize>(fields[physical_end], line.number, "entity boundary count")?;
                let record_end = physical_end + 1 + boundary_count;
                if fields.len() != record_end {
                    return err(
                        line.number,
                        "entity boundary tag count does not match record".to_string(),
                    );
                }
                for field in &fields[physical_end + 1..record_end] {
                    number::<i64>(field, line.number, "signed boundary tag")?;
                }
            }

            if !physical_tags.is_empty() {
                groups.insert((dimension, tag), physical_tags);
            }
        }
    }
    Ok(groups)
}

fn physical_element_sets(
    elements: &[RawElement],
    physical_names: &[PhysicalName],
    entity_groups: &HashMap<(usize, u32), Vec<u32>>,
    report: &mut ParseReport,
) -> Vec<ElementSet> {
    let names: HashMap<(usize, u32), usize> = physical_names
        .iter()
        .enumerate()
        .map(|(index, group)| ((group.dimension, group.tag), index))
        .collect();
    let mut memberships = vec![Vec::new(); physical_names.len()];
    let mut membership_sets = vec![HashSet::new(); physical_names.len()];

    for element in elements {
        let tags = if element.physical_tags.is_empty() {
            entity_groups
                .get(&(element.entity_dimension, element.entity_tag))
                .cloned()
                .unwrap_or_default()
        } else {
            element.physical_tags.clone()
        };

        for tag in tags {
            let Some(&index) = names.get(&(element.entity_dimension, tag)) else {
                report.warnings.push(format!(
                    "line {}: element {} references undefined physical tag {tag}",
                    element.line, element.id
                ));
                continue;
            };
            if membership_sets[index].insert(element.id) {
                memberships[index].push(element.id);
            }
        }
    }

    physical_names
        .iter()
        .enumerate()
        .filter(|(index, _)| !memberships[*index].is_empty())
        .map(|(index, group)| ElementSet {
            name: group.name.clone(),
            element_ids: memberships[index].clone(),
        })
        .collect()
}

fn push_node(
    nodes: &mut Vec<Node>,
    tag: u32,
    coords: [f64; 3],
    line: usize,
) -> Result<(), ParseError> {
    if !coords[0].is_finite() || !coords[1].is_finite() || !coords[2].is_finite() {
        return err(line, format!("node {tag} has a non-finite coordinate"));
    }
    nodes.push(Node { id: tag, coords });
    Ok(())
}

fn supported_kind(element_type: u32) -> Option<ElemKind> {
    match element_type {
        2 => Some(ElemKind::Cps3),
        3 => Some(ElemKind::Cps4),
        4 => Some(ElemKind::C3d4),
        5 => Some(ElemKind::C3d8),
        _ => None,
    }
}

fn element_dimension(kind: ElemKind) -> usize {
    match kind {
        ElemKind::Cps3 | ElemKind::Cps4 | ElemKind::Cps8 | ElemKind::Cpe4 | ElemKind::Cpe8 => 2,
        ElemKind::C3d4 | ElemKind::C3d8 => 3,
        ElemKind::B31 => 1,
    }
}

fn nonempty_lines(text: &str) -> Vec<Line<'_>> {
    text.lines()
        .enumerate()
        .filter_map(|(index, raw)| {
            let trimmed = raw.trim();
            (!trimmed.is_empty()).then_some(Line {
                number: index + 1,
                text: trimmed,
            })
        })
        .collect()
}

fn number<T: FromStr>(field: &str, line: usize, what: &str) -> Result<T, ParseError> {
    field.parse::<T>().map_err(|_| ParseError {
        line,
        message: format!("malformed {what} '{field}'"),
    })
}

fn err<T>(line: usize, message: String) -> Result<T, ParseError> {
    Err(ParseError { line, message })
}
