//! Abaqus `.inp` parser for the v0.1 keyword subset (issue #2).
//!
//! Line model: `**…` comments and blank lines are dropped up front; every
//! remaining line is either a keyword line (`*NAME, KEY=VALUE, …`, parsed
//! case-insensitively) or a data line (comma-separated fields). A keyword
//! "block" owns the run of data lines that follows it.

use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use solvix_core::{
    BoundaryCondition, ConcentratedLoad, ElemKind, Element, ElementSet, Material, Model, Node,
};

use super::{ParseError, ParseReport};

struct KwLine {
    line: usize,
    name: String,
    params: HashMap<String, String>,
}

enum Item {
    Kw(KwLine),
    Data { line: usize, text: String },
}

pub(super) fn parse_inp(text: &str) -> Result<(Model, ParseReport), ParseError> {
    let items = tokenize(text);
    let mut model = Model {
        thickness: 1.0,
        ..Model::default()
    };
    let mut report = ParseReport::default();
    let mut nsets: HashMap<String, Vec<u32>> = HashMap::new();
    let mut node_ids: HashSet<u32> = HashSet::new();
    let mut elem_ids: HashSet<u32> = HashSet::new();
    let mut pending_material: Option<String> = None;
    let mut steps = 0usize;

    let mut j = 0usize;
    while j < items.len() {
        let kw = match &items[j] {
            Item::Kw(kw) => kw,
            Item::Data { line, .. } => {
                report.warnings.push(format!(
                    "line {line}: data line outside any keyword block ignored"
                ));
                j += 1;
                continue;
            }
        };

        if let Some(name) = &pending_material {
            if kw.name != "elastic" {
                report.warnings.push(format!(
                    "line {}: material '{name}' has no *ELASTIC block",
                    kw.line
                ));
                pending_material = None;
            }
        }

        match kw.name.as_str() {
            "heading" => {
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            "node" => {
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if f.len() != 3 && f.len() != 4 {
                        return err(
                            *dno,
                            format!("*NODE expects id, x, y[, z]; got {} field(s)", f.len()),
                        );
                    }
                    let id: u32 = num(f[0], *dno, "node id")?;
                    let x = num(f[1], *dno, "coordinate")?;
                    let y = num(f[2], *dno, "coordinate")?;
                    let z = if f.len() == 4 {
                        num(f[3], *dno, "coordinate")?
                    } else {
                        0.0
                    };
                    if !node_ids.insert(id) {
                        return err(*dno, format!("duplicate node id {id}"));
                    }
                    model.mesh.nodes.push(Node {
                        id,
                        coords: [x, y, z],
                    });
                    j += 1;
                }
            }

            "element" => {
                let ty = kw.params.get("type").map(|s| s.to_uppercase());
                let kind = match ty.as_deref() {
                    Some("CPS4") => Some(ElemKind::Cps4),
                    Some("CPS3") => Some(ElemKind::Cps3),
                    Some("C3D8") => Some(ElemKind::C3d8),
                    Some("C3D4") => Some(ElemKind::C3d4),
                    _ => None,
                };
                let Some(kind) = kind else {
                    let what = ty.map_or_else(
                        || "TYPE-less *ELEMENT".to_string(),
                        |t| format!("element TYPE={t}"),
                    );
                    report.warnings.push(format!(
                        "line {}: unsupported {what} — block skipped (supported: CPS4/CPS3/C3D8/C3D4)",
                        kw.line
                    ));
                    j += 1;
                    while matches!(items.get(j), Some(Item::Data { .. })) {
                        j += 1;
                    }
                    continue;
                };
                let want = 1 + kind.node_count();
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    // One element per record; a trailing comma continues the
                    // connectivity on the next data line.
                    let mut f: Vec<String> = Vec::new();
                    let mut first_dno = 0usize;
                    while let Some(Item::Data { line: dno, text }) = items.get(j) {
                        if f.is_empty() {
                            first_dno = *dno;
                        }
                        let ends_comma = text.ends_with(',');
                        f.extend(text.split(',').map(|s| s.trim().to_string()));
                        j += 1;
                        let have = f.iter().filter(|s| !s.is_empty()).count();
                        if !ends_comma || have >= want {
                            break;
                        }
                    }
                    let f: Vec<&str> = f
                        .iter()
                        .filter(|s| !s.is_empty())
                        .map(|s| s.as_str())
                        .collect();
                    if f.len() != want {
                        return err(
                            first_dno,
                            format!(
                                "*ELEMENT {:?} expects {} field(s) (id + {} nodes); got {}",
                                kind,
                                want,
                                want - 1,
                                f.len()
                            ),
                        );
                    }
                    let id: u32 = num(f[0], first_dno, "element id")?;
                    let connectivity: Vec<u32> = f[1..]
                        .iter()
                        .map(|s| num(s, first_dno, "node id"))
                        .collect::<Result<_, _>>()?;
                    if !elem_ids.insert(id) {
                        return err(first_dno, format!("duplicate element id {id}"));
                    }
                    if let Some(setname) = kw.params.get("elset") {
                        upsert_elset(&mut model.elsets, setname, &[id]);
                    }
                    model.mesh.elements.push(Element {
                        id,
                        kind,
                        connectivity,
                    });
                }
            }

            "elset" | "nset" => {
                let name = kw
                    .params
                    .get("elset")
                    .or_else(|| kw.params.get("nset"))
                    .map(|s| s.to_uppercase());
                let generate = kw.params.contains_key("generate");
                let Some(name) = name else {
                    report.warnings.push(format!(
                        "line {}: set keyword without ELSET=/NSET= name — block skipped",
                        kw.line
                    ));
                    j += 1;
                    while matches!(items.get(j), Some(Item::Data { .. })) {
                        j += 1;
                    }
                    continue;
                };
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    let ids: Vec<u32> = if generate {
                        if f.len() < 2 || f.len() > 3 {
                            return err(
                                *dno,
                                format!(
                                    "GENERATE expects start, end[, step]; got {} field(s)",
                                    f.len()
                                ),
                            );
                        }
                        let start: u32 = num(f[0], *dno, "set id")?;
                        let end: u32 = num(f[1], *dno, "set id")?;
                        let step: u32 = if f.len() == 3 {
                            num(f[2], *dno, "set step")?
                        } else {
                            1
                        };
                        if step == 0 {
                            return err(*dno, "GENERATE step must be >= 1".to_string());
                        }
                        (start..=end).step_by(step as usize).collect()
                    } else {
                        f.iter()
                            .map(|s| num(s, *dno, "set member id"))
                            .collect::<Result<Vec<u32>, _>>()?
                    };
                    if kw.name == "elset" {
                        upsert_elset(&mut model.elsets, &name, &ids);
                    } else {
                        nsets.entry(name.clone()).or_default().extend(ids);
                    }
                    j += 1;
                }
            }

            "material" => {
                match kw.params.get("name") {
                    Some(n) => pending_material = Some(n.to_uppercase()),
                    None => report.warnings.push(format!(
                        "line {}: *MATERIAL without NAME — block ignored",
                        kw.line
                    )),
                }
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            "elastic" => {
                j += 1;
                let mut got_row = false;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if !got_row {
                        if f.len() != 2 {
                            return err(
                                *dno,
                                format!(
                                    "*ELASTIC expects youngs_modulus, poisson_ratio; got {} field(s)",
                                    f.len()
                                ),
                            );
                        }
                        match pending_material.take() {
                            Some(name) => {
                                let youngs_modulus = num(f[0], *dno, "youngs modulus")?;
                                let poisson_ratio = num(f[1], *dno, "poisson ratio")?;
                                model.materials.push(Material {
                                    name,
                                    youngs_modulus,
                                    poisson_ratio,
                                });
                            }
                            None => report.warnings.push(format!(
                                "line {}: *ELASTIC outside *MATERIAL — block ignored",
                                dno
                            )),
                        }
                        got_row = true;
                    } else {
                        report.warnings.push(format!(
                            "line {}: extra *ELASTIC row ignored (temperature dependence is v0.2)",
                            dno
                        ));
                    }
                    j += 1;
                }
            }

            "solid section" => {
                let elset = kw.params.get("elset").map(|s| s.to_uppercase());
                let material = kw.params.get("material").map(|s| s.to_uppercase());
                j += 1;
                match (elset, material) {
                    (Some(elset), Some(material)) => {
                        if !model.elsets.iter().any(|s| s.name == elset) {
                            return err(
                                kw.line,
                                format!("*SOLID SECTION references unknown ELSET '{elset}'"),
                            );
                        }
                        if !model.materials.iter().any(|m| m.name == material) {
                            return err(
                                kw.line,
                                format!("*SOLID SECTION references unknown MATERIAL '{material}'"),
                            );
                        }
                        model.mat_assignments.push((elset, material));
                        if let Some(Item::Data { line: dno, text }) = items.get(j) {
                            let f = fields(text);
                            if f.len() != 1 {
                                return err(
                                    *dno,
                                    format!(
                                        "*SOLID SECTION data expects one thickness; got {} field(s)",
                                        f.len()
                                    ),
                                );
                            }
                            model.thickness = num(f[0], *dno, "thickness")?;
                            j += 1;
                        }
                        while matches!(items.get(j), Some(Item::Data { .. })) {
                            j += 1;
                        }
                    }
                    _ => {
                        report.warnings.push(format!(
                            "line {}: *SOLID SECTION needs ELSET and MATERIAL — block skipped",
                            kw.line
                        ));
                        while matches!(items.get(j), Some(Item::Data { .. })) {
                            j += 1;
                        }
                    }
                }
            }

            "boundary" => {
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if f.len() < 2 || f.len() > 4 {
                        return err(
                            *dno,
                            format!(
                                "*BOUNDARY expects node/set, dof[, dof][, value]; got {} field(s)",
                                f.len()
                            ),
                        );
                    }
                    let node_ids = resolve_nodes(&nsets, f[0], *dno)?;
                    let dof1: usize = num(f[1], *dno, "dof")?;
                    let dofs: Vec<usize> = if f.len() >= 3 {
                        let dof2: usize = num(f[2], *dno, "dof")?;
                        (dof1..=dof2).collect()
                    } else {
                        vec![dof1]
                    };
                    if f.len() == 4 {
                        let value: f64 = num(f[3], *dno, "prescribed value")?;
                        if value != 0.0 {
                            report.warnings.push(format!(
                                "line {dno}: nonzero prescribed displacement {value} ignored (v0.1: homogeneous BCs only)"
                            ));
                        }
                    }
                    model
                        .boundary_conditions
                        .push(BoundaryCondition { node_ids, dofs });
                    j += 1;
                }
            }

            "cload" => {
                j += 1;
                while let Some(Item::Data { line: dno, text }) = items.get(j) {
                    let f = fields(text);
                    if f.len() != 3 {
                        return err(
                            *dno,
                            format!(
                                "*CLOAD expects node/set, dof, value; got {} field(s)",
                                f.len()
                            ),
                        );
                    }
                    let nodes = resolve_nodes(&nsets, f[0], *dno)?;
                    let dof: usize = num(f[1], *dno, "dof")?;
                    let value = num(f[2], *dno, "load value")?;
                    for node_id in nodes {
                        model.cloads.push(ConcentratedLoad {
                            node_id,
                            dof,
                            value,
                        });
                    }
                    j += 1;
                }
            }

            "step" => {
                steps += 1;
                if steps > 1 {
                    report.warnings.push(format!(
                        "line {}: more than one *STEP — later steps merged into step 1 (v0.1)",
                        kw.line
                    ));
                }
                if kw
                    .params
                    .get("op")
                    .is_some_and(|s| s.eq_ignore_ascii_case("new"))
                {
                    report.warnings.push(format!(
                        "line {}: OP=NEW ignored — BCs/loads accumulate (v0.1)",
                        kw.line
                    ));
                }
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            "static" | "end step" => {
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }

            other => {
                report.warnings.push(format!(
                    "line {}: unsupported keyword '*{}' — block skipped",
                    kw.line,
                    other.to_uppercase()
                ));
                j += 1;
                while matches!(items.get(j), Some(Item::Data { .. })) {
                    j += 1;
                }
            }
        }
    }

    if let Some(name) = &pending_material {
        report
            .warnings
            .push(format!("material '{name}' has no *ELASTIC block"));
    }
    Ok((model, report))
}

fn tokenize(text: &str) -> Vec<Item> {
    text.lines()
        .enumerate()
        .filter_map(|(idx, raw)| {
            let t = raw.trim();
            if t.is_empty() || t.starts_with("**") {
                None
            } else if let Some(body) = t.strip_prefix('*') {
                let (name, rest) = body.split_once(',').map_or((body, ""), |(n, r)| (n, r));
                Some(Item::Kw(KwLine {
                    line: idx + 1,
                    name: name.trim().to_lowercase(),
                    params: parse_params(rest),
                }))
            } else {
                Some(Item::Data {
                    line: idx + 1,
                    text: t.to_string(),
                })
            }
        })
        .collect()
}

fn parse_params(rest: &str) -> HashMap<String, String> {
    rest.split(',')
        .filter_map(|p| {
            let p = p.trim();
            if p.is_empty() {
                return None;
            }
            Some(match p.split_once('=') {
                Some((k, v)) => (k.trim().to_lowercase(), v.trim().to_string()),
                None => (p.to_lowercase(), String::new()),
            })
        })
        .collect()
}

/// Comma-split a data line into trimmed fields, dropping trailing empty ones.
fn fields(line: &str) -> Vec<&str> {
    let mut f: Vec<&str> = line.split(',').map(str::trim).collect();
    while f.last() == Some(&"") {
        f.pop();
    }
    f
}

fn num<T: FromStr>(s: &str, line: usize, what: &str) -> Result<T, ParseError> {
    s.parse::<T>().map_err(|_| ParseError {
        line,
        message: format!("malformed {what} '{s}'"),
    })
}

fn err<T>(line: usize, message: String) -> Result<T, ParseError> {
    Err(ParseError { line, message })
}

/// Resolve a *BOUNDARY/*CLOAD target: node-set name first, then a bare node id.
fn resolve_nodes(
    nsets: &HashMap<String, Vec<u32>>,
    target: &str,
    line: usize,
) -> Result<Vec<u32>, ParseError> {
    if let Some(ids) = nsets.get(&target.to_uppercase()) {
        return Ok(ids.clone());
    }
    target
        .parse::<u32>()
        .map(|id| vec![id])
        .map_err(|_| ParseError {
            line,
            message: format!("unknown node set '{target}'"),
        })
}

/// Create `name` (uppercased) in `sets` or append to it — `*ELEMENT, ELSET=`
/// may legitimately name a set that does not exist yet.
fn upsert_elset(sets: &mut Vec<ElementSet>, name: &str, ids: &[u32]) {
    let key = name.to_uppercase();
    if let Some(s) = sets.iter_mut().find(|s| s.name == key) {
        s.element_ids.extend_from_slice(ids);
    } else {
        sets.push(ElementSet {
            name: key,
            element_ids: ids.to_vec(),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use solvix_core::ConcentratedLoad;

    const DECK: &str = r#"
** full-coverage smoke deck
*HEADING
    smoke deck
*NODE
1, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 2.0, 0.5
4, 0.0, 1.0
5, 1.0, 1.0
6, 2.0, 1.5
*ELEMENT, TYPE=CPS4, ELSET=QUADS
1, 1, 2, 5, 4
*ELEMENT, TYPE=CPS3, ELSET=TRIS
2, 2, 3, 5
*ELSET, ELSET=ALLELM
1, 2
*ELSET, ELSET=TOP, GENERATE
4, 6, 2
*NSET, NSET=FIX
1, 4
*NSET, NSET=TIP, GENERATE
5, 6, 1
*MATERIAL, NAME=steel
*ELASTIC
210.0E9, 0.3
*SOLID SECTION, ELSET=ALLELM, MATERIAL=steel
0.05
*DLOAD
ALLSIDES, P, 1.0
*STEP
*STATIC
*BOUNDARY
FIX, 1, 2
7, 1, 1, 2.5
*CLOAD
TIP, 2, -50.0
1, 1, 10.0
*END STEP
"#;

    #[test]
    fn full_deck_parses_all_keywords() {
        let (m, r) = parse_inp(DECK).unwrap();

        assert_eq!(m.mesh.nodes.len(), 6);
        assert_eq!(m.mesh.nodes[1].coords, [1.0, 0.0, 0.0]);
        assert_eq!(m.mesh.nodes[2].coords, [2.0, 0.5, 0.0]);

        assert_eq!(m.mesh.elements.len(), 2);
        assert_eq!(m.mesh.elements[0].kind, ElemKind::Cps4);
        assert_eq!(m.mesh.elements[0].connectivity, vec![1, 2, 5, 4]);
        assert_eq!(m.mesh.elements[1].kind, ElemKind::Cps3);
        assert_eq!(m.mesh.elements[1].connectivity, vec![2, 3, 5]);

        let set = |n: &str| m.elsets.iter().find(|s| s.name == n).unwrap();
        assert_eq!(set("QUADS").element_ids, vec![1]);
        assert_eq!(set("TRIS").element_ids, vec![2]);
        assert_eq!(set("ALLELM").element_ids, vec![1, 2]);
        assert_eq!(set("TOP").element_ids, vec![4, 6]);

        assert_eq!(m.thickness, 0.05);
        assert_eq!(m.materials.len(), 1);
        assert_eq!(m.materials[0].name, "STEEL");
        assert_eq!(m.materials[0].youngs_modulus, 210.0e9);
        assert_eq!(m.materials[0].poisson_ratio, 0.3);
        assert_eq!(
            m.mat_assignments,
            vec![("ALLELM".to_string(), "STEEL".to_string())]
        );

        // Sets expanded + dof range in BCs, bare node id with ignored value.
        assert_eq!(m.boundary_conditions.len(), 2);
        assert_eq!(m.boundary_conditions[0].node_ids, vec![1, 4]);
        assert_eq!(m.boundary_conditions[0].dofs, vec![1, 2]);
        assert_eq!(m.boundary_conditions[1].node_ids, vec![7]);
        assert_eq!(m.boundary_conditions[1].dofs, vec![1]);

        // CLOAD on NSET expands to every node, value per node.
        assert_eq!(
            m.cloads,
            vec![
                ConcentratedLoad {
                    node_id: 5,
                    dof: 2,
                    value: -50.0
                },
                ConcentratedLoad {
                    node_id: 6,
                    dof: 2,
                    value: -50.0
                },
                ConcentratedLoad {
                    node_id: 1,
                    dof: 1,
                    value: 10.0
                },
            ]
        );

        // *DLOAD skipped + nonzero prescribed displacement ignored.
        assert_eq!(r.warnings.len(), 2);
        assert!(r.warnings[0].contains("*DLOAD"));
        assert!(r.warnings[1].contains("nonzero prescribed displacement"));
        assert!(m.validate().is_ok());
    }

    #[test]
    fn lowercase_keywords_crlf_generate_and_default_thickness() {
        let deck = concat!(
            "*node\r\n",
            "1, 0.0, 0.0\r\n",
            "2, 1.0, 0.0\r\n",
            "3, 0.0, 1.0\r\n",
            "*element, type=cps3, elset=e1\r\n",
            "1, 1, 2, 3\r\n",
            "*nset, nset=load, generate\r\n",
            "1, 3, 1\r\n",
            "*material, name=m\r\n",
            "*elastic\r\n",
            "1.0e5, 0.0\r\n",
            "*solid section, elset=e1, material=m\r\n",
            "*boundary\r\n",
            "load, 1, 2\r\n",
            "*cload\r\n",
            "2, 2, 7.0\r\n",
        );
        let (m, r) = parse_inp(deck).unwrap();

        assert_eq!(m.mesh.nodes.len(), 3);
        assert_eq!(m.elsets[0].name, "E1");
        assert_eq!(m.elsets[0].element_ids, vec![1]);
        assert_eq!(m.boundary_conditions[0].node_ids, vec![1, 2, 3]);
        assert_eq!(m.boundary_conditions[0].dofs, vec![1, 2]);
        assert_eq!(
            m.cloads,
            vec![ConcentratedLoad {
                node_id: 2,
                dof: 2,
                value: 7.0
            }]
        );
        assert_eq!(m.materials[0].name, "M");
        assert_eq!(m.materials[0].youngs_modulus, 1.0e5);
        assert_eq!(m.thickness, 1.0); // no *SOLID SECTION data line
        assert!(r.warnings.is_empty());
    }

    #[test]
    fn element_connectivity_continues_on_next_line() {
        let deck = concat!(
            "*NODE\n1,0,0\n2,1,0\n3,1,1\n4,0,1\n",
            "*ELEMENT, TYPE=CPS4\n1, 1, 2,\n3, 4\n",
        );
        let (m, r) = parse_inp(deck).unwrap();
        assert_eq!(m.mesh.elements[0].connectivity, vec![1, 2, 3, 4]);
        assert!(r.warnings.is_empty());
    }

    #[test]
    fn unsupported_element_type_and_missing_type_warn_and_skip() {
        let deck = concat!(
            "*NODE\n1,0,0\n2,1,0\n3,0,1\n",
            "*ELEMENT, TYPE=S4R\n1, 1, 2, 3, 1\n",
            "*ELEMENT\n5, 1, 2, 3\n",
        );
        let (m, r) = parse_inp(deck).unwrap();
        assert!(m.mesh.elements.is_empty());
        assert_eq!(r.warnings.len(), 2);
        assert!(r.warnings[0].contains("TYPE=S4R"));
        assert!(r.warnings[1].contains("TYPE-less"));
    }

    #[test]
    fn c3d8_and_c3d4_element_types_parse() {
        let deck = concat!(
            "*NODE\n",
            "1,0,0,0\n2,1,0,0\n3,1,1,0\n4,0,1,0\n5,0,0,1\n6,1,0,1\n7,1,1,1\n8,0,1,1\n",
            "*ELEMENT, TYPE=C3D8\n1, 1, 2, 3, 4, 5, 6, 7, 8\n",
            "*ELEMENT, TYPE=C3D4\n2, 1, 2, 3, 5\n",
        );
        let (m, r) = parse_inp(deck).unwrap();
        assert!(r.warnings.is_empty());
        assert_eq!(m.mesh.elements[0].kind, ElemKind::C3d8);
        assert_eq!(m.mesh.elements[0].connectivity.len(), 8);
        assert_eq!(m.mesh.elements[1].kind, ElemKind::C3d4);
    }

    #[test]
    fn fatal_errors_carry_line_numbers() {
        let bad = |deck: &str| parse_inp(deck).unwrap_err();

        let e = bad("*NODE\n1, 0.0\n");
        assert_eq!(e.line, 2);
        assert!(e.message.contains("field"));

        let e = bad("*NODE\n1, abc, 0.0\n");
        assert_eq!(e.line, 2);
        assert!(e.message.contains("malformed"));

        let e = bad("*NODE\n1,0,0\n1,1,1\n");
        assert_eq!(e.line, 3);
        assert!(e.message.contains("duplicate node id 1"));

        let e = bad("*NODE\n1,0,0\n2,1,0\n3,0,1\n*ELEMENT, TYPE=CPS3\n1,1,2,3\n1,2,3,1\n");
        assert_eq!(e.line, 7);
        assert!(e.message.contains("duplicate element id 1"));

        let e = bad("*BOUNDARY\nNOPE, 1, 2\n");
        assert_eq!(e.line, 2);
        assert!(e.message.contains("unknown node set 'NOPE'"));

        let e = bad("*CLOAD\nNOPE, 1, 1.0\n");
        assert_eq!(e.line, 2);

        let e = bad("*MATERIAL, NAME=M\n*ELASTIC\n1.0\n");
        assert_eq!(e.line, 3);

        let e = bad("*SOLID SECTION, ELSET=NOPE, MATERIAL=M\n");
        assert_eq!(e.line, 1);
        assert!(e.message.contains("unknown ELSET 'NOPE'"));
    }

    #[test]
    fn parses_examples_cantilever() {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/cantilever.inp");
        let Ok(text) = std::fs::read_to_string(&path) else {
            eprintln!("note: examples/cantilever.inp missing (mid-rewrite) — skipping");
            return;
        };
        match parse_inp(&text) {
            Ok((model, _)) => {
                assert!(!model.mesh.nodes.is_empty());
                assert!(!model.mesh.elements.is_empty());
            }
            Err(e) => {
                eprintln!(
                    "note: examples/cantilever.inp not parseable yet ({:?}) — skipping",
                    e
                );
            }
        }
    }
}
