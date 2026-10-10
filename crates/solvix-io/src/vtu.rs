//! ASCII and appended-binary inline-XML `.vtu` writers.

use std::collections::HashMap;
use std::path::Path;

use solvix_core::{ElemKind, Model, Solution};

pub(super) fn write_vtu(path: &Path, model: &Model, solution: &Solution) -> std::io::Result<()> {
    let nodes = &model.mesh.nodes;
    let elements = &model.mesh.elements;
    let index: HashMap<u32, usize> = nodes.iter().enumerate().map(|(i, n)| (n.id, i)).collect();

    let mut s = String::new();
    s.push_str("<?xml version=\"1.0\"?>\n");
    s.push_str("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_format=\"ascii\">\n");
    s.push_str("  <UnstructuredGrid>\n");
    s.push_str(&format!(
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
        nodes.len(),
        elements.len()
    ));

    s.push_str("      <Points>\n");
    s.push_str("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for n in nodes {
        s.push_str(&format!(
            "          {} {} {}\n",
            n.coords[0], n.coords[1], n.coords[2]
        ));
    }
    s.push_str("        </DataArray>\n      </Points>\n");

    s.push_str("      <Cells>\n");
    s.push_str(
        "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ",
    );
    for e in elements {
        for &nid in &e.connectivity {
            // ponytail: unresolved node id maps to index 0; Model::validate() upstream prevents this
            s.push_str(&format!("{} ", index.get(&nid).copied().unwrap_or(0)));
        }
    }
    s.push_str("\n        </DataArray>\n");
    s.push_str("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ");
    let mut offset = 0usize;
    for e in elements {
        offset += e.connectivity.len();
        s.push_str(&format!("{offset} "));
    }
    s.push_str("\n        </DataArray>\n");
    s.push_str("        <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n          ");
    for e in elements {
        let vtk_type = match e.kind {
            ElemKind::Cps3 => 5,                   // VTK_TRIANGLE
            ElemKind::Cps4 => 9,                   // VTK_QUAD
            ElemKind::Cps8 | ElemKind::Cpe8 => 23, // VTK_QUADRATIC_QUAD
            ElemKind::Cpe4 => 9,                   // VTK_QUAD
            ElemKind::C3d8 => 12,                  // VTK_HEXAHEDRON
            ElemKind::C3d4 => 10,                  // VTK_TETRA
            ElemKind::B31 => 3,                    // VTK_LINE
        };
        s.push_str(&format!("{vtk_type} "));
    }
    s.push_str("\n        </DataArray>\n      </Cells>\n");

    s.push_str(
        "      <PointData>\n        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n",
    );
    for i in 0..nodes.len() {
        let d = solution
            .displacements
            .get(i)
            .copied()
            .unwrap_or([0.0, 0.0, 0.0]);
        s.push_str(&format!("          {} {} {}\n", d[0], d[1], d[2]));
    }
    s.push_str("        </DataArray>\n      </PointData>\n");

    s.push_str("      <CellData>\n        <DataArray type=\"Float64\" Name=\"VonMises\" format=\"ascii\">\n");
    for i in 0..elements.len() {
        s.push_str(&format!(
            "          {}\n",
            solution.von_mises.get(i).copied().unwrap_or(0.0)
        ));
    }
    s.push_str("        </DataArray>\n      </CellData>\n");

    s.push_str("    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n");
    std::fs::write(path, s)
}

struct AppendedData {
    bytes: Vec<u8>,
    offsets: Vec<usize>,
}

impl AppendedData {
    fn new() -> Self {
        Self {
            bytes: Vec::new(),
            offsets: Vec::new(),
        }
    }

    fn append(&mut self, block: &[u8]) -> std::io::Result<()> {
        let byte_count = i32::try_from(block.len()).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "VTK appended-data block exceeds the Int32 byte-count header",
            )
        })?;
        self.offsets.push(self.bytes.len());
        self.bytes.extend_from_slice(&byte_count.to_le_bytes());
        self.bytes.extend_from_slice(block);
        Ok(())
    }

    fn offset(&self, array: usize) -> usize {
        self.offsets[array]
    }
}

pub(super) fn write_vtu_binary(
    path: &Path,
    model: &Model,
    solution: &Solution,
) -> std::io::Result<()> {
    let nodes = &model.mesh.nodes;
    let elements = &model.mesh.elements;
    let index: HashMap<u32, usize> = nodes.iter().enumerate().map(|(i, n)| (n.id, i)).collect();

    let mut points = Vec::with_capacity(nodes.len() * 3 * 8);
    for node in nodes {
        for coordinate in node.coords {
            points.extend_from_slice(&coordinate.to_le_bytes());
        }
    }

    let mut connectivity = Vec::new();
    for element in elements {
        for &node_id in &element.connectivity {
            // ponytail: unresolved node id maps to index 0; Model::validate() upstream prevents this
            let node_index = index.get(&node_id).copied().unwrap_or(0);
            connectivity.extend_from_slice(&(node_index as i64).to_le_bytes());
        }
    }

    let mut offsets = Vec::with_capacity(elements.len() * 8);
    let mut connectivity_offset = 0_usize;
    for element in elements {
        connectivity_offset += element.connectivity.len();
        offsets.extend_from_slice(&(connectivity_offset as i64).to_le_bytes());
    }

    let mut types = Vec::with_capacity(elements.len() * 8);
    for element in elements {
        let vtk_type: i64 = match element.kind {
            ElemKind::Cps3 => 5,                   // VTK_TRIANGLE
            ElemKind::Cps4 => 9,                   // VTK_QUAD
            ElemKind::Cps8 | ElemKind::Cpe8 => 23, // VTK_QUADRATIC_QUAD
            ElemKind::Cpe4 => 9,                   // VTK_QUAD
            ElemKind::C3d8 => 12,                  // VTK_HEXAHEDRON
            ElemKind::C3d4 => 10,                  // VTK_TETRA
            ElemKind::B31 => 3,                    // VTK_LINE
        };
        types.extend_from_slice(&vtk_type.to_le_bytes());
    }

    let mut displacements = Vec::with_capacity(nodes.len() * 3 * 8);
    for node_index in 0..nodes.len() {
        let displacement = solution
            .displacements
            .get(node_index)
            .copied()
            .unwrap_or([0.0, 0.0, 0.0]);
        for component in displacement {
            displacements.extend_from_slice(&component.to_le_bytes());
        }
    }

    let mut von_mises = Vec::with_capacity(elements.len() * 8);
    for element_index in 0..elements.len() {
        let stress = solution
            .von_mises
            .get(element_index)
            .copied()
            .unwrap_or(0.0);
        von_mises.extend_from_slice(&stress.to_le_bytes());
    }

    let mut appended = AppendedData::new();
    appended.append(&points)?;
    appended.append(&connectivity)?;
    appended.append(&offsets)?;
    appended.append(&types)?;
    appended.append(&displacements)?;
    appended.append(&von_mises)?;

    let mut file = Vec::new();
    file.extend_from_slice(
        br#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
"#,
    );
    file.extend_from_slice(
        format!(
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
            nodes.len(),
            elements.len()
        )
        .as_bytes(),
    );
    file.extend_from_slice(
        format!(
            "      <Points>\n        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"{}\"/>\n      </Points>\n",
            appended.offset(0)
        )
        .as_bytes(),
    );
    file.extend_from_slice(
        format!(
            "      <Cells>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"{}\"/>\n",
            appended.offset(1)
        )
        .as_bytes(),
    );
    file.extend_from_slice(
        format!(
            "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"{}\"/>\n",
            appended.offset(2)
        )
        .as_bytes(),
    );
    file.extend_from_slice(
        format!(
            "        <DataArray type=\"Int64\" Name=\"types\" format=\"appended\" offset=\"{}\"/>\n      </Cells>\n",
            appended.offset(3)
        )
        .as_bytes(),
    );
    file.extend_from_slice(
        format!(
            "      <PointData>\n        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"appended\" offset=\"{}\"/>\n      </PointData>\n",
            appended.offset(4)
        )
        .as_bytes(),
    );
    file.extend_from_slice(
        format!(
            "      <CellData>\n        <DataArray type=\"Float64\" Name=\"VonMises\" format=\"appended\" offset=\"{}\"/>\n      </CellData>\n",
            appended.offset(5)
        )
        .as_bytes(),
    );
    file.extend_from_slice(
        b"    </Piece>\n  </UnstructuredGrid>\n  <AppendedData encoding=\"raw\">_",
    );
    file.extend_from_slice(&appended.bytes);
    file.extend_from_slice(b"\n  </AppendedData>\n</VTKFile>\n");

    std::fs::write(path, file)
}

#[cfg(test)]
mod tests {
    use super::*;
    use solvix_core::{Element, Mesh, Node};

    const APPENDED_MARKER: &[u8] = b"<AppendedData encoding=\"raw\">_";

    fn xml_and_blob_start(file: &[u8]) -> (&str, usize) {
        let marker = file
            .windows(APPENDED_MARKER.len())
            .position(|window| window == APPENDED_MARKER)
            .expect("raw AppendedData marker");
        let blob_start = marker + APPENDED_MARKER.len();
        (
            std::str::from_utf8(&file[..blob_start]).expect("UTF-8 XML header"),
            blob_start,
        )
    }

    fn section<'a>(xml: &'a str, name: &str) -> &'a str {
        let start = xml.find(&format!("<{name}>")).expect("section start");
        let end = xml.find(&format!("</{name}>")).expect("section end");
        &xml[start..end]
    }

    fn data_array<'a>(section: &'a str, name: Option<&str>) -> &'a str {
        let wanted_name = name.map(|name| format!("Name=\"{name}\""));
        section
            .match_indices("<DataArray")
            .find_map(|(start, _)| {
                let relative_end = section[start..].find("/>")?;
                let array = &section[start..start + relative_end + 2];
                wanted_name
                    .as_ref()
                    .is_none_or(|wanted_name| array.contains(wanted_name))
                    .then_some(array)
            })
            .unwrap_or_else(|| panic!("missing DataArray for {name:?}"))
    }

    fn attribute<'a>(element: &'a str, name: &str) -> &'a str {
        let prefix = format!("{name}=\"");
        let start = element.find(&prefix).unwrap_or_else(|| {
            panic!("missing {name} attribute in {element}");
        }) + prefix.len();
        let end = element[start..].find('"').unwrap_or_else(|| {
            panic!("unterminated {name} attribute in {element}");
        }) + start;
        &element[start..end]
    }

    fn appended_block(file: &[u8], blob_start: usize, offset: usize) -> &[u8] {
        let header_start = blob_start + offset;
        let header: [u8; 4] = file[header_start..header_start + 4]
            .try_into()
            .expect("Int32 byte-count header");
        let byte_count = i32::from_le_bytes(header);
        assert!(byte_count >= 0, "negative byte-count header: {byte_count}");
        let payload_start = header_start + 4;
        let payload_end = payload_start + usize::try_from(byte_count).expect("usize length");
        &file[payload_start..payload_end]
    }

    fn f64_values(block: &[u8]) -> Vec<f64> {
        block
            .as_chunks::<8>()
            .0
            .iter()
            .map(|bytes| f64::from_le_bytes(*bytes))
            .collect()
    }

    fn i64_values(block: &[u8]) -> Vec<i64> {
        block
            .as_chunks::<8>()
            .0
            .iter()
            .map(|bytes| i64::from_le_bytes(*bytes))
            .collect()
    }

    fn assert_exact_f64(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "f64 mismatch at index {index}"
            );
        }
    }

    #[test]
    fn vtu_is_well_formed_and_round_trips() {
        let model = Model {
            mesh: Mesh {
                nodes: vec![
                    Node {
                        id: 1,
                        coords: [0.0, 0.0, 0.0],
                    },
                    Node {
                        id: 2,
                        coords: [1.0, 0.0, 0.0],
                    },
                    Node {
                        id: 3,
                        coords: [1.0, 1.0, 0.0],
                    },
                    Node {
                        id: 4,
                        coords: [0.0, 1.0, 0.0],
                    },
                ],
                elements: vec![
                    Element {
                        id: 1,
                        kind: ElemKind::Cps4,
                        connectivity: vec![1, 2, 3, 4],
                    },
                    Element {
                        id: 2,
                        kind: ElemKind::Cps3,
                        connectivity: vec![2, 4, 3],
                    },
                ],
            },
            thickness: 1.0,
            ..Model::default()
        };
        let solution = Solution {
            displacements: vec![
                [0.1, 0.2, 0.0],
                [0.0, -0.5, 0.0],
                [-0.3, 0.0, 0.5],
                [0.0, 0.0, 0.0],
            ],
            von_mises: vec![42.5, 7.25],
            ..Solution::default()
        };

        let path =
            std::env::temp_dir().join(format!("solvix_io_vtu_test_{}.vtu", std::process::id()));
        write_vtu(&path, &model, &solution).unwrap();
        let xml = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();

        assert!(xml.contains(r#"type="UnstructuredGrid""#));
        assert!(xml.contains(r#"byte_format="ascii""#));
        assert!(xml.contains(r#"NumberOfPoints="4""#));
        assert!(xml.contains(r#"NumberOfCells="2""#));
        // Quad (indices 0 1 2 3) then tri (1 3 2); offsets 4 7; types 9 5.
        assert!(xml.contains("0 1 2 3 1 3 2"));
        assert!(xml.contains("4 7"));
        assert!(xml.contains("9 5"));
        assert!(xml.contains("Name=\"Displacement\""));
        assert!(xml.contains("0.1 0.2 0"));
        assert!(xml.contains("-0.3 0 0.5"));
        assert!(xml.contains("Name=\"VonMises\""));
        assert!(xml.contains("42.5"));
        assert!(xml.contains("7.25"));
        assert!(xml.trim_end().ends_with("</VTKFile>"));
    }

    #[test]
    fn vtu_types_for_new_kinds() {
        // CPS8/CPE8 -> VTK_QUADRATIC_QUAD (23), B31 -> VTK_LINE (3).
        let model = Model {
            mesh: Mesh {
                nodes: (1..=11)
                    .map(|i| Node {
                        id: i,
                        coords: [i as f64, 0.0, 0.0],
                    })
                    .collect(),
                elements: vec![
                    Element {
                        id: 1,
                        kind: ElemKind::Cps8,
                        connectivity: vec![1, 2, 3, 4, 5, 6, 7, 8],
                    },
                    Element {
                        id: 2,
                        kind: ElemKind::Cpe8,
                        connectivity: vec![1, 2, 3, 4, 5, 6, 7, 8],
                    },
                    Element {
                        id: 3,
                        kind: ElemKind::Cpe4,
                        connectivity: vec![1, 2, 3, 4],
                    },
                    Element {
                        id: 4,
                        kind: ElemKind::B31,
                        connectivity: vec![9, 10],
                    },
                ],
            },
            thickness: 1.0,
            ..Model::default()
        };
        let path =
            std::env::temp_dir().join(format!("solvix_io_vtu_kinds_{}.vtu", std::process::id()));
        write_vtu(&path, &model, &Solution::default()).unwrap();
        let xml = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert!(xml.contains("23 23 9 3"), "types line: {xml}");
    }

    #[test]
    fn binary_vtu_round_trips_exact_appended_values() {
        let coordinate = f64::from_bits(0x3ff3_c0d6_e2f4_b0d1);
        let displacement = f64::from_bits(0x3fdd_5555_5555_5555);
        let stress = f64::from_bits(0x404a_cccc_cccc_cccc);
        let model = Model {
            mesh: Mesh {
                nodes: vec![
                    Node {
                        id: 10,
                        coords: [f64::MIN_POSITIVE, -0.0, coordinate],
                    },
                    Node {
                        id: 20,
                        coords: [f64::MAX, coordinate, 0.0],
                    },
                    Node {
                        id: 30,
                        coords: [-coordinate, f64::MIN_POSITIVE, -0.0],
                    },
                    Node {
                        id: 40,
                        coords: [0.0, -coordinate, f64::MAX],
                    },
                ],
                elements: vec![
                    Element {
                        id: 1,
                        kind: ElemKind::Cps4,
                        connectivity: vec![20, 40, 10, 30],
                    },
                    Element {
                        id: 2,
                        kind: ElemKind::Cps3,
                        connectivity: vec![20, 30, 40],
                    },
                ],
            },
            thickness: 1.0,
            ..Model::default()
        };
        let solution = Solution {
            displacements: vec![
                [displacement, f64::MAX, -f64::MIN_POSITIVE],
                [-0.0, coordinate, -displacement],
                [f64::MIN_POSITIVE, -coordinate, 0.0],
                [-f64::MAX, displacement, coordinate],
            ],
            von_mises: vec![stress, -stress],
            ..Solution::default()
        };
        let path = std::env::temp_dir().join(format!(
            "solvix_io_vtu_binary_exact_{}.vtu",
            std::process::id()
        ));
        write_vtu_binary(&path, &model, &solution).unwrap();
        let file = std::fs::read(&path).unwrap();
        std::fs::remove_file(&path).unwrap();

        let (xml, blob_start) = xml_and_blob_start(&file);
        assert!(xml.contains(
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
        ));
        assert!(xml.contains("<Piece NumberOfPoints=\"4\" NumberOfCells=\"2\">"));

        let points = data_array(section(xml, "Points"), None);
        assert_eq!(attribute(points, "format"), "appended");
        let points_offset = attribute(points, "offset").parse::<usize>().unwrap();
        let points_block = appended_block(&file, blob_start, points_offset);
        assert_eq!(points_block.len(), 4 * 3 * 8);
        let expected_points: Vec<f64> = model.mesh.nodes.iter().flat_map(|n| n.coords).collect();
        assert_exact_f64(&f64_values(points_block), &expected_points);

        let cells = section(xml, "Cells");
        let connectivity = data_array(cells, Some("connectivity"));
        assert_eq!(attribute(connectivity, "format"), "appended");
        let connectivity_offset = attribute(connectivity, "offset").parse::<usize>().unwrap();
        let connectivity_block = appended_block(&file, blob_start, connectivity_offset);
        assert_eq!(connectivity_block.len(), 7 * 8);
        assert_eq!(i64_values(connectivity_block), vec![1, 3, 0, 2, 1, 2, 3]);

        let offsets = data_array(cells, Some("offsets"));
        let offsets_offset = attribute(offsets, "offset").parse::<usize>().unwrap();
        let offsets_block = appended_block(&file, blob_start, offsets_offset);
        assert_eq!(offsets_block.len(), 2 * 8);
        assert_eq!(i64_values(offsets_block), vec![4, 7]);

        let types = data_array(cells, Some("types"));
        let types_offset = attribute(types, "offset").parse::<usize>().unwrap();
        let types_block = appended_block(&file, blob_start, types_offset);
        assert_eq!(types_block.len(), 2 * 8);
        assert_eq!(i64_values(types_block), vec![9, 5]);

        let displacement_array = data_array(section(xml, "PointData"), Some("Displacement"));
        assert_eq!(attribute(displacement_array, "format"), "appended");
        let displacement_offset = attribute(displacement_array, "offset")
            .parse::<usize>()
            .unwrap();
        let displacement_block = appended_block(&file, blob_start, displacement_offset);
        assert_eq!(displacement_block.len(), 4 * 3 * 8);
        let expected_displacements: Vec<f64> = solution
            .displacements
            .iter()
            .flat_map(|values| *values)
            .collect();
        assert_exact_f64(&f64_values(displacement_block), &expected_displacements);

        let von_mises = data_array(section(xml, "CellData"), Some("VonMises"));
        assert_eq!(attribute(von_mises, "format"), "appended");
        let von_mises_offset = attribute(von_mises, "offset").parse::<usize>().unwrap();
        let von_mises_block = appended_block(&file, blob_start, von_mises_offset);
        assert_eq!(von_mises_block.len(), 2 * 8);
        assert_eq!(f64_values(von_mises_block), solution.von_mises);
    }

    #[test]
    fn binary_vtu_is_smaller_than_ascii_for_a_large_mesh() {
        let node_count = 300_usize;
        let coordinate = f64::from_bits(0x3ff3_c0d6_e2f4_b0d1);
        let displacement = f64::from_bits(0x3fdd_5555_5555_5555);
        let stress = f64::from_bits(0x404a_cccc_cccc_cccc);
        let model = Model {
            mesh: Mesh {
                nodes: (1..=node_count)
                    .map(|id| Node {
                        id: id as u32,
                        coords: [
                            id as f64 * coordinate,
                            -(id as f64) * coordinate,
                            coordinate / f64::from(id as u32),
                        ],
                    })
                    .collect(),
                elements: (1..node_count)
                    .map(|id| Element {
                        id: id as u32,
                        kind: ElemKind::B31,
                        connectivity: vec![id as u32, id as u32 + 1],
                    })
                    .collect(),
            },
            thickness: 1.0,
            ..Model::default()
        };
        let solution = Solution {
            displacements: (1..=node_count)
                .map(|id| {
                    [
                        displacement,
                        -displacement,
                        displacement / f64::from(id as u32),
                    ]
                })
                .collect(),
            von_mises: (1..node_count)
                .map(|id| stress * f64::from(id as u32))
                .collect(),
            ..Solution::default()
        };
        let binary_path = std::env::temp_dir().join(format!(
            "solvix_io_vtu_binary_size_{}.vtu",
            std::process::id()
        ));
        let ascii_path = std::env::temp_dir().join(format!(
            "solvix_io_vtu_ascii_size_{}.vtu",
            std::process::id()
        ));
        write_vtu_binary(&binary_path, &model, &solution).unwrap();
        write_vtu(&ascii_path, &model, &solution).unwrap();
        let binary_size = std::fs::metadata(&binary_path).unwrap().len();
        let ascii_size = std::fs::metadata(&ascii_path).unwrap().len();
        std::fs::remove_file(&binary_path).unwrap();
        std::fs::remove_file(&ascii_path).unwrap();

        assert!(
            binary_size < ascii_size,
            "binary size {binary_size} should be smaller than ASCII size {ascii_size}"
        );
    }
}
