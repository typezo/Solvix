//! ASCII inline-XML `.vtu` writer (VTK UnstructuredGrid, issue #2).

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
            ElemKind::Cps3 => 5,  // VTK_TRIANGLE
            ElemKind::Cps4 => 9,  // VTK_QUAD
            ElemKind::C3d8 => 12, // VTK_HEXAHEDRON
            ElemKind::C3d4 => 10, // VTK_TETRA
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

#[cfg(test)]
mod tests {
    use super::*;
    use solvix_core::{Element, Mesh, Node};

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
}
