use solvix_core::{BoundaryCondition, ConcentratedLoad, ElemKind, ModelError, SolveError};
use solvix_io::parse_msh;

fn fixture(name: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    std::fs::read_to_string(path).unwrap()
}

fn element_set<'a>(model: &'a solvix_core::Model, name: &str) -> &'a solvix_core::ElementSet {
    model
        .elsets
        .iter()
        .find(|set| set.name == name)
        .unwrap_or_else(|| panic!("missing element set {name}"))
}

#[test]
fn parses_legacy_quad_mesh_and_physical_groups() {
    let (model, report) = parse_msh(&fixture("square4.msh")).unwrap();

    assert_eq!(model.mesh.nodes.len(), 6);
    assert_eq!(model.mesh.elements.len(), 2);
    assert_eq!(model.mesh.elements[0].kind, ElemKind::Cps4);
    assert_eq!(model.mesh.elements[0].connectivity, vec![1, 2, 5, 4]);
    assert_eq!(model.mesh.elements[1].kind, ElemKind::Cps4);
    assert_eq!(model.mesh.elements[1].connectivity, vec![2, 3, 6, 5]);

    assert_eq!(element_set(&model, "LEFT_BLOCK").element_ids, vec![1]);
    assert_eq!(element_set(&model, "RIGHT_BLOCK").element_ids, vec![2]);
    assert_eq!(element_set(&model, "ALL").element_ids, vec![1, 2]);
    assert_eq!(
        model.mat_assignments,
        vec![("ALL".to_string(), "DEFAULT".to_string())]
    );
    assert_eq!(model.materials.len(), 1);
    assert_eq!(model.materials[0].name, "DEFAULT");
    assert_eq!(model.materials[0].youngs_modulus, 210.0e9);
    assert_eq!(model.materials[0].poisson_ratio, 0.3);
    assert_eq!(model.thickness, 1.0);
    assert!(model.boundary_conditions.is_empty());
    assert!(model.cloads.is_empty());
    assert!(report
        .warnings
        .iter()
        .any(|warning| warning.contains("created DEFAULT")));
    assert!(matches!(
        model.validate(),
        Err(ModelError::NoBoundaryConditions)
    ));
}

#[test]
fn parses_legacy_triangle_mesh() {
    let (model, _) = parse_msh(&fixture("square_tri.msh")).unwrap();

    assert_eq!(model.mesh.nodes.len(), 6);
    assert_eq!(model.mesh.elements.len(), 4);
    assert!(model
        .mesh
        .elements
        .iter()
        .all(|element| element.kind == ElemKind::Cps3));
    assert_eq!(element_set(&model, "DOMAIN").element_ids, vec![1, 2, 3, 4]);
    assert_eq!(model.mesh.elements[3].connectivity, vec![2, 6, 5]);
}

#[test]
fn parses_msh41_hexa_with_nonsequential_tags() {
    let (model, _) = parse_msh(&fixture("cube_hex.msh")).unwrap();

    assert_eq!(model.mesh.nodes.len(), 8);
    assert_eq!(model.mesh.elements.len(), 1);
    assert_eq!(model.mesh.elements[0].kind, ElemKind::C3d8);
    assert_eq!(
        model.mesh.elements[0].connectivity,
        vec![2, 4, 88, 6, 55, 31, 8, 1]
    );
    assert_eq!(element_set(&model, "CUBE").element_ids, vec![1]);
    assert_eq!(element_set(&model, "STEEL").element_ids, vec![1]);

    let node = model.mesh.nodes.iter().find(|node| node.id == 88).unwrap();
    assert_eq!(node.coords, [1.0, 1.0, 0.0]);
}

#[test]
fn parses_legacy_tetra_mesh() {
    let (model, _) = parse_msh(&fixture("tet.msh")).unwrap();

    assert_eq!(model.mesh.nodes.len(), 4);
    assert_eq!(model.mesh.elements.len(), 1);
    assert_eq!(model.mesh.elements[0].kind, ElemKind::C3d4);
    assert_eq!(model.mesh.elements[0].connectivity, vec![1, 2, 3, 4]);
    assert_eq!(element_set(&model, "TETRA").element_ids, vec![1]);
}

#[test]
fn unsupported_element_type_warns_and_skips() {
    let text = concat!(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n",
        "$Nodes\n3\n1 0 0 0\n2 1 0 0\n3 0 1 0\n$EndNodes\n",
        "$Elements\n2\n",
        "1 1 0 1 2\n",
        "2 2 2 1 1 1 2 3\n",
        "$EndElements\n"
    );
    let (model, report) = parse_msh(text).unwrap();

    assert_eq!(model.mesh.elements.len(), 1);
    assert_eq!(model.mesh.elements[0].kind, ElemKind::Cps3);
    assert!(report
        .warnings
        .iter()
        .any(|warning| warning.contains("unsupported MSH element type 1")));
}

#[test]
fn malformed_number_carries_line_number() {
    let text = concat!(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n",
        "$Nodes\n1\n1 0 abc 0\n$EndNodes\n",
        "$Elements\n1\n1 2 2 1 1 1 1 2\n$EndElements\n"
    );
    let error = parse_msh(text).unwrap_err();

    assert_eq!(error.line, 6);
    assert!(error.message.contains("malformed node coordinate 'abc'"));
}

#[test]
fn imported_square_passes_uniaxial_patch_test() {
    let (mut model, _) = parse_msh(&fixture("square4.msh")).unwrap();
    model.boundary_conditions = vec![
        BoundaryCondition {
            node_ids: vec![1, 4],
            dofs: vec![1],
        },
        BoundaryCondition {
            node_ids: vec![1],
            dofs: vec![2],
        },
    ];
    model.cloads = vec![
        ConcentratedLoad {
            node_id: 3,
            dof: 1,
            value: 500.0e3,
        },
        ConcentratedLoad {
            node_id: 6,
            dof: 1,
            value: 500.0e3,
        },
    ];

    let solution = solvix_core::solve(&model).unwrap();
    let displacement = solution.displacements[2][0];
    let expected = 1.0e6 / 210.0e9;
    assert!((displacement - expected).abs() <= expected * 1.0e-8);
    assert_eq!(solution.displacements[0][0], 0.0);
    for stress in &solution.von_mises {
        assert!((*stress - 1.0e6).abs() <= 1.0e6 * 1.0e-7);
    }

    let bare = parse_msh(&fixture("square4.msh")).unwrap().0;
    assert!(matches!(
        solvix_core::solve(&bare),
        Err(SolveError::NoConstraints)
    ));
}
