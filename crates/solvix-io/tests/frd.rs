use solvix_io::{parse_frd, FrdResult};

fn fixture(name: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name);
    std::fs::read_to_string(path).unwrap()
}

#[test]
fn parses_fixed_column_fortran_records() {
    let result = parse_frd(&fixture("beam_fixed.frd")).unwrap();

    assert_eq!(
        result.nodes,
        vec![
            (1, [0.0, 0.0, 0.0]),
            (2, [1.0, 0.0, 0.0]),
            (10, [-1.25, 0.25, 3.75]),
        ]
    );
    assert_eq!(
        result.displacements,
        vec![
            (1, [1.23456e-3, -2.34567e-4, 3_456.78]),
            (2, [-1.0, 0.0, 9.87654e-5]),
            (10, [0.5, -0.6, 0.7]),
        ]
    );
    assert_eq!(
        result.stresses6,
        vec![
            (1, [100.0, 200.0, 300.0, 400.0, -500.0, 600.0]),
            (2, [-7.0, 8.0, -9.0, 10.0, 11.0, 12.0]),
        ]
    );

    // The ELEM block is unsupported and reported once; its data lines are
    // not individually reported. Everything after `-5` is ignored.
    assert_eq!(
        result.report.warnings.len(),
        1,
        "{:?}",
        result.report.warnings
    );
    assert!(result.report.warnings[0].contains("line 10"));
    assert!(result.report.warnings[0].contains("'ELEM'"));
}

#[test]
fn parses_comma_records_and_missing_stress() {
    let result = parse_frd(&fixture("plate_comma.frd")).unwrap();

    assert_eq!(
        result.nodes,
        vec![
            (1, [0.0, 0.0, 0.0]),
            (2, [2.0, 0.0, 0.0]),
            (3, [0.0, 1.5, 0.0])
        ]
    );
    assert_eq!(
        result.displacements,
        vec![
            (1, [1.0e-3, -2.0e-3, 0.0]),
            (2, [0.0, 0.0, 350.0]),
            (3, [-4.25, 5.5, -6.75]),
        ]
    );
    assert!(result.stresses6.is_empty());
    assert!(result.report.warnings.is_empty());
}

#[test]
fn parses_stress_continuations_and_abaqus_result_marker() {
    let frd = concat!(
        "    1UDISP\n",
        "-1          5 1.00000E+00 2.00000E+00 3.00000E+00\n",
        "    1USTRESS\n",
        "7, 1.0, -2.0, 3.0\n",
        "-2, 4.0, -5.0, 6.0\n",
        "1, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5\n",
        "-3\n"
    );
    let result = parse_frd(frd).unwrap();

    assert_eq!(result.displacements, vec![(5, [1.0, 2.0, 3.0])]);
    assert_eq!(
        result.stresses6,
        vec![
            (7, [1.0, -2.0, 3.0, 4.0, -5.0, 6.0]),
            (1, [7.5, 8.5, 9.5, 10.5, 11.5, 12.5]),
        ]
    );
    assert!(result.report.warnings.is_empty());
}

#[test]
fn malformed_numbers_carry_exact_line_numbers() {
    let comma = concat!("1CNODE\n", "1, 0.0, 0.0, 0.0\n", "2, 0.0, oops, 0.0\n");
    let error = parse_frd(comma).unwrap_err();
    assert_eq!(error.line, 3);
    assert!(
        error.message.contains("malformed node coordinate 'oops'"),
        "{}",
        error.message
    );

    // Concatenated negative values leave no separator, so the bad field
    // must be found through the fixed columns rather than whitespace.
    let fixed = concat!(
        "    1CNODE\n",
        "-1          1 0.00000E+00 0.00000E+00 0.00000E+00\n",
        "-1          2 0.00000E+0O 0.00000E+00 0.00000E+00\n",
        "-4\n"
    );
    let error = parse_frd(fixed).unwrap_err();
    assert_eq!(error.line, 3);
    assert!(
        error.message.contains("malformed node coordinate"),
        "{}",
        error.message
    );

    let disp = concat!("1CNODE\n", "1, 0, 0, 0\n", "1CDISP\n", "1, 1.0, bad, 2.0\n");
    let error = parse_frd(disp).unwrap_err();
    assert_eq!(error.line, 4);
    assert!(
        error
            .message
            .contains("malformed displacement component 'bad'"),
        "{}",
        error.message
    );

    // Non-ASCII input must be rejected, never sliced as if it were ASCII.
    let unicode = concat!("1CNODE\n", "1 1.0 ó 3.0\n");
    let error = parse_frd(unicode).unwrap_err();
    assert_eq!(error.line, 2);
    assert!(
        error.message.contains("malformed node coordinate 'ó'"),
        "{}",
        error.message
    );

    let stress = concat!(
        "1CSTRESS\n",
        "1, 10.0, 20.0, 30.0\n",
        "-2 40.0, bad, 60.0\n"
    );
    let error = parse_frd(stress).unwrap_err();
    assert_eq!(error.line, 3);
    assert!(
        error.message.contains("malformed stress component 'bad'"),
        "{}",
        error.message
    );
}

#[test]
fn stress_record_shape_errors_are_reported() {
    let error = parse_frd("1CSTRESS\n7, 1.0, 2.0\n-1\n").unwrap_err();
    assert_eq!(error.line, 2);
    assert!(
        error
            .message
            .contains("incomplete STRESS record for node 7 — got 2 of six components"),
        "{}",
        error.message
    );

    let error = parse_frd("1CSTRESS\n-2 1.0 2.0 3.0\n").unwrap_err();
    assert_eq!(error.line, 2);
    assert_eq!(error.message, "STRESS continuation has no matching record");

    let error = parse_frd("1CSTRESS\n7, 1.0, 2.0, 3.0\n4.0, 5.0, 6.0, 7.0\n").unwrap_err();
    assert_eq!(error.line, 3);
    assert_eq!(
        error.message,
        "STRESS continuation overflows six components (7)"
    );

    let error = parse_frd("1CDISP\n1, 1.0, 2.0\n").unwrap_err();
    assert_eq!(error.line, 2);
    assert_eq!(
        error.message,
        "DISP record expects id, ux, uy, uz — got 2 component(s)"
    );
}

#[test]
fn tolerates_headers_unsupported_blocks_and_stray_records() {
    let frd = concat!(
        "    1C\n",
        "    1PSTEP                       1\n",
        "  100N     1           1    1         1\n",
        "-4\n",
        "    1CELEM\n",
        "         1         2         3         4\n",
        "-4\n",
        "    1CRF\n",
        "1, 1.0\n",
        "-3\n",
        "    1CNODE\n",
        "1 1.5 2.5 3.5\n",
        "-5\n",
        "1, 9.0, 9.0, 9.0\n"
    );
    let result = parse_frd(frd).unwrap();

    assert_eq!(result.nodes, vec![(1, [1.5, 2.5, 3.5])]);
    assert_eq!(
        result.report.warnings.len(),
        2,
        "{:?}",
        result.report.warnings
    );
    assert!(result.report.warnings[0].contains("'ELEM'"));
    assert!(result.report.warnings[1].contains("'RF'"));

    // Records outside a known block are skipped, but visible in the report.
    let stray = parse_frd("1CNODE\n1, 0, 0, 0\n9999\nstray junk\n").unwrap();
    assert_eq!(stray.nodes, vec![(1, [0.0, 0.0, 0.0])]);
    assert_eq!(stray.report.warnings.len(), 1);
    assert!(stray.report.warnings[0].contains("line 4"));
}

#[test]
fn parses_abaqus_written_block_headers() {
    let frd = concat!(
        "    1C\n",
        "    1PSTEP                       1\n",
        "  100N     2           1    1         1\n",
        "-4\n",
        "-1          1 0.00000E+00 0.00000E+00 0.00000E+00\n",
        "-1          2 1.00000E+00 0.00000E+00 0.00000E+00\n",
        "-3\n",
        "  100C        4011                       1    1         1\n",
        "    1UDISP\n",
        "-4\n",
        "-1          1 1.00000E-03 0.00000E+00 0.00000E+00\n",
        "-1          2-1.00000E-03 2.00000E-03 0.00000E+00\n",
        "-3\n",
        "-5\n"
    );
    let result = parse_frd(frd).unwrap();

    assert_eq!(
        result.nodes,
        vec![(1, [0.0, 0.0, 0.0]), (2, [1.0, 0.0, 0.0])]
    );
    assert_eq!(
        result.displacements,
        vec![(1, [1.0e-3, 0.0, 0.0]), (2, [-1.0e-3, 2.0e-3, 0.0])]
    );
    assert!(result.report.warnings.is_empty());
}

#[test]
fn numeric_data_is_not_mistaken_for_a_block_header() {
    let result = parse_frd("1CNODE\n100 1.0 2.0 3.0\n-1\n").unwrap();

    assert_eq!(result.nodes, vec![(100, [1.0, 2.0, 3.0])]);
    assert!(result.report.warnings.is_empty());
}

#[test]
fn result_type_is_model_independent_and_defaultable() {
    let empty: FrdResult = parse_frd("    1C\n").unwrap();
    assert!(empty.nodes.is_empty());
    assert!(empty.displacements.is_empty());
    assert!(empty.stresses6.is_empty());
    assert!(empty.report.warnings.is_empty());
}
