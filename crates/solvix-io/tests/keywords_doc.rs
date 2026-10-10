//! Sync guard: every keyword documented in docs/KEYWORDS.md must have a
//! match arm in the parser, and vice versa (issue #22).

/// Keyword names documented in docs/KEYWORDS.md: rows whose first table
/// cell is `*NAME` (uppercase, spaces allowed).
fn doc_keywords(doc: &str) -> Vec<String> {
    doc.lines()
        .filter_map(|l| {
            let l = l.trim();
            if !l.starts_with("| *") {
                return None;
            }
            let cell = l.split('|').nth(1)?.trim();
            let name = cell.strip_prefix('*')?;
            (!name.is_empty() && name.chars().all(|c| c.is_ascii_uppercase() || c == ' '))
                .then(|| name.to_lowercase())
        })
        .collect()
}

/// Keyword names handled by the parser: lines of the big match in
/// parser.rs, which all have the shape `"name" | "name" => {`.
fn parser_keywords(src: &str) -> Vec<String> {
    let mut out = Vec::new();
    for line in src.lines() {
        let t = line.trim();
        if t.starts_with('"') && t.ends_with("=> {") {
            let arms = &t[..t.len() - 4];
            let parts: Vec<&str> = arms.split('"').collect();
            for name in parts.iter().skip(1).step_by(2) {
                out.push((*name).to_string());
            }
        }
    }
    out
}

#[test]
fn keyword_docs_in_sync() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let doc = std::fs::read_to_string(root.join("docs/KEYWORDS.md")).unwrap();
    let src = std::fs::read_to_string(root.join("crates/solvix-io/src/parser.rs")).unwrap();

    let doc_kw = doc_keywords(&doc);
    let src_kw = parser_keywords(&src);

    // Sanity: the curated core is definitely present in both.
    let core = [
        "node",
        "element",
        "nset",
        "elset",
        "material",
        "elastic",
        "density",
        "plastic",
        "solid section",
        "boundary",
        "cload",
        "dload",
        "step",
        "static",
        "end step",
    ];
    for k in core {
        assert!(
            doc_kw.contains(&k.to_string()),
            "doc missing core keyword {k}"
        );
        assert!(
            src_kw.contains(&k.to_string()),
            "parser missing core keyword {k}"
        );
    }
    assert!(doc_kw.len() >= 20, "doc table suspiciously small");

    let mut doc_sorted = doc_kw.clone();
    doc_sorted.sort();
    let mut src_sorted = src_kw.clone();
    src_sorted.sort();
    assert_eq!(
        doc_sorted, src_sorted,
        "docs/KEYWORDS.md and parser.rs keyword sets diverged"
    );
}
