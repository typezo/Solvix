//! Structured mesh generation for translated Abaqus parts (issue #35).
//!
//! v1 shape ceiling (documented in the issue): a part sketch is reduced
//! to its line-segment bounding box plus circular holes. The mesh is a
//! rectangle CPS4 grid; elements whose centroid lies inside a hole
//! circle are dropped and kept nodes within half a cell of the hole
//! radius are snapped onto the circle — the classic "poor-man's hole".
//! A real geometry kernel is v0.7 work.

use std::collections::HashSet;

use solvix_core::{ElemKind, Element, Mesh, Node};

/// Circular hole: center (cx, cy) and radius r.
pub struct Hole {
    pub cx: f64,
    pub cy: f64,
    pub r: f64,
}

/// Build the structured grid over `[x0, x0+w] x [y0, y0+h]` with
/// `nx * ny` CPS4 elements, minus hole interiors. Returns the mesh and
/// human-readable notes for the translation report.
pub fn grid(
    x0: f64,
    y0: f64,
    w: f64,
    h: f64,
    nx: usize,
    ny: usize,
    holes: &[Hole],
) -> (Mesh, Vec<String>) {
    let dx = w / nx as f64;
    let dy = h / ny as f64;
    let id = |i: usize, j: usize| (j * (nx + 1) + i + 1) as u32;

    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push(Node {
                id: id(i, j),
                coords: [x0 + i as f64 * dx, y0 + j as f64 * dy, 0.0],
            });
        }
    }

    let mut elements = Vec::with_capacity(nx * ny);
    let mut dropped = 0usize;
    for j in 0..ny {
        for i in 0..nx {
            let (cx, cy) = (x0 + (i as f64 + 0.5) * dx, y0 + (j as f64 + 0.5) * dy);
            let inside = holes
                .iter()
                .any(|c| ((cx - c.cx).powi(2) + (cy - c.cy).powi(2)).sqrt() < c.r);
            if inside {
                dropped += 1;
                continue;
            }
            elements.push(Element {
                id: (j * nx + i + 1) as u32,
                kind: ElemKind::Cps4,
                connectivity: vec![id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)],
            });
        }
    }

    // Snap nodes in the half-cell band around each hole onto the circle.
    let band = 0.5 * dx.max(dy);
    let mut snapped = 0usize;
    for n in &mut nodes {
        for c in holes {
            let d = ((n.coords[0] - c.cx).powi(2) + (n.coords[1] - c.cy).powi(2)).sqrt();
            if (d - c.r).abs() <= band {
                let s = c.r / d.max(1e-300);
                n.coords[0] = c.cx + s * (n.coords[0] - c.cx);
                n.coords[1] = c.cy + s * (n.coords[1] - c.cy);
                snapped += 1;
                break;
            }
        }
    }

    // Prune unreferenced (hole-interior) nodes so the stiffness stays
    // non-singular.
    let referenced: HashSet<u32> = elements
        .iter()
        .flat_map(|e| e.connectivity.iter().copied())
        .collect();
    let kept = nodes.len() - referenced.len();
    nodes.retain(|n| referenced.contains(&n.id));

    let mut notes = Vec::new();
    if dropped > 0 {
        notes.push(format!(
            "structured mesh: {dropped} element(s) dropped inside circular hole(s)"
        ));
    }
    if snapped > 0 {
        notes.push(format!(
            "structured mesh: {snapped} node(s) snapped onto hole boundary (faceted hole)"
        ));
    }
    if kept > 0 {
        notes.push(format!("structured mesh: {kept} interior node(s) pruned"));
    }

    (Mesh { nodes, elements }, notes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_grid_is_full() {
        let (m, notes) = grid(0.0, 0.0, 1.0, 0.1, 10, 2, &[]);
        assert_eq!(m.nodes.len(), 11 * 3);
        assert_eq!(m.elements.len(), 20);
        assert!(notes.is_empty());
        // CCW ordering: first element's signed area is positive.
        let e = &m.elements[0];
        let c: Vec<[f64; 2]> = e
            .connectivity
            .iter()
            .map(|&id| {
                let n = m.nodes.iter().find(|n| n.id == id).unwrap();
                [n.coords[0], n.coords[1]]
            })
            .collect();
        let area = (c[1][0] - c[0][0]) * (c[1][1] + c[0][1])
            + (c[2][0] - c[1][0]) * (c[2][1] + c[1][1])
            + (c[3][0] - c[2][0]) * (c[3][1] + c[2][1])
            + (c[0][0] - c[3][0]) * (c[0][1] + c[3][1]);
        assert!(area < 0.0, "shoelace (negative = CCW in xy): {area}");
    }

    #[test]
    fn hole_drops_elements_and_snaps_boundary() {
        let hole = [Hole {
            cx: 0.5,
            cy: 0.5,
            r: 0.2,
        }];
        let (m, notes) = grid(0.0, 0.0, 1.0, 1.0, 20, 20, &hole);
        assert!(m.elements.len() < 400);
        // Every node of every element is on or outside the circle.
        for n in &m.nodes {
            let d = ((n.coords[0] - 0.5).powi(2) + (n.coords[1] - 0.5).powi(2)).sqrt();
            assert!(d > 0.2 - 0.06, "node {} deep inside the hole: d={d}", n.id);
        }
        // Some node was snapped onto the radius exactly.
        assert!(
            m.nodes.iter().any(|n| {
                (((n.coords[0] - 0.5).powi(2) + (n.coords[1] - 0.5).powi(2)).sqrt() - 0.2).abs()
                    < 1e-12
            }),
            "no node on the hole radius"
        );
        assert!(!notes.is_empty());
    }
}
