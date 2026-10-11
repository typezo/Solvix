//! 2D geometry kernel (issue #53, v0.8): sketch primitives, boolean
//! circle subtraction, and tessellated boundary regions ready for
//! meshing.
//!
//! Tessellation ceiling (documented): every arc/circle is polygonized
//! before any boolean or containment test. A circle becomes
//! `clamp(round(2*pi*r/target), 16, 128)` chord segments, so hole
//! geometry is exact only to the chord sagitta — worst case
//! `r*(1 - cos(pi/16))` ≈ 1.9% of `r` at the coarse 16-segment cap,
//! ≈ 0.04% at the 128-segment cap. Only circles strictly inside the
//! closed profile can be subtracted; a circle crossing the profile
//! boundary is an error (no general boolean ops in v1). Hole loops
//! must be pairwise disjoint.

/// One sketch primitive. Coordinates in the sketch plane.
///
/// [`Primitive::Arc`] sweeps from `a0` to `a1` (radians); the sweep
/// `a1 - a0` may be negative (clockwise) and `|a1 - a0| <= 2*pi`.
/// [`Primitive::Circle`] is a standalone closed hole candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Primitive {
    Line {
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
    },
    Arc {
        cx: f64,
        cy: f64,
        r: f64,
        a0: f64,
        a1: f64,
    },
    Circle {
        cx: f64,
        cy: f64,
        r: f64,
    },
}

/// A sketch: open profile primitives (lines/arcs) plus circles.
pub type Sketch = Vec<Primitive>;

fn circle_segments(r: f64, target: f64) -> usize {
    (std::f64::consts::TAU * r / target)
        .round()
        .clamp(16.0, 128.0) as usize
}

/// Polygonize `p` into points along its path. Lines are split into
/// `ceil(len/target)` collinear segments; arcs/circles are
/// chord-tessellated within the segment ceiling. A circle yields `n`
/// points (closed loop, not repeated).
pub(crate) fn tessellate(p: &Primitive, target: f64) -> Vec<[f64; 2]> {
    match *p {
        Primitive::Line { x1, y1, x2, y2 } => {
            let len = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt();
            let n = (len / target).ceil().max(1.0) as usize;
            (0..=n)
                .map(|i| {
                    let t = i as f64 / n as f64;
                    [x1 + t * (x2 - x1), y1 + t * (y2 - y1)]
                })
                .collect()
        }
        Primitive::Arc { cx, cy, r, a0, a1 } => {
            let sweep = a1 - a0;
            let max_segs = ((sweep.abs() / std::f64::consts::TAU) * 128.0)
                .ceil()
                .max(1.0);
            let n = (sweep.abs() * r / target).round().clamp(1.0, max_segs) as usize;
            (0..=n)
                .map(|i| {
                    let a = a0 + sweep * i as f64 / n as f64;
                    [cx + r * a.cos(), cy + r * a.sin()]
                })
                .collect()
        }
        Primitive::Circle { cx, cy, r } => {
            let n = circle_segments(r, target);
            (0..n)
                .map(|i| {
                    let a = std::f64::consts::TAU * i as f64 / n as f64;
                    [cx + r * a.cos(), cy + r * a.sin()]
                })
                .collect()
        }
    }
}

/// A planar region: one outer boundary loop plus hole loops, all
/// polygonal (arcs already tessellated). Loop orientation is
/// normalized to counter-clockwise; use [`Region::contains`].
#[derive(Debug, Clone, PartialEq)]
pub struct Region {
    pub outer: Vec<[f64; 2]>,
    pub holes: Vec<Vec<[f64; 2]>>,
}

impl Region {
    /// Even-odd ray-cast containment for a single closed loop.
    pub fn inside_loop(pts: &[[f64; 2]], p: [f64; 2]) -> bool {
        let n = pts.len();
        if n < 3 {
            return false;
        }
        let mut inside = false;
        for i in 0..n {
            let (a, b) = (pts[i], pts[(i + 1) % n]);
            if (a[1] > p[1]) != (b[1] > p[1]) {
                let t = (p[1] - a[1]) / (b[1] - a[1]);
                if p[0] < a[0] + t * (b[0] - a[0]) {
                    inside = !inside;
                }
            }
        }
        inside
    }

    /// Inside the outer loop and outside every hole.
    pub fn contains(&self, p: [f64; 2]) -> bool {
        Self::inside_loop(&self.outer, p) && !self.holes.iter().any(|h| Self::inside_loop(h, p))
    }

    /// Net enclosed area (holes subtracted).
    pub fn area(&self) -> f64 {
        loop_area(&self.outer) - self.holes.iter().map(|h| loop_area(h)).sum::<f64>()
    }
}

fn loop_area(pts: &[[f64; 2]]) -> f64 {
    shoelace(pts).abs()
}

fn shoelace(pts: &[[f64; 2]]) -> f64 {
    let n = pts.len();
    let mut s = 0.0;
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        s += a[0] * b[1] - b[0] * a[1];
    }
    0.5 * s
}

fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

/// Boolean subtract (v1: circle-in-polygon, per issue #53): chain the
/// sketch's lines/arcs into closed loops, take the largest as the
/// outer boundary, treat the remaining loops as holes, and subtract
/// every circle that lies strictly inside the outer loop as a
/// tessellated hole loop.
pub fn region_from_sketch(sketch: &Sketch, target: f64) -> Result<Region, String> {
    let tol = 1e-9
        * sketch
            .iter()
            .flat_map(|p| tessellate(p, target))
            .map(|p| p[0].abs().max(p[1].abs()))
            .fold(1.0f64, f64::max)
            .max(1.0);

    let mut loops: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut circles: Vec<Vec<[f64; 2]>> = Vec::new();
    for p in sketch {
        let pts = tessellate(p, target);
        match p {
            Primitive::Circle { .. } => circles.push(pts),
            _ => loops.push(pts),
        }
    }
    let loops = chain_loops(&loops, tol)?;
    if loops.is_empty() {
        return Err("sketch has no closed line/arc profile to bound a region".into());
    }

    let outer_idx = loops
        .iter()
        .enumerate()
        .max_by(|a, b| loop_area(a.1).total_cmp(&loop_area(b.1)))
        .unwrap()
        .0;
    let outer = oriented(&loops[outer_idx]);
    let mut holes: Vec<Vec<[f64; 2]>> = Vec::new();
    for (i, l) in loops.iter().enumerate() {
        if i == outer_idx {
            continue;
        }
        if !l.iter().all(|&p| Region::inside_loop(&outer, p)) {
            return Err("interior profile loop touches or crosses the outer profile (v1)".into());
        }
        holes.push(oriented(l));
    }

    for c in &circles {
        if c.iter().any(|&p| !Region::inside_loop(&outer, p))
            || c.iter()
                .any(|&p| holes.iter().any(|h| Region::inside_loop(h, p)))
        {
            return Err(format!(
                "circle at ({:.4}, {:.4}) crosses the profile or another hole (v1: strictly-inside holes only)",
                c[0][0], c[0][1]
            ));
        }
        holes.push(c.clone());
    }
    Ok(Region { outer, holes })
}

/// Normalize a loop to counter-clockwise and weld near-duplicate
/// consecutive points.
fn oriented(pts: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut v: Vec<[f64; 2]> = Vec::with_capacity(pts.len());
    for &p in pts {
        if v.last().is_none_or(|&q| dist(p, q) > 1e-12) {
            v.push(p);
        }
    }
    while v.len() > 1 && dist(v[0], v[v.len() - 1]) <= 1e-12 {
        v.pop();
    }
    if shoelace(&v) < 0.0 {
        v.reverse();
    }
    v
}

/// Weld open paths' endpoints (within `tol`) and chain them into
/// closed loops. Linear scan; path counts are tens at most.
fn chain_loops(paths: &[Vec<[f64; 2]>], tol: f64) -> Result<Vec<Vec<[f64; 2]>>, String> {
    let mut used = vec![false; paths.len()];
    let mut loops = Vec::new();
    for (i, p) in paths.iter().enumerate() {
        if !used[i] && p.len() > 2 && dist(p[0], p[p.len() - 1]) <= tol {
            loops.push(p.clone());
            used[i] = true;
        }
    }
    for start in 0..paths.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain: Vec<[f64; 2]> = paths[start].clone();
        loop {
            let tail = chain[chain.len() - 1];
            if chain.len() > 2 && dist(tail, chain[0]) <= tol {
                break; // closed
            }
            let next = (0..paths.len()).find(|&i| {
                !used[i]
                    && (dist(paths[i][0], tail) <= tol
                        || dist(paths[i][paths[i].len() - 1], tail) <= tol)
            });
            let Some(j) = next else {
                return Err(format!(
                    "profile chain dead-ends at ({:.4}, {:.4}) — sketch not closed",
                    tail[0], tail[1]
                ));
            };
            used[j] = true;
            if dist(paths[j][0], tail) <= tol {
                for &p in &paths[j][1..] {
                    chain.push(p);
                }
            } else {
                for &p in paths[j].iter().rev().skip(1) {
                    chain.push(p);
                }
            }
        }
        loops.push(chain);
    }
    Ok(loops)
}

/// Squared distance from `p` to the segment `ab`.
pub(crate) fn point_segment_dist2(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let len2 = ab[0] * ab[0] + ab[1] * ab[1];
    if len2 <= 0.0 {
        return dist(p, a).powi(2);
    }
    let t = (((p[0] - a[0]) * ab[0] + (p[1] - a[1]) * ab[1]) / len2).clamp(0.0, 1.0);
    let q = [a[0] + t * ab[0], a[1] + t * ab[1]];
    dist(p, q).powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(x1: f64, y1: f64, x2: f64, y2: f64) -> Primitive {
        Primitive::Line { x1, y1, x2, y2 }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Sketch {
        vec![
            line(x0, y0, x1, y0),
            line(x1, y0, x1, y1),
            line(x1, y1, x0, y1),
            line(x0, y1, x0, y0),
        ]
    }

    #[test]
    fn rectangle_region_area_and_contains() {
        let r = region_from_sketch(&rect(0.0, 0.0, 2.0, 1.0), 0.1).unwrap();
        assert_eq!(r.holes.len(), 0);
        assert!((r.area() - 2.0).abs() < 1e-9);
        assert!(r.contains([1.0, 0.5]));
        assert!(!r.contains([-1e-6, 0.5]));
        assert!(!r.contains([3.0, 0.5]));
    }

    #[test]
    fn circle_subtract_makes_hole_loop() {
        let mut sk = rect(0.0, 0.0, 1.0, 1.0);
        sk.push(Primitive::Circle {
            cx: 0.5,
            cy: 0.5,
            r: 0.2,
        });
        let r = region_from_sketch(&sk, 0.05).unwrap();
        assert_eq!(r.holes.len(), 1);
        // Inscribed polygon: slightly below pi*r^2, within the ceiling.
        let hole_area = loop_area(&r.holes[0]);
        assert!(
            hole_area > 0.97 * std::f64::consts::PI * 0.04
                && hole_area < std::f64::consts::PI * 0.04,
            "hole area {hole_area}"
        );
        assert!((r.area() - (1.0 - hole_area)).abs() < 1e-12);
        assert!(!r.contains([0.5, 0.5]));
        assert!(r.contains([0.05, 0.05]));
    }

    #[test]
    fn circle_crossing_profile_is_rejected() {
        let mut sk = rect(0.0, 0.0, 1.0, 1.0);
        sk.push(Primitive::Circle {
            cx: 0.95,
            cy: 0.5,
            r: 0.2,
        });
        assert!(region_from_sketch(&sk, 0.05).is_err());
    }

    #[test]
    fn open_profile_is_rejected() {
        let sk = vec![line(0.0, 0.0, 1.0, 0.0)];
        assert!(region_from_sketch(&sk, 0.1).is_err());
    }

    #[test]
    fn arc_tessellation_runs_correct_direction() {
        // Quarter arc from 90 deg down to 0 deg (clockwise sweep).
        let pts = tessellate(
            &Primitive::Arc {
                cx: 0.0,
                cy: 0.0,
                r: 1.0,
                a0: std::f64::consts::FRAC_PI_2,
                a1: 0.0,
            },
            0.1,
        );
        assert!((pts[0][0] - 0.0).abs() < 1e-9 && (pts[0][1] - 1.0).abs() < 1e-9);
        let last = pts[pts.len() - 1];
        assert!((last[0] - 1.0).abs() < 1e-9 && last[1].abs() < 1e-9);
        for p in &pts {
            assert!((dist(*p, [0.0, 0.0]) - 1.0).abs() < 1e-12);
        }
    }
}
