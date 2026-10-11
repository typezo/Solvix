//! Conformal triangular meshing inside polygons-with-holes
//! (issue #53, v0.8). Zero dependencies.
//!
//! Pipeline (constrained Delaunay v1):
//! 1. tessellate the region boundary loops (arcs at the documented
//!    [`crate::geom`] segment ceiling, lines split to ~`target`);
//! 2. interior points on a spaced grid accepted by the size field
//!    `h(p) = min(dist(p, boundary))/2 + target`: a candidate is kept
//!    only when no existing point lies within `0.7 * h(p)`, which
//!    grades the mesh from `target` at the boundary to coarse far
//!    field;
//! 3. Bowyer-Watson incremental Delaunay over all points (super
//!    triangle, walking point location, cavity BFS);
//! 4. boundary recovery: every boundary edge missing from the
//!    triangulation is restored by flipping the nearest crossing edge
//!    (convex quad only), falling back to splitting the constraint at
//!    its midpoint; bounded iteration budget, hard error on overrun;
//! 5. cleanup: 5 Laplacian-smoothing passes (interior nodes only,
//!    moves reverted if they invert a triangle) plus sliver deletion
//!    by shortest-edge collapse for aspect > 20.
//!
//! Ceilings (documented): cocircular degeneracies are broken with an
//! epsilon tie-break (slightly non-Delaunay but planar in practice);
//! containment trimming is by triangle centroid (exact once boundary
//! edges are recovered); no guaranteed angle bounds — quality is
//! measured per mesh and asserted in tests.

use std::collections::{HashMap, HashSet};

use solvix_core::{ElemKind, Element, Mesh, Node};

use crate::geom::{self, Region, Sketch};

const NIL: usize = usize::MAX;
const SMOOTH_PASSES: usize = 5;
const SLIVER_ASPECT: f64 = 20.0;
const NODE_CAP: usize = 80_000;

/// Mesh quality summary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quality {
    pub min_angle_deg: f64,
    pub max_aspect: f64,
}

/// Result of [`mesh_sketch`]: nodes, CCW triangles (node indices into
/// `nodes`), boundary loops as node-index cycles (outer loop first),
/// and measured quality.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshResult {
    pub nodes: Vec<[f64; 2]>,
    pub tris: Vec<[usize; 3]>,
    pub boundary: Vec<Vec<usize>>,
    pub quality: Quality,
}

/// Conformally mesh `sketch` with nominal boundary edge size
/// `target_size`. Errors on unclosed profiles, holes crossing the
/// profile, or (rare) failed boundary recovery.
pub fn mesh_sketch(sketch: &Sketch, target_size: f64) -> Result<MeshResult, String> {
    if !target_size.is_finite() || target_size <= 0.0 {
        return Err(format!("invalid target size {target_size}"));
    }
    let region = geom::region_from_sketch(sketch, target_size)?;
    let mut loops: Vec<Vec<[f64; 2]>> = Vec::with_capacity(1 + region.holes.len());
    loops.push(region.outer.clone());
    loops.extend(region.holes.iter().cloned());

    let mut data: Vec<[f64; 2]> = Vec::new();
    let mut constraints: Vec<(usize, usize)> = Vec::new();
    for l in &loops {
        let cleaned = clean_loop(l);
        if cleaned.len() < 3 {
            return Err("degenerate boundary loop (fewer than 3 points)".into());
        }
        let start = data.len();
        for p in cleaned {
            data.push(p);
        }
        // +3: Delaunay::new prepends the super-triangle vertices.
        for i in start..data.len() - 1 {
            constraints.push((i + 3, i + 4));
        }
        constraints.push((data.len() - 1 + 3, start + 3));
    }

    fill_interior(&region, &loops, target_size, &mut data)?;

    let mut d = Delaunay::new(data)?;
    for pi in 3..d.pts.len() {
        d.insert(pi)?;
    }
    d.recover_edges(&constraints)?;
    d.trim(&region);

    let boundary = mesh_boundary_cycles(&d.pts, &d.tris);
    let mut boundary_nodes = vec![false; d.pts.len()];
    let mut boundary_edges: HashSet<(usize, usize)> = HashSet::new();
    for cyc in &boundary {
        for i in 0..cyc.len() {
            boundary_nodes[cyc[i]] = true;
            let (a, b) = (cyc[i], cyc[(i + 1) % cyc.len()]);
            boundary_edges.insert((a.min(b), a.max(b)));
        }
    }
    d.smooth(SMOOTH_PASSES, &boundary_nodes, &region);
    d.collapse_slivers(&boundary_edges);

    // Compact to the used vertices.
    let mut remap = vec![NIL; d.pts.len()];
    let mut nodes: Vec<[f64; 2]> = Vec::new();
    let mut tris: Vec<[usize; 3]> = Vec::new();
    for t in d.tris.iter().filter(|t| t.alive) {
        let v = t.v.map(|v| match remap[v] {
            NIL => {
                remap[v] = nodes.len();
                nodes.push(d.pts[v]);
                nodes.len() - 1
            }
            n => n,
        });
        tris.push(v);
    }
    let boundary: Vec<Vec<usize>> = boundary
        .iter()
        .filter(|c| c.iter().all(|&v| remap[v] != NIL))
        .map(|c| c.iter().map(|&v| remap[v]).collect())
        .collect();

    if tris.is_empty() {
        return Err("triangulation came out empty".into());
    }
    let mut quality = Quality {
        min_angle_deg: f64::MAX,
        max_aspect: 0.0,
    };
    for t in &tris {
        let q = tri_quality(&nodes, t);
        quality.min_angle_deg = quality.min_angle_deg.min(q.0);
        quality.max_aspect = quality.max_aspect.max(q.1);
    }

    Ok(MeshResult {
        nodes,
        tris,
        boundary,
        quality,
    })
}

/// Build a `solvix_core` CPS3 mesh (node/element ids from 1).
pub fn to_cps3_mesh(r: &MeshResult) -> Mesh {
    Mesh {
        nodes: r
            .nodes
            .iter()
            .enumerate()
            .map(|(i, &p)| Node {
                id: i as u32 + 1,
                coords: [p[0], p[1], 0.0],
            })
            .collect(),
        elements: r
            .tris
            .iter()
            .enumerate()
            .map(|(i, t)| Element {
                id: i as u32 + 1,
                kind: ElemKind::Cps3,
                connectivity: vec![t[0] as u32 + 1, t[1] as u32 + 1, t[2] as u32 + 1],
            })
            .collect(),
    }
}

fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}

fn clean_loop(l: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut v: Vec<[f64; 2]> = Vec::with_capacity(l.len());
    for &p in l {
        if v.last().is_none_or(|&q| dist(p, q) > 1e-12) {
            v.push(p);
        }
    }
    while v.len() > 1 && dist(v[0], v[v.len() - 1]) <= 1e-12 {
        v.pop();
    }
    v
}

/// Interior fill: spaced grid candidates over the region bbox,
/// accepted when no existing point lies within `0.7 * h(p)` with
/// `h(p) = dist(p, boundary)/2 + target`. Spatial-hash neighbor
/// lookup; distance to boundary over the tessellated segments.
fn fill_interior(
    region: &Region,
    loops: &[Vec<[f64; 2]>],
    target: f64,
    data: &mut Vec<[f64; 2]>,
) -> Result<(), String> {
    let nboundary = data.len();
    let mut min = [f64::MAX, f64::MAX];
    let mut max = [f64::MIN, f64::MIN];
    for l in loops {
        for p in l.iter() {
            for k in 0..2 {
                min[k] = min[k].min(p[k]);
                max[k] = max[k].max(p[k]);
            }
        }
    }
    let key = |p: [f64; 2]| {
        (
            (p[0] / target).floor() as i64,
            (p[1] / target).floor() as i64,
        )
    };
    let mut bins: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, p) in data.iter().take(nboundary).enumerate() {
        bins.entry(key(*p)).or_default().push(i);
    }

    let nx = ((max[0] - min[0]) / target) as usize + 1;
    let ny = ((max[1] - min[1]) / target) as usize + 1;
    for j in 0..ny {
        for i in 0..nx {
            let p = [
                min[0] + (i as f64 + 0.5) * target,
                min[1] + (j as f64 + 0.5) * target,
            ];
            if !region.contains(p) {
                continue;
            }
            let mut d = f64::MAX;
            for l in loops {
                let n = l.len();
                for s in 0..n {
                    d = d.min(geom::point_segment_dist2(p, l[s], l[(s + 1) % n]).sqrt());
                }
            }
            let r = 0.7 * (d / 2.0 + target);
            let rings = (r / target).ceil() as i64;
            let (bx, by) = key(p);
            let mut blocked = false;
            'outer: for dy in -rings..=rings {
                for dx in -rings..=rings {
                    if let Some(idx) = bins.get(&(bx + dx, by + dy)) {
                        for &v in idx {
                            if dist(p, data[v]) < r {
                                blocked = true;
                                break 'outer;
                            }
                        }
                    }
                }
            }
            if blocked {
                continue;
            }
            data.push(p);
            bins.entry(key(p)).or_default().push(data.len() - 1);
            if data.len() > NODE_CAP {
                return Err(format!(
                    "interior fill exceeded the {NODE_CAP}-node cap (target too small)"
                ));
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------
// Delaunay core
// ---------------------------------------------------------------------

#[derive(Clone)]
struct Tri {
    v: [usize; 3],
    nb: [usize; 3],
    alive: bool,
}

struct Delaunay {
    pts: Vec<[f64; 2]>,
    tris: Vec<Tri>,
    mark: Vec<u32>,
    stamp: u32,
    eps_edge: f64,
    eps_in: f64,
    pert: f64,
    last: usize,
}

fn orient2d(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
}

fn incircle(a: [f64; 2], b: [f64; 2], c: [f64; 2], p: [f64; 2]) -> f64 {
    let (adx, ady) = (a[0] - p[0], a[1] - p[1]);
    let (bdx, bdy) = (b[0] - p[0], b[1] - p[1]);
    let (cdx, cdy) = (c[0] - p[0], c[1] - p[1]);
    let (aa, bb, cc) = (
        adx * adx + ady * ady,
        bdx * bdx + bdy * bdy,
        cdx * cdx + cdy * cdy,
    );
    adx * (bdy * cc - cdy * bb) - ady * (bdx * cc - cdx * bb) + aa * (bdx * cdy - cdx * bdy)
}

/// Proper crossing of the two segments (strict sign change on both).
fn segments_cross(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2], eps: f64) -> bool {
    let d1 = orient2d(c, d, a);
    let d2 = orient2d(c, d, b);
    let d3 = orient2d(a, b, c);
    let d4 = orient2d(a, b, d);
    ((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps))
        && ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps))
}

/// Deterministic equidistributed offset in [-0.5, 0.5) (golden-ratio
/// low-discrepancy sequence).
fn fract_golden(i: usize) -> f64 {
    (i as f64 * 0.6180339887498949).fract() - 0.5
}

/// Perturb `p` by <= `pert` per axis (breaks exact collinearity and
/// cocircularity of lattice-like input; the perturbation is ~1000x
/// the orientation epsilons and 1e-9 of the bbox, far below every
/// geometric tolerance downstream).
fn perturbed(p: [f64; 2], idx: usize, pert: f64) -> [f64; 2] {
    [
        p[0] + pert * fract_golden(2 * idx + 1),
        p[1] + pert * fract_golden(2 * idx + 2),
    ]
}

impl Delaunay {
    /// `data` points; three super-triangle vertices are prepended
    /// (indices 0..3).
    fn new(data: Vec<[f64; 2]>) -> Result<Self, String> {
        if data.len() < 3 {
            return Err("fewer than three points to triangulate".into());
        }
        let mut min = [f64::MAX, f64::MAX];
        let mut max = [f64::MIN, f64::MIN];
        for p in &data {
            for k in 0..2 {
                min[k] = min[k].min(p[k]);
                max[k] = max[k].max(p[k]);
            }
        }
        let (cx, cy) = (0.5 * (min[0] + max[0]), 0.5 * (min[1] + max[1]));
        let diag = (max[0] - min[0]).hypot(max[1] - min[1]);
        let m = 10.0 * (diag + 1.0);
        let pert = 1e-9 * (diag + 1.0);
        let mut pts = vec![[cx - m, cy - 0.9876543 * m], [cx + m, cy - m], [cx, cy + m]];
        for (i, p) in data.into_iter().enumerate() {
            pts.push(perturbed(p, i, pert));
        }
        let eps_edge = 1e-12 * (diag + 1.0).powi(2);
        let eps_in = 1e-11 * (diag + 1.0).powi(4);
        Ok(Delaunay {
            pts,
            tris: vec![Tri {
                v: [0, 1, 2],
                nb: [NIL; 3],
                alive: true,
            }],
            mark: vec![0],
            stamp: 0,
            eps_edge,
            eps_in,
            pert,
            last: 0,
        })
    }

    fn incircle_tri(&self, t: usize, p: [f64; 2]) -> f64 {
        let v = self.tris[t].v;
        incircle(self.pts[v[0]], self.pts[v[1]], self.pts[v[2]], p)
    }

    /// Walk to a triangle whose closed interior contains `p`
    /// (straight walk: follow the segment from the current
    /// triangle's centroid to `p`).
    fn locate(&mut self, p: [f64; 2]) -> Result<usize, String> {
        if !self.tris[self.last].alive {
            self.last = self.tris.iter().position(|t| t.alive).unwrap_or(0);
        }
        let mut t = self.last;
        for _ in 0..4 * self.tris.len() + 64 {
            let v = self.tris[t].v;
            let o = [
                orient2d(self.pts[v[0]], self.pts[v[1]], p),
                orient2d(self.pts[v[1]], self.pts[v[2]], p),
                orient2d(self.pts[v[2]], self.pts[v[0]], p),
            ];
            if o.iter().all(|&x| x > -self.eps_edge) {
                self.last = t;
                return Ok(t);
            }
            // Cross the edge strictly intersected by centroid -> p.
            let q = [
                (self.pts[v[0]][0] + self.pts[v[1]][0] + self.pts[v[2]][0]) / 3.0,
                (self.pts[v[0]][1] + self.pts[v[1]][1] + self.pts[v[2]][1]) / 3.0,
            ];
            let mut step = None;
            for i in 0..3 {
                let (a, b) = (self.pts[v[(i + 1) % 3]], self.pts[v[(i + 2) % 3]]);
                if segments_cross(q, p, a, b, self.eps_edge) {
                    step = Some(i);
                    break;
                }
            }
            // Fallback for degenerate centroid lines: most-negative
            // orientation.
            let s = step.unwrap_or_else(|| {
                o.iter()
                    .enumerate()
                    .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .unwrap()
                    .0
            });
            let u = self.tris[t].nb[s];
            if u == NIL {
                self.last = t;
                return Ok(t); // numeric escape at the hull
            }
            t = u;
        }
        Err("point location walk failed (degenerate triangulation)".into())
    }

    /// Bowyer-Watson insertion of existing vertex `pi`.
    fn insert(&mut self, pi: usize) -> Result<(), String> {
        let p = self.pts[pi];
        let t0 = self.locate(p)?;
        self.stamp += 1;
        let s = self.stamp;
        self.mark[t0] = s;
        let mut queue = vec![t0];
        let mut cavity = Vec::new();
        while let Some(t) = queue.pop() {
            cavity.push(t);
            for i in 0..3 {
                let u = self.tris[t].nb[i];
                if u != NIL && self.mark[u] != s && self.incircle_tri(u, p) > self.eps_in {
                    self.mark[u] = s;
                    queue.push(u);
                }
            }
        }
        // Boundary edges of the cavity: (a, b, outer neighbor, owner).
        let mut bnd: Vec<(usize, usize, usize, usize)> = Vec::new();
        for &t in &cavity {
            for i in 0..3 {
                let u = self.tris[t].nb[i];
                if u == NIL || self.mark[u] != s {
                    bnd.push((
                        self.tris[t].v[(i + 1) % 3],
                        self.tris[t].v[(i + 2) % 3],
                        u,
                        t,
                    ));
                }
            }
        }
        // Order the boundary edges into one cycle around the cavity.
        let mut start_of: HashMap<usize, usize> = HashMap::with_capacity(bnd.len());
        for (i, e) in bnd.iter().enumerate() {
            if start_of.insert(e.0, i).is_some() {
                return Err("degenerate cavity (pinched vertex)".into());
            }
        }
        let mut order = vec![0usize];
        let mut cur = bnd[0].1;
        while cur != bnd[0].0 {
            let Some(&i) = start_of.get(&cur) else {
                return Err("degenerate cavity (broken boundary)".into());
            };
            order.push(i);
            cur = bnd[i].1;
        }
        if order.len() != bnd.len() {
            return Err("degenerate cavity (multi-cycle)".into());
        }
        for &t in &cavity {
            self.tris[t].alive = false;
        }
        let base = self.tris.len();
        let nfan = order.len();
        for &i in &order {
            let (a, b, outer, _) = bnd[i];
            self.tris.push(Tri {
                v: [a, b, pi],
                nb: [NIL, NIL, outer],
                alive: true,
            });
            self.mark.push(0);
        }
        for k in 0..nfan {
            let (t, tn) = (base + k, base + (k + 1) % nfan);
            self.tris[t].nb[0] = tn;
            self.tris[tn].nb[1] = t;
            let (_, _, outer, owner) = bnd[order[k]];
            if outer != NIL {
                for j in 0..3 {
                    if self.tris[outer].nb[j] == owner {
                        self.tris[outer].nb[j] = t;
                        break;
                    }
                }
            }
        }
        self.last = base;
        Ok(())
    }

    fn find_edge(&self, a: usize, b: usize) -> Option<(usize, usize)> {
        self.tris.iter().enumerate().find_map(|(t, tri)| {
            if !tri.alive {
                return None;
            }
            (0..3).find_map(|i| {
                let (u, v) = (tri.v[(i + 1) % 3], tri.v[(i + 2) % 3]);
                (u == a && v == b || u == b && v == a).then_some((t, i))
            })
        })
    }

    fn has_edge(&self, a: usize, b: usize) -> bool {
        self.find_edge(a, b).is_some()
    }

    /// Flip edge `slot` of triangle `t` when the adjoining quad is
    /// strictly convex. Returns whether the flip happened.
    fn try_flip(&mut self, t: usize, slot: usize) -> bool {
        let u = self.tris[t].nb[slot];
        if u == NIL {
            return false;
        }
        let (a, b, c) = (
            self.tris[t].v[(slot + 1) % 3],
            self.tris[t].v[(slot + 2) % 3],
            self.tris[t].v[slot],
        );
        let su = match (0..3).find(|&i| self.tris[u].nb[i] == t) {
            Some(su) => su,
            None => return false,
        };
        let x = self.tris[u].v[su];
        let (pa, pb, pc, px) = (self.pts[a], self.pts[b], self.pts[c], self.pts[x]);
        // Quad c -> a -> x -> b must be convex (orient(b,c,a) > 0
        // already holds because t is CCW).
        if orient2d(pc, pa, px) <= self.eps_edge
            || orient2d(pa, px, pb) <= self.eps_edge
            || orient2d(px, pb, pc) <= self.eps_edge
        {
            return false;
        }
        // Old t = [C, A, B] (slots 0/1/2 = edges (A,B), (B,C), (C,A));
        // old u = [X, B, A] (slot su = (B,A), then (A,X), (X,B)).
        let t_nb1 = self.tris[t].nb[(slot + 1) % 3];
        let t_nb2 = self.tris[t].nb[(slot + 2) % 3];
        let u_nb1 = self.tris[u].nb[(su + 1) % 3];
        let u_nb2 = self.tris[u].nb[(su + 2) % 3];
        self.tris[t] = Tri {
            v: [c, a, x],
            nb: [u_nb1, u, t_nb2],
            alive: true,
        };
        self.tris[u] = Tri {
            v: [c, x, b],
            nb: [u_nb2, t_nb1, t],
            alive: true,
        };
        // Edge (B,C) moved from t to u; edge (A,X) moved from u to t.
        if t_nb1 != NIL {
            for j in 0..3 {
                if self.tris[t_nb1].nb[j] == t {
                    self.tris[t_nb1].nb[j] = u;
                }
            }
        }
        if u_nb1 != NIL {
            for j in 0..3 {
                if self.tris[u_nb1].nb[j] == u {
                    self.tris[u_nb1].nb[j] = t;
                }
            }
        }
        true
    }

    /// Constrained-edge recovery by flipping the nearest crossing
    /// edge; falls back to splitting the constraint at its midpoint.
    fn recover_edges(&mut self, constraints: &[(usize, usize)]) -> Result<(), String> {
        let mut work: Vec<(usize, usize)> = constraints.to_vec();
        let mut budget = 20 * constraints.len() + 10_000;
        while let Some((a, b)) = work.pop() {
            if self.has_edge(a, b) {
                continue;
            }
            budget -= 1;
            if budget == 0 {
                return Err("boundary edge recovery exceeded its iteration budget".into());
            }
            let (pa, pb) = (self.pts[a], self.pts[b]);
            let mut best: Option<(usize, usize, f64)> = None;
            for t in 0..self.tris.len() {
                if !self.tris[t].alive {
                    continue;
                }
                for i in 0..3 {
                    let (u, v) = (self.tris[t].v[(i + 1) % 3], self.tris[t].v[(i + 2) % 3]);
                    if u == a || u == b || v == a || v == b {
                        continue;
                    }
                    if segments_cross(pa, pb, self.pts[u], self.pts[v], self.eps_edge) {
                        // Parameter of the crossing along a -> b.
                        let ab = [pb[0] - pa[0], pb[1] - pa[1]];
                        let cd = [
                            self.pts[v][0] - self.pts[u][0],
                            self.pts[v][1] - self.pts[u][1],
                        ];
                        let denom = ab[0] * cd[1] - ab[1] * cd[0];
                        if denom.abs() <= self.eps_edge {
                            continue;
                        }
                        let s = ((self.pts[u][0] - pa[0]) * cd[1]
                            - (self.pts[u][1] - pa[1]) * cd[0])
                            / denom;
                        if best.is_none_or(|(_, _, bs)| s < bs) {
                            best = Some((t, i, s));
                        }
                    }
                }
            }
            match best {
                None => {
                    if !self.has_edge(a, b) {
                        return Err("boundary edge lost (degenerate triangulation)".into());
                    }
                }
                Some((t, i, _)) => {
                    if self.try_flip(t, i) {
                        work.push((a, b));
                    } else {
                        let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1])];
                        let m = perturbed(mid, self.pts.len(), self.pert);
                        let mi = match self.near_vertex(m) {
                            Some(i) => i,
                            None => {
                                self.pts.push(m);
                                let mi = self.pts.len() - 1;
                                self.insert(mi)?;
                                mi
                            }
                        };
                        work.push((a, mi));
                        work.push((mi, b));
                    }
                }
            }
        }
        Ok(())
    }

    fn near_vertex(&self, p: [f64; 2]) -> Option<usize> {
        let tol = self.eps_edge.sqrt();
        (3..self.pts.len()).find(|&i| dist(p, self.pts[i]) <= tol)
    }

    /// Debug-only topology validation: alive CCW triangles, every
    /// internal edge shared by exactly two triangles, neighbor links
    /// symmetric, no vertex strictly inside another triangle.
    #[cfg(test)]
    fn debug_validate(&self, just: usize) {
        let alive: Vec<usize> = self
            .tris
            .iter()
            .enumerate()
            .filter(|(_, t)| t.alive)
            .map(|(i, _)| i)
            .collect();
        for &t in &alive {
            let tri = &self.tris[t];
            let a2 = orient2d(self.pts[tri.v[0]], self.pts[tri.v[1]], self.pts[tri.v[2]]);
            assert!(
                a2 > 0.0,
                "insert {just}: tri {t} not CCW (area2={a2}): {:?}",
                tri.v.map(|v| self.pts[v])
            );
            for i in 0..3 {
                let u = tri.nb[i];
                if u != NIL {
                    assert!(
                        self.tris[u].alive,
                        "insert {just}: tri {t} dead neighbor {u}"
                    );
                    let back = self.tris[u].nb.contains(&t);
                    assert!(back, "insert {just}: nb link {t}->{u} not symmetric");
                }
            }
        }
        let mut count: HashMap<(usize, usize), usize> = HashMap::new();
        for &t in &alive {
            for i in 0..3 {
                let (a, b) = (self.tris[t].v[(i + 1) % 3], self.tris[t].v[(i + 2) % 3]);
                *count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        assert!(
            count.values().all(|&c| c <= 2),
            "insert {just}: edge over-shared"
        );
        // Every inserted data vertex (3..=just) must appear in some
        // alive triangle.
        for v in 3..=just.min(self.pts.len() - 1) {
            assert!(
                alive.iter().any(|&t| self.tris[t].v.contains(&v)),
                "insert {just}: vertex {v} {:?} lost from triangulation",
                self.pts[v]
            );
        }
    }

    /// Drop super triangles and everything outside the region
    /// (centroid containment — exact once boundary edges recovered).
    fn trim(&mut self, region: &Region) {
        for t in self.tris.iter_mut() {
            if !t.alive {
                continue;
            }
            if t.v.iter().any(|&v| v < 3) {
                t.alive = false;
                continue;
            }
            let p = [
                (self.pts[t.v[0]][0] + self.pts[t.v[1]][0] + self.pts[t.v[2]][0]) / 3.0,
                (self.pts[t.v[0]][1] + self.pts[t.v[1]][1] + self.pts[t.v[2]][1]) / 3.0,
            ];
            if !region.contains(p) {
                t.alive = false;
            }
        }
    }

    /// Laplacian smoothing: interior nodes relax toward their
    /// neighbors' centroid; moves that would invert a triangle or
    /// leave the region are reverted.
    fn smooth(&mut self, passes: usize, boundary: &[bool], region: &Region) {
        for _ in 0..passes {
            let mut acc = vec![[0.0f64; 2]; self.pts.len()];
            let mut cnt = vec![0usize; self.pts.len()];
            for t in self.tris.iter().filter(|t| t.alive) {
                for i in 0..3 {
                    let (a, b) = (t.v[i], t.v[(i + 1) % 3]);
                    acc[a][0] += self.pts[b][0];
                    acc[a][1] += self.pts[b][1];
                    cnt[a] += 1;
                }
            }
            let mut moved = vec![None; self.pts.len()];
            for v in 3..self.pts.len() {
                if boundary.get(v).copied().unwrap_or(true) || cnt[v] == 0 {
                    continue;
                }
                moved[v] = Some([
                    0.5 * self.pts[v][0] + 0.5 * acc[v][0] / cnt[v] as f64,
                    0.5 * self.pts[v][1] + 0.5 * acc[v][1] / cnt[v] as f64,
                ]);
            }
            for (v, mv) in moved.iter().enumerate().skip(3) {
                let Some(q) = *mv else { continue };
                let old = self.pts[v];
                self.pts[v] = q;
                let bad = !region.contains(q)
                    || self.tris.iter().any(|t| {
                        t.alive
                            && t.v.contains(&v)
                            && orient2d(self.pts[t.v[0]], self.pts[t.v[1]], self.pts[t.v[2]])
                                <= self.eps_edge
                    });
                if bad {
                    self.pts[v] = old;
                }
            }
        }
    }

    /// Delete slivers (aspect > [`SLIVER_ASPECT`]) by collapsing
    /// their shortest edge into an interior endpoint. Skips
    /// collapses that would invert a neighbor or orphan a boundary
    /// edge.
    fn collapse_slivers(&mut self, boundary_edges: &HashSet<(usize, usize)>) {
        for _ in 0..3 {
            let mut slivers: Vec<(f64, usize)> = self
                .tris
                .iter()
                .enumerate()
                .filter(|(_, t)| t.alive)
                .filter_map(|(i, t)| {
                    let (_, aspect) = tri_quality(&self.pts, &t.v);
                    (aspect > SLIVER_ASPECT).then_some((aspect, i))
                })
                .collect();
            if slivers.is_empty() {
                return;
            }
            slivers.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
            let mut any = false;
            for (_, t) in slivers {
                if !self.tris[t].alive {
                    continue;
                }
                if self.collapse_one(t, boundary_edges) {
                    any = true;
                }
            }
            if !any {
                return;
            }
        }
    }

    fn collapse_one(&mut self, t: usize, boundary_edges: &HashSet<(usize, usize)>) -> bool {
        let v = self.tris[t].v;
        let lens = [
            dist(self.pts[v[0]], self.pts[v[1]]),
            dist(self.pts[v[1]], self.pts[v[2]]),
            dist(self.pts[v[2]], self.pts[v[0]]),
        ];
        let s = (0..3)
            .min_by(|&a, &b| lens[a].partial_cmp(&lens[b]).unwrap())
            .unwrap();
        // Shortest edge (i, j), apex k; collapse the interior endpoint.
        let (i, j, k) = (v[s], v[(s + 1) % 3], v[(s + 2) % 3]);
        let bnd = |x: usize| boundary_edges.iter().any(|&(p, q)| p == x || q == x);
        let (d, keep) = if !bnd(i) {
            (i, j)
        } else if !bnd(j) {
            (j, i)
        } else if !bnd(k) {
            (k, i)
        } else {
            return false; // all vertices on the boundary — leave it
        };
        let newpos = if bnd(keep) {
            self.pts[keep]
        } else {
            [
                0.5 * (self.pts[d][0] + self.pts[keep][0]),
                0.5 * (self.pts[d][1] + self.pts[keep][1]),
            ]
        };
        // Simulate: replace d by keep everywhere, drop degenerates,
        // keep only if every affected triangle stays CCW and no
        // boundary edge is orphaned.
        let old = self.pts[keep];
        self.pts[keep] = newpos;
        let mut ok = true;
        'chk: for tt in self.tris.iter() {
            if !tt.alive {
                continue;
            }
            let w = tt.v.map(|x| if x == d { keep } else { x });
            if w[0] == w[1] || w[1] == w[2] || w[2] == w[0] {
                // This triangle dies; it must not carry a boundary
                // edge other than through keep (which survives).
                for x in 0..3 {
                    let (p, q) = (tt.v[x], tt.v[(x + 1) % 3]);
                    if p != d && q != d && boundary_edges.contains(&(p.min(q), p.max(q))) {
                        ok = false;
                        break 'chk;
                    }
                }
                continue;
            }
            if (tt.v.contains(&d) || tt.v.contains(&keep))
                && orient2d(self.pts[w[0]], self.pts[w[1]], self.pts[w[2]]) <= self.eps_edge
            {
                ok = false;
                break 'chk;
            }
        }
        if !ok {
            self.pts[keep] = old;
            return false;
        }
        for tt in self.tris.iter_mut() {
            if !tt.alive {
                continue;
            }
            for x in tt.v.iter_mut() {
                if *x == d {
                    *x = keep;
                }
            }
            if tt.v[0] == tt.v[1] || tt.v[1] == tt.v[2] || tt.v[2] == tt.v[0] {
                tt.alive = false;
            }
        }
        true
    }
}

/// Boundary node indices of the alive triangulation, as cycles
/// (largest-|area| loop first). Directed edges come from the CCW
/// triangles, so each boundary vertex has at most one outgoing edge.
fn mesh_boundary_cycles(pts: &[[f64; 2]], tris: &[Tri]) -> Vec<Vec<usize>> {
    let mut count: HashMap<(usize, usize), usize> = HashMap::new();
    for t in tris.iter().filter(|t| t.alive) {
        for i in 0..3 {
            let (a, b) = (t.v[(i + 1) % 3], t.v[(i + 2) % 3]);
            *count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let mut next: HashMap<usize, usize> = HashMap::new();
    for t in tris.iter().filter(|t| t.alive) {
        for i in 0..3 {
            let (a, b) = (t.v[(i + 1) % 3], t.v[(i + 2) % 3]);
            if count[&(a.min(b), a.max(b))] == 1 {
                next.insert(a, b);
            }
        }
    }
    let mut starts: Vec<usize> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut cycles = Vec::new();
    for st in starts {
        if !next.contains_key(&st) {
            continue;
        }
        let mut cyc = vec![st];
        let mut cur = st;
        while let Some(&n) = next.get(&cur) {
            next.remove(&cur);
            if n == st {
                break;
            }
            cyc.push(n);
            cur = n;
        }
        if cyc.len() >= 3 {
            cycles.push(cyc);
        }
    }
    cycles.sort_by(|a, b| {
        cycle_area(pts, a)
            .abs()
            .partial_cmp(&cycle_area(pts, b).abs())
            .unwrap()
            .reverse()
    });
    cycles
}

fn cycle_area(pts: &[[f64; 2]], cyc: &[usize]) -> f64 {
    let n = cyc.len();
    let mut s = 0.0;
    for i in 0..n {
        let a = pts[cyc[i]];
        let b = pts[cyc[(i + 1) % n]];
        s += a[0] * b[1] - b[0] * a[1];
    }
    0.5 * s
}

/// (min angle in degrees, max edge / min edge)
fn tri_quality(pts: &[[f64; 2]], t: &[usize; 3]) -> (f64, f64) {
    let (a, b, c) = (pts[t[0]], pts[t[1]], pts[t[2]]);
    let ab = dist(a, b);
    let bc = dist(b, c);
    let ca = dist(c, a);
    let area2 = orient2d(a, b, c);
    let (lmin, lmax) = (ab.min(bc).min(ca), ab.max(bc).max(ca));
    let angle = |x: f64, y: f64, z: f64| {
        ((x * x + y * y - z * z) / (2.0 * x * y))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    };
    let mut mn = angle(ab, ca, bc)
        .min(angle(ab, bc, ca))
        .min(angle(bc, ca, ab));
    if area2 <= 0.0 {
        mn = 0.0; // inverted triangles count as zero quality
    }
    (mn, lmax / lmin.max(1e-300))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Primitive;

    fn square(x0: f64, y0: f64, s: f64) -> Sketch {
        vec![
            Primitive::Line {
                x1: x0,
                y1: y0,
                x2: x0 + s,
                y2: y0,
            },
            Primitive::Line {
                x1: x0 + s,
                y1: y0,
                x2: x0 + s,
                y2: y0 + s,
            },
            Primitive::Line {
                x1: x0 + s,
                y1: y0 + s,
                x2: x0,
                y2: y0 + s,
            },
            Primitive::Line {
                x1: x0,
                y1: y0 + s,
                x2: x0,
                y2: y0,
            },
        ]
    }

    #[test]
    fn delaunay_unit_square_two_triangles() {
        let mut d = Delaunay::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]).unwrap();
        for pi in 3..d.pts.len() {
            d.insert(pi).unwrap();
            d.debug_validate(pi);
        }
        let alive: Vec<&Tri> = d
            .tris
            .iter()
            .filter(|t| t.alive && t.v.iter().all(|&v| v >= 3))
            .collect();
        assert_eq!(alive.len(), 2);
        let area: f64 = alive
            .iter()
            .map(|t| 0.5 * orient2d(d.pts[t.v[0]], d.pts[t.v[1]], d.pts[t.v[2]]))
            .sum();
        assert!((area - 1.0).abs() < 1e-7, "hull area {area}");
    }

    #[test]
    fn cocircular_grid_stays_planar() {
        // 4x4 lattice: exactly cocircular quadruples everywhere.
        let mut pts = Vec::new();
        for j in 0..4 {
            for i in 0..4 {
                pts.push([i as f64 * 0.25, j as f64 * 0.25]);
            }
        }
        let mut d = Delaunay::new(pts).unwrap();
        for pi in 3..d.pts.len() {
            d.insert(pi).unwrap();
        }
        let alive: Vec<Tri> = d
            .tris
            .iter()
            .filter(|t| t.alive && t.v.iter().all(|&v| v >= 3))
            .cloned()
            .collect();
        assert_eq!(alive.len(), 18); // 3x3 cells x 2
        for t in &alive {
            assert!(orient2d(d.pts[t.v[0]], d.pts[t.v[1]], d.pts[t.v[2]]) > 0.0);
        }
        let mut count: HashMap<(usize, usize), usize> = HashMap::new();
        for t in &alive {
            for i in 0..3 {
                let (a, b) = (t.v[(i + 1) % 3], t.v[(i + 2) % 3]);
                *count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        assert!(count.values().all(|&c| c <= 2));
    }

    #[test]
    fn square_with_hole_recovers_boundary() {
        let mut sk = square(0.0, 0.0, 1.0);
        sk.push(Primitive::Circle {
            cx: 0.5,
            cy: 0.5,
            r: 0.2,
        });
        let r = mesh_sketch(&sk, 0.1).unwrap();
        assert!(r.quality.min_angle_deg > 15.0, "{:?}", r.quality);
        assert_eq!(r.boundary.len(), 2);
        let mesh_area: f64 = r
            .tris
            .iter()
            .map(|t| 0.5 * orient2d(r.nodes[t[0]], r.nodes[t[1]], r.nodes[t[2]]))
            .sum();
        // Inscribed hole polygon: area up to ~1% above 1 - pi r^2.
        assert!(
            (mesh_area - (1.0 - std::f64::consts::PI * 0.04)).abs() < 0.02,
            "area {mesh_area}"
        );
    }
}
