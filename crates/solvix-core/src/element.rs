//! Plane-stress element stiffness: CPS3 (closed form) and CPS4
//! (2x2 Gauss isoparametric bilinear quad).

use crate::SolveError;

/// 3x3 plane-stress constitutive matrix.
pub(crate) fn d_matrix(e: f64, nu: f64) -> [[f64; 3]; 3] {
    let c = e / (1.0 - nu * nu);
    [
        [c, nu * c, 0.0],
        [nu * c, c, 0.0],
        [0.0, 0.0, c * (1.0 - nu) / 2.0],
    ]
}

/// Strain-displacement matrices kept per element for stress recovery:
/// the centroid matrix, plus the 2x2 Gauss-point matrices for quads
/// (ordering matches the stiffness loop: xi-major over (-,+), then (+,+)).
pub(crate) enum CentroidB {
    Tri([[f64; 6]; 3]),
    Quad {
        center: [[f64; 8]; 3],
        gauss: Box<[[[f64; 8]; 3]; 4]>,
    },
}

/// CPS3 stiffness. Returns `None` for a zero-area (degenerate) triangle.
/// Works for either node ordering: the barycentric identity behind `B`
/// holds with the signed area, and the integral uses `|A|`.
pub(crate) fn cps3_stiffness(
    coords: &[[f64; 2]; 3],
    d: &[[f64; 3]; 3],
    t: f64,
) -> Option<([[f64; 6]; 6], CentroidB)> {
    let b = [
        coords[1][1] - coords[2][1],
        coords[2][1] - coords[0][1],
        coords[0][1] - coords[1][1],
    ];
    let c = [
        coords[2][0] - coords[1][0],
        coords[0][0] - coords[2][0],
        coords[1][0] - coords[0][0],
    ];
    let two_a = b[0] * c[1] - b[1] * c[0]; // signed 2*area
    if two_a == 0.0 {
        return None;
    }
    let inv = 1.0 / two_a;
    let mut bm = [[0.0; 6]; 3];
    for i in 0..3 {
        bm[0][2 * i] = b[i] * inv;
        bm[1][2 * i + 1] = c[i] * inv;
        bm[2][2 * i] = c[i] * inv;
        bm[2][2 * i + 1] = b[i] * inv;
    }
    let mut ke = [[0.0; 6]; 6];
    btdb(&mut ke, &bm, d, t * two_a.abs() / 2.0);
    Some((ke, CentroidB::Tri(bm)))
}

const QUAD_CORNERS: [(f64, f64); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];

/// Shape function derivatives at (xi, eta): (dN/dxi, dN/deta), each [f64; 4].
fn quad_shape_derivs(xi: f64, eta: f64) -> ([f64; 4], [f64; 4]) {
    let mut dn_dxi = [0.0; 4];
    let mut dn_deta = [0.0; 4];
    for i in 0..4 {
        let (xi_i, eta_i) = QUAD_CORNERS[i];
        dn_dxi[i] = 0.25 * xi_i * (1.0 + eta * eta_i);
        dn_deta[i] = 0.25 * eta_i * (1.0 + xi * xi_i);
    }
    (dn_dxi, dn_deta)
}

/// B matrix at (xi, eta), or `None` if det(J) is not positive there.
fn quad_b(coords: &[[f64; 2]; 4], xi: f64, eta: f64) -> Option<([[f64; 8]; 3], f64)> {
    let (dn_dxi, dn_deta) = quad_shape_derivs(xi, eta);
    let mut j = [[0.0; 2]; 2];
    for i in 0..4 {
        j[0][0] += dn_dxi[i] * coords[i][0];
        j[0][1] += dn_dxi[i] * coords[i][1];
        j[1][0] += dn_deta[i] * coords[i][0];
        j[1][1] += dn_deta[i] * coords[i][1];
    }
    let detj = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    if detj <= 0.0 {
        return None;
    }
    let inv = 1.0 / detj;
    let mut bm = [[0.0; 8]; 3];
    for i in 0..4 {
        let dnx = inv * (j[1][1] * dn_dxi[i] - j[0][1] * dn_deta[i]);
        let dny = inv * (j[0][0] * dn_deta[i] - j[1][0] * dn_dxi[i]);
        bm[0][2 * i] = dnx;
        bm[1][2 * i + 1] = dny;
        bm[2][2 * i] = dny;
        bm[2][2 * i + 1] = dnx;
    }
    Some((bm, detj))
}

/// CPS4 stiffness with 2x2 Gauss. Returns `None` if det(J) <= 0 at any
/// Gauss point (distorted/inverted element).
pub(crate) fn cps4_stiffness(
    coords: &[[f64; 2]; 4],
    d: &[[f64; 3]; 3],
    t: f64,
) -> Option<([[f64; 8]; 8], CentroidB)> {
    let g = 1.0f64 / (3.0f64).sqrt();
    let mut ke = [[0.0; 8]; 8];
    let mut gauss = [[[0.0; 8]; 3]; 4];
    let mut k = 0;
    for &xi in &[-g, g] {
        for &eta in &[-g, g] {
            let (bm, detj) = quad_b(coords, xi, eta)?;
            btdb(&mut ke, &bm, d, t * detj); // unit Gauss weights
            gauss[k] = bm;
            k += 1;
        }
    }
    let (center, _) = quad_b(coords, 0.0, 0.0)?; // centroid is interior if GPs are
    Some((
        ke,
        CentroidB::Quad {
            center,
            gauss: Box::new(gauss),
        },
    ))
}

/// ke += w * B^T D B for a 3xN strain-displacement matrix.
fn btdb<const N: usize>(ke: &mut [[f64; N]; N], b: &[[f64; N]; 3], d: &[[f64; 3]; 3], w: f64) {
    let mut db = [[0.0; N]; 3];
    for r in 0..3 {
        for c in 0..N {
            db[r][c] = d[r][0] * b[0][c] + d[r][1] * b[1][c] + d[r][2] * b[2][c];
        }
    }
    for a in 0..N {
        for c in a..N {
            let v = w * (b[0][a] * db[0][c] + b[1][a] * db[1][c] + b[2][a] * db[2][c]);
            ke[a][c] += v;
            if a != c {
                ke[c][a] += v;
            }
        }
    }
}

/// Stress components (sxx, syy, sxy) from a strain-displacement matrix,
/// the constitutive matrix and the element nodal displacements.
pub(crate) fn stress_components<const N: usize>(
    b: &[[f64; N]; 3],
    d: &[[f64; 3]; 3],
    ue: &[f64],
) -> [f64; 3] {
    let (exx, eyy, gxy) = strain(b, ue);
    [
        d[0][0] * exx + d[0][1] * eyy,
        d[1][0] * exx + d[1][1] * eyy,
        d[2][2] * gxy,
    ]
}

/// Von Mises from plane-stress components (sxx, syy, sxy).
pub(crate) fn von_mises(s: &[f64; 3]) -> f64 {
    (s[0] * s[0] - s[0] * s[1] + s[1] * s[1] + 3.0 * s[2] * s[2]).sqrt()
}

/// Von Mises stress from a centroid strain-displacement matrix, the
/// constitutive matrix and the element nodal displacements.
pub(crate) fn centroid_von_mises(b: &CentroidB, d: &[[f64; 3]; 3], ue: &[f64]) -> f64 {
    match b {
        CentroidB::Tri(bm) => von_mises(&stress_components(bm, d, ue)),
        CentroidB::Quad { center, .. } => von_mises(&stress_components(center, d, ue)),
    }
}

/// Gauss-point sign pairs in stiffness-loop order (xi-major), matching
/// the `gauss` slot ordering of [`CentroidB::Quad`].
const GP_SIGNS: [(f64, f64); 4] = [(-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), (1.0, 1.0)];

/// Extrapolate 2x2 Gauss-point values to the 4 corner nodes: bilinear
/// shape functions inverted at the Gauss points (Abaqus-style). The
/// equivalent 4x4 matrix is [[a,b,c,b],[b,a,b,c],[c,b,a,b],[b,c,b,a]]
/// with a = 1 + sqrt(3)/2, b = -1/2, c = 1 - sqrt(3)/2 (rows sum to 1;
/// the b = -1/2 sign is what makes constant and bilinear fields exact).
pub(crate) fn extrapolate_to_corners<const M: usize>(gp: &[[f64; M]; 4]) -> [[f64; M]; 4] {
    let sqrt3 = (3.0f64).sqrt();
    let h = |s: f64| 0.5 * (1.0 + s * sqrt3);
    let mut out = [[0.0; M]; 4];
    for (i, &(cx, cy)) in QUAD_CORNERS.iter().enumerate() {
        for (j, &(gx, gy)) in GP_SIGNS.iter().enumerate() {
            let w = h(cx * gx) * h(cy * gy);
            for k in 0..M {
                out[i][k] += w * gp[j][k];
            }
        }
    }
    out
}

fn strain<const N: usize>(b: &[[f64; N]; 3], ue: &[f64]) -> (f64, f64, f64) {
    let mut e = [0.0; 3];
    for r in 0..3 {
        for (c, &u) in ue.iter().enumerate() {
            e[r] += b[r][c] * u;
        }
    }
    (e[0], e[1], e[2])
}

/// Convenience: map a failed element computation to `InvalidModel` with the id.
pub(crate) fn invalid(elem_id: u32, what: &str) -> SolveError {
    SolveError::InvalidModel(format!("element {elem_id}: {what}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_square() -> [[f64; 2]; 4] {
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
    }

    #[test]
    fn cps4_stiffness_is_symmetric_and_kills_rigid_translation() {
        let d = d_matrix(1.0, 0.25);
        let (ke, _) = cps4_stiffness(&unit_square(), &d, 1.0).unwrap();
        for (a, row) in ke.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                assert!((v - ke[c][a]).abs() < 1e-12, "K[{a},{c}] != K[{c},{a}]");
            }
        }
        let mut u = [0.0; 8];
        for i in 0..4 {
            u[2 * i] = 0.7;
            u[2 * i + 1] = -1.3;
        }
        for (a, row) in ke.iter().enumerate() {
            let ku: f64 = row.iter().zip(u).map(|(k, ui)| k * ui).sum();
            assert!(ku.abs() < 1e-12, "row {a} not translation-invariant");
        }
    }

    #[test]
    fn extrapolation_to_corners_is_exact_for_bilinear_fields() {
        let f = |xi: f64, eta: f64| 2.0 + 3.0 * xi - 5.0 * eta + 1.5 * xi * eta;
        let g = 1.0f64 / (3.0f64).sqrt();
        let gp = [[f(-g, -g)], [f(-g, g)], [f(g, -g)], [f(g, g)]];
        let corners = extrapolate_to_corners(&gp);
        for (i, &(cx, cy)) in QUAD_CORNERS.iter().enumerate() {
            assert!(
                (corners[i][0] - f(cx, cy)).abs() < 1e-12,
                "corner {i}: {} != {}",
                corners[i][0],
                f(cx, cy)
            );
        }
        // Constant field: every row of the extrapolation matrix sums to 1.
        let c = extrapolate_to_corners(&[[7.0]; 4]);
        assert!(c.iter().all(|&v| (v[0] - 7.0).abs() < 1e-12));
    }

    #[test]
    fn inverted_quad_is_rejected() {
        let d = d_matrix(1.0, 0.25);
        let mut coords = unit_square();
        coords.swap(2, 3); // clockwise ordering inverts J
        assert!(cps4_stiffness(&coords, &d, 1.0).is_none());
    }

    #[test]
    fn degenerate_triangle_is_rejected() {
        let d = d_matrix(1.0, 0.25);
        let coords = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]];
        assert!(cps3_stiffness(&coords, &d, 1.0).is_none());
    }
}
