//! Element stiffness: CPS3/CPS4 plane stress (2 dofs/node) and
//! C3D8/C3D4 linear solids (3 dofs/node), plus stress recovery helpers.

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

/// 6x6 isotropic 3D Hooke matrix. Component order throughout the crate:
/// (xx, yy, zz, xy, yz, xz). Requires nu < 0.5 (caller checks).
pub(crate) fn d3_matrix(e: f64, nu: f64) -> [[f64; 6]; 6] {
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let diag = lam + 2.0 * mu;
    let mut d = [[0.0; 6]; 6];
    for (r, row) in d.iter_mut().enumerate().take(3) {
        for (c, v) in row.iter_mut().enumerate().take(3) {
            *v = if r == c { diag } else { lam };
        }
    }
    d[3][3] = mu;
    d[4][4] = mu;
    d[5][5] = mu;
    d
}

/// Constitutive law stored per element, dispatched by element kind:
/// plane elements get the plane-stress 3x3, solids the 3D 6x6.
pub(crate) enum Constitutive {
    Plane([[f64; 3]; 3]),
    Solid(Box<[[f64; 6]; 6]>),
}

/// Strain-displacement matrices kept per element for stress recovery:
/// the centroid matrix, plus the Gauss-point matrices for quads and hexes
/// (ordering matches the stiffness loops; see GP_SIGNS / GP3_SIGNS).
pub(crate) enum CentroidB {
    Tri([[f64; 6]; 3]),
    Quad {
        center: [[f64; 8]; 3],
        gauss: Box<[[[f64; 8]; 3]; 4]>,
    },
    Tet(Box<[[f64; 12]; 6]>),
    Hex {
        center: Box<[[f64; 24]; 6]>,
        gauss: Box<[[[f64; 24]; 6]; 8]>,
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
/// Also the det(J) source for the CPS4 consistent mass (modal).
pub(crate) fn quad_b(coords: &[[f64; 2]; 4], xi: f64, eta: f64) -> Option<([[f64; 8]; 3], f64)> {
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

/// Bilinear shape functions at (xi, eta), node order matching
/// [`QUAD_CORNERS`] (used by the consistent mass matrix).
pub(crate) fn quad_shapes(xi: f64, eta: f64) -> [f64; 4] {
    let mut n = [0.0; 4];
    for i in 0..4 {
        let (xi_i, eta_i) = QUAD_CORNERS[i];
        n[i] = 0.25 * (1.0 + xi * xi_i) * (1.0 + eta * eta_i);
    }
    n
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

/// C3D8 corner signs in natural coordinates, Abaqus ordering: nodes 1-4
/// bottom face counter-clockwise, 5-8 top face the same way.
const HEX_CORNERS: [(f64, f64, f64); 8] = [
    (-1.0, -1.0, -1.0),
    (1.0, -1.0, -1.0),
    (1.0, 1.0, -1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
    (1.0, -1.0, 1.0),
    (1.0, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
];

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Inverse of a 3x3 matrix with known non-zero determinant (adjugate/det).
fn inv3(m: &[[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let cof = |r: usize, c: usize| {
        let (a, b) = ((r + 1) % 3, (r + 2) % 3);
        let (d, e) = ((c + 1) % 3, (c + 2) % 3);
        let sign = if (r + c).is_multiple_of(2) { 1.0 } else { -1.0 };
        sign * (m[a][d] * m[b][e] - m[a][e] * m[b][d])
    };
    let id = 1.0 / det;
    let mut out = [[0.0; 3]; 3];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v = cof(j, i) * id; // transposed cofactors
        }
    }
    out
}

/// Trilinear shape function derivatives at (xi, eta, zeta):
/// dn[a][i] = dN_i/d(xi_a) with xi_a in (xi, eta, zeta).
fn hex_shape_derivs(xi: f64, eta: f64, zeta: f64) -> [[f64; 8]; 3] {
    let mut dn = [[0.0; 8]; 3];
    for (i, &(a, b, c)) in HEX_CORNERS.iter().enumerate() {
        dn[0][i] = 0.125 * a * (1.0 + eta * b) * (1.0 + zeta * c);
        dn[1][i] = 0.125 * b * (1.0 + xi * a) * (1.0 + zeta * c);
        dn[2][i] = 0.125 * c * (1.0 + xi * a) * (1.0 + eta * b);
    }
    dn
}

/// B matrix (6x24) at a natural point, or `None` if det(J) is not
/// positive there. J[a][b] = dx_b/d(xi_a); dN/dx = J^-1 * dN/d(xi).
/// Also the det(J) source for the C3D8 consistent mass (modal).
pub(crate) fn hex_b(
    coords: &[[f64; 3]; 8],
    xi: f64,
    eta: f64,
    zeta: f64,
) -> Option<([[f64; 24]; 6], f64)> {
    let dn = hex_shape_derivs(xi, eta, zeta);
    let mut j = [[0.0; 3]; 3];
    for i in 0..8 {
        for a in 0..3 {
            for b in 0..3 {
                j[a][b] += dn[a][i] * coords[i][b];
            }
        }
    }
    let detj = det3(&j);
    if detj <= 0.0 {
        return None;
    }
    let inv = inv3(&j, detj);
    let mut bm = [[0.0; 24]; 6];
    let cols = dn[0].iter().zip(&dn[1]).zip(&dn[2]);
    for (i, ((&d0, &d1), &d2)) in cols.enumerate() {
        let dnv = [d0, d1, d2];
        let mut g = [0.0; 3]; // shape-function gradient in xyz
        for b in 0..3 {
            for a in 0..3 {
                g[b] += inv[b][a] * dnv[a];
            }
        }
        let c = 3 * i;
        bm[0][c] = g[0];
        bm[1][c + 1] = g[1];
        bm[2][c + 2] = g[2];
        bm[3][c] = g[1]; // gxy
        bm[3][c + 1] = g[0];
        bm[4][c + 1] = g[2]; // gyz
        bm[4][c + 2] = g[1];
        bm[5][c] = g[2]; // gxz
        bm[5][c + 2] = g[0];
    }
    Some((bm, detj))
}

/// Trilinear shape functions at (xi, eta, zeta), node order matching
/// [`HEX_CORNERS`] (used by the consistent mass matrix).
pub(crate) fn hex_shapes(xi: f64, eta: f64, zeta: f64) -> [f64; 8] {
    let mut n = [0.0; 8];
    for (i, &(a, b, c)) in HEX_CORNERS.iter().enumerate() {
        n[i] = 0.125 * (1.0 + xi * a) * (1.0 + eta * b) * (1.0 + zeta * c);
    }
    n
}

/// C3D8 stiffness with 2x2x2 Gauss. Returns `None` if det(J) <= 0 at any
/// Gauss point (distorted/inverted element).
pub(crate) fn c3d8_stiffness(
    coords: &[[f64; 3]; 8],
    d: &[[f64; 6]; 6],
) -> Option<([[f64; 24]; 24], CentroidB)> {
    let g = 1.0f64 / (3.0f64).sqrt();
    let mut ke = [[0.0; 24]; 24];
    let mut gauss = [[[0.0; 24]; 6]; 8];
    let mut k = 0;
    for &xi in &[-g, g] {
        for &eta in &[-g, g] {
            for &zeta in &[-g, g] {
                let (bm, detj) = hex_b(coords, xi, eta, zeta)?;
                btdb(&mut ke, &bm, d, detj); // unit Gauss weights
                gauss[k] = bm;
                k += 1;
            }
        }
    }
    let (center, _) = hex_b(coords, 0.0, 0.0, 0.0)?; // centroid is interior if GPs are
    Some((
        ke,
        CentroidB::Hex {
            center: Box::new(center),
            gauss: Box::new(gauss),
        },
    ))
}

/// C3D4 stiffness, constant-strain tetrahedron. Returns `None` for zero
/// volume. Works for either node ordering: shape-function gradients are
/// built from the signed 6V determinant, the integral uses |V|.
pub(crate) fn c3d4_stiffness(
    coords: &[[f64; 3]; 4],
    d: &[[f64; 6]; 6],
) -> Option<([[f64; 12]; 12], CentroidB)> {
    let sub = |i: usize| {
        [
            coords[i][0] - coords[3][0],
            coords[i][1] - coords[3][1],
            coords[i][2] - coords[3][2],
        ]
    };
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let (p1, p2, p3) = (sub(0), sub(1), sub(2));
    let det = p1[0] * (p2[1] * p3[2] - p2[2] * p3[1]) - p1[1] * (p2[0] * p3[2] - p2[2] * p3[0])
        + p1[2] * (p2[0] * p3[1] - p2[1] * p3[0]); // signed 6V
    if det == 0.0 {
        return None;
    }
    let inv = 1.0 / det;
    // grad L_i = (p_j x p_k) / det (cyclic), grad L_4 = -sum of the rest.
    let g1 = cross(p2, p3);
    let g2 = cross(p3, p1);
    let g3 = cross(p1, p2);
    let g4 = [
        -(g1[0] + g2[0] + g3[0]),
        -(g1[1] + g2[1] + g3[1]),
        -(g1[2] + g2[2] + g3[2]),
    ];
    let mut bm = [[0.0; 12]; 6];
    for (i, gi) in [g1, g2, g3, g4].iter().enumerate() {
        let c = 3 * i;
        bm[0][c] = gi[0] * inv;
        bm[1][c + 1] = gi[1] * inv;
        bm[2][c + 2] = gi[2] * inv;
        bm[3][c] = gi[1] * inv; // gxy
        bm[3][c + 1] = gi[0] * inv;
        bm[4][c + 1] = gi[2] * inv; // gyz
        bm[4][c + 2] = gi[1] * inv;
        bm[5][c] = gi[2] * inv; // gxz
        bm[5][c + 2] = gi[0] * inv;
    }
    let mut ke = [[0.0; 12]; 12];
    btdb(&mut ke, &bm, d, det.abs() / 6.0);
    Some((ke, CentroidB::Tet(Box::new(bm))))
}

/// ke += w * B^T D B for an SxN strain-displacement matrix (S = 3 or 6).
fn btdb<const N: usize, const S: usize>(
    ke: &mut [[f64; N]; N],
    b: &[[f64; N]; S],
    d: &[[f64; S]; S],
    w: f64,
) {
    let mut db = [[0.0; N]; S];
    for r in 0..S {
        for c in 0..N {
            db[r][c] = (0..S).map(|k| d[r][k] * b[k][c]).sum();
        }
    }
    for a in 0..N {
        for c in a..N {
            let v = w * (0..S).map(|r| b[r][a] * db[r][c]).sum::<f64>();
            ke[a][c] += v;
            if a != c {
                ke[c][a] += v;
            }
        }
    }
}

/// Strain vector from a strain-displacement matrix and nodal displacements.
fn strain<const N: usize, const S: usize>(b: &[[f64; N]; S], ue: &[f64]) -> [f64; S] {
    let mut e = [0.0; S];
    for r in 0..S {
        for (c, &u) in ue.iter().enumerate() {
            e[r] += b[r][c] * u;
        }
    }
    e
}

/// Stress components (xx, yy, zz, xy, yz, xz) from a plane B/D pair;
/// out-of-plane components are zero by the plane-stress assumption.
pub(crate) fn plane_stress<const N: usize>(
    b: &[[f64; N]; 3],
    d: &[[f64; 3]; 3],
    ue: &[f64],
) -> [f64; 6] {
    let e = strain(b, ue);
    [
        d[0][0] * e[0] + d[0][1] * e[1],
        d[1][0] * e[0] + d[1][1] * e[1],
        0.0,
        d[2][2] * e[2],
        0.0,
        0.0,
    ]
}

/// Stress components (xx, yy, zz, xy, yz, xz) from a solid B/D pair.
pub(crate) fn solid_stress<const N: usize>(
    b: &[[f64; N]; 6],
    d: &[[f64; 6]; 6],
    ue: &[f64],
) -> [f64; 6] {
    let e = strain(b, ue);
    let mut s = [0.0; 6];
    for r in 0..6 {
        s[r] = (0..6).map(|k| d[r][k] * e[k]).sum();
    }
    s
}

/// Von Mises from full 3D components (xx, yy, zz, xy, yz, xz). Reduces to
/// the plane-stress formula when only in-plane components are non-zero.
pub(crate) fn von_mises(s: &[f64; 6]) -> f64 {
    (0.5 * ((s[0] - s[1]) * (s[0] - s[1])
        + (s[1] - s[2]) * (s[1] - s[2])
        + (s[2] - s[0]) * (s[2] - s[0]))
        + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
        .sqrt()
}

/// Centroid stress components for an element.
pub(crate) fn centroid_stress(b: &CentroidB, d: &Constitutive, ue: &[f64]) -> [f64; 6] {
    match (b, d) {
        (CentroidB::Tri(bm), Constitutive::Plane(dp)) => plane_stress(bm, dp, ue),
        (CentroidB::Quad { center, .. }, Constitutive::Plane(dp)) => plane_stress(center, dp, ue),
        (CentroidB::Tet(bm), Constitutive::Solid(ds)) => solid_stress(bm, ds, ue),
        (CentroidB::Hex { center, .. }, Constitutive::Solid(ds)) => solid_stress(center, ds, ue),
        (CentroidB::Tri(_) | CentroidB::Quad { .. }, Constitutive::Solid(_))
        | (CentroidB::Tet(_) | CentroidB::Hex { .. }, Constitutive::Plane(_)) => {
            unreachable!("element kind and constitutive law out of sync")
        }
    }
}

/// Stress components at each element node (extrapolated from the Gauss
/// points for quads/hexes, constant for triangles/tets), parallel to the
/// element's connectivity.
pub(crate) fn nodal_stresses(b: &CentroidB, d: &Constitutive, ue: &[f64]) -> Vec<[f64; 6]> {
    match (b, d) {
        (CentroidB::Tri(bm), Constitutive::Plane(dp)) => vec![plane_stress(bm, dp, ue); 3],
        (CentroidB::Quad { gauss, .. }, Constitutive::Plane(dp)) => {
            let mut gp = [[0.0; 6]; 4];
            for (j, bm) in gauss.iter().enumerate() {
                gp[j] = plane_stress(bm, dp, ue);
            }
            extrapolate_to_corners(&gp).into()
        }
        (CentroidB::Tet(bm), Constitutive::Solid(ds)) => vec![solid_stress(bm, ds, ue); 4],
        (CentroidB::Hex { gauss, .. }, Constitutive::Solid(ds)) => {
            let mut gp = [[0.0; 6]; 8];
            for (j, bm) in gauss.iter().enumerate() {
                gp[j] = solid_stress(bm, ds, ue);
            }
            extrapolate_to_corners3(&gp).into()
        }
        (CentroidB::Tri(_) | CentroidB::Quad { .. }, Constitutive::Solid(_))
        | (CentroidB::Tet(_) | CentroidB::Hex { .. }, Constitutive::Plane(_)) => {
            unreachable!("element kind and constitutive law out of sync")
        }
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

/// 2x2x2 Gauss-point sign triples in stiffness-loop order (xi-major,
/// then eta, then zeta), matching the `gauss` slot ordering of
/// [`CentroidB::Hex`].
const GP3_SIGNS: [(f64, f64, f64); 8] = [
    (-1.0, -1.0, -1.0),
    (-1.0, -1.0, 1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, 1.0, 1.0),
    (1.0, -1.0, -1.0),
    (1.0, -1.0, 1.0),
    (1.0, 1.0, -1.0),
    (1.0, 1.0, 1.0),
];

/// Extrapolate 2x2x2 Gauss-point values to the 8 corner nodes: the 3D
/// analog of the quad scheme. The 8 GP values fix a trilinear
/// interpolant in GP coordinates (the GPs sit at +-1/sqrt(3), so the
/// corners map to +-sqrt(3) there); extrapolation is that interpolant
/// evaluated at the corners, i.e. W[i][j] = N_j(sqrt(3)*xi_i, ...). No
/// matrix inversion needed and exact for trilinear fields; rows sum to 1.
fn extrapolate_to_corners3<const M: usize>(gp: &[[f64; M]; 8]) -> [[f64; M]; 8] {
    let sqrt3 = (3.0f64).sqrt();
    let h = |s: f64| 0.5 * (1.0 + s * sqrt3);
    let mut out = [[0.0; M]; 8];
    for (i, &(cx, cy, cz)) in HEX_CORNERS.iter().enumerate() {
        for (j, &(gx, gy, gz)) in GP3_SIGNS.iter().enumerate() {
            let w = h(cx * gx) * h(cy * gy) * h(cz * gz);
            for k in 0..M {
                out[i][k] += w * gp[j][k];
            }
        }
    }
    out
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

    fn unit_cube() -> [[f64; 3]; 8] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
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

    #[test]
    fn c3d8_stiffness_is_symmetric_and_kills_rigid_translations() {
        let d = d3_matrix(1.0, 0.25);
        let (ke, _) = c3d8_stiffness(&unit_cube(), &d).unwrap();
        for (a, row) in ke.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                assert!((v - ke[c][a]).abs() < 1e-12, "K[{a},{c}] != K[{c},{a}]");
            }
        }
        for t in 0..3 {
            let u: Vec<f64> = (0..8)
                .flat_map(|_| {
                    let mut v = [0.0; 3];
                    v[t] = 1.0;
                    v
                })
                .collect();
            for (a, row) in ke.iter().enumerate() {
                let ku: f64 = row.iter().zip(&u).map(|(&k, &ui)| k * ui).sum();
                assert!(ku.abs() < 1e-10, "row {a} not invariant in axis {t}");
            }
        }
    }

    #[test]
    fn c3d4_stiffness_kills_rigid_translations() {
        let d = d3_matrix(1.0, 0.25);
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let (ke, _) = c3d4_stiffness(&coords, &d).unwrap();
        for t in 0..3 {
            let u: Vec<f64> = (0..4)
                .flat_map(|_| {
                    let mut v = [0.0; 3];
                    v[t] = 1.0;
                    v
                })
                .collect();
            for (a, row) in ke.iter().enumerate() {
                let ku: f64 = row.iter().zip(&u).map(|(&k, &ui)| k * ui).sum();
                assert!(ku.abs() < 1e-12, "row {a} not invariant in axis {t}");
            }
        }
    }

    #[test]
    fn hex_extrapolation_is_exact_for_trilinear_fields() {
        // Patch test on a full trilinear field (all 8 monomials).
        let f = |xi: f64, eta: f64, zeta: f64| {
            2.0 + 3.0 * xi - 5.0 * eta + 1.5 * zeta + 0.5 * xi * eta + 0.3 * eta * zeta
                - 0.2 * xi * zeta
                + 1.1 * xi * eta * zeta
        };
        let g = 1.0f64 / (3.0f64).sqrt();
        let mut gp = [[0.0f64; 1]; 8];
        for (j, &(gx, gy, gz)) in GP3_SIGNS.iter().enumerate() {
            gp[j][0] = f(g * gx, g * gy, g * gz);
        }
        let corners = extrapolate_to_corners3(&gp);
        for (i, &(cx, cy, cz)) in HEX_CORNERS.iter().enumerate() {
            assert!(
                (corners[i][0] - f(cx, cy, cz)).abs() < 1e-12,
                "corner {i}: {} != {}",
                corners[i][0],
                f(cx, cy, cz)
            );
        }
        // Constant field: every row of the extrapolation matrix sums to 1.
        let c = extrapolate_to_corners3(&[[7.0; 1]; 8]);
        assert!(c.iter().all(|&v| (v[0] - 7.0).abs() < 1e-12));
    }

    #[test]
    fn inverted_hex_is_rejected() {
        let d = d3_matrix(1.0, 0.25);
        let mut coords = unit_cube();
        coords.swap(1, 2); // clockwise bottom face inverts J
        assert!(c3d8_stiffness(&coords, &d).is_none());
    }

    #[test]
    fn degenerate_tet_is_rejected() {
        let d = d3_matrix(1.0, 0.25);
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0], // coplanar: zero volume
        ];
        assert!(c3d4_stiffness(&coords, &d).is_none());
    }
}
