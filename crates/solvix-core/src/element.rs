//! Element stiffness: CPS3/CPS4/CPS8 plane stress, CPE4/CPE8 plane
//! strain (2 dofs/node), C3D8/C3D4 linear solids and B31 Hermite
//! beams (3 dofs/node), plus stress recovery helpers.

use std::sync::OnceLock;

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

/// 3x3 plane-strain constitutive matrix (in-plane block; eps_zz = 0).
/// Caller checks nu < 0.5. Entry [0][1] equals lambda, which the stress
/// recovery uses for the non-zero sigma_zz = lambda*(eps_xx + eps_yy).
pub(crate) fn dpe_matrix(e: f64, nu: f64) -> [[f64; 3]; 3] {
    let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    [
        [c * (1.0 - nu), c * nu, 0.0],
        [c * nu, c * (1.0 - nu), 0.0],
        [0.0, 0.0, c * (1.0 - 2.0 * nu) / 2.0],
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
/// plane elements get the plane-stress/plane-strain 3x3, solids the 3D
/// 6x6, beams carry their section constants (E from the material, A and
/// I from `*BEAM SECTION`).
pub(crate) enum Constitutive {
    Plane([[f64; 3]; 3]),
    PlaneStrain([[f64; 3]; 3]),
    Solid(Box<[[f64; 6]; 6]>),
    Beam { e: f64, a: f64, i: f64 },
}

/// Strain-displacement matrices kept per element for stress recovery:
/// the centroid matrix, plus the Gauss-point matrices for quads and hexes
/// (ordering matches the stiffness loops; see GP_SIGNS / GP9_SIGNS).
pub(crate) enum CentroidB {
    Tri([[f64; 6]; 3]),
    Quad {
        center: [[f64; 8]; 3],
        gauss: Box<[[[f64; 8]; 3]; 4]>,
    },
    Quad8 {
        center: [[f64; 16]; 3],
        gauss: Box<[[[f64; 16]; 3]; 9]>,
    },
    /// B31: (axial-strain row, curvature row) in GLOBAL dofs, at
    /// mid-span (centroid stress) and at both end nodes (nodal
    /// recovery — Hermite strain varies along the axis, unlike quads).
    Beam {
        center: Box<([f64; 6], [f64; 6])>,
        ends: [([f64; 6], [f64; 6]); 2],
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

/// CPS8 node natural coordinates, Abaqus serendipity ordering: corners
/// 1-4 CCW, then midsides 5(1-2), 6(2-3), 7(3-4), 8(4-1).
const QUAD8_NODES: [(f64, f64); 8] = [
    (-1.0, -1.0),
    (1.0, -1.0),
    (1.0, 1.0),
    (-1.0, 1.0),
    (0.0, -1.0),
    (1.0, 0.0),
    (0.0, 1.0),
    (-1.0, 0.0),
];

/// 8-node serendipity shape functions at (xi, eta), order matching
/// [`QUAD8_NODES`]. Span: {1, xi, eta, xi^2, xi*eta, eta^2, xi^2*eta,
/// xi*eta^2}.
pub(crate) fn quad8_shapes(xi: f64, eta: f64) -> [f64; 8] {
    let mut n = [0.0; 8];
    for i in 0..4 {
        let (xi_i, eta_i) = QUAD8_NODES[i];
        n[i] = 0.25 * (1.0 + xi * xi_i) * (1.0 + eta * eta_i) * (xi * xi_i + eta * eta_i - 1.0);
    }
    n[4] = 0.5 * (1.0 - xi * xi) * (1.0 - eta);
    n[5] = 0.5 * (1.0 + xi) * (1.0 - eta * eta);
    n[6] = 0.5 * (1.0 - xi * xi) * (1.0 + eta);
    n[7] = 0.5 * (1.0 - xi) * (1.0 - eta * eta);
    n
}

/// Shape function derivatives at (xi, eta): (dN/dxi, dN/deta), each
/// aligned with [`QUAD8_NODES`].
fn quad8_shape_derivs(xi: f64, eta: f64) -> ([f64; 8], [f64; 8]) {
    let mut dn_dxi = [0.0; 8];
    let mut dn_deta = [0.0; 8];
    for i in 0..4 {
        let (xi_i, eta_i) = QUAD8_NODES[i];
        dn_dxi[i] = 0.25 * xi_i * (1.0 + eta * eta_i) * (2.0 * xi * xi_i + eta * eta_i);
        dn_deta[i] = 0.25 * eta_i * (1.0 + xi * xi_i) * (xi * xi_i + 2.0 * eta * eta_i);
    }
    dn_dxi[4] = -xi * (1.0 - eta);
    dn_deta[4] = -0.5 * (1.0 - xi * xi);
    dn_dxi[5] = 0.5 * (1.0 - eta * eta);
    dn_deta[5] = -eta * (1.0 + xi);
    dn_dxi[6] = -xi * (1.0 + eta);
    dn_deta[6] = 0.5 * (1.0 - xi * xi);
    dn_dxi[7] = -0.5 * (1.0 - eta * eta);
    dn_deta[7] = -eta * (1.0 - xi);
    (dn_dxi, dn_deta)
}

/// B matrix (3x16) at a natural point, or `None` if det(J) is not
/// positive there. Also the det(J) source for the CPS8 consistent mass.
pub(crate) fn quad8_b(coords: &[[f64; 2]; 8], xi: f64, eta: f64) -> Option<([[f64; 16]; 3], f64)> {
    let (dn_dxi, dn_deta) = quad8_shape_derivs(xi, eta);
    let mut j = [[0.0; 2]; 2];
    for i in 0..8 {
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
    let mut bm = [[0.0; 16]; 3];
    for i in 0..8 {
        let dnx = inv * (j[1][1] * dn_dxi[i] - j[0][1] * dn_deta[i]);
        let dny = inv * (j[0][0] * dn_deta[i] - j[1][0] * dn_dxi[i]);
        bm[0][2 * i] = dnx;
        bm[1][2 * i + 1] = dny;
        bm[2][2 * i] = dny;
        bm[2][2 * i + 1] = dnx;
    }
    Some((bm, detj))
}

/// 3x3 Gauss abscissae in stiffness-loop order (xi-major, then eta),
/// matching the `gauss` slot ordering of [`CentroidB::Quad8`].
const GP9_SIGNS: [(f64, f64); 9] = [
    (-1.0, -1.0),
    (-1.0, 0.0),
    (-1.0, 1.0),
    (0.0, -1.0),
    (0.0, 0.0),
    (0.0, 1.0),
    (1.0, -1.0),
    (1.0, 0.0),
    (1.0, 1.0),
];

pub(crate) const G3: f64 = 0.7745_9666_9241_4834; // sqrt(3/5)
pub(crate) const W3: [f64; 3] = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];

/// CPS8/CPE8 stiffness with 3x3 Gauss (shared by both plane families —
/// the geometry is identical, only D differs). Returns `None` if
/// det(J) <= 0 at any Gauss point (distorted/inverted element).
pub(crate) fn cps8_stiffness(
    coords: &[[f64; 2]; 8],
    d: &[[f64; 3]; 3],
    t: f64,
) -> Option<([[f64; 16]; 16], CentroidB)> {
    let mut ke = [[0.0; 16]; 16];
    let mut gauss = [[[0.0; 16]; 3]; 9];
    let mut k = 0;
    for (i, &xi) in [-G3, 0.0, G3].iter().enumerate() {
        for (j, &eta) in [-G3, 0.0, G3].iter().enumerate() {
            let (bm, detj) = quad8_b(coords, xi, eta)?;
            btdb(&mut ke, &bm, d, t * detj * W3[i] * W3[j]);
            gauss[k] = bm;
            k += 1;
        }
    }
    let (center, _) = quad8_b(coords, 0.0, 0.0)?; // centroid is interior if GPs are
    Some((
        ke,
        CentroidB::Quad8 {
            center,
            gauss: Box::new(gauss),
        },
    ))
}

/// In-place Gaussian elimination with partial pivoting: solve
/// `a * x = b` for square `a` (N x N) and row-major RHS `b` (N rows of
/// M columns), std only. Used once per process for the CPS8
/// extrapolation matrix.
fn gauss_solve<const N: usize, const M: usize>(a: &mut [[f64; N]; N], b: &mut [[f64; M]; N]) {
    for col in 0..N {
        let piv = (col..N)
            .max_by(|&r, &s| a[r][col].abs().total_cmp(&a[s][col].abs()))
            .unwrap_or(col);
        if piv != col {
            a.swap(col, piv);
            b.swap(col, piv);
        }
        let d = a[col][col];
        let piv_row = b[col];
        for r in col + 1..N {
            let f = a[r][col] / d;
            if f == 0.0 {
                continue;
            }
            let piv_a = a[col];
            for (c, av) in a[r].iter_mut().enumerate().take(N).skip(col) {
                *av -= f * piv_a[c];
            }
            for (br, bcol) in b[r].iter_mut().zip(&piv_row) {
                *br -= f * bcol;
            }
        }
    }
    // Back-substitution.
    for r in (0..N).rev() {
        for c in r + 1..N {
            let f = a[r][c];
            let row = b[c];
            for (br, bcol) in b[r].iter_mut().zip(&row) {
                *br -= f * bcol;
            }
        }
        let d = a[r][r];
        for br in b[r].iter_mut() {
            *br /= d;
        }
    }
}

/// Extrapolation matrix W (8x9) mapping the 9 Gauss-point stress values
/// to the 8 nodes, computed once: the 8 serendipity shape functions are
/// least-squares fitted to the 9 GP values (normal equations
/// (A^T A) W = A^T with A[k][j] = N_j(GP_k), both fixed in natural
/// coordinates). Because the serendipity span contains all quadratics,
/// the fit is exact for stress fields the element represents exactly —
/// the patch-test guarantee. Nodal values equal the fitted coefficients
/// (shape functions are cardinal at the nodes).
fn cps8_extrapolation() -> &'static [[f64; 9]; 8] {
    static W: OnceLock<[[f64; 9]; 8]> = OnceLock::new();
    W.get_or_init(|| {
        let mut a = [[0.0f64; 8]; 9]; // A[k][j] = N_j(GP_k)
        for (k, &(sx, sy)) in GP9_SIGNS.iter().enumerate() {
            a[k] = quad8_shapes(sx * G3, sy * G3); // sign 0 lands on 0
        }
        let mut ata = [[0.0f64; 8]; 8];
        let mut at = [[0.0f64; 9]; 8]; // A^T, RHS of the normal equations
        for r in 0..8 {
            for c in 0..8 {
                ata[r][c] = (0..9).map(|k| a[k][r] * a[k][c]).sum();
            }
            for k in 0..9 {
                at[r][k] = a[k][r];
            }
        }
        gauss_solve(&mut ata, &mut at);
        at
    })
}

/// Extrapolate 3x3 Gauss-point values to the 8 nodes (see
/// [`cps8_extrapolation`]).
fn extrapolate_to_nodes8<const M: usize>(gp: &[[f64; M]; 9]) -> [[f64; M]; 8] {
    let w = cps8_extrapolation();
    let mut out = [[0.0; M]; 8];
    for i in 0..8 {
        for (j, wv) in w[i].iter().enumerate() {
            for (k, o) in out[i].iter_mut().enumerate() {
                *o += wv * gp[j][k];
            }
        }
    }
    out
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

/// Rotate a 6x6 element matrix from local to global dofs:
/// `K_g = T^T K_l T`, T block-diagonal with per-node
/// [[c, s, 0], [-s, c, 0], [0, 0, 1]] (all 2x2 node-block pairs,
/// including the coupling blocks).
fn rotate_beam(kl: &[[f64; 6]; 6], c: f64, s: f64) -> [[f64; 6]; 6] {
    let r = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
    let mut out = [[0.0; 6]; 6];
    for bi in 0..2 {
        for bj in 0..2 {
            for i in 0..3 {
                for j in 0..3 {
                    out[3 * bi + i][3 * bj + j] = (0..3)
                        .map(|m| {
                            (0..3)
                                .map(|n| r[m][i] * kl[3 * bi + m][3 * bj + n] * r[n][j])
                                .sum::<f64>()
                        })
                        .sum::<f64>();
                }
            }
        }
    }
    out
}

/// B31 stiffness: 2-node Euler-Bernoulli beam with cubic (Hermite)
/// transverse interpolation, in-plane dofs [u1, v1, rotz1, u2, v2,
/// rotz2]. Closed forms of EA*int(B_a^T B_a) + EI*int(B_b^T B_b) (the
/// integral is exact for the Hermite basis). The local x-axis runs
/// node 1 -> node 2.
///
/// ponytail: shear-free Euler-Bernoulli — exact for point-load linear
/// static; Timoshenko (shear-flexible) upgrade needed once squat beams
/// matter; add GA_s/L terms + reduced integration here.
/// Returns `None` for a zero-length beam.
pub(crate) fn b31_stiffness(
    coords: &[[f64; 2]; 2],
    e: f64,
    area: f64,
    inertia: f64,
) -> Option<([[f64; 6]; 6], CentroidB)> {
    let (dx, dy) = (coords[1][0] - coords[0][0], coords[1][1] - coords[0][1]);
    let l = (dx * dx + dy * dy).sqrt();
    if l == 0.0 {
        return None;
    }
    let (c, s) = (dx / l, dy / l);
    let ea_l = e * area / l;
    let w = e * inertia / (l * l * l);
    let mut kl = [[0.0; 6]; 6];
    kl[0][0] = ea_l;
    kl[0][3] = -ea_l;
    kl[3][0] = -ea_l;
    kl[3][3] = ea_l;
    let her = [
        [12.0, 6.0 * l, -12.0, 6.0 * l],
        [6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l],
        [-12.0, -6.0 * l, 12.0, -6.0 * l],
        [6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l],
    ];
    // Bending dofs are (v1, th1, v2, th2) = dofs 1, 2, 4, 5 — not a
    // contiguous block (dof 0/3 are axial).
    const BDOF: [usize; 4] = [1, 2, 4, 5];
    for (i, &di) in BDOF.iter().enumerate() {
        for (j, &dj) in BDOF.iter().enumerate() {
            kl[di][dj] = w * her[i][j];
        }
    }

    // Strain-displacement rows at axial position x (local dofs),
    // rotated into global dofs: row_global = row_local * T (T
    // block-diagonal), so eps/chi come straight off global ue.
    let rows = |x: f64| -> ([f64; 6], [f64; 6]) {
        let ba_l = [-1.0 / l, 0.0, 0.0, 1.0 / l, 0.0, 0.0];
        let t = x / l;
        // Curvature row on the four bending dofs (v1, th1, v2, th2) =
        // local dofs 1, 2, 4, 5.
        let bb4 = [
            (-6.0 + 12.0 * t) / (l * l),
            (-4.0 + 6.0 * t) / l,
            (6.0 - 12.0 * t) / (l * l),
            (-2.0 + 6.0 * t) / l,
        ];
        let mut bb_l = [0.0; 6];
        for (loc, &v) in [1usize, 2, 4, 5].iter().zip(&bb4) {
            bb_l[*loc] = v;
        }
        let r = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]];
        let mut ba = [0.0; 6];
        let mut bb = [0.0; 6];
        for blk in 0..2 {
            for j in 0..3 {
                for m in 0..3 {
                    ba[3 * blk + j] += ba_l[3 * blk + m] * r[m][j];
                    bb[3 * blk + j] += bb_l[3 * blk + m] * r[m][j];
                }
            }
        }
        (ba, bb)
    };
    Some((
        rotate_beam(&kl, c, s),
        CentroidB::Beam {
            center: Box::new(rows(l / 2.0)),
            ends: [rows(0.0), rows(l)],
        },
    ))
}

/// B31 consistent mass: translational only (axial + Hermite bending,
/// no rotary inertia — the Euler-Bernoulli-consistent choice), rotated
/// into global dofs. `None` for a zero-length beam.
pub(crate) fn b31_mass(coords: &[[f64; 2]; 2], rho: f64, area: f64) -> Option<[[f64; 6]; 6]> {
    let (dx, dy) = (coords[1][0] - coords[0][0], coords[1][1] - coords[0][1]);
    let l = (dx * dx + dy * dy).sqrt();
    if l == 0.0 {
        return None;
    }
    let (c, s) = (dx / l, dy / l);
    let m = rho * area * l;
    let mut ml = [[0.0; 6]; 6];
    ml[0][0] = m / 3.0;
    ml[0][3] = m / 6.0;
    ml[3][0] = m / 6.0;
    ml[3][3] = m / 3.0;
    let w = m / 420.0;
    let her = [
        [156.0, 22.0 * l, 54.0, -13.0 * l],
        [22.0 * l, 4.0 * l * l, 13.0 * l, -3.0 * l * l],
        [54.0, 13.0 * l, 156.0, -22.0 * l],
        [-13.0 * l, -3.0 * l * l, -22.0 * l, 4.0 * l * l],
    ];
    const BDOF: [usize; 4] = [1, 2, 4, 5];
    for (i, &di) in BDOF.iter().enumerate() {
        for (j, &dj) in BDOF.iter().enumerate() {
            ml[di][dj] = w * her[i][j];
        }
    }
    Some(rotate_beam(&ml, c, s))
}

/// B31 stress: uniaxial sigma = E*(eps_axial ± c*curvature), reported
/// as the larger-magnitude fiber value (max at +-c).
///
/// ponytail: fiber distance c = sqrt(3I/A) assumes an equivalent
/// rectangular section; an explicit c in `*BEAM SECTION` (v0.4) for
/// non-rectangular shapes.
fn beam_stress(
    e: f64,
    area: f64,
    inertia: f64,
    ba: &[f64; 6],
    bb: &[f64; 6],
    ue: &[f64],
) -> [f64; 6] {
    let eps = (0..6).map(|c| ba[c] * ue[c]).sum::<f64>();
    let chi = (0..6).map(|c| bb[c] * ue[c]).sum::<f64>();
    let c = (3.0 * inertia / area).sqrt();
    let s = (e * (eps + c * chi)).abs().max((e * (eps - c * chi)).abs());
    [s, 0.0, 0.0, 0.0, 0.0, 0.0]
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

/// Stress components from a plane-strain B/D pair: same in-plane
/// machinery, but sigma_zz = lambda*(eps_xx + eps_yy) with lambda =
/// d[0][1] (the plane-strain D carries it directly).
pub(crate) fn plane_strain_stress<const N: usize>(
    b: &[[f64; N]; 3],
    d: &[[f64; 3]; 3],
    ue: &[f64],
) -> [f64; 6] {
    let e = strain(b, ue);
    [
        d[0][0] * e[0] + d[0][1] * e[1],
        d[1][0] * e[0] + d[1][1] * e[1],
        d[0][1] * (e[0] + e[1]),
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
    use Constitutive as C;
    match (b, d) {
        (CentroidB::Tri(bm), C::Plane(dp)) => plane_stress(bm, dp, ue),
        (CentroidB::Quad { center, .. }, C::Plane(dp)) => plane_stress(center, dp, ue),
        (CentroidB::Quad { center, .. }, C::PlaneStrain(dp)) => plane_strain_stress(center, dp, ue),
        (CentroidB::Quad8 { center, .. }, C::Plane(dp)) => plane_stress(center, dp, ue),
        (CentroidB::Quad8 { center, .. }, C::PlaneStrain(dp)) => {
            plane_strain_stress(center, dp, ue)
        }
        (CentroidB::Beam { center, .. }, C::Beam { e, a, i }) => {
            beam_stress(*e, *a, *i, &center.0, &center.1, ue)
        }
        (CentroidB::Tet(bm), C::Solid(ds)) => solid_stress(bm, ds, ue),
        (CentroidB::Hex { center, .. }, C::Solid(ds)) => solid_stress(center, ds, ue),
        _ => unreachable!("element kind and constitutive law out of sync"),
    }
}

/// Stress components at each element node (extrapolated from the Gauss
/// points for quads/hexes, constant for triangles/tets, end-section
/// values for beams), parallel to the element's connectivity.
pub(crate) fn nodal_stresses(b: &CentroidB, d: &Constitutive, ue: &[f64]) -> Vec<[f64; 6]> {
    use Constitutive as C;
    match (b, d) {
        (CentroidB::Tri(bm), C::Plane(dp)) => vec![plane_stress(bm, dp, ue); 3],
        (CentroidB::Quad { gauss, .. }, C::Plane(dp)) => {
            let mut gp = [[0.0; 6]; 4];
            for (j, bm) in gauss.iter().enumerate() {
                gp[j] = plane_stress(bm, dp, ue);
            }
            extrapolate_to_corners(&gp).into()
        }
        (CentroidB::Quad { gauss, .. }, C::PlaneStrain(dp)) => {
            let mut gp = [[0.0; 6]; 4];
            for (j, bm) in gauss.iter().enumerate() {
                gp[j] = plane_strain_stress(bm, dp, ue);
            }
            extrapolate_to_corners(&gp).into()
        }
        (CentroidB::Quad8 { gauss, .. }, C::Plane(dp)) => {
            let mut gp = [[0.0; 6]; 9];
            for (j, bm) in gauss.iter().enumerate() {
                gp[j] = plane_stress(bm, dp, ue);
            }
            extrapolate_to_nodes8(&gp).into()
        }
        (CentroidB::Quad8 { gauss, .. }, C::PlaneStrain(dp)) => {
            let mut gp = [[0.0; 6]; 9];
            for (j, bm) in gauss.iter().enumerate() {
                gp[j] = plane_strain_stress(bm, dp, ue);
            }
            extrapolate_to_nodes8(&gp).into()
        }
        (CentroidB::Beam { ends, .. }, C::Beam { e, a, i }) => ends
            .iter()
            .map(|(ba, bb)| beam_stress(*e, *a, *i, ba, bb, ue))
            .collect(),
        (CentroidB::Tet(bm), C::Solid(ds)) => vec![solid_stress(bm, ds, ue); 4],
        (CentroidB::Hex { gauss, .. }, C::Solid(ds)) => {
            let mut gp = [[0.0; 6]; 8];
            for (j, bm) in gauss.iter().enumerate() {
                gp[j] = solid_stress(bm, ds, ue);
            }
            extrapolate_to_corners3(&gp).into()
        }
        _ => unreachable!("element kind and constitutive law out of sync"),
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

    fn unit_square8() -> [[f64; 2]; 8] {
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.0],
            [1.0, 0.5],
            [0.5, 1.0],
            [0.0, 0.5],
        ]
    }

    #[test]
    fn quad8_shapes_partition_unity_and_are_cardinal_at_nodes() {
        for &(xi, eta) in &[(0.0, 0.0), (0.3, -0.7), (-0.9, 0.2), (G3, -G3)] {
            let n = quad8_shapes(xi, eta);
            assert!(
                (n.iter().sum::<f64>() - 1.0).abs() < 1e-12,
                "at ({xi},{eta})"
            );
        }
        for (k, &(nx, ny)) in QUAD8_NODES.iter().enumerate() {
            let n = quad8_shapes(nx, ny);
            for (j, &nj) in n.iter().enumerate() {
                let want = if j == k { 1.0 } else { 0.0 };
                assert!((nj - want).abs() < 1e-12, "node {k}, shape {j}: {nj}");
            }
        }
    }

    #[test]
    fn cps8_stiffness_is_symmetric_and_kills_rigid_translation() {
        let d = d_matrix(1.0, 0.25);
        let (ke, _) = cps8_stiffness(&unit_square8(), &d, 1.0).unwrap();
        for (a, row) in ke.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                assert!((v - ke[c][a]).abs() < 1e-12, "K[{a},{c}] != K[{c},{a}]");
            }
        }
        let mut u = [0.0; 16];
        for i in 0..8 {
            u[2 * i] = 0.7;
            u[2 * i + 1] = -1.3;
        }
        for (a, row) in ke.iter().enumerate() {
            let ku: f64 = row.iter().zip(u).map(|(k, ui)| k * ui).sum();
            assert!(ku.abs() < 1e-12, "row {a} not translation-invariant");
        }
    }

    #[test]
    fn inverted_quad8_is_rejected() {
        let d = d_matrix(1.0, 0.25);
        let mut coords = unit_square8();
        coords.swap(1, 3); // clockwise corners invert J
        assert!(cps8_stiffness(&coords, &d, 1.0).is_none());
    }

    #[test]
    fn cps8_extrapolation_is_exact_through_the_serendipity_span() {
        // Constant + full quadratic: the serendipity span, so the
        // least-squares fit recovers nodal values exactly.
        let f = |xi: f64, eta: f64| 2.0 + 3.0 * xi - 5.0 * eta + 1.5 * xi * xi - 0.8 * eta * eta;
        let mut gp = [[0.0f64; 1]; 9];
        for (k, &(sx, sy)) in GP9_SIGNS.iter().enumerate() {
            gp[k][0] = f(sx * G3, sy * G3);
        }
        let nodes = extrapolate_to_nodes8(&gp);
        for (i, &(nx, ny)) in QUAD8_NODES.iter().enumerate() {
            assert!(
                (nodes[i][0] - f(nx, ny)).abs() < 1e-10,
                "node {i}: {} != {}",
                nodes[i][0],
                f(nx, ny)
            );
        }
        // Constant field: every row of the extrapolation matrix sums to 1.
        let c = extrapolate_to_nodes8(&[[7.0; 1]; 9]);
        assert!(c.iter().all(|&v| (v[0] - 7.0).abs() < 1e-12));
    }

    #[test]
    fn gauss_solve_matches_known_2x2() {
        let mut a = [[2.0, 1.0], [1.0, 3.0]];
        let mut b = [[3.0, 1.0], [5.0, 0.0]];
        gauss_solve(&mut a, &mut b);
        // [2 1;1 3]^-1 = (1/5)[3 -1; -1 2]: x1 = (0.8, 1.4), x2 = (0.6, -0.2)
        assert!((b[0][0] - 0.8).abs() < 1e-12);
        assert!((b[1][0] - 1.4).abs() < 1e-12);
        assert!((b[0][1] - 0.6).abs() < 1e-12);
        assert!((b[1][1] - (-0.2)).abs() < 1e-12);
    }

    #[test]
    fn b31_stiffness_kills_rigid_translation_and_is_symmetric() {
        let coords = [[0.0, 0.0], [2.0, 1.0]]; // angled beam
        let (ke, _) = b31_stiffness(&coords, 1.0e6, 0.05, 1.0e-4).unwrap();
        for (a, row) in ke.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                assert!((v - ke[c][a]).abs() < 1e-9 * ke[c][a].abs().max(1.0));
            }
        }
        let u = [0.4, -0.7, 0.0, 0.4, -0.7, 0.0];
        for (a, row) in ke.iter().enumerate() {
            let ku: f64 = row.iter().zip(u).map(|(k, ui)| k * ui).sum();
            assert!(ku.abs() < 1e-9, "row {a} not translation-invariant");
        }
    }

    #[test]
    fn b31_stiffness_kills_rigid_rotation_for_vertical_beam() {
        // Beam along +y: a rigid rotation about the origin (ux = -y,
        // uy = x, rotz = 1) must produce zero nodal forces.
        let (ke, _) = b31_stiffness(&[[0.0, 0.0], [0.0, 3.0]], 1.0e6, 0.05, 1.0e-4).unwrap();
        let u = [0.0, 0.0, 1.0, -3.0, 0.0, 1.0];
        for (a, row) in ke.iter().enumerate() {
            let ku: f64 = row.iter().zip(u).map(|(k, ui)| k * ui).sum();
            assert!(ku.abs() < 1e-6, "row {a}: {ku}");
        }
    }

    #[test]
    fn b31_mass_axial_columns_carry_half_each() {
        let m = b31_mass(&[[0.0, 0.0], [2.5, 0.0]], 7800.0, 0.02).unwrap();
        let total = 7800.0 * 0.02 * 2.5;
        // Axial columns: m/3 + m/6 each — together the full mass for a
        // rigid axial velocity. (Transverse columns also carry Hermite
        // coupling into rotz, so only the axial pair has the clean sum.)
        for d in [0, 3] {
            let col: f64 = (0..6).map(|r| m[r][d]).sum();
            assert!((col - total / 2.0).abs() < 1e-9, "dof {d}: {col}");
        }
        for (a, row) in m.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                assert!((v - m[c][a]).abs() < 1e-9 * v.abs().max(1e-9));
            }
        }
    }
}
