/*!
Defines the interpolation schemes used in the numerically-integrated elements
*/

/// Gauss point, includes natural coordinates (2D) and weight
#[derive(Clone, Copy)]
pub struct Gp {
    /// Natural coordinates
    pub coords: [f64; 2],
    /// Integration weight
    pub wt: f64,
}
/// Natural coordinates and weights at the triangle nodes
#[allow(dead_code)]
pub static TRINODES: [Gp; 3] = [
    Gp {
        coords: [0.0, 0.0],
        wt: 1. / 3.,
    },
    Gp {
        coords: [1.0, 0.0],
        wt: 1. / 3.,
    },
    Gp {
        coords: [0.0, 1.0],
        wt: 1. / 3.,
    },
];
/// Integration scheme for 1 point triangle
pub static TRIGP1: [Gp; 1] = [Gp {
    coords: [1. / 3., 1. / 3.],
    wt: 1.0,
}];
/// Integration scheme for 3 point triangle
pub static TRIGP3: [Gp; 3] = [
    Gp {
        coords: [1. / 6., 1. / 6.],
        wt: 1. / 3.,
    },
    Gp {
        coords: [1. / 6., 2. / 3.],
        wt: 1. / 3.,
    },
    Gp {
        coords: [2. / 3., 1. / 6.],
        wt: 1. / 3.,
    },
];
/// Integration scheme for 6 point triangle
#[allow(dead_code)]
pub static TRIGP6: [Gp; 6] = [
    Gp {
        coords: [0.091576213509771, 0.091576213509771],
        wt: 0.109951743655322,
    },
    Gp {
        coords: [0.816847572980459, 0.091576213509771],
        wt: 0.109951743655322,
    },
    Gp {
        coords: [0.091576213509771, 0.816847572980459],
        wt: 0.109951743655322,
    },
    Gp {
        coords: [0.445948490915965, 0.445948490915965],
        wt: 0.223381589678011,
    },
    Gp {
        coords: [0.445948490915965, 0.108103018168070],
        wt: 0.223381589678011,
    },
    Gp {
        coords: [0.108103018168070, 0.445948490915965],
        wt: 0.223381589678011,
    },
];
static ONEOVERSQRT3: f64 = 0.57735026919;
/// Natural coordinates and weights at the quadrilateral nodes
#[allow(dead_code)]
pub static QUADNODES: [Gp; 4] = [
    Gp {
        coords: [-1.0, -1.0],
        wt: 1.0,
    },
    Gp {
        coords: [1.0, -1.0],
        wt: 1.0,
    },
    Gp {
        coords: [1.0, 1.0],
        wt: 1.0,
    },
    Gp {
        coords: [-1.0, 1.0],
        wt: 1.0,
    },
];
/// Integration scheme for 1 point quad
pub static QUADGP1: [Gp; 1] = [Gp {
    coords: [0., 0.],
    wt: 4.0,
}];
/// Integration scheme for 4 point quad
pub static QUADGP4: [Gp; 4] = [
    Gp {
        coords: [ONEOVERSQRT3, ONEOVERSQRT3],
        wt: 1.0,
    },
    Gp {
        coords: [-ONEOVERSQRT3, ONEOVERSQRT3],
        wt: 1.0,
    },
    Gp {
        coords: [ONEOVERSQRT3, -ONEOVERSQRT3],
        wt: 1.0,
    },
    Gp {
        coords: [-ONEOVERSQRT3, -ONEOVERSQRT3],
        wt: 1.0,
    },
];
