pub mod interpolation {
    #[derive(Clone, Copy)]
    pub struct Gp {
        pub coords: [f64; 2],
        pub wt: f64,
    }
    pub static TRINODES: [Gp; 3] = [Gp{coords: [0.0, 0.0], wt: 1./3.}, 
                                    Gp{coords: [1.0, 0.0], wt: 1./3.}, 
                                    Gp{coords: [0.0, 1.0], wt: 1./3.}];
    pub static TRIGP3: [Gp; 3] = [Gp{coords: [1./6., 1./6.], wt: 1./3.}, 
                                  Gp{coords: [1./6., 2./3.], wt: 1./3.}, 
                                  Gp{coords: [2./3., 1./6.], wt: 1./3.}];
    #[allow(dead_code)]
    pub static TRIGP6: [Gp; 6] = [Gp{coords: [0.091576213509771, 0.091576213509771], wt: 0.109951743655322},
                                  Gp{coords: [0.816847572980459, 0.091576213509771], wt: 0.109951743655322},
                                  Gp{coords: [0.091576213509771, 0.816847572980459], wt: 0.109951743655322},
                                  Gp{coords: [0.445948490915965, 0.445948490915965], wt: 0.223381589678011},
                                  Gp{coords: [0.445948490915965, 0.108103018168070], wt: 0.223381589678011},
                                  Gp{coords: [0.108103018168070, 0.445948490915965], wt: 0.223381589678011}];
    static ONEOVERSQRT3: f64 = 0.57735026919;
    pub static QUADNODES: [Gp; 4] = [Gp{coords: [-1.0, -1.0], wt: 1.0}, 
                                     Gp{coords: [ 1.0, -1.0], wt: 1.0}, 
                                     Gp{coords: [ 1.0,  1.0], wt: 1.0},
                                     Gp{coords: [-1.0,  1.0], wt: 1.0}];
    pub static QUADGP4: [Gp; 4] = [Gp{coords: [ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                   Gp{coords: [-ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                   Gp{coords: [ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0},
                                   Gp{coords: [-ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0}];
}   

use interpolation::*;
use na::{DMatrix, Vector3};
use crate::Cplx;
use crate::preprocess::mesh_data::{CollocationPoint, Coords, ElementType, Mesh};

/// Calculate classical and 'hypersingular' Green's function (dg) and its derivative (dh) for the given origin (x), destination (y),
/// normal vector at y (non-normalized), and wavenumber (k)
pub fn get_greens_functions(k: f64, x: &Coords, n_x: &Vector3<f64>, 
                           y: &Coords, n_y: &Vector3<f64>, use_hypersingular: bool) -> (Cplx, Cplx) {
    let e_nx = n_x / n_x.norm();
    let e_ny = n_y / n_y.norm();
    let r = x - y;
    let rdist = r.norm();
    let e_r = r / rdist;
    // g = e^(ikr) / (4 pi r)
    let mut g = Cplx::new(0.0, k*rdist).exp() / (4.0 * std::f64::consts::PI * rdist);
    // h = (1/r - ik) * g * (r dot n_y), where r and n are unit vectors -> (r dot n) is a direction cosine
    let f1 = Cplx::new(1.0 / rdist, -k);
    let rdoty = e_r.dot(&e_ny);
    let mut h = g * f1 * rdoty;
    if use_hypersingular {
        // dg = -(1/r - ik) * g * (r dot n_x), where r and n are unit vectors -> (r dot n) is a direction cosine
        let rdotx = e_r.dot(&e_nx);
        let dg = -g * f1 * rdotx;
        // dh = {(1/r - ik) * (n_y dot n_x) + (k^2 * r - 3 * (1/r - ik)) * (r dot n_y) * (r dot n_x)} * g / r
        let dh = 1.0 / rdist * g * (f1 * e_ny.dot(&e_nx) + (k*k*rdist - 3.0*f1) * rdotx * rdoty);
        // coupling parameter gamma = i/k
        let beta = Cplx::new(0.0, 1.0 / k);
        g = Cplx::new(0.0, 0.0);
        h = Cplx::new(0.0, 0.0);
        g += beta * dg;
        h += beta * dh;
    }
    return (g, h)
}

// Methods for numerically-integrated elements
pub struct NIElement <'a> {
    integration: Vec<Gp>,
    mesh: &'a Mesh,
    pub element_id: usize,
    element_type: ElementType
}
impl NIElement <'_> {
    pub fn new<'a>(meshdata: &'a Mesh, element: usize) -> NIElement <'a> {
        let etype = &meshdata.elements[element].etype;
        let int = match etype.clone() {
            ElementType::Tri => TRIGP3.to_vec(),
            ElementType::Quad => QUADGP4.to_vec(),
            _ => {error!("Invalid numerically integrated element"); Vec::new()}
        };
        NIElement{integration: int, 
                  mesh: &meshdata, 
                  element_id: element,
                  element_type: etype.clone()}
    }
    #[inline]
    pub fn get_num_nodes(&self) -> usize {
        match self.element_type {
            ElementType::Tri => 3,
            ElementType::Quad => 4,
            _ => 0
        }
    }
    #[inline]
    pub fn get_num_cpts(&self) -> usize {
        return self.integration.len();
    }
    /// Get shape functions for this element at the given natural coordinates
    fn shape_functions_at(&self, gp: &Gp) -> Vec<f64> {
        let xi = &gp.coords[0];
        let eta = &gp.coords[1];
        match self.element_type {
            ElementType::Tri => {
                vec![1.0 - *xi - *eta,
                *xi,
                *eta]
            },
            ElementType::Quad => {
                vec![0.25*(1.-*xi)*(1.-*eta),
                     0.25*(1.+*xi)*(1.-*eta),
                     0.25*(1.+*xi)*(1.+*eta),
                     0.25*(1.-*xi)*(1.+*eta)]
            },
            _ => vec![0.0]
        }
    }
    /// Get shape function derivatives for this element at the given natural coordinates
    fn shape_derivatives_at(&self, _gp: &Gp) -> DMatrix<f64> {
        match self.element_type {
            ElementType::Tri => {
                DMatrix::from_row_slice(3, 2,
                    &[-1.0, -1.0,
                            1.0, 0.0,
                            0.0, 1.0])
            },
            ElementType::Quad => {
                DMatrix::from_row_slice(4, 2,
                    &[-0.25, -0.25,
                            0.25, -0.25,
                            0.25, 0.25,
                           -0.25, 0.25])
            },
            _ => DMatrix::from_element(1, 1, 0.0)
        }
    }
    /// Get physical coordinates at the given natural coordinates, using shape functions
    fn coordinates_at(&self, gp: &Gp) -> Coords {
        let n = self.shape_functions_at(gp);
        let mut x = Coords::from_element(0.0);
        let element = &self.mesh.elements[self.element_id];
        for ni in 0..n.len() {
            let node_index = &element.node_ids[ni];
            let icoord = &self.mesh.nodes[*node_index].coords;
            x += n[ni] * icoord;
        }
        return x;
    }
    /// Get natural coordinates corresponding to node index
    pub fn natural_coordinates_at_node(&self, i: usize) -> Gp {
        match self.element_type {
            ElementType::Tri => {
                TRINODES[i]
            },
            ElementType::Quad => {
                QUADNODES[i]
            },
            _ => Gp {coords: [0.0, 0.0], wt: 0.0}
        }
    }
    /// Get normal vector (non-normalized) at the given natural coordinates
    pub fn normal_vector_at_gp(&self, gp: &Gp) -> Vector3<f64> {
        match self.element_type {
            ElementType::Tri => {
                let enodes = &self.mesh.elements[self.element_id].node_ids;
                let e0 = &self.mesh.nodes[enodes[0]].coords;
                let e1 = &self.mesh.nodes[enodes[1]].coords;
                let e2 = &self.mesh.nodes[enodes[2]].coords;
                let a = e1 - e0;
                let b = e2 - e0;
                a.cross(&b)
            },
            // for quads or higher, use shape functions to calculate the derivative dx / dxi,
            // i.e., change in physical coordinates w.r.t. natural coordinates
            _ => {
                let dn = self.shape_derivatives_at(gp);
                let mut dndxi = Vector3::from_element(0.0);
                let mut dndeta = dndxi.clone();
                let element = &self.mesh.elements[self.element_id];
                for i in 0..self.get_num_nodes() {
                    let node_index = &element.node_ids[i];
                    let icoord = &self.mesh.nodes[*node_index].coords;
                    dndxi += dn[(i,0)] * icoord;
                    dndeta += dn[(i,1)] * icoord;
                }
                // normal vector is then the cross product of dx/dxi and dx/deta
                dndxi.cross(&dndeta)
            }
        }
    }
    /// Get the determinant of the transformation from natural to physical coordinates at the given 
    /// natural coordinates
    fn detj_at(&self, gp: &Gp) -> f64 {
        // Calculate using the normal vector, since it is return un-scaled
        match self.element_type {
            ElementType::Tri => 0.5 * self.normal_vector_at_gp(gp).norm(),
            ElementType::Quad => 0.25 * self.normal_vector_at_gp(gp).norm(),
            _ => 0.0
        }
    }
    /// Numerically integrate the influence matrices from this element to a given origin point (x)
    pub fn influence_matrices_at(&self, k: f64, x: &Coords, n_x: &Vector3<f64>, use_hypersingular: bool) 
        -> (Vec::<Cplx>, Vec::<Cplx>) {
        let mut h = vec![Cplx::new(0.0, 0.0); self.get_num_nodes()];
        let mut g = h.clone();
        for gp in &self.integration {
            let n_y = self.normal_vector_at_gp(gp);
            let detj = self.detj_at(gp);
            let y = self.coordinates_at(gp);
            let (g_gp, h_gp) = get_greens_functions(k, x, n_x, &y, &n_y, use_hypersingular);
            let n = self.shape_functions_at(gp);
            for i in 0..n.len() {
                h[i] += h_gp * n[i] * detj * gp.wt;
                g[i] += g_gp * n[i] * detj * gp.wt;
            }
        }
        return (h, g)
    }
    pub fn get_integration_points_and_normals(&self) -> Vec<CollocationPoint> {
        let mut result: Vec<CollocationPoint> = Vec::new();
        for gp in &self.integration {
            let y = self.coordinates_at(gp);
            let n_y = self.normal_vector_at_gp(gp);
            let dw = self.detj_at(gp) * gp.wt;
            result.push(CollocationPoint { id: 0, coords: y, normal: n_y, dw: dw })
        }
        return result
    }
}