pub mod interpolation {
    #[derive(Clone)]
    pub struct Gp {
        pub coords: [f64; 2],
        pub wt: f64,
    }
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
    pub static QUADGP4: [Gp; 4] = [Gp{coords: [ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                   Gp{coords: [-ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                   Gp{coords: [ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0},
                                   Gp{coords: [-ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0}];
}   

use interpolation::*;
use na::{DMatrix, Vector3};
use crate::Cplx;
use crate::preprocess::mesh_data::{Coords, ElementType, Mesh};

fn get_greens_functions(k: f64, origin: &Coords, x: &Coords, normal: &Vector3<f64>) -> (Cplx, Cplx) {
    let n = normal / normal.norm();
    let r = origin - x;
    let rdist = r.norm();
    let runit = r / rdist;
    let g = Cplx::new(0.0, k*rdist).exp() / (4.0 * std::f64::consts::PI * rdist);
    let h = g * Cplx::new(1.0 / rdist, -k) * runit.dot(&n);
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
    fn get_num_nodes(&self) -> usize {
        match self.element_type {
            ElementType::Tri => 3,
            ElementType::Quad => 4,
            _ => 0
        }
    }
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
    pub fn coordinates_at(&self, gp: &Gp) -> Coords {
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
    pub fn normal_vector_at(&self, gp: &Gp) -> Vector3<f64> {
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
            // for quads or higher
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
                dndxi.cross(&dndeta)
            }
        }
    }
    pub fn detj_at(&self, gp: &Gp) -> f64 {
        match self.element_type {
            ElementType::Tri => 0.5 * self.normal_vector_at(gp).norm(),
            ElementType::Quad => 0.25 * self.normal_vector_at(gp).norm(),
            _ => 0.0
        }
    }
    pub fn influence_matrices_at(&self, k: f64, origin: &Coords) -> (Vec::<Cplx>, Vec::<Cplx>) {
        let mut h = vec![Cplx::new(0.0, 0.0); self.get_num_nodes()];
        let mut g = h.clone();
        for gp in &self.integration {
            let normal = self.normal_vector_at(gp);
            let detj = self.detj_at(gp);
            let x = self.coordinates_at(gp);
            let (g_gp, h_gp) = get_greens_functions(k, origin, &x, &normal);
            let n = self.shape_functions_at(gp);
            for i in 0..n.len() {
                h[i] += h_gp * n[i] * detj * gp.wt;
                g[i] += g_gp * n[i] * detj * gp.wt;
            }
        }
        return (h, g)
    }
}