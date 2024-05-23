pub mod interpolation {
    pub struct Gp {
        pub coords: [f64; 2],
        pub wt: f64,
    }
    pub static TRIGP3: [Gp; 3] = [Gp{coords: [1./6., 1./6.], wt: 1./6.}, 
                                    Gp{coords: [1./6., 2./3.], wt: 1./6.}, 
                                    Gp{coords: [2./3., 1./6.], wt: 1./6.}];
    static ONEOVERSQRT3: f64 = 0.57735026919;
    pub static QUADGP4: [Gp; 4] = [Gp{coords: [ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                    Gp{coords: [-ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                    Gp{coords: [ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0},
                                    Gp{coords: [-ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0}];
}   

use interpolation::*;
use na::{Vector3, ComplexField, Complex};
type Cplx = Complex<f64>;
use crate::preprocess::mesh_data::{Coords, Mesh};

// Methods for numerically-integrated elements
pub trait NumIntElement {
    fn shape_functions_at(gp: &Gp) -> Vec<f64>;
    fn coordinates_at(&self, gp: &Gp) -> Coords;
    fn normal_vector_at(&self) -> Vector3<f64>;
    fn influence_matrices_at(&self, k: f64, origin: &Coords) -> (Vec::<Cplx>, Vec::<Cplx>);
}

pub struct Triangle <'a> {
    pub integration: &'a [Gp; 3],
    pub mesh: &'a Mesh,
    pub element_id: usize,
}
impl Triangle <'_> {
    pub fn new<'a>(meshdata: &'a Mesh, element: usize) -> Triangle <'a> {
        Triangle{integration: &TRIGP3, mesh: &meshdata, element_id: element}
    }
}
impl NumIntElement for Triangle <'_> {
    fn shape_functions_at(gp: &Gp) -> Vec<f64> {
        let xi = gp.coords[0];
        let eta = gp.coords[1];
        let n = [1.0 - xi - eta,
                            xi,
                            eta];
        return n.to_vec();
    }
    fn coordinates_at(&self, gp: &Gp) -> Coords {
        let n = Triangle::shape_functions_at(gp);
        let mut x = Coords::from_element(0.0);
        let element = &self.mesh.elements[self.element_id-1];
        for ni in 0..n.len() {
            let node_index = &element.node_ids[ni];
            let icoord = &self.mesh.nodes[*node_index].coords;
            x += n[ni] * icoord;
        }
        return x;
    }
    fn normal_vector_at(&self) -> Vector3<f64> {
        let enodes = &self.mesh.elements[self.element_id-1].node_ids;
        let e0 = &self.mesh.nodes[enodes[0]].coords;
        let e1 = &self.mesh.nodes[enodes[1]].coords;
        let e2 = &self.mesh.nodes[enodes[2]].coords;
        let a = e1 - e0;
        let b = e2 - e0;
        b.cross(&a)
    }
    fn influence_matrices_at(&self, k: f64, origin: &Coords) -> (Vec::<Cplx>, Vec::<Cplx>) {

        let mut h = vec![Cplx::new(0.0, 0.0); 3];
        let mut g = h.clone();

        let normal = self.normal_vector_at();
        let detj = normal.norm();
        for gp in self.integration {
            let x = self.coordinates_at(gp);
            let r = origin - x;
            let rdist = r.norm();
            let runit = r / rdist;
            let g_gp = Cplx::new(0.0, k*rdist).exp() / (4.0 * std::f64::consts::PI * rdist);
            let h_gp = g_gp * Cplx::new(1.0 / rdist, -k) * runit.dot(&normal);
            let n = Triangle::shape_functions_at(gp);
            for i in 0..3 {
                h[i] += h_gp * n[i] * detj * gp.wt;
                g[i] += g_gp * n[i] * detj * gp.wt;
            }
        }
        return (h, g);
    }
}

pub struct Quad <'a> {
    pub integration: &'a [Gp; 4],
    pub mesh: &'a Mesh,
    pub element_id: usize,
}
impl Quad <'_> {
    pub fn new<'a>(meshdata: &'a Mesh, element: usize) -> Quad <'a> {
        Quad{integration: &QUADGP4, mesh: &meshdata, element_id: element}
    }
}
impl NumIntElement for Quad <'_> {
    fn shape_functions_at(gp: &Gp) -> Vec<f64> {
        let xi = gp.coords[0];
        let eta = gp.coords[1];
        let n = [0.25*(1.-xi)*(1.-eta),
                            0.25*(1.+xi)*(1.-eta),
                            0.25*(1.+xi)*(1.+eta),
                            0.25*(1.-xi)*(1.+eta)];
        return n.to_vec();
    }
    fn coordinates_at(&self, gp: &Gp) -> Coords {
        let n = Quad::shape_functions_at(gp);
        let mut x = Coords::from_element(0.0);
        let element = &self.mesh.elements[self.element_id-1];
        for ni in 0..n.len() {
            let node_index = &element.node_ids[ni];
            let icoord = &self.mesh.nodes[*node_index].coords;
            x += n[ni] * icoord;
        }
        return x;
    }
    fn normal_vector_at(&self) -> Vector3<f64> {
        let enodes = &self.mesh.elements[self.element_id-1].node_ids;
        let e0 = &self.mesh.nodes[enodes[0]].coords;
        let e1 = &self.mesh.nodes[enodes[1]].coords;
        let e2 = &self.mesh.nodes[enodes[2]].coords;
        let a = e1 - e0;
        let b = e2 - e0;
        b.cross(&a)
    }
    fn influence_matrices_at(&self, k: f64, origin: &Coords) -> (Vec::<Cplx>, Vec::<Cplx>) {

        let mut h = vec![Cplx::new(0.0, 0.0); 4];
        let mut g = h.clone();

        let normal = self.normal_vector_at();
        let detj = normal.norm();
        for gp in self.integration {
            let x = self.coordinates_at(gp);
            let r = origin - x;
            let rdist = r.norm();
            let runit = r / rdist;
            let g_gp = Cplx::new(0.0, k*rdist).exp() / (4.0 * std::f64::consts::PI * rdist);
            let h_gp = g_gp * Cplx::new(1.0 / rdist, -k) * runit.dot(&normal);
            let n = Quad::shape_functions_at(gp);
            for i in 0..4 {
                h[i] += h_gp * n[i] * detj * gp.wt;
                g[i] += g_gp * n[i] * detj * gp.wt;
            }
        }
        return (h, g);
    }
}