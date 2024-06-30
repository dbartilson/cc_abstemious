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
use na::{Complex, DMatrix, Vector3};
type Cplx = Complex<f64>;
use crate::preprocess::mesh_data::{Coords, Mesh};

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
pub trait NumIntElement {
    fn shape_functions_at(gp: &Gp) -> Vec<f64>;
    fn shape_derivatives_at(gp: &Gp) -> DMatrix<f64>;
    fn coordinates_at(&self, gp: &Gp) -> Coords;
    fn normal_vector_at(&self, gp: &Gp) -> Vector3<f64>;
    fn detj_at(&self, gp: &Gp) -> f64;
    fn get_integration(&self) -> &Vec<Gp>;
    #[allow(dead_code)]
    fn get_num_nodes(&self) -> usize;
    fn influence_matrices_at(&self, k: f64, origin: &Coords, h: &mut Vec::<Cplx>, g: &mut Vec::<Cplx>) {
        h.fill(Cplx::new(0.0, 0.0));
        g.fill(Cplx::new(0.0, 0.0));
        for gp in self.get_integration() {
            let normal = self.normal_vector_at(gp);
            let detj = self.detj_at(gp);
            let x = self.coordinates_at(gp);
            let (g_gp, h_gp) = get_greens_functions(k, origin, &x, &normal);
            let n = Self::shape_functions_at(gp);
            for i in 0..n.len() {
                h[i] += h_gp * n[i] * detj * gp.wt;
                g[i] += g_gp * n[i] * detj * gp.wt;
            }
        }
    }
}

pub struct Triangle <'a> {
    pub integration: Vec<Gp>,
    pub mesh: &'a Mesh,
    pub element_id: usize,
}
impl Triangle <'_> {
    pub fn new<'a>(meshdata: &'a Mesh, element: usize) -> Triangle <'a> {
        Triangle{integration: TRIGP3.to_vec(), mesh: &meshdata, element_id: element}
    }
}
impl NumIntElement for Triangle <'_> {
    fn shape_functions_at(gp: &Gp) -> Vec<f64> {
        let xi = gp.coords[0];
        let eta = gp.coords[1];
        let n = vec![1.0 - xi - eta,
                               xi,
                               eta];
        return n
    }
    fn shape_derivatives_at(_gp: &Gp) -> DMatrix<f64> {
        let dn = DMatrix::from_row_slice(3, 2,
            &[-1.0, -1.0,
                    1.0, 0.0,
                    0.0, 1.0]);
        return dn
    }
    fn coordinates_at(&self, gp: &Gp) -> Coords {
        let n = Triangle::shape_functions_at(gp);
        let mut x = Coords::from_element(0.0);
        let element = &self.mesh.elements[self.element_id];
        for ni in 0..n.len() {
            let node_index = &element.node_ids[ni];
            let icoord = &self.mesh.nodes[*node_index].coords;
            x += n[ni] * icoord;
        }
        return x;
    }
    fn normal_vector_at(&self, gp: &Gp) -> Vector3<f64> {
        let dn = Self::shape_derivatives_at(gp);
        let mut dndxi = Vector3::from_element(0.0);
        let mut dndeta = dndxi.clone();
        let element = &self.mesh.elements[self.element_id];
        for i in 0..3 {
            let node_index = &element.node_ids[i];
            let icoord = &self.mesh.nodes[*node_index].coords;
            dndxi += dn[(i,0)] * icoord;
            dndeta += dn[(i,1)] * icoord;
        }
        return dndxi.cross(&dndeta)
    }
    fn detj_at(&self, gp: &Gp) -> f64 {
        return self.normal_vector_at(gp).norm() * 0.5
    }
    fn get_integration(&self) -> &Vec<Gp> {
        &self.integration
    }
    fn get_num_nodes(&self) -> usize {3}
}

pub struct Quad <'a> {
    pub integration: Vec<Gp>,
    pub mesh: &'a Mesh,
    pub element_id: usize,
}
impl Quad <'_> {
    pub fn new<'a>(meshdata: &'a Mesh, element: usize) -> Quad <'a> {
        Quad{integration: QUADGP4.to_vec(), mesh: &meshdata, element_id: element}
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
    fn shape_derivatives_at(_gp: &Gp) -> DMatrix<f64> {
        let dn = DMatrix::from_row_slice(4, 2,
            &[-0.25, -0.25,
                    0.25, -0.25,
                    0.25, 0.25,
                   -0.25, 0.25]);
        return dn
    }
    fn coordinates_at(&self, gp: &Gp) -> Coords {
        let n = Quad::shape_functions_at(gp);
        let mut x = Coords::from_element(0.0);
        let element = &self.mesh.elements[self.element_id];
        for ni in 0..n.len() {
            let node_index = &element.node_ids[ni];
            let icoord = &self.mesh.nodes[*node_index].coords;
            x += n[ni] * icoord;
        }
        return x;
    }
    fn normal_vector_at(&self, gp: &Gp) -> Vector3<f64> {
        let dn = Self::shape_derivatives_at(gp);
        let mut dndxi = Vector3::from_element(0.0);
        let mut dndeta = dndxi.clone();
        let element = &self.mesh.elements[self.element_id];
        for i in 0..4 {
            let node_index = &element.node_ids[i];
            let icoord = &self.mesh.nodes[*node_index].coords;
            dndxi += dn[(i,0)] * icoord;
            dndeta += dn[(i,1)] * icoord;
        }
        return dndxi.cross(&dndeta)
    }
    fn detj_at(&self, gp: &Gp) -> f64 {
        return 0.25 * self.normal_vector_at(gp).norm()
    }
    fn get_integration(&self) -> &Vec<Gp> {
        &self.integration
    }
    fn get_num_nodes(&self) -> usize {4}
}