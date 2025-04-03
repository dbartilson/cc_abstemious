pub mod interpolation {
    #[derive(Clone, Copy)]
    pub struct Gp {
        pub coords: [f64; 2],
        pub wt: f64,
    }
    #[allow(dead_code)]
    pub static TRINODES: [Gp; 3] = [Gp{coords: [0.0, 0.0], wt: 1./3.}, 
                                    Gp{coords: [1.0, 0.0], wt: 1./3.}, 
                                    Gp{coords: [0.0, 1.0], wt: 1./3.}];
    pub static TRIGP1: [Gp; 1] = [Gp{coords: [1./3., 1./3.], wt: 1.0}];
    #[allow(dead_code)]
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
    #[allow(dead_code)]
    pub static QUADNODES: [Gp; 4] = [Gp{coords: [-1.0, -1.0], wt: 1.0}, 
                                     Gp{coords: [ 1.0, -1.0], wt: 1.0}, 
                                     Gp{coords: [ 1.0,  1.0], wt: 1.0},
                                     Gp{coords: [-1.0,  1.0], wt: 1.0}];
    pub static QUADGP1: [Gp; 1] = [Gp{coords: [0., 0.], wt: 4.0}];
    #[allow(dead_code)]
    pub static QUADGP4: [Gp; 4] = [Gp{coords: [ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                   Gp{coords: [-ONEOVERSQRT3, ONEOVERSQRT3], wt: 1.0}, 
                                   Gp{coords: [ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0},
                                   Gp{coords: [-ONEOVERSQRT3, -ONEOVERSQRT3], wt: 1.0}];
}   

use interpolation::*;
use na::{DMatrix, Vector3};
use crate::preprocess::mesh_data::{CollocationPoint, Coords, ElementType, Mesh};

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
            ElementType::Tri => TRIGP1.to_vec(),
            ElementType::Quad => QUADGP1.to_vec(),
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
    #[allow(dead_code)]
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
    pub fn get_integration_points_and_normals(&self) -> Vec<CollocationPoint> {
        let mut result: Vec<CollocationPoint> = Vec::new();
        for gp in &self.integration {
            let y = self.coordinates_at(gp);
            let n_y = self.normal_vector_at_gp(gp);
            let area = self.detj_at(gp);
            result.push(CollocationPoint { id: 0, coords: y, normal: n_y.normalize(), area: area, wt: gp.wt })
        }
        return result
    }
}