/*!
Defines the surface elements
*/

pub mod interpolation;

use crate::preprocess::mesh::{CollocationPoint, Coords, ElementType, Mesh};
use interpolation::*;
use na::{DMatrix, Vector3};

/// Numerically integrated elements
///
/// Contains reference to the mesh and a copy of the integration scheme
pub struct NIElement<'a> {
    /// Vector of Gauss Points (integration scheme)
    integration: Vec<Gp>,
    mesh: &'a Mesh,
    pub element_id: usize,
    element_type: ElementType,
}
impl NIElement<'_> {
    pub fn new(meshdata: &Mesh, element: usize) -> NIElement<'_> {
        let etype = &meshdata.elements[element].etype;
        let int = match etype.clone() {
            ElementType::Tri => TRIGP3.to_vec(),
            ElementType::Quad => QUADGP4.to_vec(),
            _ => {
                error!("Invalid numerically integrated element");
                Vec::new()
            }
        };
        NIElement {
            integration: int,
            mesh: meshdata,
            element_id: element,
            element_type: etype.clone(),
        }
    }
    /// Set up element with one integration point (for the collocation point setup)
    pub fn new_1_point_integration(meshdata: &Mesh, element: usize) -> NIElement<'_> {
        let etype = &meshdata.elements[element].etype;
        let int = match etype.clone() {
            ElementType::Tri => TRIGP1.to_vec(),
            ElementType::Quad => QUADGP1.to_vec(),
            _ => {
                error!("Invalid numerically integrated element");
                Vec::new()
            }
        };
        NIElement {
            integration: int,
            mesh: meshdata,
            element_id: element,
            element_type: etype.clone(),
        }
    }
    /// Get number of nodes for element
    #[inline]
    pub fn get_num_nodes(&self) -> usize {
        match self.element_type {
            ElementType::Tri => 3,
            ElementType::Quad => 4,
            _ => 0,
        }
    }
    /// Get shape functions for this element at the given natural coordinates
    fn shape_functions_at(&self, gp: &Gp) -> Vec<f64> {
        let xi = &gp.coords[0];
        let eta = &gp.coords[1];
        match self.element_type {
            ElementType::Tri => {
                vec![1.0 - *xi - *eta, *xi, *eta]
            }
            ElementType::Quad => {
                vec![
                    0.25 * (1. - *xi) * (1. - *eta),
                    0.25 * (1. + *xi) * (1. - *eta),
                    0.25 * (1. + *xi) * (1. + *eta),
                    0.25 * (1. - *xi) * (1. + *eta),
                ]
            }
            _ => vec![0.0],
        }
    }
    /// Get shape function derivatives for this element at the given natural coordinates
    fn shape_derivatives_at(&self, _gp: &Gp) -> DMatrix<f64> {
        match self.element_type {
            ElementType::Tri => DMatrix::from_row_slice(3, 2, &[-1.0, -1.0, 1.0, 0.0, 0.0, 1.0]),
            ElementType::Quad => {
                DMatrix::from_row_slice(4, 2, &[-0.25, -0.25, 0.25, -0.25, 0.25, 0.25, -0.25, 0.25])
            }
            _ => DMatrix::from_element(1, 1, 0.0),
        }
    }
    /// Get physical coordinates at the given natural coordinates, using shape functions
    fn coordinates_at(&self, gp: &Gp) -> Coords {
        let n = self.shape_functions_at(gp);
        let mut x = Coords::from_element(0.0);
        let element = &self.mesh.elements[self.element_id];
        for (i, ni) in n.iter().enumerate() {
            let node_index = &element.node_ids[i];
            let icoord = &self.mesh.nodes[*node_index].coords;
            x += *ni * icoord;
        }
        x
    }
    /// Get natural coordinates corresponding to node index
    #[allow(dead_code)]
    pub fn natural_coordinates_at_node(&self, i: usize) -> Gp {
        match self.element_type {
            ElementType::Tri => TRINODES[i],
            ElementType::Quad => QUADNODES[i],
            _ => Gp {
                coords: [0.0, 0.0],
                wt: 0.0,
            },
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
            }
            // for quads or higher, use shape functions to calculate the derivative dx / dxi,
            // i.e., change in physical coordinates w.r.t. natural coordinates
            _ => {
                let dn = self.shape_derivatives_at(gp);
                let mut dndxi = Vector3::from_element(0.0);
                let mut dndeta = dndxi;
                let element = &self.mesh.elements[self.element_id];
                for i in 0..self.get_num_nodes() {
                    let node_index = &element.node_ids[i];
                    let icoord = &self.mesh.nodes[*node_index].coords;
                    dndxi += dn[(i, 0)] * icoord;
                    dndeta += dn[(i, 1)] * icoord;
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
            ElementType::Quad => self.normal_vector_at_gp(gp).norm(),
            _ => 0.0,
        }
    }
    /// Define a set of collocation points corresponding to the integration scheme of the element
    ///
    /// Sets up physical coordinates, normal vector, area (detj), and weight from the natural coordinates
    pub fn get_integration_points_and_normals(&self) -> Vec<CollocationPoint> {
        let mut result: Vec<CollocationPoint> = Vec::new();
        for gp in &self.integration {
            let y = self.coordinates_at(gp);
            let n_y = self.normal_vector_at_gp(gp);
            let area = self.detj_at(gp);
            result.push(CollocationPoint {
                id: 0,
                coords: y,
                normal: n_y.normalize(),
                area,
                wt: gp.wt,
            })
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::{
        NIElement,
        interpolation::{Gp, TRIGP1},
    };
    use crate::{
        elements::interpolation::QUADGP1,
        preprocess::mesh::{Body, CollocationPoint, Coords, Element, ElementType, Mesh, Node},
    };
    use approx::assert_relative_eq;

    fn setup_dummy_mesh() -> Mesh {
        let nodes = vec![
            Node {
                id: 1,
                coords: Coords::new(0.0, 0.0, 0.0),
            },
            Node {
                id: 2,
                coords: Coords::new(1.0, 0.0, 0.0),
            },
            Node {
                id: 3,
                coords: Coords::new(1.0, 1.0, 0.0),
            },
            Node {
                id: 4,
                coords: Coords::new(0.0, 1.0, 0.0),
            },
        ];
        let elements = vec![
            Element {
                id: 1,
                body_id: 1,
                etype: ElementType::Tri,
                node_ids: vec![0, 1, 2],
                intpts: Vec::new(),
            },
            Element {
                id: 2,
                body_id: 1,
                etype: ElementType::Quad,
                node_ids: vec![0, 1, 2, 3],
                intpts: Vec::new(),
            },
        ];
        let bodies = vec![Body {
            id: 1,
            element_ids: vec![1, 2],
        }];
        let cpts: Vec<CollocationPoint> = Vec::new();
        Mesh {
            nodes,
            elements,
            bodies,
            cpts,
        }
    }

    #[test]
    fn tri_shape_function() {
        let mesh = setup_dummy_mesh();
        let element = NIElement::new(&mesh, 0);
        let gp = Gp {
            coords: [0.25, 0.25],
            wt: 0.0,
        };
        let sf = element.shape_functions_at(&gp);
        assert_relative_eq!(sf[0], 0.5, epsilon = 1e-8);
        assert_relative_eq!(sf[1], 0.25, epsilon = 1e-8);
        assert_relative_eq!(sf[2], 0.25, epsilon = 1e-8);
    }

    #[test]
    fn quad_shape_function() {
        let mesh = setup_dummy_mesh();
        let element = NIElement::new(&mesh, 1);
        let gp = Gp {
            coords: [0.25, 0.25],
            wt: 0.0,
        };
        let sf = element.shape_functions_at(&gp);
        assert_relative_eq!(sf[0], 0.140625, epsilon = 1e-8);
        assert_relative_eq!(sf[1], 0.234375, epsilon = 1e-8);
        assert_relative_eq!(sf[2], 0.390625, epsilon = 1e-8);
        assert_relative_eq!(sf[3], 0.234375, epsilon = 1e-8);
    }

    #[test]
    fn tri_area() {
        let mesh = setup_dummy_mesh();
        let element = NIElement::new(&mesh, 0);
        let gp = TRIGP1[0];
        let detj = element.detj_at(&gp);
        assert_relative_eq!(detj * gp.wt, 0.5, epsilon = 1e-9)
    }

    #[test]
    fn quad_area() {
        let mesh = setup_dummy_mesh();
        let element = NIElement::new(&mesh, 1);
        let gp = QUADGP1[0];
        let detj = element.detj_at(&gp);
        assert_relative_eq!(detj * gp.wt, 1.0, epsilon = 1e-9)
    }
}
