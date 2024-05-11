use std::path::Path;

//#[macro_use]
extern crate nalgebra as na;

pub mod model_data;
use model_data::model_data::MeshData;

pub mod elements;
use elements::elements::interpolation::TRIGP;

use na::{ComplexField, DMatrix, Vector3, Complex};

use crate::elements::elements::Triangle;

fn main() -> std::io::Result<()> {

    let path = Path::new("D:/Downloads/folder/docs/untitled.vtk");
    let mut meshdata: MeshData = Default::default();
    let result = meshdata.read_from_vtk(path);

    let ibody = 3;
    let nelem = &meshdata.bodies[ibody-1].element_ids.len();
    let nnode = &meshdata.nodes.len();

    type Cplx = Complex<f64>;
    type Coords = Vector3<f64>;
    let mut h = DMatrix::<Cplx>::from_element(*nnode, *nnode, Cplx::new(0.,0.));
    let mut g = h.clone();

    let k = 1.0;

    let ibody = &meshdata.bodies[ibody-1];
    for i in 0..*nnode {
        h[(i, i)].re = 0.5;
        let o = Vector3::from_column_slice(&meshdata.nodes[i].coords);
        for e in 0..*nelem {
            let e_id = ibody.element_ids[e] - 1;
            let enodes = &meshdata.elements[e_id].node_ids;
            let mut ecoords: Vec<Coords> = Vec::new();
            for enode in enodes {
                ecoords.push(Vector3::from_column_slice(&meshdata.nodes[*enode].coords));
            }
            let a = ecoords[1] - ecoords[0];
            let b = ecoords[2] - ecoords[0];
            let normal = b.cross(&a);
            let detj = normal.norm();
            let mut x = Coords::from_element(0.0);
            for gp in &TRIGP {
                let n = Triangle::shape_functions_at(gp);
                for ni in 0..n.len() {
                    x += n[ni] * ecoords[ni];
                }
                let r = o - x;
                let rdist = r.norm();
                let runit = r / rdist;
                let g_gp = Cplx::new(0.0, k*rdist).exp() / (4.0 * std::f64::consts::PI * rdist);
                let h_gp = g_gp * Cplx::new(1.0 / rdist, -k) * runit.dot(&normal);
                for enode in enodes {
                    h[(i, *enode)] += h_gp * detj * gp.wt;
                    g[(i, *enode)] += g_gp * detj * gp.wt;
                }
            }
        }
    }
    return result;
}
