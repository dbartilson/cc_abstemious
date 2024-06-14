use std::collections::HashMap;
use na::{DMatrix, Complex, Vector3};
use crate::preprocess::input_data as id;
use crate::preprocess::mesh_data as mesh;
use crate::elements::*;
use std::f64::consts::PI;
type Cplx = Complex<f64>;

pub fn get_surface_influence_matrices(user_input: &id::UserInput, mesh: &mesh::Mesh, 
    eqn_map: &HashMap<usize, usize>) -> (DMatrix::<Cplx>, DMatrix::<Cplx>) {

    let omega = 2.0 * PI * user_input.frequency;
    let c = &user_input.sound_speed;
    let k = omega / c;

    let num_eqn = eqn_map.len();
    let body_id = &user_input.body_index;
    let nelem = &mesh.bodies[body_id-1].element_ids.len();

    let hdiag = match user_input.problem_type {
        id::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
        id::ProblemType::Interior => Cplx::new(0.0, 0.0)
    };
    let mut h = DMatrix::<Cplx>::from_diagonal_element(num_eqn, num_eqn, hdiag);
    let mut g = DMatrix::<Cplx>::from_element(num_eqn, num_eqn, Cplx::new(0.,0.));
    for (inode, ieqn) in eqn_map {
        let o = &mesh.nodes[*inode].coords;
        for e in 0..*nelem {
            let e_id = &mesh.bodies[body_id-1].element_ids[e];

            let enodes = &mesh.elements[*e_id-1].node_ids;
            let mut e_eqns = Vec::<usize>::new();
            for enode in enodes {
                match eqn_map.get(enode) {
                    Some(eqn) => e_eqns.push(*eqn),
                    None => println!("Eqn not found for element node {}", enode)
                }
            }
            let mut he = Vec::new();
            let mut ge = Vec::new();
            match &mesh.elements[*e_id-1].etype {
                mesh::ElementType::Tri => {
                    let tri = Triangle::new(&mesh, *e_id);
                    (he, ge) = tri.influence_matrices_at(k, o);
                },
                mesh::ElementType::Quad => {
                    let quad = Quad::new(&mesh, *e_id);
                    (he, ge) = quad.influence_matrices_at(k, o);
                }
                _ => {
                    println!("Invalid element!");
                }
            }

            for j in 0..e_eqns.len() {
                h[(*ieqn, e_eqns[j])] += he[j];
                g[(*ieqn, e_eqns[j])] += ge[j];
            }
        }
    }
    return (h, g)

}

pub fn get_field_influence_matrices(user_input: &id::UserInput, mesh: &mesh::Mesh, 
    eqn_map: &HashMap<usize, usize>) -> (DMatrix::<Cplx>, DMatrix::<Cplx>) {

    let omega = 2.0 * PI * user_input.frequency;
    let c = &user_input.sound_speed;
    let k = omega / c;

    let nfp = user_input.field_points.len();
    let body_id = &user_input.body_index;
    let nelem = &mesh.bodies[body_id-1].element_ids.len();
    let num_eqn = eqn_map.len();

    let mut m = DMatrix::<Cplx>::from_element(nfp, num_eqn,  Cplx::new(0.,0.));
    let mut l = m.clone();
    for i in 0..nfp {
        let fieldpt = &user_input.field_points[i];
        let coord = Vector3::from_column_slice(fieldpt);
        for e in 0..*nelem {
            let e_id = &mesh.bodies[body_id-1].element_ids[e];
            let enodes = &mesh.elements[*e_id-1].node_ids;
            let mut e_eqns = Vec::<usize>::new();
            for enode in enodes {
                match eqn_map.get(enode) {
                    Some(eqn) => e_eqns.push(*eqn),
                    None => println!("Eqn not found for element node {}", enode)
                }
            }
            let mut me = Vec::new();
            let mut le = Vec::new();
            match &mesh.elements[*e_id-1].etype {
                mesh::ElementType::Tri => {
                    let tri = Triangle::new(&mesh, *e_id);
                    (me, le) = tri.influence_matrices_at(k, &coord);
                },
                mesh::ElementType::Quad => {
                    let quad = Quad::new(&mesh, *e_id);
                    (me, le) = quad.influence_matrices_at(k, &coord);
                }
                _ => {
                    println!("Invalid element!");
                }
            }
            for j in 0..e_eqns.len() {
                m[(i, e_eqns[j])] += me[j];
                l[(i, e_eqns[j])] += le[j];
            }
        }
    }
    return (m, l);
}