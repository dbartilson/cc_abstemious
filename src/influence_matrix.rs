use na::{DMatrix, Vector3};
use crate::preprocess;
use crate::preprocess::input_data as id;
use crate::preprocess::mesh_data as mesh;
use crate::elements::*;
use crate::Cplx;

pub fn get_surface_influence_matrices(predata: &preprocess::PreData) 
    -> (DMatrix::<Cplx>, DMatrix::<Cplx>) {
    // evaluate the surface BEM influence matrices. These matrices are complex-valued,
    // square, and non-symmetric in general

    info!(" Assembling surface BEM influence matrices...");

    let mesh = predata.get_mesh();
    let eqn_map = predata.get_eqn_map();
    let k = predata.get_wavenumber();

    let num_eqn = eqn_map.len();
    let mesh_body = predata.get_mesh_body();
    let nelem = mesh_body.element_ids.len();

    let hdiag = match predata.get_problem_type() {
        // the H matrix has -1/2 added along the diagonal for exterior problems
        id::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
        id::ProblemType::Interior => Cplx::new(0.0, 0.0)
    };
    let mut h = DMatrix::<Cplx>::from_diagonal_element(num_eqn, num_eqn, hdiag);
    let mut g = DMatrix::<Cplx>::from_element(num_eqn, num_eqn, Cplx::new(0.,0.));
    for e in 0..nelem {
        let e_id = &mesh_body.element_ids[e];
        let enodes = &mesh.elements[*e_id].node_ids;
        let mut e_eqns = Vec::<usize>::new();
        let mut he = Vec::new();
        let mut ge = Vec::new();
        for enode in enodes {
            match eqn_map.get(enode) {
                Some(eqn) => e_eqns.push(*eqn),
                None => error!("Eqn not found for element node {}", enode)
            }
        }
        match &mesh.elements[*e_id].etype {
            mesh::ElementType::Tri => {
                let tri = Triangle::new(&mesh, *e_id);
                he = vec![Cplx::new(0.0, 0.0); 3];
                ge = he.clone();
                for (inode, ieqn) in eqn_map {
                    let o = &mesh.nodes[*inode].coords;
                    tri.influence_matrices_at(k, o, &mut he, &mut ge);
                    for j in 0..e_eqns.len() {
                        h[(*ieqn, e_eqns[j])] += he[j];
                        g[(*ieqn, e_eqns[j])] += ge[j];
                    }
                }
            },
            mesh::ElementType::Quad => {
                let quad = Quad::new(&mesh, *e_id);
                he = vec![Cplx::new(0.0, 0.0); 4];
                ge = he.clone();
                for (inode, ieqn) in eqn_map {
                    let o = &mesh.nodes[*inode].coords;
                    quad.influence_matrices_at(k, o, &mut he, &mut ge);
                    for j in 0..e_eqns.len() {
                        h[(*ieqn, e_eqns[j])] += he[j];
                        g[(*ieqn, e_eqns[j])] += ge[j];
                    }
                }
            }
            _ => {
                error!("Invalid element!");
            }
        }
    }
    return (h, g)
}

pub fn get_field_influence_matrices(predata: &preprocess::PreData) -> (DMatrix::<Cplx>, DMatrix::<Cplx>) {
    // evaluate the field BEM influence matrices. These matrices are complex-valued, and typically rectangular

    info!(" Calculating field results...");

    let mesh = predata.get_mesh();
    let eqn_map = predata.get_eqn_map();
    let k = predata.get_wavenumber();

    let field_points = predata.get_field_points();
    let nfp = field_points.len();
    let mesh_body = predata.get_mesh_body();
    let nelem = mesh_body.element_ids.len();
    let num_eqn = eqn_map.len();

    let mut m = DMatrix::<Cplx>::from_element(nfp, num_eqn,  Cplx::new(0.,0.));
    let mut l = m.clone();
    for (i, fieldpt) in field_points.iter().enumerate() {
        let coord = Vector3::from_column_slice(fieldpt);
        for e in 0..nelem {
            let e_id = &mesh_body.element_ids[e];
            let enodes = &mesh.elements[*e_id].node_ids;
            let mut e_eqns = Vec::<usize>::new();
            for enode in enodes {
                match eqn_map.get(enode) {
                    Some(eqn) => e_eqns.push(*eqn),
                    None => error!("Eqn not found for element node {}", enode)
                }
            }
            let mut me = Vec::new();
            let mut le = Vec::new();
            match &mesh.elements[*e_id].etype {
                mesh::ElementType::Tri => {
                    let tri = Triangle::new(&mesh, *e_id);
                    (me, le) = tri.influence_matrices_at(k, &coord);
                },
                mesh::ElementType::Quad => {
                    let quad = Quad::new(&mesh, *e_id);
                    (me, le) = quad.influence_matrices_at(k, &coord);
                }
                _ => {
                    error!("Invalid element!");
                }
            }
            // assemble
            for j in 0..e_eqns.len() {
                m[(i, e_eqns[j])] += me[j];
                l[(i, e_eqns[j])] += le[j];
            }
        }
    }
    return (m, l);
}