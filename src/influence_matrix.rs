use std::sync::{Arc, Mutex};
use scoped_threadpool::Pool;
use na::{DMatrix, Vector3};
use serde_json::error;

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
    let num_threads = match std::thread::available_parallelism() {
        Ok(result) => std::cmp::max(result.get() / 2, 2),
        Err(_) => 2
    };
    info!(" Using {} threads...", num_threads);
    let mut pool = Pool::new(num_threads as u32);
    let h_share = Arc::new(Mutex::new(DMatrix::<Cplx>::from_diagonal_element(num_eqn, num_eqn, hdiag)));
    let g_share = Arc::new(Mutex::new(DMatrix::<Cplx>::from_element(num_eqn, num_eqn, Cplx::new(0.,0.))));
    pool.scoped(|scope| {
        for e in 0..nelem {
            let h_share = h_share.clone();
            let g_share = g_share.clone();
            scope.execute(move|| {
                let e_id = &mesh_body.element_ids[e];
                let enodes = &mesh.elements[*e_id].node_ids;
                let mut e_eqns = Vec::<usize>::new();
                for enode in enodes {
                    match eqn_map.get(enode) {
                        Some(eqn) => e_eqns.push(*eqn),
                        None => error!("Eqn not found for element node {}", enode)
                    }
                }
                let element = NIElement::new(&mesh, *e_id);
                for (inode, ieqn) in eqn_map {
                    let o = &mesh.nodes[*inode].coords;
                    let (he, ge) = element.influence_matrices_at(k, o);
                    let mut hi = h_share.lock().unwrap();
                    let mut gi = g_share.lock().unwrap();
                    for j in 0..e_eqns.len() {
                        hi[(*ieqn, e_eqns[j])] += he[j];
                        gi[(*ieqn, e_eqns[j])] += ge[j];
                    }
                }
            });
        }
    });
    let h = Arc::try_unwrap(h_share).unwrap().into_inner().unwrap();
    let g = Arc::try_unwrap(g_share).unwrap().into_inner().unwrap();
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
        let element = NIElement::new(&mesh, *e_id);
        for (i, fieldpt) in field_points.iter().enumerate()  {
            let coord = Vector3::from_column_slice(fieldpt);
            let (me, le) = element.influence_matrices_at(k, &coord);
            for j in 0..e_eqns.len() {
                m[(i, e_eqns[j])] += me[j];
                l[(i, e_eqns[j])] += le[j];
            }
        }
    }
    return (m, l);
}

pub fn get_surface_matrices_row(predata: &preprocess::PreData, i: &usize) {

    let mesh = predata.get_mesh();
    let eqn_map = predata.get_eqn_map();
    let k = predata.get_wavenumber();

    let num_eqn = eqn_map.len();
    let mesh_body = predata.get_mesh_body();
    let nelem = mesh_body.element_ids.len();

    // find node index corresponding to equation j
    let i_node = match predata.get_node_map().get(i) {
        Some(node) => node,
        None => {error!("Node index not found for equation {}", i); &0}
    };
}