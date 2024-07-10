use std::sync::{Arc, Mutex};
use scoped_threadpool::Pool;
use na::{DMatrix, Vector3, DVector};

use crate::preprocess;
use crate::preprocess::input_data as id;
use crate::elements::*;
use crate::Cplx;

/// return alpha and beta for system matrix [alpha*H + beta*G]
pub fn get_lhs_factors(predata: &preprocess::PreData) -> (Cplx, Cplx) {
    let sbc = predata.get_surface_bc();
    match sbc.bc_type {
        preprocess::input_data::BCType::Pressure => {
            (Cplx::new(0.0, 0.0), Cplx::new(1.0, 0.0))
        }   
        preprocess::input_data::BCType::NormalVelocity => {
            (Cplx::new(1.0, 0.0), Cplx::new(0.0, 0.0))
        }
        preprocess::input_data::BCType::Impedance => {
            let impedance_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            let omega = predata.get_angular_frequency();
            let rho = predata.get_mass_density();
            (Cplx::new(1.0, 0.0), Cplx::new(0.0, -omega * rho) / impedance_bc)
        }
    }
}

/// return gamma and delta for rhs [gamma*H + delta*G] * vn or phi
pub fn get_rhs_factors(predata: &preprocess::PreData) -> (Cplx, Cplx) {
    let sbc = predata.get_surface_bc();
    match sbc.bc_type {
        preprocess::input_data::BCType::Pressure => {
            (Cplx::new(-1.0, 0.0), Cplx::new(0.0, 0.0))
        }   
        preprocess::input_data::BCType::NormalVelocity => {
            (Cplx::new(0.0, 0.0), Cplx::new(1.0, 0.0))
        }
        preprocess::input_data::BCType::Impedance => {
            (Cplx::new(0.0, 0.0), Cplx::new(0.0, 0.0))
        }
    }
}

/// evaluate the surface BEM influence matrices. These matrices are complex-valued,
/// square, and non-symmetric in general
pub fn get_dense_surface_matrices(predata: &preprocess::PreData) 
    -> (DMatrix::<Cplx>, DMatrix::<Cplx>) {

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
    // use a parallel pool of threads
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

/// evaluate the field BEM influence matrices. These matrices are complex-valued, and typically rectangular
pub fn get_dense_field_matrices(predata: &preprocess::PreData) -> (DMatrix::<Cplx>, DMatrix::<Cplx>) {

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

/// Get the row of the LHS (for ACA)
pub fn get_lhs_row(predata: &preprocess::PreData, i: &usize) -> Vec<Cplx> {
    let (h, g) = get_surface_matrices_row(predata, &i);
    let (alpha, beta) = get_lhs_factors(predata);
    let mut r = h;
    // r = alpha * h + beta * g
    r.axpy(beta, &g, alpha);
    let rr: Vec<Cplx> = r.data.into();
    return rr
}
/// Get the column of the LHS (for ACA)
pub fn get_lhs_column(predata: &preprocess::PreData, j: &usize) -> Vec<Cplx> {
    let (h, g) = get_surface_matrices_column(predata, &j);
    let (alpha, beta) = get_lhs_factors(predata);
    let mut r = h;
    // r = alpha * h + beta * g
    r.axpy(beta, &g, alpha);
    let rr: Vec<Cplx> = r.data.into();
    return rr
}
/// Get the row of the RHS (for ACA)
pub fn get_rhs_row(predata: &preprocess::PreData, i: &usize) -> Vec<Cplx> {
    let (h, g) = get_surface_matrices_row(predata, &i);
    let (alpha, beta) = get_rhs_factors(predata);
    let mut r = h;
    // r = alpha * h + beta * g
    r.axpy(beta, &g, alpha);
    let rr: Vec<Cplx> = r.data.into();
    return rr
}
/// Get the column of the RHS (for ACA)
pub fn get_rhs_column(predata: &preprocess::PreData, j: &usize) -> Vec<Cplx> {
    let (h, g) = get_surface_matrices_column(predata, &j);
    let (alpha, beta) = get_rhs_factors(predata);
    let mut r = h;
    // r = alpha * h + beta * g
    r.axpy(beta, &g, alpha);
    let rr: Vec<Cplx> = r.data.into();
    return rr
}

fn get_surface_matrices_row(predata: &preprocess::PreData, i: &usize) 
        -> (DVector::<Cplx>, DVector::<Cplx>) {

    let mesh = predata.get_mesh();
    let eqn_map = predata.get_eqn_map();
    let k = predata.get_wavenumber();

    let num_eqn = eqn_map.len();
    let mesh_body = predata.get_mesh_body();
    let nelem = mesh_body.element_ids.len();

    // find node index corresponding to equation j
    let inode = match predata.get_node_map().get(i) {
        Some(node) => node,
        None => {error!("Node index not found for equation {}", i); &0}
    };
    let o = &mesh.nodes[*inode].coords;
    let mut h = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0.0, 0.0));
    let mut g = h.clone();
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
        let (he, ge) = element.influence_matrices_at(k, o);
        for j in 0..e_eqns.len() {
            h[e_eqns[j]] += he[j];
            g[e_eqns[j]] += ge[j];
        }
    }
    let hdiag = match predata.get_problem_type() {
        // the H matrix has -1/2 added along the diagonal for exterior problems
        id::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
        id::ProblemType::Interior => Cplx::new(0.0, 0.0)
    };
    h[*i] += hdiag;

    return (h, g)
}

fn get_surface_matrices_column(predata: &preprocess::PreData, j: &usize) 
        -> (DVector::<Cplx>, DVector::<Cplx>) {

    let mesh = predata.get_mesh();
    let eqn_map = predata.get_eqn_map();
    let k = predata.get_wavenumber();

    let num_eqn = eqn_map.len();
    let mesh_body = predata.get_mesh_body();

    // find node index corresponding to equation j
    let jnode = match predata.get_node_map().get(j) {
        Some(node) => node,
        None => {error!("Node index not found for equation {}", j); &0}
    };
    // get list of elements at this node
    let el_list = &predata.get_revcon()[*jnode];
    if el_list.is_empty() {error!("Element list is empty for node {}",j);}
    let mut h = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0.0, 0.0));
    let mut g = h.clone();
    for e in el_list {
        let e_id = &mesh_body.element_ids[*e];
        let enodes = &mesh.elements[*e_id].node_ids;
        // find which node index of this element corresponds to this column
        let index = match enodes.iter().position(|&r| r == *jnode) {
            Some(found) => found,
            None => {error!("Node not found for element"); 0}
        };
        let element = NIElement::new(&mesh, *e_id);
        for (inode, ieqn) in eqn_map {
            let o = &mesh.nodes[*inode].coords;
            let (he, ge) = element.influence_matrices_at(k, o);
            h[*ieqn] += he[index];
            g[*ieqn] += ge[index];
        }
    }
    let hdiag = match predata.get_problem_type() {
        // the H matrix has -1/2 added along the diagonal for exterior problems
        id::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
        id::ProblemType::Interior => Cplx::new(0.0, 0.0)
    };
    h[*j] += hdiag;

    return (h, g)
}