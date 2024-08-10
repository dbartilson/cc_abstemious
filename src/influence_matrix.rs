use std::collections::HashSet;
use std::sync::{Arc, Mutex};
use scoped_threadpool::Pool;
use na::{DMatrix, Vector3, DVector};

use crate::preprocess::{self, input_data};
use crate::elements::*;
use crate::Cplx;

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

    let hdiag = predata.get_hdiag();
    let gdiag = predata.get_gdiag();
    let num_threads = preprocess::get_num_threads();
    let use_hypersingular = *predata.get_method_type() == input_data::MethodType::BurtonMiller;
    // use a parallel pool of threads
    info!(" Using {} threads...", num_threads);
    let mut pool = Pool::new(num_threads as u32);
    let h_share = Arc::new(Mutex::new(DMatrix::<Cplx>::from_diagonal_element(num_eqn, num_eqn, hdiag)));
    let g_share = Arc::new(Mutex::new(DMatrix::<Cplx>::from_diagonal_element(num_eqn, num_eqn, gdiag)));
    pool.scoped(|scope| {
        for e in 0..nelem {
            let h_share = h_share.clone();
            let g_share = g_share.clone();
            scope.execute(move|| {
                let e_id = &mesh_body.element_ids[e];
                let e_eqns = &mesh.elements[*e_id].eqn_idx;
                let element = NIElement::new(&mesh, *e_id);
                for (inode, ieqn) in eqn_map {
                    let x = &mesh.nodes[*inode].coords;
                    let n_x = &mesh.nodes[*inode].normal;
                    let (he, ge) = element.influence_matrices_at(k, x, n_x, use_hypersingular);
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
    let n_x = Vector3::new(0.0, 0.0, 0.0); // dummy
    for e in 0..nelem {
        let e_id = &mesh_body.element_ids[e];
        let e_eqns = &mesh.elements[*e_id].eqn_idx;
        let element = NIElement::new(&mesh, *e_id);
        for (i, fieldpt) in field_points.iter().enumerate()  {
            let x = Vector3::from_column_slice(fieldpt);
            let (me, le) = element.influence_matrices_at(k, &x, &n_x, false);
            for j in 0..e_eqns.len() {
                m[(i, e_eqns[j])] += me[j];
                l[(i, e_eqns[j])] += le[j];
            }
        }
    }
    return (m, l);
}

/// return alpha and beta for system matrix [alpha*H + beta*G]
fn get_lhs_factors(predata: &preprocess::PreData) -> (Cplx, Cplx) {
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
fn get_rhs_factors(predata: &preprocess::PreData) -> (Cplx, Cplx) {
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

pub enum EqnSide {
    RHS,
    LHS
}

/// Get row or column of the LHS or RHS, i or j must be singleton vector
pub fn get_surface_row_or_column(predata: &preprocess::PreData, 
                                 i: Vec<usize>, 
                                 j: Vec<usize>, 
                                 side: EqnSide) -> Vec<Cplx> {
    if !((i.len() == 1)^(j.len() == 1)) {
        error!("Invalid call to get_surface_row_or_column: i or j must be singleton");
    }
    let (alpha, beta) = match side {
        EqnSide::LHS => get_lhs_factors(predata),
        EqnSide::RHS => get_rhs_factors(predata),
    };
    if i.len() == 1 {
        let (mut h, g) = get_surface_matrices_row(predata, &i[0], &j);
        h.axpy(beta, &g, alpha);
        let rr: Vec<Cplx> = h.data.into();
        return rr      
    }
    else {
        let (mut h, g) = get_surface_matrices_column(predata, &i, &j[0]);
        h.axpy(beta, &g, alpha);
        let rr: Vec<Cplx> = h.data.into();
        return rr
    };  
}

fn get_surface_matrices_row(predata: &preprocess::PreData, i: &usize, j: &Vec<usize>) 
        -> (DVector::<Cplx>, DVector::<Cplx>) {

    let mesh = predata.get_mesh();
    let k = predata.get_wavenumber();
    let num_column = j.len();
    let use_hypersingular = *predata.get_method_type() == input_data::MethodType::BurtonMiller;
    let mut el_list = HashSet::<usize>::new();
    for jeqn in j {
        if let Some(node_index) = predata.get_node_map().get(jeqn) {
            let els = &predata.get_revcon()[*node_index];
            for el in els {
                el_list.insert(*el);
            }
        }
    }
    // find node index corresponding to equation j
    let inode = match predata.get_node_map().get(i) {
        Some(node) => node,
        None => {error!("Node index not found for equation {}", i); &0}
    };
    let x = &mesh.nodes[*inode].coords;
    let n_x = &mesh.nodes[*inode].normal;
    let mut h = DVector::<Cplx>::from_element(num_column, Cplx::new(0.0, 0.0));
    let mut g = h.clone();
    for e in el_list {
        let e_eqns = &mesh.elements[e].eqn_idx;
        let element = NIElement::new(&mesh, e);
        let (he, ge) = element.influence_matrices_at(k, x, n_x, use_hypersingular);
        for k in 0..e_eqns.len() {
            if let Some(index) = j.iter().position(|&r| r == e_eqns[k]) {
                h[index] += he[k];
                g[index] += ge[k];
            }
        }
    }
    if let Some(diag) = j.iter().position(|&r| r == *i) {
        h[diag] += predata.get_hdiag();
        g[diag] += predata.get_gdiag();
    }

    return (h, g)
}

fn get_surface_matrices_column(predata: &preprocess::PreData, i: &Vec<usize>, j: &usize) 
        -> (DVector::<Cplx>, DVector::<Cplx>) {

    let mesh = predata.get_mesh();
    let k = predata.get_wavenumber();
    let use_hypersingular = *predata.get_method_type() == input_data::MethodType::BurtonMiller;
    let num_eqn = i.len();

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
    for e_id in el_list {
        let enodes = &mesh.elements[*e_id].node_ids;
        // find which node index of this element corresponds to this column
        let index = match enodes.iter().position(|&r| r == *jnode) {
            Some(found) => found,
            None => {error!("Node not found for element"); 0}
        };
        let element = NIElement::new(&mesh, *e_id);
        for (eqn_index, ieqn) in i.iter().enumerate() {
            let inode = match predata.get_node_map().get(ieqn) {
                Some(node) => node,
                None => {error!("Node index not found for equation {}", ieqn); &0}
            };
            let x = &mesh.nodes[*inode].coords;
            let n_x = &mesh.nodes[*inode].normal;
            let (he, ge) = element.influence_matrices_at(k, x, n_x, use_hypersingular);
            h[eqn_index] += he[index];
            g[eqn_index] += ge[index];
        }
    }
    if let Some(diag) = i.iter().position(|&r| r == *j) {
        h[diag] += predata.get_hdiag();
        g[diag] += predata.get_gdiag();
    }

    return (h, g)
}