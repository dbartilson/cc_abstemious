use std::sync::{Arc, Mutex};
use scoped_threadpool::Pool;
use na::{DMatrix, Vector3, DVector};

use crate::preprocess::{self, input_data};
use crate::elements::*;
use crate::Cplx;

fn get_gh_functions(predata: &preprocess::PreData, i: usize, j: usize) -> (Cplx, Cplx) {
    if i == j { return (Cplx::new(0.0, 0.0), Cplx::new(0.0, 0.0)) }
    
    let cpti = &predata.get_cpts()[i];
    let y = &cpti.coords;
    let n_y = &cpti.normal;
    let cptj = &predata.get_cpts()[j];
    let x = &cptj.coords;
    let n_x = &cptj.normal;
    let (g, h) = get_greens_functions(predata.get_wavenumber(), x, n_x, y, n_y, predata.use_hypersingular());
    return (g * cptj.dw, h * cptj.dw)
}

/// evaluate the surface BEM influence matrices. These matrices are complex-valued,
/// square, and non-symmetric in general
pub fn get_dense_surface_matrices(predata: &preprocess::PreData) 
    -> (DMatrix::<Cplx>, DMatrix::<Cplx>) {

    info!(" Assembling surface BEM influence matrices...");
    let k = predata.get_wavenumber();

    let num_eqn = predata.get_num_eqn();
    let ncpts = predata.get_cpts().len();

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
        for i in 0..ncpts {
            let h_share = h_share.clone();
            let g_share = g_share.clone();
            scope.execute(move|| {
                for j in 0..ncpts {
                    let (g_j, h_j) = get_gh_functions(predata, i, j);
                    let mut hi = h_share.lock().unwrap();
                    let mut gi = g_share.lock().unwrap();
                    hi[(i, j)] += h_j;
                    gi[(i, j)] += g_j;
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
    let k = predata.get_wavenumber();

    let field_points = predata.get_field_points();
    let nfp = field_points.len();
    let ncpts = mesh.cpts.len();
    let num_eqn = predata.get_num_eqn();

    let mut m = DMatrix::<Cplx>::from_element(nfp, num_eqn,  Cplx::new(0.,0.));
    let mut l = m.clone();
    let n_x = Vector3::new(0.0, 0.0, 0.0); // dummy
    for j in 0..ncpts {
        let cptj = &mesh.cpts[j];
        let y = &cptj.coords;
        let n_y = &cptj.normal;
        for (i, fieldpt) in field_points.iter().enumerate()  {
            let x = Vector3::from_column_slice(fieldpt);
            let (g_j, h_j) = get_greens_functions(k, &x, &n_x, y, n_y, false);
            m[(i, j)] += h_j * cptj.dw;
            l[(i, j)] += g_j * cptj.dw;
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

    let num_column = j.len();
    let mut h = DVector::<Cplx>::from_element(num_column, Cplx::new(0.0, 0.0));
    let mut g = h.clone();
    for (index, jeqn) in j.iter().enumerate() {
        let (g_j, h_j) = get_gh_functions(&predata, *i, *jeqn);
        h[index] += h_j;
        g[index] += g_j;
    }
    if let Some(diag) = j.iter().position(|&r| r == *i) {
        h[diag] += predata.get_hdiag();
        g[diag] += predata.get_gdiag();
    }

    return (h, g)
}

fn get_surface_matrices_column(predata: &preprocess::PreData, i: &Vec<usize>, j: &usize) 
        -> (DVector::<Cplx>, DVector::<Cplx>) {

    let num_row = i.len();
    let mut h = DVector::<Cplx>::from_element(num_row, Cplx::new(0.0, 0.0));
    let mut g = h.clone();
    for (index, ieqn) in i.iter().enumerate() {
        let (g_j, h_j) = get_gh_functions(&predata, *ieqn, *j);
        h[index] += h_j;
        g[index] += g_j;
    }
    if let Some(diag) = i.iter().position(|&r| r == *j) {
        h[diag] += predata.get_hdiag();
        g[diag] += predata.get_gdiag();
    }

    return (h, g)
}