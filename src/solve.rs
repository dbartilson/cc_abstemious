pub mod gmres;
pub mod aca;
pub mod h_matrix;

use na::{DMatrix, DVector};
use crate::influence_matrix;
use crate::influence_matrix::{get_surface_row_or_column, EqnSide};
use crate::preprocess::{self, input_data};
use crate::Cplx;

pub fn solve_for_surface<'a>(predata: &'a preprocess::PreData, rhs_inc: &DVector::<Cplx>) 
    -> (DVector::<Cplx>, DVector::<Cplx>) {
    match predata.get_solver() {
        preprocess::input_data::Solver::Direct {  } |
        preprocess::input_data::Solver::Iterative { .. } => {
            get_surface_dense(predata, rhs_inc)
        }
        preprocess::input_data::Solver::Hierarchical { .. } => {
            get_surface_hmatrix(predata, rhs_inc)
        }
    }  
}

fn get_surface_dense(predata: &preprocess::PreData, rhs_inc: &DVector::<Cplx>) 
    -> (DVector::<Cplx>, DVector::<Cplx>) {

    let (h, g) = influence_matrix::get_dense_surface_matrices(predata);
    info!(" Solving system of equations...");
    let omega = predata.get_angular_frequency();
    let rho = predata.get_mass_density();
    let num_eqn = predata.get_num_eqn();

    let mut phi = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    let mut vn = phi.clone();
    let sbc = predata.get_surface_bc();
    match sbc.bc_type {
        preprocess::input_data::BCType::Pressure => {
            let pressure_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            phi.fill(pressure_bc);
            vn = rhs_inc.clone();
            vn.gemv(Cplx::new(-1.0,0.0), &h, &phi, Cplx::new(1.0,0.0));
            //vn = phi_inc - h * phi; // rhs
            let lhs = g.clone();
            solve_dense(predata, &lhs, &mut vn);
        }   
        preprocess::input_data::BCType::NormalVelocity => {
            let velocity_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            vn.fill(velocity_bc);
            phi = rhs_inc.clone();
            phi.gemv(Cplx::new(1.0,0.0), &g, &vn, Cplx::new(-1.0,0.0));
            //phi = g * vn - phi_inc; // rhs
            let lhs = h.clone();
            solve_dense(predata, &lhs, &mut phi);
        }
        preprocess::input_data::BCType::Impedance => {
            let impedance_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            let factor = Cplx::new(0., -omega * rho) / impedance_bc;
            let mut lhs = h.clone();
            for i in 0..lhs.nrows() {
                for j in 0..lhs.ncols() {
                    lhs[(i,j)] += factor * g[(i,j)];
                }
            }
            phi = -rhs_inc.clone(); // rhs
            solve_dense(predata, &lhs, &mut phi);
            // back-compute for surface velocities
            vn.axpy(factor, &phi, Cplx::new(0.0,0.0));
        }
    }
    return (phi, vn);
}


fn get_surface_hmatrix(predata: &preprocess::PreData, rhs_inc: &DVector::<Cplx>) 
    -> (DVector::<Cplx>, DVector::<Cplx>) {
    let num_eqn = predata.get_num_eqn();
    let mut phi = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    let mut vn = phi.clone();
    let mut rhs = rhs_inc.clone();
    {
        info!(" Calculating RHS...");
        let get_row_or_column = |i: Vec<usize>, j: Vec<usize>| get_surface_row_or_column(predata, i, j, EqnSide::RHS);
        let sbc = predata.get_surface_bc();
        let hmatrix = if sbc.value[0] == 0.0 && sbc.value[1] == 0.0 {
            h_matrix::HMatrix::new()
        } else {
            h_matrix::HMatrix::new_from(num_eqn, 
                &get_row_or_column, 
                &predata.get_cpts(), 
                predata.get_eqn_map(),
                32,
                1e-4)
        };
        match sbc.bc_type {
            preprocess::input_data::BCType::Pressure => {
                let pressure_bc = Cplx::new(sbc.value[0], sbc.value[1]);
                phi.fill(pressure_bc);
                // rhs = phi_inc - H * phi
                hmatrix.gemv(Cplx::new(1.0, 0.0), &phi, Cplx::new(1.0, 0.0), &mut rhs);
            }   
            preprocess::input_data::BCType::NormalVelocity => {
                let velocity_bc = Cplx::new(sbc.value[0], sbc.value[1]);
                vn.fill(velocity_bc);
                // rhs = G*vn - phi_inc
                hmatrix.gemv(Cplx::new(1.0, 0.0), &vn, Cplx::new(-1.0, 0.0), &mut rhs);
            }
            preprocess::input_data::BCType::Impedance => {
                // solve for phi, but need to post-process for vn later
                rhs.axpy(Cplx::new(0.0, 0.0), &rhs_inc, Cplx::new(-1.0, 0.0));
            }
        }
    }
    info!(" Calculating LHS...");
    let get_row_or_column = |i, j| get_surface_row_or_column(predata, i, j, EqnSide::LHS);
    let hmatrix = h_matrix::HMatrix::new_from(num_eqn, 
                                                       &get_row_or_column, 
                                                       &predata.get_cpts(), 
                                                       predata.get_eqn_map(),
                                                       32,
                                                       1e-4);
    if let input_data::Solver::Hierarchical { tolerance, max_iterations } = predata.get_solver() {
        let max_it = max_iterations;
        let tol = tolerance;
        let mut gm = gmres::GMRES::new(*max_it, *tol);
        gm.hmatrix = Some(hmatrix);
        info!(" Solving system of equations...");
        gm.solve(&mut rhs);
    }
    let sbc = predata.get_surface_bc();
    match sbc.bc_type {
        preprocess::input_data::BCType::Pressure => {
            // rhs contains solved-for velocity
            vn = rhs;
        }   
        preprocess::input_data::BCType::NormalVelocity => {
            // rhs contains solved-for phi
            phi = rhs;
        }
        preprocess::input_data::BCType::Impedance => {
            phi = rhs;
            // back-compute for surface velocities
            let omega = predata.get_angular_frequency();
            let rho = predata.get_mass_density();
            let impedance_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            let factor = Cplx::new(0., -omega * rho) / impedance_bc;
            vn.axpy(factor, &phi, Cplx::new(0.0,0.0));
        }
    }
    return (phi, vn);
}

/// Solve the system of equations for dense cases (direct or iterative)
fn solve_dense(predata: &preprocess::PreData, a: &DMatrix<Cplx>, x: &mut DVector<Cplx>) {
    match predata.get_solver() {
        preprocess::input_data::Solver::Direct { } => solve_lu(a, x),
        preprocess::input_data::Solver::Iterative { tolerance, max_iterations } => {
            let mut gm = gmres::GMRES::new(*max_iterations, *tolerance);
            gm.a = Some(a.clone());
            gm.solve(x);
        },
        _ => {}
    }
}

/// Solve the dense system of equations using direct (LU) approach
pub fn solve_lu(a: &DMatrix<Cplx>, x: &mut DVector<Cplx>) {
    // solve A*x = b using direct LU, where the input vector b is overwritten by the solution
    info!(" Using direct (LU) solver...");
    let a_lu = a.clone().lu();
    a_lu.solve_mut(x);
}

pub fn get_field(predata: &preprocess::PreData, m: &DMatrix::<Cplx>, l: &DMatrix::<Cplx>, 
    phi: &DVector::<Cplx>, vn: &DVector::<Cplx>, phi_inc_fp: &DVector::<Cplx>) -> DVector::<Cplx> {

    let want_total = *predata.get_output_type() == preprocess::input_data::OutputType::Total;
    let mut phi_fp = phi_inc_fp.clone();
    if !want_total {phi_fp.fill(Cplx::new(0.0,0.0));}

    phi_fp.gemv(-Cplx::new(1.0, 0.0), &m, &phi, Cplx::new(1.0, 0.0));
    phi_fp.gemv(-Cplx::new(-1.0, 0.0), &l, &vn, Cplx::new(1.0, 0.0));
    //let phi_fp = m * phi - l * vn + phi_inc_fp;

    return phi_fp;
}

#[cfg(test)]
mod tests {
    extern crate approx;
    use crate::solve::gmres;
    use crate::Cplx;
    use crate::solve;

    pub fn generate_random_ab(m: usize, n: usize, seed: u64) -> (na::DMatrix::<Cplx>, na::DVector::<Cplx>) {
        use rand::{Rng, SeedableRng, rngs::StdRng};
        let mut a = na::DMatrix::<Cplx>::from_diagonal_element(m, n, Cplx::new(n as f64, 0.0));

        let mut rng = StdRng::seed_from_u64(seed);
        for ai in a.iter_mut() {
            ai.re += rng.random::<f64>();
            ai.im += rng.random::<f64>();
        }
        let mut b = na::DVector::<Cplx>::from_element(n, Cplx::new(0.0, 0.0));
        for bi in b.iter_mut() {
            bi.re = rng.random::<f64>();
            bi.im = rng.random::<f64>();
        }
        return (a, b)
    }

    #[test]
    fn solve_lu() {
        let n = 100;
        let (a, mut b) = generate_random_ab(n, n, 10);
        solve::solve_lu(&a, &mut b);
        approx::assert_relative_eq!(b[0].re, 0.0035386336544537098, epsilon = 1.0e-10);
        approx::assert_relative_eq!(b[0].im, 0.0017794290945383035, epsilon = 1.0e-10);
    }

    #[test]
    fn solve_gmres() {
        let n = 100;
        let (a, mut b) = generate_random_ab(n, n, 10);
        let mut gm = gmres::GMRES::new(100, 1.0e-8);
        gm.a = Some(a);
        gm.solve(&mut b);
        approx::assert_relative_eq!(b[0].re, 0.0035386336544537098, epsilon = 1.0e-6);
        approx::assert_relative_eq!(b[0].im, 0.0017794290945383035, epsilon = 1.0e-6);
    }
}
