pub mod gmres;

use na::{DMatrix, DVector};
use crate::preprocess;
use crate::Cplx;

pub fn get_surface(predata: &preprocess::PreData, h: &DMatrix::<Cplx>, g: &DMatrix::<Cplx>, phi_inc: &DVector::<Cplx>) -> (DVector::<Cplx>, DVector::<Cplx>) {

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
            vn = phi_inc.clone();
            vn.gemv(Cplx::new(-1.0,0.0), &h, &phi, Cplx::new(1.0,0.0));
            //vn = phi_inc - h * phi; // rhs
            let lhs = g.clone();
            solve(predata, lhs, &mut vn);
        }   
        preprocess::input_data::BCType::NormalVelocity => {
            let velocity_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            vn.fill(velocity_bc);
            phi = phi_inc.clone();
            phi.gemv(Cplx::new(1.0,0.0), &g, &vn, Cplx::new(-1.0,0.0));
            //phi = g * vn - phi_inc; // rhs
            let lhs = h.clone();
            solve(predata, lhs, &mut phi);
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
            phi = -phi_inc.clone(); // rhs
            solve(predata, lhs, &mut phi);
        }
    }
    return (phi, vn);
}

fn solve(predata: &preprocess::PreData, a: DMatrix<Cplx>, x: &mut DVector<Cplx>) {
    match predata.get_solver_type() {
        preprocess::input_data::SolverType::Direct => solve_lu(a, x),
        preprocess::input_data::SolverType::Iterative => gmres::solve_gmresk(&a, x, predata.get_solver_max_it(), predata.get_solver_tolerance())
    }
}

pub fn solve_lu(a: DMatrix<Cplx>, x: &mut DVector<Cplx>) {
    // solve A*x = b using direct LU, where the input vector b is overwritten by the solution
    info!(" Using direct (LU) solver...");
    let a_lu = a.lu();
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
    use crate::Cplx;
    use crate::solve;

    fn generate_random_ab(n: usize, seed: u64) -> (na::DMatrix::<Cplx>, na::DVector::<Cplx>) {
        use rand::{Rng, SeedableRng, rngs::StdRng};
        let mut a = na::DMatrix::<Cplx>::from_diagonal_element(n, n, Cplx::new(n as f64, 0.0));

        let mut rng = StdRng::seed_from_u64(seed);
        for ai in a.iter_mut() {
            ai.re += rng.gen::<f64>();
            ai.im += rng.gen::<f64>();
        }
        let mut b = na::DVector::<Cplx>::from_element(n, Cplx::new(0.0, 0.0));
        for bi in b.iter_mut() {
            bi.re = rng.gen::<f64>();
            bi.im = rng.gen::<f64>();
        }
        return (a, b)
    }

    #[test]
    fn random_matrix_lu() {
        let n = 100;
        let (a, mut b) = generate_random_ab(n, 10);
        solve::solve_lu(a, &mut b);
        approx::assert_relative_eq!(b[0].re, 0.0035386336544537098, epsilon = 1.0e-10, max_relative = 1.0);
        approx::assert_relative_eq!(b[0].im, 0.0017794290945383035, epsilon = 1.0e-10, max_relative = 1.0);
    }

    #[test]
    fn random_matrix_gmres() {
        let n = 100;
        let (a, mut b) = generate_random_ab(n, 10);
        solve::gmres::solve_gmresk(&a, &mut b, 100, 1.0e-8);
        approx::assert_relative_eq!(b[0].re, 0.0035386336544537098, epsilon = 1.0e-10, max_relative = 1.0);
        approx::assert_relative_eq!(b[0].im, 0.0017794290945383035, epsilon = 1.0e-10, max_relative = 1.0);
    }
}
