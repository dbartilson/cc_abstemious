use na::{DMatrix, DVector, Complex};
use std::f64::consts::PI;
use crate::preprocess::input_data::*;
type Cplx = Complex<f64>;   

pub fn solve_lu(user_input: &UserInput, h: &DMatrix::<Cplx>, g: &DMatrix::<Cplx>, phi_inc: &DVector::<Cplx>,
num_eqn: &usize) -> (DVector::<Cplx>, DVector::<Cplx>) {

    let omega = 2.0 * PI * user_input.frequency;
    let rho = &user_input.mass_density;

    let mut phi = DVector::<Cplx>::from_element(*num_eqn, Cplx::new(0., 0.));
    let mut vn = phi.clone();
    let sbc = &user_input.surface_bc;
    match sbc.bc_type {
        BCType::Pressure => {
            let pressure_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            phi.fill(pressure_bc);
            vn = phi_inc.clone();
            vn.gemv(Cplx::new(-1.0,0.0), &h, &phi, Cplx::new(1.0,0.0));
            //vn = phi_inc - h * phi; // rhs
            let lhs = g.clone();
            let glu = lhs.lu();
            glu.solve_mut(&mut vn);
        }   
        BCType::NormalVelocity => {
            let velocity_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            vn.fill(velocity_bc);
            phi = phi_inc.clone();
            phi.gemv(Cplx::new(1.0,0.0), &g, &vn, Cplx::new(-1.0,0.0));
            //phi = g * vn - phi_inc; // rhs
            let lhs = g.clone();
            let hlu = lhs.lu();
            hlu.solve_mut(&mut phi);
        }
        BCType::Impedance => {
            let impedance_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            let factor = Cplx::new(0., -omega * rho) / impedance_bc;
            let mut lhs = h.clone();
            for i in 0..lhs.nrows() {
                for j in 0..lhs.ncols() {
                    lhs[(i,j)] += factor * g[(i,j)];
                }
            }
            phi = -phi_inc.clone(); // rhs
            let hlu = lhs.lu();
            hlu.solve_mut(&mut phi);
        }
    }
    return (phi, vn);
}

pub fn get_fp(m: &DMatrix::<Cplx>, l: &DMatrix::<Cplx>, phi: &DVector::<Cplx>, 
    vn: &DVector::<Cplx>, phi_inc_fp: &DVector::<Cplx>) -> DVector::<Cplx> {

    let mut phi_fp = phi_inc_fp.clone();

    phi_fp.gemv(-Cplx::new(-1.0, 0.0), &m, &phi, Cplx::new(1.0, 0.0));
    phi_fp.gemv(-Cplx::new(1.0, 0.0), &l, &vn, Cplx::new(1.0, 0.0));
    //let phi_fp = - m * phi + l * vn + phi_inc_fp;
    return phi_fp;
}