use na::{DMatrix, DVector, Complex};
use crate::preprocess;
type Cplx = Complex<f64>;   

pub fn solve_lu(predata: &preprocess::PreData, h: &DMatrix::<Cplx>, g: &DMatrix::<Cplx>, phi_inc: &DVector::<Cplx>) -> (DVector::<Cplx>, DVector::<Cplx>) {

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
            let glu = lhs.lu();
            glu.solve_mut(&mut vn);
        }   
        preprocess::input_data::BCType::NormalVelocity => {
            let velocity_bc = Cplx::new(sbc.value[0], sbc.value[1]);
            vn.fill(velocity_bc);
            phi = phi_inc.clone();
            phi.gemv(Cplx::new(1.0,0.0), &g, &vn, Cplx::new(-1.0,0.0));
            //phi = g * vn - phi_inc; // rhs
            let lhs = g.clone();
            let hlu = lhs.lu();
            hlu.solve_mut(&mut phi);
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