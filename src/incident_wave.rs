use crate::preprocess;
use na::{Complex, DVector, Vector3};
use std::f64::consts::PI;
type Cplx = Complex<f64>;

/// return the incident velocity potential and normal velocity on the surface
pub fn get_incident_surface(predata: &preprocess::PreData) -> DVector::<Cplx> {
    let cpts = predata.get_cpts();
    let k = predata.get_wavenumber();
    let hypersingular = predata.get_hypersingular();

    let num_eqn = predata.get_num_eqn();

    // incident rhs can be 
    // phi_inc                  for classical method
    // phi_inc + beta * vn_inc  for Burton-Miller method
    let mut rhs_inc = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    
    match predata.get_incident_wave() {
        &preprocess::input_data::IncidentWaveInput::PlaneWave { direction, amplitude } => {
            // amplitude in pressure units
            let p_amp = Cplx::new(amplitude[0], amplitude[1]);
            // amplitude in velocity potential units (phi = p / (i omega rho))
            let amp = p_amp / Cplx::new(0.0, predata.get_angular_frequency() * predata.get_mass_density());
            let dir3 = Vector3::from_column_slice(&direction).normalize();
            for cpt in cpts {
                let coord = &cpt.coords;
                // phi_inc = A * exp(ik (x dot d))
                let phi_inc = amp * Cplx::new(0., k * dir3.dot(coord)).exp();
                rhs_inc[cpt.id] = phi_inc;
                if hypersingular.is {
                    // vn_inc = phi_inc * ik * (e_n dot d)
                    let normal = &cpt.normal;
                    let vn_inc = phi_inc * Cplx::new(0.0, k) * dir3.dot(normal);
                    rhs_inc[cpt.id] += hypersingular.factor * vn_inc;
                }
            }
        }
        &preprocess::input_data::IncidentWaveInput::SphericalWave { origin, amplitude } => {
            // amplitude in pressure units
            let p_amp = Cplx::new(amplitude[0], amplitude[1]);
            // amplitude in velocity potential units (phi = p / (i omega rho))
            let amp = p_amp / Cplx::new(0.0, predata.get_angular_frequency() * predata.get_mass_density());
            let origin3 = Vector3::from_column_slice(&origin).normalize();
            for cpt in cpts {
                let coord = &cpt.coords;
                let rvec = coord - origin3;
                let r = rvec.magnitude();
                let e_r = rvec / r;
                // phi_inc = A / (4 pi r) * exp(ikr)
                let phi_inc = amp * Cplx::new(0., k * r).exp() / (4.0 * PI * r);
                rhs_inc[cpt.id] = phi_inc;
                if hypersingular.is {
                    // vn_inc = phi_inc * (ik - 1/r) * (e_n dot e_r)
                    let normal = &cpt.normal;
                    let vn_inc = phi_inc * Cplx::new(-1.0 / r, k) * e_r.dot(normal);
                    rhs_inc[cpt.id] += hypersingular.factor * vn_inc;
                }
            }
        }
    }
    return rhs_inc;
}

/// return the incident velocity potential at field points
pub fn get_incident_field(predata: &preprocess::PreData) -> DVector::<Cplx> {
    let k = predata.get_wavenumber();

    let field_points = predata.get_field_points();
    let num_fp = field_points.len();

    let mut phi_inc_fp = DVector::<Cplx>::from_element(num_fp, Cplx::new(0., 0.));
    
    match predata.get_incident_wave() {
        &preprocess::input_data::IncidentWaveInput::PlaneWave { direction, amplitude } => {
            // amplitude in pressure units
            let p_amp = Cplx::new(amplitude[0], amplitude[1]);
            // amplitude in velocity potential units (phi = p / (i omega rho))
            let amp = p_amp / Cplx::new(0.0, predata.get_angular_frequency() * predata.get_mass_density());
            let dir3 = Vector3::from_column_slice(&direction).normalize();
            for (i, fp) in field_points.iter().enumerate() {
                let coord = Vector3::from_column_slice(fp);
                phi_inc_fp[i] = amp * Cplx::new(0., k * dir3.dot(&coord)).exp();
            }
        }
        &preprocess::input_data::IncidentWaveInput::SphericalWave { origin, amplitude } => {
            // amplitude in pressure units
            let p_amp = Cplx::new(amplitude[0], amplitude[1]);
            // amplitude in velocity potential units (phi = p / (i omega rho))
            let amp = p_amp / Cplx::new(0.0, predata.get_angular_frequency() * predata.get_mass_density());
            let origin3 = Vector3::from_column_slice(&origin).normalize();
            for (i, fp) in field_points.iter().enumerate() {
                let coord = Vector3::from_column_slice(fp);
                let r = (coord - origin3).magnitude();
                phi_inc_fp[i] = amp * Cplx::new(0., k * r).exp() / (4.0 * PI * r);
            }
        }
    }
    return phi_inc_fp;
}