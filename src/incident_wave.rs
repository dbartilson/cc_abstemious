use crate::preprocess::{self, input_data};
use na::{Complex, DVector, Vector3};
use std::f64::consts::PI;
type Cplx = Complex<f64>;

/// return the incident velocity potential and normal velocity on the surface
pub fn get_incident_surface(predata: &preprocess::PreData) -> DVector::<Cplx> {
    let mesh = predata.get_mesh();
    let eqn_map = predata.get_eqn_map();
    let k = predata.get_wavenumber();
    let is_burton_miller = *predata.get_method_type() == input_data::MethodType::BurtonMiller;
    let beta = Cplx::new(0.0, 1.0 / k);

    let num_eqn = eqn_map.len();

    // incident rhs can be 
    // phi_inc                  for classical method
    // phi_inc + beta * vn_inc  for Burton-Miller method
    let mut rhs_inc = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    
    let inc_wave = predata.get_incident_wave();
    // amplitude in pressure units
    let p_amp = Cplx::new(inc_wave.amplitude[0], inc_wave.amplitude[1]);
    // amplitude in velocity potential units (phi = p / (i omega rho))
    let amp = p_amp / Cplx::new(0.0, predata.get_angular_frequency() * predata.get_mass_density());
    let vector = &inc_wave.origin;
    let mut vec3 = Vector3::from_column_slice(vector);
    vec3.normalize_mut();
    match inc_wave.wave_type {
        preprocess::input_data::WaveType::PlaneWave => {
            for (inode, ieqn) in eqn_map {
                let coord = &mesh.nodes[*inode].coords;
                // phi_inc = A * exp(ik (x dot d))
                let phi_inc = amp * Cplx::new(0., k * vec3.dot(coord)).exp();
                rhs_inc[*ieqn] = phi_inc;
                if is_burton_miller {
                    // vn_inc = phi_inc * ik * (e_n dot d)
                    let normal = &mesh.nodes[*inode].normal;
                    let vn_inc = phi_inc * Cplx::new(0.0, k) * vec3.dot(normal);
                    rhs_inc[*ieqn] += beta * vn_inc
                }
            }
        }
        preprocess::input_data::WaveType::SphericalWave => {
            for (inode, ieqn) in eqn_map {
                let coord = &mesh.nodes[*inode].coords;
                let rvec = coord - vec3;
                let r = rvec.magnitude();
                let e_r = rvec / r;
                // phi_inc = A / (4 pi r) * exp(ikr)
                let phi_inc = amp * Cplx::new(0., k * r).exp() / (4.0 * PI * r);
                rhs_inc[*ieqn] = phi_inc;
                if is_burton_miller {
                    // vn_inc = phi_inc * (ik - 1/r) * (e_n dot e_r)
                    let normal = &mesh.nodes[*inode].normal;
                    let vn_inc = phi_inc * Cplx::new(-1.0 / r, k) * e_r.dot(normal);
                    rhs_inc[*ieqn] += beta * vn_inc
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
    
    let inc_wave = predata.get_incident_wave();
    // amplitude in pressure units
    let p_amp = Cplx::new(inc_wave.amplitude[0], inc_wave.amplitude[1]);
    // amplitude in velocity potential units (phi = p / (i omega rho))
    let amp = p_amp / Cplx::new(0.0, predata.get_angular_frequency() * predata.get_mass_density());
    let vector = &inc_wave.origin;
    let mut vec3 = Vector3::from_column_slice(vector);
    vec3 /= vec3.magnitude();
    match inc_wave.wave_type {
        preprocess::input_data::WaveType::PlaneWave => {
            for (i, fp) in field_points.iter().enumerate() {
                let coord = Vector3::from_column_slice(fp);
                phi_inc_fp[i] = amp * Cplx::new(0., k * vec3.dot(&coord)).exp();
            }
        }
        preprocess::input_data::WaveType::SphericalWave => {
            for (i, fp) in field_points.iter().enumerate() {
                let coord = Vector3::from_column_slice(fp);
                let r = (coord - vec3).magnitude();
                phi_inc_fp[i] = amp * Cplx::new(0., k * r).exp() / (4.0 * PI * r);
            }
        }
    }
    return phi_inc_fp;
}