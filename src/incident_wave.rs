use std::collections::HashMap;
use crate::preprocess::{input_data::*, mesh_data::Mesh};
use na::{DVector, Complex, Vector3};
use std::f64::consts::PI;
type Cplx = Complex<f64>;

pub fn get_incident_wave(user_input: &UserInput, mesh: &Mesh, 
        eqn_map: &HashMap<usize, usize>) -> (DVector::<Cplx>, DVector::<Cplx>) {

    let omega = 2.0 * PI * user_input.frequency;
    let c = &user_input.sound_speed;
    let k = omega / c;

    let num_eqn = eqn_map.len();
    let num_fp = user_input.field_points.len();

    let mut phi_inc = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    let mut phi_inc_fp = DVector::<Cplx>::from_element(num_fp, Cplx::new(0., 0.));
    
    let fp = &user_input.field_points;
    let inc_wave = &user_input.incident_wave;
    let amp = Cplx::new(inc_wave.amplitude[0], inc_wave.amplitude[1]);
    let vector = &inc_wave.origin;
    let vec3 = Vector3::from_column_slice(vector);
    match inc_wave.wave_type {
        WaveType::PlaneWave => {
            for (inode, ieqn) in eqn_map {
                let coord = &mesh.nodes[*inode].coords;
                phi_inc[*ieqn] = amp * Cplx::new(0., k * vec3.dot(coord));
            }
            for i in 0..num_fp {
                let coord = Vector3::from_column_slice(&fp[i]);
                phi_inc_fp[i] = amp * Cplx::new(0., k * vec3.dot(&coord));
            }
        }
        WaveType::SphericalWave => {
            for (inode, ieqn) in eqn_map {
                let coord = &mesh.nodes[*inode].coords;
                let r = (coord - vec3).magnitude();
                phi_inc[*ieqn] = amp * Cplx::new(0., k * r) / (4.0 * PI * r);
            }
            for i in 0..num_fp {
                let coord = Vector3::from_column_slice(&fp[i]);
                let r = (coord - vec3).magnitude();
                phi_inc_fp[i] = amp * Cplx::new(0., k * r) / (4.0 * PI * r);
            }
        }
    }
    return (phi_inc, phi_inc_fp);
}