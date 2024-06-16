use crate::preprocess;
use na::{Complex, ComplexField, DVector, Vector3};
use std::f64::consts::PI;
type Cplx = Complex<f64>;

pub fn get_incident_wave(predata: &preprocess::PreData) -> (DVector::<Cplx>, DVector::<Cplx>) {

    let mesh = predata.get_mesh();
    let eqn_map = predata.get_eqn_map();
    let k = predata.get_wavenumber();

    let num_eqn = eqn_map.len();
    let field_points = predata.get_field_points();
    let num_fp = field_points.len();

    let mut phi_inc = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    let mut phi_inc_fp = DVector::<Cplx>::from_element(num_fp, Cplx::new(0., 0.));
    
    let inc_wave = predata.get_incident_wave();
    let amp = Cplx::new(inc_wave.amplitude[0], inc_wave.amplitude[1]);
    let vector = &inc_wave.origin;
    let mut vec3 = Vector3::from_column_slice(vector);
    vec3 /= vec3.magnitude();
    match inc_wave.wave_type {
        preprocess::input_data::WaveType::PlaneWave => {
            for (inode, ieqn) in eqn_map {
                let coord = &mesh.nodes[*inode].coords;
                phi_inc[*ieqn] = amp * Cplx::new(0., k * vec3.dot(coord)).exp();
            }
            for (i, fp) in field_points.iter().enumerate() {
                let coord = Vector3::from_column_slice(fp);
                phi_inc_fp[i] = amp * Cplx::new(0., k * vec3.dot(&coord)).exp();
            }
        }
        preprocess::input_data::WaveType::SphericalWave => {
            for (inode, ieqn) in eqn_map {
                let coord = &mesh.nodes[*inode].coords;
                let r = (coord - vec3).magnitude();
                phi_inc[*ieqn] = amp * Cplx::new(0., k * r).exp() / (4.0 * PI * r);
            }
            for (i, fp) in field_points.iter().enumerate() {
                let coord = Vector3::from_column_slice(fp);
                let r = (coord - vec3).magnitude();
                phi_inc_fp[i] = amp * Cplx::new(0., k * r).exp() / (4.0 * PI * r);
            }
        }
    }
    return (phi_inc, phi_inc_fp);
}