use std::{io::Write, path::Path};

extern crate nalgebra as na;

pub mod input_data;
pub mod model_data;
pub mod preprocess;
pub mod elements;

use model_data::model_data::MeshData;
use input_data::input_data::read_input_json;
use crate::{elements::elements::*, input_data::input_data::{BCType, ProblemType, WaveType}};
use na::{Complex, DMatrix, DVector, Vector3};
use std::f64::consts::PI;

type Cplx = Complex<f64>;

fn main() -> std::io::Result<()> {

    println!("=== CC-ABSTEMIOUS <=> BEM-ACOUSTICS ===");
    println!("Ver. 0.0");
    println!("");

    let input_path_str = "D:/Repositories/cc-abstemious/docs/input_1.json";
    let user_input = read_input_json(Path::new(input_path_str)).unwrap();

    let mut mesh: MeshData = Default::default();
    let result = mesh.read_from_vtk(Path::new(&user_input.mesh_file));

    let body_id = &user_input.body_index;
    let nelem = &mesh.bodies[body_id-1].element_ids.len();

    let omega = 2.0 * PI * user_input.frequency;
    let c = &user_input.sound_speed;
    let rho = &user_input.mass_density;
    let k = omega / c;

    // preprocess to get node to eqn map
    print!(" Preprocessing...");
    std::io::stdout().flush().unwrap();
    let (num_eqn, eqn_map) = preprocess::preprocess::get_eqn_map(&mesh, *body_id);
    let num_fp = user_input.field_points.len();
    println!(" Complete!");
    
    print!(" Assembling surface BEM influence matrices...");
    std::io::stdout().flush().unwrap();
    let hdiag = match user_input.problem_type {
        ProblemType::Exterior => Cplx::new(-0.5, 0.0),
        ProblemType::Interior => Cplx::new(0.0, 0.0)
    };
    let mut h = DMatrix::<Cplx>::from_diagonal_element(num_eqn, num_eqn, hdiag);
    let mut g = DMatrix::<Cplx>::from_element(num_eqn, num_eqn, Cplx::new(0.,0.));
    for (inode, ieqn) in &eqn_map {
        let o = &mesh.nodes[*inode].coords;
        for e in 0..*nelem {
            let e_id = &mesh.bodies[body_id-1].element_ids[e];
            let tri = Triangle::new(&mesh, *e_id);
            let enodes = &mesh.elements[*e_id-1].node_ids;
            let mut e_eqns = Vec::<usize>::new();
            for enode in enodes {
                match eqn_map.get(enode) {
                    Some(eqn) => e_eqns.push(*eqn),
                    None => println!("Eqn not found for element node {}", enode)
                }
            }
            let (he, ge) = tri.influence_matrices_at(k, o);
            for j in 0..e_eqns.len() {
                h[(*ieqn, e_eqns[j])] += he[j];
                g[(*ieqn, e_eqns[j])] += ge[j];
            }
        }
    }
    println!(" Complete!");

    let mut phi_inc = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    let mut phi_inc_fp = DVector::<Cplx>::from_element(num_fp, Cplx::new(0., 0.));
    {
        let fp = &user_input.field_points;
        let inc_wave = &user_input.incident_wave;
        let amp = Cplx::new(inc_wave.amplitude[0], inc_wave.amplitude[1]);
        let vector = &inc_wave.origin;
        let vec3 = Vector3::from_column_slice(vector);
        match inc_wave.wave_type {
            WaveType::PlaneWave => {
                for (inode, ieqn) in &eqn_map {
                    let coord = &mesh.nodes[*inode].coords;
                    phi_inc[*ieqn] = amp * Cplx::new(0., k * vec3.dot(coord));
                }
                for i in 0..num_fp {
                    let coord = Vector3::from_column_slice(&fp[i]);
                    phi_inc_fp[i] = amp * Cplx::new(0., k * vec3.dot(&coord));
                }
            }
            WaveType::SphericalWave => {
                for (inode, ieqn) in &eqn_map {
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
    }

    print!(" Solving system (direct LU)...");
    std::io::stdout().flush().unwrap();
    let mut phi = DVector::<Cplx>::from_element(num_eqn, Cplx::new(0., 0.));
    let mut vn = phi.clone();
    {
        let sbc = &user_input.surface_bc;
        match sbc.bc_type {
            BCType::Pressure => {
                let pressure_bc = Cplx::new(sbc.value[0], sbc.value[1]);
                phi.fill(pressure_bc);
                vn = phi_inc.clone();
                vn.gemv(Cplx::new(-1.0,0.0), &h, &phi, Cplx::new(1.0,0.0));
                //vn = phi_inc - h * phi; // rhs
                let glu = g.lu();
                glu.solve_mut(&mut vn);
            }   
            BCType::NormalVelocity => {
                let velocity_bc = Cplx::new(sbc.value[0], sbc.value[1]);
                vn.fill(velocity_bc);
                phi = phi_inc.clone();
                phi.gemv(Cplx::new(1.0,0.0), &g, &vn, Cplx::new(-1.0,0.0));
                //phi = g * vn - phi_inc; // rhs
                let hlu = h.lu();
                hlu.solve_mut(&mut phi);
            }
            BCType::Impedance => {
                let impedance_bc = Cplx::new(sbc.value[0], sbc.value[1]);
                let factor = Cplx::new(0., -omega * rho) / impedance_bc;
                for i in 0..h.nrows() {
                    for j in 0..h.ncols() {
                        h[(i,j)] += factor * g[(i,j)];
                    }
                }
                phi = -phi_inc.clone(); // rhs
                let hlu = h.lu();
                hlu.solve_mut(&mut phi);
            }
        }
    }
    println!(" Complete!");

    print!(" Post-processing...");
    std::io::stdout().flush().unwrap();
    let mut phi_fp = phi_inc_fp.clone();
    {
        let nfp = user_input.field_points.len();
        let mut m   = DMatrix::<Cplx>::from_diagonal_element(nfp, num_eqn,  Cplx::new(-1.,0.));
        let mut l = DMatrix::<Cplx>::from_element(nfp, num_eqn,  Cplx::new(0.,0.));
        for i in 0..nfp {
            let fieldpt = &user_input.field_points[i];
            let coord = Vector3::from_column_slice(fieldpt);
            for e in 0..*nelem {
                let e_id = &mesh.bodies[body_id-1].element_ids[e];
                let tri = Triangle::new(&mesh, *e_id);
                let enodes = &mesh.elements[*e_id-1].node_ids;
                let mut e_eqns = Vec::<usize>::new();
                for enode in enodes {
                    match eqn_map.get(enode) {
                        Some(eqn) => e_eqns.push(*eqn),
                        None => println!("Eqn not found for element node {}", enode)
                    }
                }
                let (me, le) = tri.influence_matrices_at(k, &coord);
                for j in 0..e_eqns.len() {
                    m[(i, e_eqns[j])] += me[j];
                    l[(i, e_eqns[j])] += le[j];
                }
            }
        }
        phi_fp.gemv(-Cplx::new(-1.0, 0.0), &m, &phi, Cplx::new(1.0, 0.0));
        phi_fp.gemv(-Cplx::new(1.0, 0.0), &l, &vn, Cplx::new(1.0, 0.0));
        //let phi_fp = - m * phi + l * vn + phi_inc_fp;
    }
    println!(" Complete!");

    println!(" Exiting...");
    return result;
}
