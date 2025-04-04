#[macro_use]
extern crate approx;

use std::f64::consts::PI;

use cc_abstemious::preprocess::input_data::*;

fn default_input() -> UserInput {
    UserInput {
        mesh_file: "./src/tests/sphere.vtk".to_string(),
        body_index: 3,
        frequency: FrequencyInput::List { values: vec![100.0] },
        sound_speed: 1.0,
        mass_density: 1.0,
        problem_type: ProblemType::Exterior,
        method_type: MethodType::Classical,
        solver: Solver::Direct {  },
        incident_wave: IncidentWaveInput {
            origin: [0.0, 0.0, 0.0],
            wave_type: WaveType::SphericalWave,
            amplitude: [0.0, 0.0]
        },
        surface_bc: SurfaceBoundaryCondition {
            bc_type: BCType::NormalVelocity,
            value: [0.0, 0.0]
        },
        output: Output {
            o_type: OutputType::Scattered,
            field: OutputField::Pressure,
            field_points: Vec::new(),
            file: "".to_string()
        }
    }
}

#[allow(dead_code)]
fn rigid_sphere_plane_wave_ring() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![10.0] };
    // water
    input.sound_speed = 1500.0;
    input.mass_density = 1000.0;
    // incident wave
    input.incident_wave.origin = [1.0, 0.0, 0.0];
    input.incident_wave.wave_type = WaveType::PlaneWave;
    input.incident_wave.amplitude = [1.0, 0.0];
    // output fule
    input.output.file = "./src/tests/rigid_sphere_plane_wave_bem.csv".to_string();
    // set up field points (ring in XY plane)
    let num_fp = 100;
    let radius = 10.0;
    for i in 0..num_fp {
        let theta = 2.0 * PI * (i as f64) / (num_fp as f64);
        let x = radius * f64::cos(theta);
        let y = radius * f64::sin(theta);
        input.output.field_points.push([x, y, 0.0]);
    }

    analysis.set_input(input);
    analysis.run();
    analysis.write_results_at_frequency(0);
}

//#[allow(dead_code)]
#[test]
fn rigid_sphere_plane_wave_sweep() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    //input.mesh_file = "./src/tests/refined_sphere.vtk".to_string();
    input.method_type = MethodType::BurtonMiller;
    input.frequency = FrequencyInput::LinearSpaced { start: 10.0, end: 1000.0, number: 50 };
    // water
    input.sound_speed = 1500.0;
    input.mass_density = 1000.0;
    // incident wave
    input.incident_wave.origin = [1.0, 0.0, 0.0];
    input.incident_wave.wave_type = WaveType::PlaneWave;
    input.incident_wave.amplitude = [1.0, 0.0];
    // output fule
    input.output.file = "./src/tests/rigid_sphere_plane_wave_bem.csv".to_string();
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    analysis.write_results_at_point(0);
    let _fp = analysis.get_result();
}

#[test]
fn rigid_sphere_plane_wave() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    // water
    input.sound_speed = 1500.0;
    input.mass_density = 1000.0;
    // incident wave
    input.incident_wave.origin = [1.0, 0.0, 0.0];
    input.incident_wave.wave_type = WaveType::PlaneWave;
    input.incident_wave.amplitude = [1.0, 0.0];
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_result();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.re, -0.000039380091293979903, epsilon = 1.0e-10);
    assert_relative_eq!(fpi.im, 0.00047605686656319906, epsilon = 1.0e-10);
}

#[test]
fn rigid_sphere_plane_wave_iterative() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    // water
    input.sound_speed = 1500.0;
    input.mass_density = 1000.0;
    // incident wave
    input.incident_wave.origin = [1.0, 0.0, 0.0];
    input.incident_wave.wave_type = WaveType::PlaneWave;
    input.incident_wave.amplitude = [1.0, 0.0];
    input.solver = Solver::Iterative { max_iterations: 1000, tolerance: 1.0e-5 };
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_result();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.re, -0.000039380091293979903, epsilon = 1.0e-8);
    assert_relative_eq!(fpi.im, 0.00047605686656319906, epsilon = 1.0e-8);
}

#[test]
fn rigid_sphere_plane_wave_hmatrix() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    // water
    input.sound_speed = 1500.0;
    input.mass_density = 1000.0;
    // incident wave
    input.incident_wave.origin = [1.0, 0.0, 0.0];
    input.incident_wave.wave_type = WaveType::PlaneWave;
    input.incident_wave.amplitude = [1.0, 0.0];
    input.solver = Solver::Hierarchical { max_iterations: 1000, tolerance: 1.0e-5 };
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_result();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.re, -3.9380308176032198e-05, epsilon = 1.0e-10);
    assert_relative_eq!(fpi.im, 0.00047605563272935873, epsilon = 1.0e-10);
}

#[test]
fn rigid_sphere_plane_wave_burton_miller() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.method_type = MethodType::BurtonMiller;
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    // water
    input.sound_speed = 1500.0;
    input.mass_density = 1000.0;
    // incident wave
    input.incident_wave.origin = [1.0, 0.0, 0.0];
    input.incident_wave.wave_type = WaveType::PlaneWave;
    input.incident_wave.amplitude = [1.0, 0.0];
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_result();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.re, -0.000039380091293979903, epsilon = 1.0e-10);
    assert_relative_eq!(fpi.im, 0.00047605686656319906, epsilon = 1.0e-10);
}