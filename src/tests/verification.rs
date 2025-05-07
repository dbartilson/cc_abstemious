#[macro_use]
extern crate approx;

use std::{f64::consts::PI, path::Path};

use cc_abstemious::preprocess::input::*;

fn default_input() -> UserInput {
    UserInput {
        mesh_file: "./src/tests/sphere.vtk".to_string(),
        body_index: 3,
        frequency: FrequencyInput::List { values: vec![100.0] },
        // water
        sound_speed: 1500.0,
        mass_density: 1000.0,
        problem_type: ProblemType::Exterior,
        method_type: MethodType::Classical,
        solver: Solver::Direct {  },
        incident_wave: Vec::new(),
        surface_bc: SurfaceBoundaryCondition {
            bc_type: BCType::NormalVelocity,
            value: [0.0, 0.0]
        },
        output: Output {
            o_type: OutputType::Scattered,
            field: OutputField::Pressure,
            field_points: Vec::new(),
            file: "".to_string(),
            request_power: false
        }
    }
}

#[allow(dead_code)]
fn rigid_sphere_plane_wave_ring() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![10.0] };
    // incident wave
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0]
    }];
    // output file
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
    analysis.write_results();
}

#[allow(dead_code)]
//#[test]
fn rigid_sphere_plane_wave_sweep() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    //input.mesh_file = "./src/tests/refined_sphere.vtk".to_string();
    input.method_type = MethodType::Classical;
    input.frequency = FrequencyInput::LinearSpaced { start: 10.0, end: 1000.0, number: 50 };
    // incident wave
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0]
    }];
    // output file
    input.output.file = "./src/tests/rigid_sphere_plane_wave_bem.csv".to_string();
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    analysis.write_results();
    let _fp = analysis.get_results();
}

#[test]
fn rigid_sphere_plane_wave() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0]
    }];
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_results();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.norm(), 0.00045892241574703284, epsilon = 1.0e-10);
    assert_relative_eq!(fpi.re, -0.000038592340803747317, epsilon = 1.0e-10);
    assert_relative_eq!(fpi.im, 0.00045729685643614465, epsilon = 1.0e-10);
}

#[test]
fn rigid_sphere_plane_wave_iterative() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    // incident wave
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0]
    }];
    input.solver = Solver::Iterative { max_iterations: 1000, tolerance: 1.0e-5 };
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_results();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.re, -0.000038592340803747317, epsilon = 1.0e-8);
    assert_relative_eq!(fpi.im, 0.00045729685643614465, epsilon = 1.0e-8);
}

#[test]
fn rigid_sphere_plane_wave_hmatrix() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    // incident wave
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0]
    }];
    input.solver = Solver::Hierarchical { max_iterations: 1000, tolerance: 1.0e-5 };
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_results();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.re, -0.000038592340803747317, epsilon = 1.0e-8);
    assert_relative_eq!(fpi.im, 0.00045729685643614465, epsilon = 1.0e-8);
}

#[test]
fn rigid_sphere_plane_wave_burton_miller() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.method_type = MethodType::BurtonMiller;
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    // incident wave
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0]
    }];
    let radius = 10.0;
    let theta = 0.0;
    let x = radius * f64::cos(theta);
    let y = radius * f64::sin(theta);
    input.output.field_points.push([x, y, 0.0]);

    analysis.set_input(input);
    analysis.run();
    let fp = analysis.get_results();
    let fpi = fp[0].scattered.as_ref().unwrap()[0];
    assert_relative_eq!(fpi.norm(),  0.0005407404080455677, epsilon = 1.0e-10);
    assert_relative_eq!(fpi.re, -0.00021387485159944078, epsilon = 1.0e-10);
    assert_relative_eq!(fpi.im, 0.000_496_646_490_722_127_4, epsilon = 1.0e-10);
}

#[test]
fn monopole_power() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    input.surface_bc = SurfaceBoundaryCondition {
        bc_type: BCType::NormalVelocity,
        value: [1.0, 0.0]
    };

    analysis.set_input(input);
    analysis.run();
    let power = analysis.get_results()[0].power.scattered;
    assert_relative_eq!(power, 95617.4138092842, epsilon = 1.0e-10);
}

#[allow(dead_code)]
//#[test]
fn monopole_power_sweep() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    // output file
    input.output.file = "./src/tests/pulsating_sphere.csv".to_string();
    input.frequency = FrequencyInput::LinearSpaced { start: 10.0, end: 1000.0, number: 50 };
    input.surface_bc = SurfaceBoundaryCondition {
        bc_type: BCType::NormalVelocity,
        value: [1.0, 0.0]
    };

    analysis.set_input(input);
    analysis.run();
    analysis.write_results();
}

#[test]
fn output_at_one_point() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::LinearSpaced { start: 10.0, end: 1000.0, number: 2 };
    input.surface_bc = SurfaceBoundaryCondition {
        bc_type: BCType::NormalVelocity,
        value: [1.0, 0.0]
    };
    input.solver = Solver::Iterative { max_iterations: 1000, tolerance: 1.0e-5 };
    input.output.field_points.push([10.0, 0.0, 0.0]);
    let ofile = "./src/tests/output_at_one_point.csv".to_string();
    input.output.file = ofile.clone();

    analysis.set_input(input);
    analysis.run();
    analysis.write_results();

    let path = Path::new(&ofile);
    assert!(path.exists());
    let _ = std::fs::remove_file(path);
}

#[test]
fn output_power() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![100.0] };
    input.surface_bc = SurfaceBoundaryCondition {
        bc_type: BCType::NormalVelocity,
        value: [1.0, 0.0]
    };
    input.solver = Solver::Iterative { max_iterations: 1000, tolerance: 1.0e-5 };
    input.output.request_power = true;
    input.output.file = "./src/tests/output_power.csv".to_string();

    analysis.set_input(input);
    analysis.run();
    analysis.write_results();

    let ofiles = ["./src/tests/output_power_power.csv".to_string(),
                               "./src/tests/output_power_index.csv".to_string(),
                               "./src/tests/output_power_1.csv".to_string()];
    for ofile in ofiles {
        let path = Path::new(&ofile);
        assert!(path.exists());
        let _ = std::fs::remove_file(path);
    }
}

#[test]
fn output() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::LinearSpaced { start: 10.0, end: 1000.0, number: 2 };
    input.surface_bc = SurfaceBoundaryCondition {
        bc_type: BCType::NormalVelocity,
        value: [1.0, 0.0]
    };
    input.solver = Solver::Iterative { max_iterations: 1000, tolerance: 1.0e-5 };
    input.output.field_points.push([10.0, 0.0, 0.0]);
    input.output.field_points.push([-10.0, 0.0, 0.0]);
    input.output.request_power = true;
    input.output.file = "./src/tests/output.csv".to_string();

    analysis.set_input(input);
    analysis.run();
    analysis.write_results();

    let ofiles = ["./src/tests/output_power.csv".to_string(),
                               "./src/tests/output_index.csv".to_string(),
                               "./src/tests/output_1.csv".to_string(),
                               "./src/tests/output_2.csv".to_string()];
    for ofile in ofiles {
        let path = Path::new(&ofile);
        assert!(path.exists());
        let _ = std::fs::remove_file(path);
    }
}