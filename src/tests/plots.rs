use cc_abstemious::preprocess::input::*;
use std::f64::consts::PI;

use crate::default_input;

#[allow(dead_code)]
//#[test]
fn rigid_sphere_plane_wave_ring() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = FrequencyInput::List { values: vec![10.0] };
    // incident wave
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0],
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
    input.method_type = MethodType::Classical;
    input.frequency = FrequencyInput::LinearSpaced {
        start: 10.0,
        end: 1000.0,
        number: 50,
    };
    // incident wave
    input.incident_wave = vec![IncidentWaveInput::PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude: [1.0, 0.0],
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

#[allow(dead_code)]
//#[test]
fn monopole_power_sweep() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    // output file
    input.output.file = "./src/tests/pulsating_sphere.csv".to_string();
    input.frequency = FrequencyInput::LinearSpaced {
        start: 10.0,
        end: 1000.0,
        number: 50,
    };
    input.surface_bc = SurfaceBoundaryCondition {
        bc_type: BCType::NormalVelocity,
        value: [1.0, 0.0],
    };

    analysis.set_input(input);
    analysis.run();
    analysis.write_results();
}

#[test]
fn speed_test() {
    for i in 8..9 {
        let mut analysis = cc_abstemious::Analysis::new();
        let mut input = default_input();
        input.mesh_file = format!("./src/tests/sphere_{}.vtk", i+1);
        input.solver = Solver::Hierarchical {
            max_iterations: 1000,
            tolerance: 1.0e-5,
        };
        input.incident_wave = vec![IncidentWaveInput::PlaneWave {
            direction: [1.0, 0.0, 0.0],
            amplitude: [1.0, 0.0],
        }];
        analysis.set_input(input);
        analysis.run();
    }
}
