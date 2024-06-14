use std::f64::consts::PI;

use cc_abstemious::preprocess::input_data::*;

fn default_input() -> UserInput {
    UserInput {
        mesh_file: "./src/tests/refined_sphere.vtk".to_string(),
        body_index: 3,
        frequency: 1.0,
        sound_speed: 1.0,
        mass_density: 1.0,
        problem_type: ProblemType::Exterior,
        field_points: Vec::new(),
        incident_wave: IncidentWaveInput {
            origin: [0.0, 0.0, 0.0],
            wave_type: WaveType::SphericalWave,
            amplitude: [0.0, 0.0]
        },
        surface_bc: SurfaceBoundaryCondition {
            bc_type: BCType::NormalVelocity,
            value: [0.0, 0.0]
        },
        output_file: "".to_string()
    }
}

#[test]
fn rigid_sphere_plane_wave() {
    let mut analysis = cc_abstemious::Analysis::new();
    let mut input = default_input();
    input.frequency = 100.0;
    // water
    input.sound_speed = 1500.0;
    input.mass_density = 1000.0;
    // incident wave
    input.incident_wave.origin = [1.0, 0.0, 0.0];
    input.incident_wave.wave_type = WaveType::PlaneWave;
    input.incident_wave.amplitude = [1.0, 0.0];
    // output fule
    input.output_file = "./src/tests/rigid_sphere_plane_wave_bem.csv".to_string();
    // set up field points (ring in XY plane)
    let num_fp = 100;
    let radius = 10.0;
    for i in 0..num_fp {
        let theta = 2.0 * PI * (i as f64) / (num_fp as f64);
        let x = radius * f64::cos(theta);
        let y = radius * f64::sin(theta);
        input.field_points.push([x, y, 0.0]);
    }

    analysis.set_input(input);
    analysis.run();
    analysis.write_fp_result();
    let _fp = analysis.get_fp_result();
}