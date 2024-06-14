use std::{io::Write, path::Path};

use preprocess::input_data::{self, UserInput};

extern crate nalgebra as na;
type Cplx = na::Complex<f64>;

pub mod preprocess;
pub mod elements;
pub mod incident_wave;
pub mod influence_matrix;
pub mod solve;
pub mod postprocess;

enum AnalysisState {
    Input,
    Solve,
    Null
}

pub struct Analysis {
    input: Option<preprocess::input_data::UserInput>,
    analysis_state: AnalysisState,
    phi_fp: Option<na::DVector<Cplx>>,
    phi_fp_inc: Option<na::DVector<Cplx>>
}

impl Analysis {
    pub fn new() -> Analysis {
        println!("=== cc_abstemious <=> BEM-ACOUSTICS ===");
        println!("Ver. 0.0");
        println!("");
    
        println!(" Current directory: {}", std::env::current_dir().unwrap().display());
        Analysis {
            input: None,
            analysis_state: AnalysisState::Null,
            phi_fp: None,
            phi_fp_inc: None
        }
    }
    pub fn input_from_file(&mut self, input_path_str: &String) {
        println!(" Attempting to read input file: {}", input_path_str);
        self.input = Some(input_data::read_input_json(input_path_str).unwrap());
        self.analysis_state = AnalysisState::Input;
    }
    pub fn input_from_string(&mut self, input_str: &str) {
        self.input = Some(input_data::read_input_string(input_str).unwrap());
        self.analysis_state = AnalysisState::Input;
    }
    pub fn set_input(&mut self, input: UserInput) {
        self.input = Some(input);
        self.analysis_state = AnalysisState::Input;
    }
    pub fn run(&mut self) {

        if self.input.is_none() {
            panic!("No input found");
        }
        let input = self.input.as_ref().unwrap();
        // read mesh VTK
        let mut mesh: preprocess::mesh_data::Mesh = Default::default();
        let _result = mesh.read_from_vtk(Path::new(&input.mesh_file));

        let body_id = &input.body_index;

        // preprocess to get node to eqn map
        print!(" Preprocessing...");
        std::io::stdout().flush().unwrap();
        let eqn_map = preprocess::get_eqn_map(&mesh, *body_id);
        println!(" Complete!");

        print!(" Assembling surface BEM influence matrices...");
        std::io::stdout().flush().unwrap();
        let (h, g) = influence_matrix::get_surface_influence_matrices(&input, &mesh, &eqn_map);
        println!(" Complete!");
    
        let (phi_inc, phi_inc_fp) = incident_wave::get_incident_wave(&input, &mesh, &eqn_map);

        print!(" Solving system (direct LU)...");
        std::io::stdout().flush().unwrap();
        let (phi, vn) = solve::solve_lu(&input, &h, &g, &phi_inc, &eqn_map.len());
        println!(" Complete!");
    
        print!(" Post-processing...");
        std::io::stdout().flush().unwrap();
        let (m, l) = influence_matrix::get_field_influence_matrices(&input, &mesh, &eqn_map);
        self.phi_fp = Some(solve::get_fp(&m, &l, &phi, &vn, &phi_inc_fp));
        println!(" Complete!");

        self.phi_fp_inc = Some(phi_inc_fp);

        self.analysis_state = AnalysisState::Solve;
    }
    pub fn get_fp_result(&self) -> na::DVector<Cplx> {
        match &self.phi_fp {
            Some(result) => {
                result.clone()
            },
            None => panic!("No FP")
        }
    }
    pub fn write_fp_result(&self) {
        let _result = postprocess::write_fp_csv(&self.input.as_ref().unwrap().output_file, 
            self.phi_fp.as_ref().unwrap(), 
            self.phi_fp_inc.as_ref().unwrap(),
            &self.input.as_ref().unwrap().field_points);
    }

}