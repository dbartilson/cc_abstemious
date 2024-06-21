use std::io::Write;

use crate::preprocess;
use crate::incident_wave;
use crate::influence_matrix;
use crate::solve;
use crate::postprocess;

enum AnalysisState {
    Input,
    Solve,
    Null
}

pub struct Analysis {
    temp_input: Option<preprocess::input_data::UserInput>,
    predata: Option<preprocess::PreData>,
    analysis_state: AnalysisState,
    freq_index: usize,
    results: Vec<postprocess::FPResult>
}

impl Analysis {
    pub fn new() -> Analysis {
        println!("=== cc_abstemious <=> BEM-ACOUSTICS ===");
        println!("Ver. 0.0");
        println!("");
    
        println!(" Current directory: {}", std::env::current_dir().unwrap().display());
        Analysis {
            temp_input: None,
            predata: None,
            analysis_state: AnalysisState::Null,
            freq_index: 0,
            results: Vec::new()
        }
    }
    pub fn input_from_file(&mut self, input_path_str: &String) {
        println!(" Attempting to read input file: {}", input_path_str);
        self.temp_input = Some(preprocess::input_data::read_input_json(input_path_str).unwrap());
        self.analysis_state = AnalysisState::Input;
    }
    pub fn input_from_string(&mut self, input_str: &str) {
        self.temp_input = Some(preprocess::input_data::read_input_string(input_str).unwrap());
        self.analysis_state = AnalysisState::Input;
    }
    pub fn set_input(&mut self, input: preprocess::input_data::UserInput) {
        self.temp_input = Some(input);
        self.analysis_state = AnalysisState::Input;
    }
    pub fn run(&mut self) {

        if self.temp_input.is_none() {
            panic!("No input found");
        }

        // preprocess to get node to eqn map
        print!(" Preprocessing...");
        std::io::stdout().flush().unwrap();
        self.predata = Some(preprocess::preprocess(self.temp_input.take().unwrap()));
        let predata = self.predata.as_mut().unwrap();
        println!(" Complete!");

        let nfreq = predata.get_frequencies().len();
        let single_freq_flag = nfreq == 1;
        for i in 0..predata.get_frequencies().len() {
            predata.set_frequency_index(&i);
            let freq = predata.get_frequency();
            if !single_freq_flag {print!(" Analyzing frequency: {} ({} of {})...", freq, i+1, nfreq);std::io::stdout().flush().unwrap();}
            self.freq_index = i;

            if single_freq_flag {print!(" Assembling surface BEM influence matrices...");std::io::stdout().flush().unwrap();}
            let (h, g) = influence_matrix::get_surface_influence_matrices(predata);
            if single_freq_flag {println!(" Complete!");}
        
            let (phi_inc, phi_inc_fp) = incident_wave::get_incident_wave(predata);
    
            if single_freq_flag {print!(" Solving system (direct LU)...");std::io::stdout().flush().unwrap();}
            let (phi, vn) = solve::solve_lu(predata, &h, &g, &phi_inc);
            if single_freq_flag {println!(" Complete!");}
        
            if single_freq_flag {print!(" Post-processing...");std::io::stdout().flush().unwrap();}
            let (m, l) = influence_matrix::get_field_influence_matrices(predata);
            
            let phi_fp = solve::get_field(predata,&m, &l, &phi, &vn, &phi_inc_fp);

            let result = postprocess::FPResult{
                frequency: freq,
                scattered: Some(phi_fp),
                incident: Some(phi_inc_fp),
                radiated_power: 0.0
            };
            self.results.push(result);
            
            println!(" Complete!");std::io::stdout().flush().unwrap();
        }

        postprocess::convert_results(predata, &mut self.results);

        self.analysis_state = AnalysisState::Solve;
    }
    pub fn get_result(&self) -> &Vec<postprocess::FPResult> {
        return &self.results;
    }
    pub fn write_results_at_frequency(&self, ifreq: usize) {
        let _u = postprocess::write_results_at_frequency(self.predata.as_ref().unwrap(), &self.results[ifreq]);
    }
    pub fn write_results_at_point(&self, index: usize) {
        let _u = postprocess::write_results_at_point(self.predata.as_ref().unwrap(), &self.results, index);
    }

}