extern crate simplelog;

use simplelog::*;
use std::fs::File;
use std::path::Path;

use crate::preprocess;
use crate::incident_wave;
use crate::influence_matrix;
use crate::solve;
use crate::postprocess;

const VER_MAJOR: usize = 1;
const VER_MINOR: usize = 0; 

enum AnalysisState {
    Input,
    Solve,
    Null
}

pub struct Analysis {
    temp_input: Option<preprocess::input_data::UserInput>,
    log_file: String,
    predata: Option<preprocess::PreData>,
    analysis_state: AnalysisState,
    freq_index: usize,
    results: Vec<postprocess::FPResult>
}

impl Analysis {
    pub fn new() -> Analysis {
        Analysis {
            temp_input: None,
            log_file: "".to_string(),
            predata: None,
            analysis_state: AnalysisState::Null,
            freq_index: 0,
            results: Vec::new()
        }
    }
    pub fn input_from_file(&mut self, input_path_str: &String) {
        println!(" Attempting to read input file: {}", input_path_str);
        let path = Path::new(input_path_str);
        self.log_file = path.file_stem().unwrap().to_str().unwrap().to_string();
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

        // set up logger, if no input file, just output to stdout
        if self.log_file.is_empty() {
            let _ = TermLogger::init(LevelFilter::Info, Config::default(), TerminalMode::Mixed, ColorChoice::Auto);
        }
        else {
            let logfile = File::create(format!("{}{}",self.log_file,".log")).unwrap();
            let _ = CombinedLogger::init(
                vec![
                    TermLogger::new(LevelFilter::Warn, Config::default(), TerminalMode::Mixed, ColorChoice::Auto),
                    WriteLogger::new(LevelFilter::Info, Config::default(), logfile),
                ]
            );
        }

        info!("=== cc_abstemious <=> BEM-ACOUSTICS ===");
        info!("Ver. {}.{}", VER_MAJOR, VER_MINOR);
        info!(" Current directory: {}", std::env::current_dir().unwrap().display());
        info!(" Starting analysis... (see log file: {}{})", self.log_file,".log");

        // preprocess
        self.predata = Some(preprocess::preprocess(self.temp_input.take().unwrap()));
        let predata = self.predata.as_mut().unwrap();

        let nfreq = predata.get_frequencies().len();
        for i in 0..nfreq {
            // set up frequency index
            predata.set_frequency_index(&i);
            let freq = predata.get_frequency();
            info!(" Analyzing frequency: {} ({} of {})...", freq, i+1, nfreq);
            self.freq_index = i;

            let (phi_inc, phi_inc_fp) = incident_wave::get_incident_wave(predata);

            let (h, g) = influence_matrix::get_surface_influence_matrices(predata);
        
            let (phi, vn) = solve::get_surface(predata, &h, &g, &phi_inc);
        
            let (m, l) = influence_matrix::get_field_influence_matrices(predata);
            
            let phi_fp = solve::get_field(predata, &m, &l, &phi, &vn, &phi_inc_fp);

            let result = postprocess::FPResult{
                frequency: freq,
                scattered: Some(phi_fp),
                incident: Some(phi_inc_fp),
                radiated_power: 0.0
            };
            self.results.push(result);
        }

        postprocess::convert_results(predata, &mut self.results);

        self.analysis_state = AnalysisState::Solve;

        info!(" Complete!");
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