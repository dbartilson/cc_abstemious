/*!
Defines the analysis object (struct) and analysis states
*/

use na::DVector;
use simplelog::*;
use std::fs::File;
use std::path::Path;

use crate::preprocess;
use crate::incident_wave;
use crate::solve;
use crate::postprocess;
use crate::Cplx;

/// Enumerate the analysis states for tracking 
enum AnalysisState {
    PostInput,
    PostSolve,
    Null
}

/// Struct containing velocity potential (phi) and normal velocity (vn), on surface
/// can be total, incident, or scattered
pub struct PrimaryVariables {
    pub phi: DVector::<Cplx>,
    pub vn: DVector::<Cplx>
}

///This contains all data related to the analysis and the main wrapper functions
/// 
/// Methods:
/// - create an Analysis
/// - input
/// - run
/// - get/write results
pub struct Analysis {
    temp_input: Option<preprocess::input::UserInput>,
    log_file: String,
    predata: Option<preprocess::PreData>,
    analysis_state: AnalysisState,
    freq_index: usize,
    results: Vec<postprocess::FPResult>
}

impl Default for Analysis {
    fn default() -> Self {
        Self::new()
    }
}

impl <'a>Analysis {
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
    /// Input from json file at path (string)
    pub fn input_from_file(&mut self, input_path_str: &String) {
        println!(" Attempting to read input file: {}", input_path_str);
        let path = Path::new(input_path_str);
        self.log_file = path.file_stem().unwrap().to_str().unwrap().to_string();
        self.temp_input = Some(preprocess::input::read_input_json(input_path_str).unwrap());
        self.analysis_state = AnalysisState::PostInput;
    }
    /// Input from json string (mostly for tests)
    pub fn input_from_string(&mut self, input_str: &str) {
        self.temp_input = Some(preprocess::input::read_input_string(input_str).unwrap());
        self.analysis_state = AnalysisState::PostInput;
    }
    /// Directly set an input_data (mostly for tests)
    pub fn set_input(&mut self, input: preprocess::input::UserInput) {
        self.temp_input = Some(input);
        self.analysis_state = AnalysisState::PostInput;
    }
    /// Run the analysis, writing out to the log file if created
    /// 
    /// Analysis stages:
    ///     - Preamble
    ///     - Loop over frequencies
    ///         - Solve for surface and field results
    ///         - Save results internally
    pub fn run(&'a mut self) {
        
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
                    TermLogger::new(LevelFilter::Info, Config::default(), TerminalMode::Mixed, ColorChoice::Auto),
                    WriteLogger::new(LevelFilter::Info, Config::default(), logfile),
                ]
            );
        }

        info!("=== cc_abstemious <=> BEM-ACOUSTICS ===");
        info!("Ver. {}.{}.{}", crate::VER_MAJOR, crate::VER_MINOR, crate::VER_SUBMINOR);
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

            let incident = incident_wave::get_incident_surface(predata);

            let total = solve::solve_for_surface(predata, &incident);

            let result = postprocess::postprocess(predata, &incident, &total);

            self.results.push(result);
        }

        self.analysis_state = AnalysisState::PostSolve;

        info!(" Complete!");
    }
    /// Return ref to field results from analysis
    pub fn get_results(&self) -> &Vec<postprocess::FPResult> { &self.results }
    /// Write results to output file at one frequency for all field points
    pub fn write_results_at_frequency(&self, ifreq: usize) {
        let _u = postprocess::write_results_at_frequency(self.predata.as_ref().unwrap(), None, &self.results[ifreq]);
    }
    /// Write results to output file for all field points at one frequency
    pub fn write_results_at_point(&self, index: usize) {
        let _u = postprocess::write_results_at_point(self.predata.as_ref().unwrap(), &self.results, index);
    }
    /// Write results to output file for all field points at one frequency
    pub fn write_power(&self) {
        let _u = postprocess::write_power(self.predata.as_ref().unwrap(), None, &self.results);
    }
}