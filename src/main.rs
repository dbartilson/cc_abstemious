fn main() {
    // exe to take in input .json file path as only argument, run, and output results at first freq
    let args: Vec<String> = std::env::args().collect();
    let input_path_str = &args[1];
    let mut analysis = cc_abstemious::Analysis::new();
    analysis.input_from_file(input_path_str);
    analysis.run();
    analysis.write_results_at_frequency(0)
}
