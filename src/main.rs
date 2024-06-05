fn main() {
    let args: Vec<String> = std::env::args().collect();
    let input_path_str = &args[1];
    let mut analysis = cc_abstemious::Analysis::new();
    analysis.input_from_file(input_path_str);
    analysis.run();
    let _fp = analysis.get_fp_result();
}
