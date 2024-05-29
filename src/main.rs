fn main() {
    let args: Vec<String> = std::env::args().collect();
    let input_path_str = &args[1];
    cc_abstemious::run(input_path_str);
}
