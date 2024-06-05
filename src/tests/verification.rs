#[test]
fn rigid_sphere_plane_wave() {
    let loc = "./src/tests/input_1.json";
    let mut analysis = cc_abstemious::Analysis::new();
    analysis.input_from_file(&loc.to_string());
    analysis.run();
    let _fp = analysis.get_fp_result();
}