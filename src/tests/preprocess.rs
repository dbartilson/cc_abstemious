#[test]
fn json_reader() {
    use cc_abstemious::preprocess::input_data::*;
    let u = read_input_json("./src/tests/input_1.json").unwrap();
    assert_eq!(u.body_index, 3);
}

#[test]
fn vtk_reader() {
    use cc_abstemious::preprocess::*;
    use std::path::Path;

    let mut mesh: mesh_data::Mesh = Default::default();
    let _result = mesh.read_from_vtk(Path::new("./src/tests/sphere.vtk"));

    assert_eq!(mesh.elements.len(), 336);
}