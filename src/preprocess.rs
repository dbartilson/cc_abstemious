pub mod input_data;
pub mod mesh_data;

use std::collections::HashMap;
use std::{io::Write, path::Path};

pub fn preprocess() -> (input_data::UserInput, mesh_data::Mesh, HashMap<usize,usize>) {
    println!("=== CC-ABSTEMIOUS <=> BEM-ACOUSTICS ===");
    println!("Ver. 0.0");
    println!("");

    println!(" Current directory: {}", std::env::current_dir().unwrap().display());

    let args: Vec<String> = std::env::args().collect();
    let input_path_str = &args[1];
    println!(" Attempting to read input file: {}", input_path_str);

    // read input json
    let user_input = input_data::read_input_json(Path::new(input_path_str)).unwrap();

    // read mesh VTK
    let mut mesh: mesh_data::Mesh = Default::default();
    let _result = mesh.read_from_vtk(Path::new(&user_input.mesh_file));

    let body_id = &user_input.body_index;

    // preprocess to get node to eqn map
    print!(" Preprocessing...");
    std::io::stdout().flush().unwrap();
    let eqn_map = get_eqn_map(&mesh, *body_id);
    println!(" Complete!");

    return (user_input, mesh, eqn_map);
}

fn get_eqn_map(meshdata: &mesh_data::Mesh, body_id: usize) -> HashMap::<usize, usize> {
    let mut eqn_map = HashMap::<usize, usize>::new();

    let nnode = meshdata.nodes.len();
    let ibody = &meshdata.bodies[body_id-1];
    let mut eqn_vector = vec![false; nnode];
    // Mark each possible node (eqn) as used or not
    for element_id in &ibody.element_ids {
        let element = &meshdata.elements[*element_id-1];
        for node_id in &element.node_ids {
            eqn_vector[*node_id] = true;
        }
    }
    // Scroll through all used eqns in order and put in map
    let mut eqn_number: usize = 0;
    for i in 0..eqn_vector.len() {
        if eqn_vector[i] {
            eqn_map.insert(i, eqn_number);
            eqn_number += 1;
        }
    }
    return eqn_map;
}