pub mod input_data;
pub mod mesh_data;

use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

pub struct PreData {
    input: input_data::UserInput,
    mesh: mesh_data::Mesh,
    eqn_map: HashMap<usize,usize>,
    ifreq: usize
}

impl PreData {
    pub fn set_frequency_index(&mut self, index: &usize) {self.ifreq = *index;}

    pub fn get_frequencies(&self) -> &Vec<f64> {return &self.input.frequency;}
    pub fn get_frequency(&self) -> f64 {return self.input.frequency[self.ifreq];}
    pub fn get_angular_frequency(&self) -> f64 {return 2.0 * PI * self.get_frequency();}
    pub fn get_wavenumber(&self) -> f64 {return self.get_angular_frequency() / self.get_sound_speed()}
    pub fn get_sound_speed(&self) -> f64 {return self.input.sound_speed;}
    pub fn get_mass_density(&self) -> f64 {return self.input.mass_density;}
    pub fn get_problem_type(&self) -> &input_data::ProblemType {return &self.input.problem_type;}
    pub fn get_mesh_body(&self) -> &mesh_data::Body {return &self.mesh.bodies[self.input.body_index - 1];}
    pub fn get_incident_wave(&self) -> &input_data::IncidentWaveInput {return &self.input.incident_wave;}
    pub fn get_surface_bc(&self) -> &input_data::SurfaceBoundaryCondition {return &self.input.surface_bc;}
    pub fn get_eqn_map(&self) -> &HashMap<usize, usize> {return &self.eqn_map;}
    pub fn get_mesh(&self) -> &mesh_data::Mesh {return &self.mesh;}
    pub fn get_num_eqn(&self) -> usize {return self.eqn_map.len();}
    pub fn get_output_filename(&self) -> &String {return &self.input.output.file;}
    pub fn get_output_field(&self) -> &input_data::OutputField {return &self.input.output.field;}
    pub fn get_output_type(&self) -> &input_data::OutputType {return &self.input.output.o_type;}
    pub fn get_field_points(&self) -> &Vec<[f64; 3]> {return &self.input.output.field_points;}
}

pub fn preprocess(input: input_data::UserInput) -> PreData {
    // read mesh VTK
    let mut mesh: mesh_data::Mesh = Default::default();
    let _result = mesh.read_from_vtk(Path::new(&input.mesh_file));

    let body_id = &input.body_index;

    // preprocess to get node to eqn map
    let eqn_map = get_eqn_map(&mesh, *body_id);

    // take ownership of input data
    return PreData{input, 
                   mesh: mesh, 
                   eqn_map: eqn_map, 
                   ifreq: 0};
}

pub fn get_eqn_map(meshdata: &mesh_data::Mesh, body_id: usize) -> HashMap::<usize, usize> {
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