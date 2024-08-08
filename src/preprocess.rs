pub mod input_data;
pub mod mesh_data;

use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

use crate::elements::NIElement;
use crate::Cplx;

pub fn get_num_threads() -> usize {
    match std::thread::available_parallelism() {
        Ok(result) => std::cmp::max(result.get() / 2, 2),
        Err(_) => 2
    }
}

pub struct PreData {
    input: input_data::UserInput,
    mesh: mesh_data::Mesh,
    eqn_map: HashMap<usize,usize>, 
    node_map: HashMap<usize,usize>,
    revcon: Vec<Vec<usize>>,
    ifreq: usize // current frequency index
}

impl PreData {
    pub fn set_frequency_index(&mut self, index: &usize) {self.ifreq = *index;}

    #[inline]
    pub fn get_frequencies(&self) -> &Vec<f64> {return &self.input.frequency;}
    #[inline]
    pub fn get_frequency(&self) -> f64 {return self.input.frequency[self.ifreq];}
    #[inline]
    pub fn get_angular_frequency(&self) -> f64 {return 2.0 * PI * self.get_frequency();}
    #[inline]
    pub fn get_wavenumber(&self) -> f64 {return self.get_angular_frequency() / self.get_sound_speed()}
    #[inline]
    pub fn get_sound_speed(&self) -> f64 {return self.input.sound_speed;}
    #[inline]
    pub fn get_mass_density(&self) -> f64 {return self.input.mass_density;}
    #[inline]
    pub fn get_problem_type(&self) -> &input_data::ProblemType {return &self.input.problem_type;}
    #[inline]
    pub fn get_hdiag(&self) -> Cplx {
        match self.get_problem_type() {
            // the H matrix has -1/2 added along the diagonal for exterior problems
            input_data::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
            input_data::ProblemType::Interior => Cplx::new(0.0, 0.0)
        }
    }
    pub fn get_mesh_body(&self) -> &mesh_data::Body {return &self.mesh.bodies[self.input.body_index - 1];}
    pub fn get_solver_type(&self) -> &input_data::SolverType {return &self.input.solver.s_type;}
    pub fn get_solver_tolerance(&self) -> f64 {return self.input.solver.tolerance;}
    pub fn get_solver_max_it(&self) -> usize {return self.input.solver.max_iterations;}
    pub fn get_incident_wave(&self) -> &input_data::IncidentWaveInput {return &self.input.incident_wave;}
    pub fn get_surface_bc(&self) -> &input_data::SurfaceBoundaryCondition {return &self.input.surface_bc;}
    /// get map from node index to equation index
    pub fn get_eqn_map(&self) -> &HashMap<usize, usize> {return &self.eqn_map;}
    /// get map from equation index to node index
    pub fn get_node_map(&self) -> &HashMap<usize, usize> {return &self.node_map;}
    /// get list of elements at each node
    pub fn get_revcon(&self) -> &Vec<Vec<usize>> {return &self.revcon;}
    pub fn get_mesh(&self) -> &mesh_data::Mesh {return &self.mesh;}
    pub fn get_num_eqn(&self) -> usize {return self.eqn_map.len();}
    pub fn get_output_filename(&self) -> &String {return &self.input.output.file;}
    pub fn get_output_field(&self) -> &input_data::OutputField {return &self.input.output.field;}
    pub fn get_output_type(&self) -> &input_data::OutputType {return &self.input.output.o_type;}
    pub fn get_field_points(&self) -> &Vec<[f64; 3]> {return &self.input.output.field_points;}
}

pub fn preprocess(input: input_data::UserInput) -> PreData {
    
    info!(" Preprocessing...");
    // read mesh VTK
    let mut mesh: mesh_data::Mesh = Default::default();
    let _result = mesh.read_from_vtk(Path::new(&input.mesh_file));

    let body_id = &input.body_index;

    process_node_normals(&mut mesh, *body_id);

    // preprocess to get node to eqn map
    let (eqn_map, node_map, revcon) = get_eqn_map(&mut mesh, *body_id);

    // take ownership of input data
    return PreData{input, 
                   mesh: mesh, 
                   eqn_map: eqn_map, 
                   node_map: node_map,
                   revcon: revcon,
                   ifreq: 0};
}

fn get_eqn_map(mesh: &mut mesh_data::Mesh, body_id: usize) 
    -> (HashMap::<usize, usize>, 
        HashMap::<usize, usize>, 
        Vec<Vec<usize>>) {
    // create a map from node index to eqn index
    let nnode = mesh.nodes.len();
    let ibody = &mesh.bodies[body_id-1];
    let mut revcon = vec![Vec::<usize>::new(); nnode];
    // Push the elements that are at each node
    for element_id in &ibody.element_ids {
        let element = &mesh.elements[*element_id];
        for node_id in &element.node_ids {
            revcon[*node_id].push(*element_id);
        }
    }
    // Scroll through all used eqns in order and put in map
    let mut eqn_map = HashMap::<usize, usize>::new();
    let mut node_map = HashMap::<usize, usize>::new();
    let mut eqn_number: usize = 0;
    for i in 0..revcon.len() {
        if !revcon[i].is_empty() {
            eqn_map.insert(i, eqn_number);
            node_map.insert(eqn_number, i);
            eqn_number += 1;
        }
    }
    // Add equation numbers to elements
    for element_id in &ibody.element_ids {
        let element = &mut mesh.elements[*element_id];
        for node_id in &element.node_ids { 
            if let Some(eqn) = eqn_map.get(node_id) {
                element.eqn_idx.push(*eqn);
            }
        }
    }
    return (eqn_map, node_map, revcon);
}

/// Calculate node normals using surface elements
fn process_node_normals(mesh: &mut mesh_data::Mesh, body_id: usize) {
    let nnode = mesh.nodes.len();
    let mut normals =  vec![na::Vector3::<f64>::new(0.0, 0.0, 0.0); nnode];
    let ibody = &mesh.bodies[body_id-1];
    for element_id in &ibody.element_ids {
        let element = NIElement::new(&mesh, *element_id);
        // get nodes at this element, sum the normal into the normal at each node
        for i in 0..element.get_num_nodes() {
            let inode = &mesh.elements[*element_id].node_ids[i];
            let gp = element.natural_coordinates_at_node(i);
            let normal = element.normal_vector_at_gp(&gp);
            // sum together non-normalized normals, this will weight the average normal towards bigger elements
            normals[*inode] += normal;
        }
    }
    // save to mesh data, normalize now
    for i in 0..nnode {
        mesh.nodes[i].normal = normals[i].normalize();
    }
}
