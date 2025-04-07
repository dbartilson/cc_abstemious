pub mod input_data;
pub mod mesh_data;

use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

use crate::elements::NIElement;
use crate::tools;
use crate::Cplx;

pub fn get_num_threads() -> usize {
    match std::thread::available_parallelism() {
        Ok(result) => std::cmp::max(result.get() / 2, 2),
        Err(_) => 2
    }
}

pub struct Hypersingular {
    pub is: bool,
    pub factor: Cplx
}

pub struct PreData {
    input: input_data::UserInput,
    mesh: mesh_data::Mesh,
    eqn_map: HashMap<usize,usize>, 
    node_map: HashMap<usize,usize>,
    frequency_list: Vec<f64>,
    ifreq: usize // current frequency index
}

impl PreData {
    pub fn set_frequency_index(&mut self, index: &usize) {self.ifreq = *index;}

    #[inline]
    pub fn get_frequencies(&self) -> &Vec<f64> {return &self.frequency_list;}
    #[inline]
    pub fn get_frequency(&self) -> f64 {return self.frequency_list[self.ifreq];}
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
    pub fn get_method_type(&self) -> &input_data::MethodType {return &self.input.method_type;}
    #[inline]
    pub fn get_hypersingular(&self) -> Hypersingular {
        return Hypersingular { is: self.use_hypersingular(), factor: self.get_burton_miller_coupling_factor()} }
    #[inline]
    fn use_hypersingular(&self) -> bool {return self.input.method_type == input_data::MethodType::BurtonMiller;}
    #[inline]
    fn get_burton_miller_coupling_factor(&self) -> Cplx {return Cplx::new(0.0, 1.0 / self.get_wavenumber())}
    #[inline]
    pub fn get_hdiag(&self) -> Cplx {
        match self.get_problem_type() {
            // the H matrix has -1/2 added along the diagonal for exterior problems
            input_data::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
            input_data::ProblemType::Interior => Cplx::new(0.0, 0.0)
        }
    }
    #[inline]
    pub fn get_gdiag(&self) -> Cplx {
        match self.get_method_type() {
            // the G matrix has 1/2 (ultimately beta/2, where beta = i/k) added along the diagonal for Burton-Miller formulation
            input_data::MethodType::Classical => Cplx::new(0.0, 0.0),
            input_data::MethodType::BurtonMiller => 0.5 * self.get_burton_miller_coupling_factor()
        }
    }
    #[inline]
    pub fn get_mesh_body(&self) -> &mesh_data::Body {return &self.mesh.bodies[self.input.body_index - 1];}
    #[inline]
    pub fn get_solver(&self) -> &input_data::Solver {return &self.input.solver}
    #[inline]
    pub fn get_incident_wave(&self) -> &input_data::IncidentWaveInput {return &self.input.incident_wave;}
    #[inline]
    pub fn get_surface_bc(&self) -> &input_data::SurfaceBoundaryCondition {return &self.input.surface_bc;}
    #[inline]
    /// get map from node index to equation index
    pub fn get_eqn_map(&self) -> &HashMap<usize, usize> {return &self.eqn_map;}
    #[inline]
    /// get map from equation index to node index
    pub fn get_node_map(&self) -> &HashMap<usize, usize> {return &self.node_map;}
    #[inline]
    pub fn get_mesh(&self) -> &mesh_data::Mesh {return &self.mesh;}
    #[inline]
    pub fn get_cpts(&self) -> &Vec<mesh_data::CollocationPoint> {return &self.mesh.cpts}
    #[inline]
    pub fn get_num_eqn(&self) -> usize {return self.mesh.cpts.len();}
    pub fn get_output_filename(&self) -> &String {return &self.input.output.file;}
    #[inline]
    pub fn get_output_field(&self) -> &input_data::OutputField {return &self.input.output.field;}
    #[inline]
    pub fn get_output_type(&self) -> &input_data::OutputType {return &self.input.output.o_type;}
    #[inline]
    pub fn get_field_points(&self) -> &Vec<[f64; 3]> {return &self.input.output.field_points;}
}

pub fn preprocess(input: input_data::UserInput) -> PreData {
    
    info!(" Preprocessing...");
    // read mesh VTK
    let mut mesh: mesh_data::Mesh = Default::default();
    let _result = mesh.read_from_vtk(Path::new(&input.mesh_file));

    let body_id = &input.body_index;

    let frequency_list = process_frequency_list(&input.frequency);

    let (eqn_map, node_map) = process_collocation_pts(&mut mesh, *body_id);

    // take ownership of input data
    return PreData{input, 
                   mesh: mesh, 
                   eqn_map: eqn_map, 
                   node_map: node_map,
                   frequency_list,
                   ifreq: 0};
}
/// Calculate the collocation points and normals using surface elements
fn process_collocation_pts(mesh: &mut mesh_data::Mesh, body_id: usize) -> 
    (HashMap::<usize, usize>, 
     HashMap::<usize, usize>) {
    let ibody = &mesh.bodies[body_id-1];
    let mut i: usize = 0;
    for element_id in &ibody.element_ids {
        let element = NIElement::new(&mesh, *element_id);
        // get collocation points for this element
        let ecpts= element.get_integration_points_and_normals();
        // dump into global data
        for mut ecpt in ecpts {
            ecpt.id = i;
            mesh.cpts.push(ecpt);
            i += 1;
        }
    }

    // Scroll through all used eqns in order and put in map
    let mut eqn_map = HashMap::<usize, usize>::new();
    let mut cpt_map = HashMap::<usize, usize>::new();
    let mut eqn_number: usize = 0;
    for i in 0..mesh.cpts.len() {
        eqn_map.insert(i, eqn_number);
        cpt_map.insert(eqn_number, i);
        eqn_number += 1;
    }
    return (eqn_map, cpt_map)
}

fn process_frequency_list(freq_input: &input_data::FrequencyInput) 
    -> Vec<f64> {
    match freq_input {
        input_data::FrequencyInput::List {values} => {
            return values.to_vec();
        },
        input_data::FrequencyInput::LinearSpaced { start, end, number } => {
            return tools::linspace(*start, *end, *number);
        },
        input_data::FrequencyInput::LogSpaced { start, end, number } => {
            return tools::logspace(*start, *end, *number);
        }
    }
}