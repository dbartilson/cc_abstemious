/*!
Preprocessing steps and input
*/

pub mod input;
pub mod mesh;

use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

use crate::elements::NIElement;
use crate::tools;
use crate::Cplx;

/// Burton-Miller method
pub struct BurtonMiller {
    /// true: use Burton-Miller method
    pub is: bool,
    /// Complex-valued coupling factor, default i/k
    pub factor: Cplx
}

/// Preprocessing data, held by analysis
pub struct PreData {
    input: input::UserInput,
    mesh: mesh::Mesh,
    eqn_map: HashMap<usize,usize>, 
    node_map: HashMap<usize,usize>,
    frequency_list: Vec<f64>,
    ifreq: usize // current frequency index
}

impl PreData {
    /// Set frequency index to input
    pub fn set_frequency_index(&mut self, index: &usize) {self.ifreq = *index;}
    /// Get list of frequencies (reference)
    #[inline]
    pub fn get_frequencies(&self) -> &Vec<f64> {return &self.frequency_list;}
    /// Get current analysis frequency (copy)
    #[inline]
    pub fn get_frequency(&self) -> f64 {return self.frequency_list[self.ifreq];}
    /// Get current analysis ANGULAR frequency (copy)
    #[inline]
    pub fn get_angular_frequency(&self) -> f64 {return 2.0 * PI * self.get_frequency();}
    /// Get current analysis wavenumber (k = omega / c)
    #[inline]
    pub fn get_wavenumber(&self) -> f64 {return self.get_angular_frequency() / self.get_sound_speed()}
    /// Get sound speed of acoustic medium (c)
    #[inline]
    pub fn get_sound_speed(&self) -> f64 {return self.input.sound_speed;}
    /// Get mass density of acoustic medium (rho)
    #[inline]
    pub fn get_mass_density(&self) -> f64 {return self.input.mass_density;}
    /// Get problem type (internal/external fluid)
    #[inline]
    pub fn get_problem_type(&self) -> &input::ProblemType {return &self.input.problem_type;}
    /// Get method type (classical or Burton-Miller)
    #[inline]
    pub fn get_method_type(&self) -> &input::MethodType {return &self.input.method_type;}
    /// Return new hypersingular struct in current state
    #[inline]
    pub fn get_hypersingular(&self) -> BurtonMiller {
        return BurtonMiller { is: self.use_hypersingular(), factor: self.get_burton_miller_coupling_factor()} }
    /// Return whether to use Burton-Miller method
    #[inline]
    fn use_hypersingular(&self) -> bool {return self.input.method_type == input::MethodType::BurtonMiller;}
    /// Calculator Burton-Miller coupling factor (beta = i/k)
    #[inline]
    fn get_burton_miller_coupling_factor(&self) -> Cplx {return Cplx::new(0.0, 1.0 / self.get_wavenumber())}
    /// Get added diagonal of H matrix. The H matrix has -1/2 added along the diagonal for exterior problems
    #[inline]
    pub fn get_hdiag(&self) -> Cplx {
        match self.get_problem_type() {
            input::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
            input::ProblemType::Interior => Cplx::new(0.0, 0.0)
        }
    }
    /// Get added diagonal of G matrix. The G matrix has 0.5 * beta added along the diagonal for exterior problems
    #[inline]
    pub fn get_gdiag(&self) -> Cplx {
        match self.get_method_type() {
            // the G matrix has 1/2 (ultimately beta/2, where beta = i/k) added along the diagonal for Burton-Miller formulation
            input::MethodType::Classical => Cplx::new(0.0, 0.0),
            input::MethodType::BurtonMiller => 0.5 * self.get_burton_miller_coupling_factor()
        }
    }
    /// Return reference to mesh body
    #[inline]
    pub fn get_mesh_body(&self) -> &mesh::Body {return &self.mesh.bodies[self.input.body_index - 1];}
    /// Return reference to solver struct
    #[inline]
    pub fn get_solver(&self) -> &input::Solver {return &self.input.solver}
    /// Return reference to incident wave input struct
    #[inline]
    pub fn get_incident_wave(&self) -> &input::IncidentWaveInput {return &self.input.incident_wave;}
    /// Return reference to surface boundary condition struct
    #[inline]
    pub fn get_surface_bc(&self) -> &input::SurfaceBoundaryCondition {return &self.input.surface_bc;}
    /// Return reference to map from node index to equation index
    #[inline]
    pub fn get_eqn_map(&self) -> &HashMap<usize, usize> {return &self.eqn_map;}
    /// Return reference to map from equation index to node index
    #[inline]
    pub fn get_node_map(&self) -> &HashMap<usize, usize> {return &self.node_map;}
    /// Return reference to mesh
    #[inline]
    pub fn get_mesh(&self) -> &mesh::Mesh {return &self.mesh;}
    /// Return refernnce to collocation point vector
    #[inline]
    pub fn get_cpts(&self) -> &Vec<mesh::CollocationPoint> {return &self.mesh.cpts}
    /// Return number of equations
    #[inline]
    pub fn get_num_eqn(&self) -> usize {return self.mesh.cpts.len();}
    /// Return reference to output file name
    #[inline]
    pub fn get_output_filename(&self) -> &String {return &self.input.output.file;}
    /// Return reference to output field (pressure or velocity potential)
    #[inline]
    pub fn get_output_field(&self) -> &input::OutputField {return &self.input.output.field;}
    /// Return reference to output type (total or scatteres)
    #[inline]
    pub fn get_output_type(&self) -> &input::OutputType {return &self.input.output.o_type;}
    /// Return reference to vector of field points
    #[inline]
    pub fn get_field_points(&self) -> &Vec<[f64; 3]> {return &self.input.output.field_points;}
}

/// Wrapper of preprocessing steps
pub fn preprocess(input: input::UserInput) -> PreData {
    info!(" Preprocessing...");
    // read mesh VTK
    let mut mesh: mesh::Mesh = Default::default();
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
fn process_collocation_pts(mesh: &mut mesh::Mesh, body_id: usize) -> 
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

/// Set up frequency vector based on input
fn process_frequency_list(freq_input: &input::FrequencyInput) -> Vec<f64> {
    match freq_input {
        input::FrequencyInput::List {values} => {
            return values.to_vec();
        },
        input::FrequencyInput::LinearSpaced { start, end, number } => {
            return tools::linspace(*start, *end, *number);
        },
        input::FrequencyInput::LogSpaced { start, end, number } => {
            return tools::logspace(*start, *end, *number);
        }
    }
}