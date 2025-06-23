/*!
Preprocessing steps and input
*/

pub mod input;
pub mod mesh;

use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;

use crate::Cplx;
use crate::elements::NIElement;
use crate::tools;

/// Burton-Miller method switch and factor (frequency dependent)
pub struct BurtonMiller {
    /// true: use Burton-Miller method
    pub is: bool,
    /// Complex-valued coupling factor, default i/k
    pub factor: Cplx,
}

/// Hash maps for indexing
///
/// from collocation point index to equation index (and vice versa), collocation point index to equation index (and vice versa)
pub struct Maps {
    // map frm cpts to elements
    cpt2el: HashMap<usize, usize>,
    // map from elements to cpts
    el2cpt: HashMap<usize, usize>,
    // map from cpts to equations
    cpt2eqn: HashMap<usize, usize>,
    // map from equations to cpts
    eqn2cpt: HashMap<usize, usize>,
}

/// Memory use and time
pub struct Usage {
    pub num_eqn: usize,
    pub max_mem: usize, // Bytes
    pub wall_time: f64, // seconds
}

impl Usage {
    /// Update memory use and wall time with new data (consumed)
    pub fn update(&mut self, new_usage: Usage) {
        self.max_mem = usize::max(self.max_mem, new_usage.max_mem);
        self.wall_time += new_usage.wall_time
    }
    pub fn print(&self) {
        let (nb, f) = tools::report_memory(self.max_mem);
        info!(
            " Duration: {:7.4} s, Peak memory used: {:4.2} {}",
            self.wall_time, nb, f
        )
    }
}

/// Preprocessing data, held by analysis
pub struct PreData {
    input: input::UserInput,
    mesh: mesh::Mesh,
    maps: Maps,
    frequency_list: Vec<f64>,
    ifreq: usize, // current frequency index
    usage: Usage,
}

impl PreData {
    /// Set frequency index to input
    pub fn set_frequency_index(&mut self, index: &usize) {
        self.ifreq = *index;
    }
    /// Get list of frequencies (reference)
    #[inline]
    pub fn get_frequencies(&self) -> &Vec<f64> {
        &self.frequency_list
    }
    /// Get current analysis frequency (copy)
    #[inline]
    pub fn get_frequency(&self) -> f64 {
        self.frequency_list[self.ifreq]
    }
    /// Get current analysis ANGULAR frequency (copy)
    #[inline]
    pub fn get_angular_frequency(&self) -> f64 {
        2.0 * PI * self.get_frequency()
    }
    /// Get current analysis wavenumber (k = omega / c)
    #[inline]
    pub fn get_wavenumber(&self) -> f64 {
        self.get_angular_frequency() / self.get_sound_speed()
    }
    /// Get sound speed of acoustic medium (c)
    #[inline]
    pub fn get_sound_speed(&self) -> f64 {
        self.input.sound_speed
    }
    /// Get mass density of acoustic medium (rho)
    #[inline]
    pub fn get_mass_density(&self) -> f64 {
        self.input.mass_density
    }
    /// Get problem type (internal/external fluid)
    #[inline]
    pub fn get_problem_type(&self) -> &input::ProblemType {
        &self.input.problem_type
    }
    /// Get method type (classical or Burton-Miller)
    #[inline]
    pub fn get_method_type(&self) -> &input::MethodType {
        &self.input.method_type
    }
    /// Return new hypersingular struct in current state
    #[inline]
    pub fn get_hypersingular(&self) -> BurtonMiller {
        BurtonMiller {
            is: self.use_hypersingular(),
            factor: self.get_burton_miller_coupling_factor(),
        }
    }
    /// Return whether to use Burton-Miller method
    #[inline]
    fn use_hypersingular(&self) -> bool {
        self.input.method_type == input::MethodType::BurtonMiller
    }
    /// Calculator Burton-Miller coupling factor (beta = i/k)
    #[inline]
    fn get_burton_miller_coupling_factor(&self) -> Cplx {
        Cplx::new(0.0, 1.0 / self.get_wavenumber())
    }
    /// Get added diagonal of H matrix. The H matrix has -1/2 added along the diagonal for exterior problems
    #[inline]
    pub fn get_hdiag(&self) -> Cplx {
        match self.get_problem_type() {
            input::ProblemType::Exterior => Cplx::new(-0.5, 0.0),
            input::ProblemType::Interior => Cplx::new(0.0, 0.0),
        }
    }
    /// Get added diagonal of G matrix. The G matrix has 0.5 * beta added along the diagonal for exterior problems
    #[inline]
    pub fn get_gdiag(&self) -> Cplx {
        match self.get_method_type() {
            // the G matrix has 1/2 (ultimately beta/2, where beta = i/k) added along the diagonal for Burton-Miller formulation
            input::MethodType::Classical => Cplx::new(0.0, 0.0),
            input::MethodType::BurtonMiller => 0.5 * self.get_burton_miller_coupling_factor(),
        }
    }
    /// Return reference to mesh body
    #[inline]
    pub fn get_mesh_body(&self) -> &mesh::Body {
        &self.mesh.bodies[self.input.body_index - 1]
    }
    /// Return reference to solver struct
    #[inline]
    pub fn get_solver(&self) -> &input::Solver {
        &self.input.solver
    }
    /// Return reference to incident wave input struct
    #[inline]
    pub fn get_incident_wave(&self) -> &Vec<input::IncidentWaveInput> {
        &self.input.incident_wave
    }
    /// Return reference to surface boundary condition struct
    #[inline]
    pub fn get_surface_bc(&self) -> &input::SurfaceBoundaryCondition {
        &self.input.surface_bc
    }
    /// Return reference to map from collocation point index to equation
    #[inline]
    pub fn get_cpt2eqn_map(&self) -> &HashMap<usize, usize> {
        &self.maps.cpt2eqn
    }
    /// Return reference to map from equation to collocation point index
    #[inline]
    pub fn get_eqn2cpt_map(&self) -> &HashMap<usize, usize> {
        &self.maps.eqn2cpt
    }
    /// Return reference to map from collocation point index to element index
    #[inline]
    pub fn get_cpt2el_map(&self) -> &HashMap<usize, usize> {
        &self.maps.cpt2el
    }
    /// Return reference to map from element number to collocation point index
    #[inline]
    pub fn get_el2cpt_map(&self) -> &HashMap<usize, usize> {
        &self.maps.el2cpt
    }
    /// Return reference to mesh
    #[inline]
    pub fn get_mesh(&self) -> &mesh::Mesh {
        &self.mesh
    }
    /// Return refernnce to collocation point vector
    #[inline]
    pub fn get_cpts(&self) -> &Vec<mesh::CollocationPoint> {
        &self.mesh.cpts
    }
    /// Return number of equations
    #[inline]
    pub fn get_num_eqn(&self) -> usize {
        self.mesh.cpts.len()
    }
    /// Return reference to output file name
    #[inline]
    pub fn get_output_filename(&self) -> &String {
        &self.input.output.file
    }
    /// Return reference to output field (pressure or velocity potential)
    #[inline]
    pub fn get_output_field(&self) -> &input::OutputField {
        &self.input.output.field
    }
    /// Return reference to output type (total or scattered)
    #[inline]
    pub fn get_output_type(&self) -> &input::OutputType {
        &self.input.output.o_type
    }
    /// Return reference to vector of field points
    #[inline]
    pub fn get_field_points(&self) -> &Vec<[f64; 3]> {
        &self.input.output.field_points
    }
    /// Return bool of whether to write output power results
    #[inline]
    pub fn get_output_power_bool(&self) -> bool {
        self.input.output.request_power
    }
    /// Return reference to resource usage data
    #[inline]
    pub fn get_usage(&self) -> &Usage {
        &self.usage
    }
    /// Return mutable reference to resource usage data
    #[inline]
    pub fn get_mut_usage(&mut self) -> &mut Usage {
        &mut self.usage
    }
}

/// Wrapper of preprocessing steps
pub fn preprocess(input: input::UserInput) -> PreData {
    info!(" Preprocessing...");
    // read mesh VTK
    let mut mesh: mesh::Mesh = Default::default();
    let _result = mesh.read_from_vtk(Path::new(&input.mesh_file));
    let body_id = &input.body_index;

    let frequency_list = process_frequency_list(&input.frequency);

    let maps = process_collocation_pts(&mut mesh, *body_id);

    process_elements(&mut mesh, *body_id);

    info!(" Using {} threads...", tools::get_num_threads());

    let usage = Usage {
        num_eqn: mesh.cpts.len(),
        max_mem: 0,
        wall_time: 0.0,
    };

    // take ownership of input data
    PreData {
        input,
        mesh,
        maps,
        frequency_list,
        ifreq: 0,
        usage,
    }
}

/// Get the element integration data so it does not need to be done repeatedly
fn process_elements(mesh: &mut mesh::Mesh, body_id: usize) {
    let ibody = &mesh.bodies[body_id - 1];
    for element_id in &ibody.element_ids {
        // get element with all integration points
        let element = NIElement::new(mesh, *element_id);
        // get integration data
        let ecpts = element.get_integration_points_and_normals();
        // assign into mesh data
        mesh.elements[*element_id].intpts = ecpts;
    }
}

/// Calculate the collocation points and normals using surface elements
fn process_collocation_pts(mesh: &mut mesh::Mesh, body_id: usize) -> Maps {
    let ibody = &mesh.bodies[body_id - 1];
    let mut cpt2el = HashMap::<usize, usize>::new();
    let mut el2cpt = HashMap::<usize, usize>::new();
    let mut i: usize = 0;
    for element_id in &ibody.element_ids {
        // get element with only central integration point
        let element = NIElement::new_1_point_integration(mesh, *element_id);
        // get collocation point for this element
        let ecpts = element.get_integration_points_and_normals();
        // dump into global data
        for mut ecpt in ecpts {
            ecpt.id = i;
            cpt2el.insert(i, *element_id);
            el2cpt.insert(*element_id, i);
            mesh.cpts.push(ecpt);
            i += 1;
        }
    }

    // Scroll through all used eqns in order and put in map
    let mut cpt2eqn = HashMap::<usize, usize>::new();
    let mut eqn2cpt = HashMap::<usize, usize>::new();
    for i in 0..mesh.cpts.len() {
        cpt2eqn.insert(i, i);
        eqn2cpt.insert(i, i);
    }
    Maps {
        cpt2el,
        el2cpt,
        cpt2eqn,
        eqn2cpt,
    }
}

/// Set up frequency vector based on input
fn process_frequency_list(freq_input: &input::FrequencyInput) -> Vec<f64> {
    match freq_input {
        input::FrequencyInput::List { values } => values.to_vec(),
        input::FrequencyInput::LinearSpaced { start, end, number } => {
            tools::linspace(*start, *end, *number)
        }
        input::FrequencyInput::LogSpaced { start, end, number } => {
            tools::logspace(*start, *end, *number)
        }
    }
}
