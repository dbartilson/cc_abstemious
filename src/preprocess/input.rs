/*!
Input processing
*/

use schemars::JsonSchema;
use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use std::error::Error;
use std::path::Path;

/// Frequency input, can be list, linear-spaced, or log-spaced
#[derive(Deserialize, JsonSchema)]
pub enum FrequencyInput {
    List { values: Vec<f64> },
    LinearSpaced { start: f64, end: f64, number: usize },
    LogSpaced { start: f64, end: f64, number: usize }
}

/// Interior or exterior
#[derive(Deserialize, JsonSchema)]
pub enum ProblemType {
    Interior,
    Exterior
}

/// Classical or Burton-Miller
#[derive(Deserialize, PartialEq, JsonSchema)]
pub enum MethodType {
    Classical,
    BurtonMiller
}

/// Solver type and fields
/// 
/// Direct = Dense matrix with LU solve
/// Iterative = Dense matrix, but GMRES iterative solve
/// Hierarchical (aka ACA) = sparse representation, GMRES solve
#[derive(Deserialize, JsonSchema)]
pub enum Solver {
    Direct {},
    Iterative { tolerance: f64, max_iterations: usize},
    Hierarchical { tolerance: f64, max_iterations: usize}
}

/// Incident wave inputs, plane wave or spherical wave
#[derive(Deserialize, JsonSchema)]
pub enum IncidentWaveInput {
    PlaneWave {
        direction: [f64;3],
        amplitude: [f64;2]    
    },
    SphericalWave {
        origin: [f64;3],
        amplitude: [f64;2]    
    }
}

/// Surface boundary condition type
#[derive(Deserialize, JsonSchema)]
pub enum BCType {
    Pressure,
    NormalVelocity,
    Impedance
}

/// Surface boundary condition (one value and type on whole surface)
#[derive(Deserialize, JsonSchema)]
pub struct SurfaceBoundaryCondition {
    pub bc_type: BCType,
    /// Complex value, but split into vec of two real numbers (a + bi)
    pub value: [f64; 2]
}

/// Output type
#[derive(Deserialize, PartialEq, JsonSchema)]
pub enum OutputType {
    Total,
    Scattered
}

/// Output field
#[derive(Deserialize, PartialEq, JsonSchema)]
pub enum OutputField {
    Pressure,
    VelocityPotential
}

/// Output struct
#[derive(Deserialize, JsonSchema)]
pub struct Output {
    /// Scattered or total
    pub o_type: OutputType,
    /// Pressure or velocity potential
    pub field: OutputField,
    /// Output file name
    pub file: String,
    /// Vector of field point coordinates
    pub field_points: Vec<[f64;3]>
}

/// Struct of user inputs
#[derive(Deserialize, JsonSchema)]
pub struct UserInput {
    /// Path to mesh file
    pub mesh_file: String,
    /// Index of analysis body in mesh file
    pub body_index: usize,
    pub frequency: FrequencyInput,
    /// Sound speed (c) of acoustic fluid
    pub sound_speed: f64,
    /// Mass density (rho) of acoustic fluid
    pub mass_density: f64,
    pub problem_type: ProblemType,
    pub method_type: MethodType,
    pub solver: Solver,
    pub incident_wave: IncidentWaveInput,
    pub surface_bc: SurfaceBoundaryCondition,
    pub output: Output
}

/// Read input from json at path
pub fn read_input_json<P: AsRef<Path>>(path: P) -> Result<UserInput, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    info!(" Reading json file '{}' ...", path.as_ref().display().to_string());
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `UserInput`
    let u = serde_json::from_reader(reader)?;
    Ok(u)
}
/// Read input from json in a string (for testing)
pub fn read_input_string(str: &str) -> Result<UserInput, Box<dyn Error>> {
    let u = serde_json::from_str(str)?;
    Ok(u)
}

#[cfg(test)]
mod tests {
    use schemars::schema_for;
    use crate::preprocess::input;
    use std::{fs::File, io::Write};

    #[test]
    fn json_reader() {
        use std::path::Path;
        // test json reader capability
        let u = input::read_input_json(Path::new("./src/tests/input_1.json")).unwrap();
        assert_eq!(u.body_index, 3);
    }

    #[test]
    fn write_schema() {
        let schema = schema_for!(input::UserInput);
        let mut file = File::create("input_schema.json").unwrap();
        let res= file.write_all(serde_json::to_string_pretty(&schema).unwrap().as_bytes());
        assert!(res.is_ok())
    }
}