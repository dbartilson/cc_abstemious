use schemars::JsonSchema;
use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use std::error::Error;
use std::path::Path;

#[derive(Deserialize, JsonSchema)]
pub enum FrequencyInput {
    List { values: Vec<f64> },
    LinearSpaced { start: f64, end: f64, number: usize },
    LogSpaced { start: f64, end: f64, number: usize }
}

#[derive(Deserialize, JsonSchema)]
pub enum ProblemType {
    Interior,
    Exterior
}

#[derive(Deserialize, PartialEq, JsonSchema)]
pub enum MethodType {
    Classical,
    BurtonMiller
}

#[derive(Deserialize, JsonSchema)]
pub enum Solver {
    Direct {},
    Iterative { tolerance: f64, max_iterations: usize},
    Hierarchical { tolerance: f64, max_iterations: usize}
}

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

#[derive(Deserialize, JsonSchema)]
pub enum BCType {
    Pressure,
    NormalVelocity,
    Impedance
}

#[derive(Deserialize, JsonSchema)]
pub struct SurfaceBoundaryCondition {
    pub bc_type: BCType,
    pub value: [f64; 2]
}

#[derive(Deserialize, PartialEq, JsonSchema)]
pub enum OutputType {
    Total,
    Scattered
}

#[derive(Deserialize, PartialEq, JsonSchema)]
pub enum OutputField {
    Pressure,
    VelocityPotential
}

#[derive(Deserialize, JsonSchema)]
pub struct Output {
    pub o_type: OutputType,
    pub field: OutputField,
    pub file: String,
    pub field_points: Vec<[f64;3]>
}

#[derive(Deserialize, JsonSchema)]
pub struct UserInput {
    pub mesh_file: String,
    pub body_index: usize,
    pub frequency: FrequencyInput,
    pub sound_speed: f64,
    pub mass_density: f64,
    pub problem_type: ProblemType,
    pub method_type: MethodType,
    pub solver: Solver,
    pub incident_wave: IncidentWaveInput,
    pub surface_bc: SurfaceBoundaryCondition,
    pub output: Output
}

pub fn read_input_json<P: AsRef<Path>>(path: P) -> Result<UserInput, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    info!(" Reading json file '{}' ...", path.as_ref().display().to_string());
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `UserInput`
    let u = serde_json::from_reader(reader)?;
    Ok(u)
}

pub fn read_input_string(str: &str) -> Result<UserInput, Box<dyn Error>> {
    let u = serde_json::from_str(str)?;
    Ok(u)
}

#[cfg(test)]
mod tests {
    use schemars::schema_for;
    use crate::preprocess::input_data;
    use std::{fs::File, io::Write};

    #[test]
    fn json_reader() {
        use std::path::Path;
        // test json reader capability
        let u = input_data::read_input_json(Path::new("./src/tests/input_1.json")).unwrap();
        assert_eq!(u.body_index, 3);
    }

    #[test]
    fn write_schema() {
        let schema = schema_for!(input_data::UserInput);
        let mut file = File::create("input_schema.json").unwrap();
        let res= file.write_all(serde_json::to_string_pretty(&schema).unwrap().as_bytes());
        assert!(res.is_ok())
    }
}