extern crate nalgebra as na;

pub mod preprocess;
pub mod analysis;
pub mod elements;
pub mod incident_wave;
pub mod influence_matrix;
pub mod solve;
pub mod postprocess;

pub use analysis::Analysis as Analysis;

pub type Cplx = na::Complex<f64>;