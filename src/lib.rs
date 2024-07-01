extern crate nalgebra as na;
#[macro_use] extern crate log;
extern crate scoped_threadpool;

pub mod preprocess;
mod elements;
mod incident_wave;
mod influence_matrix;
mod solve;
mod postprocess;

pub mod analysis;
pub use analysis::Analysis as Analysis;
pub type Cplx = na::Complex<f64>;