//! `cc_abstemious` is a numerical acoustics simulation software based on the boundary element method.

extern crate nalgebra as na;
#[macro_use] extern crate log;
extern crate scoped_threadpool;

pub mod preprocess;
pub mod elements;
pub mod incident_wave;
pub mod influence_matrix;
pub mod solve;
pub mod postprocess;
pub mod tools;

pub mod analysis;
pub use analysis::Analysis as Analysis;
pub type Cplx = na::Complex<f64>;