/*!
# cc_abstemious 

`cc_abstemious` is an anagram of BEM-Acoustics.

`cc_abstemious` is a numerical acoustics simulation software based on the boundary element method. 

*/

pub const VER_MAJOR: usize = 1;
pub const VER_MINOR: usize = 2; 

extern crate simplelog;
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