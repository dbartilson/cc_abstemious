/*!
# cc_abstemious

`cc_abstemious` is an anagram of BEM-Acoustics.

`cc_abstemious` is a numerical acoustics simulation software based on the boundary element method.
*/

pub const VER_MAJOR: usize = 1;
pub const VER_MINOR: usize = 3;
pub const VER_SUBMINOR: usize = 2;

extern crate nalgebra as na;
extern crate simplelog;
#[macro_use]
extern crate log;
extern crate scoped_threadpool;

pub mod elements;
pub mod incident_wave;
pub mod influence_matrix;
pub mod postprocess;
pub mod preprocess;
pub mod solve;
pub mod tools;

pub mod analysis;
pub use analysis::Analysis;
pub type Cplx = na::Complex<f64>;
