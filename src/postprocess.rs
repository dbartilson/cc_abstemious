use std::error::Error;

use na::{Complex, DVector};
use csv::Writer;
type Cplx = Complex<f64>;

pub struct FPResult {
    pub frequency: f64,
    pub phi: Option<na::DVector<Cplx>>,
    pub phi_inc: Option<na::DVector<Cplx>>,
    pub radiated_power: f64
}

pub fn write_fp_csv(filename: &String, fp: &DVector<Cplx>, fpi: &DVector<Cplx>, points: &Vec<[f64;3]>) 
    -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(filename)?;
    if fp.len() != points.len() {
        panic!("Error while writing results: pressure array and field point array not equal length!")
    }
    let _ = wtr.write_record(&["index","x","y","z","phi_re","phi_im","inc_re","inc_im"]);
    for i in 0..fp.len() {
        wtr.write_record(&[i.to_string(), 
                           points[i][0].to_string(), 
                           points[i][1].to_string(), 
                           points[i][2].to_string(), 
                           fp[i].re.to_string(),
                           fp[i].im.to_string(),
                           fpi[i].re.to_string(),
                           fpi[i].im.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}