use std::error::Error;

use csv::Writer;

use crate::preprocess;
use crate::Cplx;

pub struct FPResult {
    pub frequency: f64,
    pub scattered: Option<na::DVector<Cplx>>,
    pub incident: Option<na::DVector<Cplx>>,
    pub radiated_power: f64
}

pub fn write_results_at_frequency(predata: &preprocess::PreData, result: &FPResult) 
    -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(predata.get_output_filename())?;

    let fp = result.scattered.as_ref().unwrap();
    let fpi = result.incident.as_ref().unwrap();
    let points = predata.get_field_points();
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

pub fn write_results_at_point(predata: &preprocess::PreData, results: &Vec<FPResult>, index: usize) 
    -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(predata.get_output_filename())?;
    let _ = wtr.write_record(&["freq","phi_re","phi_im","inc_re","inc_im"]);
    for result in results {
        wtr.write_record(&[result.frequency.to_string(), 
                           result.scattered.as_ref().unwrap()[index].re.to_string(),
                           result.scattered.as_ref().unwrap()[index].im.to_string(),
                           result.incident.as_ref().unwrap()[index].re.to_string(),
                           result.incident.as_ref().unwrap()[index].im.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

pub fn convert_results(predata: &preprocess::PreData, results: &mut Vec<FPResult>) {
    if *predata.get_output_field() == preprocess::input_data::OutputField::Pressure {
        let rho = predata.get_mass_density();
        for result in results {
            let omega = 2.0 * std::f64::consts::PI * result.frequency;
            convert_phi_to_pressure_vec(result.scattered.as_mut().unwrap(), omega, rho);
            convert_phi_to_pressure_vec(result.incident.as_mut().unwrap(), omega, rho);
        }
    }
}

pub fn convert_phi_to_pressure_vec(phi: &mut na::DVector<Cplx>, omega: f64, rho: f64) {
    let scale = Cplx::new(0.0, omega * rho);
    for i in 0..phi.len() {
        phi[i] *= scale;
    }
}