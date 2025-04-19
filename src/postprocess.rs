/*!
Post-processing steps and writing outputs
*/

use std::error::Error;
use csv::Writer;

use crate::preprocess;
use crate::Cplx;

/// Total radiated and incident power 
pub struct Power {
    pub radiated: f64,
    pub incident: f64
}

/// Field results at a given frequency
pub struct FPResult {
    /// Drive frequency
    pub frequency: f64,
    /// default: scattered output; may be total field, if requested
    pub scattered: Option<na::DVector<Cplx>>, 
    /// incident wave field 
    pub incident: Option<na::DVector<Cplx>>,
    /// radiated and incident power
    pub power: Power
}

/// write scattered/total and incident field results to a csv file for all points at one frequency
pub fn write_results_at_frequency(predata: &preprocess::PreData, result: &FPResult) 
    -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(predata.get_output_filename())?;
    let scatt = result.scattered.as_ref().unwrap();
    let inc = result.incident.as_ref().unwrap();
    let points = predata.get_field_points();
    if scatt.len() != points.len() {
        panic!("Error while writing results: pressure array and field point array not equal length!")
    }
    let on = match predata.get_output_field() {
        preprocess::input::OutputField::Pressure => "pre",
        preprocess::input::OutputField::VelocityPotential => "vpo"
    };
    let _ = wtr.write_record(&["index","x","y","z",
                               &format!("{}{}",on,"_re"),
                               &format!("{}{}",on,"_im"),
                               "inc_re","inc_im"]);
    for i in 0..scatt.len() {
        wtr.write_record(&[i.to_string(), 
                           points[i][0].to_string(), 
                           points[i][1].to_string(), 
                           points[i][2].to_string(), 
                           scatt[i].re.to_string(),
                           scatt[i].im.to_string(),
                           inc[i].re.to_string(),
                           inc[i].im.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

/// write scattered/total and incident field results to a csv file for one point at all frequencies
pub fn write_results_at_point(predata: &preprocess::PreData, results: &Vec<FPResult>, index: usize) 
    -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(predata.get_output_filename())?;
    let on = match predata.get_output_field() {
        preprocess::input::OutputField::Pressure => "pre",
        preprocess::input::OutputField::VelocityPotential => "vpo"
    };
    let _ = wtr.write_record(&["freq",
                               &format!("{}{}",on,"_re"),
                               &format!("{}{}",on,"_im"),
                               "inc_re","inc_im"]);
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

/// convert all field results from velocity potential to pressure, if requested
pub fn convert_results(predata: &preprocess::PreData, results: &mut Vec<FPResult>) {
    if *predata.get_output_field() == preprocess::input::OutputField::Pressure {
        let rho = predata.get_mass_density();
        for result in results {
            let omega = 2.0 * std::f64::consts::PI * result.frequency;
            convert_phi_to_pressure_vec(result.scattered.as_mut().unwrap(), omega, rho);
            convert_phi_to_pressure_vec(result.incident.as_mut().unwrap(), omega, rho);
        }
    }
}

/// convert velocity potential (phi) to pressure (p) using the formula
/// p = (i omega rho) phi
fn convert_phi_to_pressure_vec(phi: &mut na::DVector<Cplx>, omega: f64, rho: f64) {
    let scale = Cplx::new(0.0, omega * rho);
    for i in 0..phi.len() {
        phi[i] *= scale;
    }
}