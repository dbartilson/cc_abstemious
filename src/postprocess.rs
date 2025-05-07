/*!
Post-processing steps and writing outputs
*/

use std::error::Error;
use std::path::Path;
use csv::Writer;
use crate::analysis::PrimaryVariables;
use crate::incident_wave;
use crate::preprocess;
use crate::solve;
use crate::Cplx;

/// Total scattered and incident power on surface
pub struct Power {
    pub scattered: f64,
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
    /// scattered and incident power
    pub power: Power
}

impl FPResult {
    /// convert all field results from velocity potential to pressure, if requested
    pub fn convert_result(&mut self, predata: &preprocess::PreData) {
        if *predata.get_output_field() == preprocess::input::OutputField::Pressure {
            let rho = predata.get_mass_density();
            let omega = 2.0 * std::f64::consts::PI * self.frequency;
            convert_phi_to_pressure_vec(self.scattered.as_mut().unwrap(), omega, rho);
            convert_phi_to_pressure_vec(self.incident.as_mut().unwrap(), omega, rho);
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

/// Post-processing steps: get incident field
pub fn postprocess(predata: &preprocess::PreData, incident: &PrimaryVariables, total: &PrimaryVariables) 
    -> FPResult
{
    let phi_inc_fp = incident_wave::get_incident_field(predata);
            
    let phi_fp = solve::get_field(predata, total, &phi_inc_fp);

    let power = calculate_surface_power(predata, incident, total);

    let mut result = FPResult{
        frequency: predata.get_frequency(),
        scattered: Some(phi_fp),
        incident: Some(phi_inc_fp),
        power 
    };

    result.convert_result(predata);

    result
}

/// Calculate power by integrating intensity over the surface
/// 
/// Intensity is equal to 0.5 * Re (p v_n*), where v_n* is the complex conjugate of the normal velocity
/// = 0.5 * omega * rho * Re (i p v_n*)
pub fn calculate_surface_power(predata: &preprocess::PreData, incident: &PrimaryVariables, total: &PrimaryVariables)
    -> Power
{
    let mut w_inc = 0.0;
    let mut w_scatt = 0.0;
    let rho = predata.get_mass_density();
    let omega = predata.get_angular_frequency();
    let iomegarho = Cplx::new(0.0, omega * rho);
    for (i, cpt) in predata.get_cpts().iter().enumerate() {
        let i_inc = 0.5 * ((incident.phi[i] * iomegarho) * incident.vn[i].conj()).re;
        let phi_scatt = total.phi[i] - incident.phi[i];
        let vn_scatt = total.vn[i] - incident.vn[i];
        let i_scatt = 0.5 * ((phi_scatt * iomegarho) * vn_scatt.conj()).re;
        w_inc += i_inc * cpt.area * cpt.wt;
        w_scatt += i_scatt * cpt.area * cpt.wt;
    }
    Power {
        incident: w_inc,
        scattered: w_scatt
    }
}

/// write scattered/total and incident field results to a csv file for all points at one frequency
/// Optional override for input file_name
pub fn write_results_at_frequency(predata: &preprocess::PreData, file_name: Option<&String>, result: &FPResult) 
    -> Result<(), Box<dyn Error>> {
    let fname = file_name.unwrap_or(predata.get_output_filename());
    if fname.is_empty() { return Ok(()) }

    let mut wtr = Writer::from_path(fname)?;
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
    let _ = wtr.write_record(["index","x","y","z",
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
    let on = match *predata.get_output_type() {
        preprocess::input::OutputType::Scattered => "scatt",
        preprocess::input::OutputType::Total => "tot"
    };
    let _ = wtr.write_record(["freq",
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

/// write scattered/radiated and incident power over the surface
pub fn write_power(predata: &preprocess::PreData, file_name: Option<&String>, results: &Vec<FPResult>) -> Result<(), Box<dyn Error>> {
    let fname = file_name.unwrap_or(predata.get_output_filename());
    if fname.is_empty() { return Ok(()) }
    
    let mut wtr = Writer::from_path(fname)?;
    let _ = wtr.write_record(["freq",
                              "w_scat",
                              "w_inc"]);
    for result in results {
    wtr.write_record(&[result.frequency.to_string(), 
                       result.power.scattered.to_string(),
                       result.power.incident.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

/// write results to file(s)
/// If no file name given, write no results
/// If only one field point and no power requested, condense it all into one csv
/// If multiple field points or power requested:
///     Write a file for each frequency, with the index appended before the file extension
///     Write a single file for the power data at all frequencies
///     Write an index file that contains data on index, file name, and frequency for the output files
pub fn write_results(predata: &preprocess::PreData, results: &Vec<FPResult>) -> Result<(), Box<dyn Error>> {
    if predata.get_output_filename().is_empty() { return Ok(()); }

    if predata.get_field_points().len() == 1 && !predata.get_output_power_bool() {
        write_results_at_point(predata, results, 0)
    }
    else {
        // Loop over results, write a file like /path/stem_{i}.ext for each frequency, along with separate file for power
        
        // Process specified file path into chunks
        let output_file = Path::new(predata.get_output_filename());
        let of_path = output_file.parent().unwrap().to_str().unwrap();
        let of_stem = output_file.file_stem().unwrap().to_str().unwrap();
        let of_ext = output_file.extension().unwrap().to_str().unwrap();

        let index_file_name = format!("{}/{}_index.{}", of_path, of_stem, of_ext);
        let mut wtr = Writer::from_path(index_file_name)?;
        wtr.write_record(["freq_index", "frequency (Hz)", "file_name"])?;
        for (i, result) in results.iter().enumerate() {
            let ofi = format!("{}/{}_{}.{}", of_path, of_stem, i+1, of_ext);
            let _ = write_results_at_frequency(predata, Some(&ofi), result);
            wtr.write_record([(i+1).to_string(), result.frequency.to_string(), ofi])?;
        }
        wtr.flush()?;
        let ofpower = format!("{}/{}_power.{}", of_path, of_stem, of_ext);
        if predata.get_output_power_bool() { let _ = write_power(predata, Some(&ofpower), results); }
        Ok(())
    }
}