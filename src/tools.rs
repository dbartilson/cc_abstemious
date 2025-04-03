pub fn linspace(start: f64, end: f64, npts: usize) -> Vec<f64> {
    // evenly-spaced npts between start and end
    let dx = (end - start) / ((npts - 1) as f64);
    let mut x = vec![start; npts];
    for i in 1..(npts-1) {
        x[i] = start + dx * (i as f64);
    }
    x[npts-1] = end;
    return x;
}

pub fn logspace(start: f64, end: f64, npts: usize) -> Vec<f64> {
    // evenly-spaced npts between start and end, evaluated 
    // such that there are equal intervals between log(start)
    // and log(end)
    let start_log = start.log10();
    let end_log = end.log10();
    let x_log = linspace(start_log, end_log, npts);
    let mut x = vec![start; npts];
    x[npts-1] = end;
    for i in 1..(npts-1) {
        x[i] = 10.0_f64.powf(x_log[i]);
    }
    return x;
}