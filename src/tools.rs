/*!
Methods with no better place
*/

pub fn linspace(start: f64, end: f64, npts: usize) -> Vec<f64> {
    //! evenly-spaced npts between start and end
    let dx = (end - start) / ((npts - 1) as f64);
    let mut x = vec![start; npts];
    for (i, xi) in x.iter_mut().enumerate().take(npts-1).skip(1) {
        *xi = start + dx * (i as f64);
    }
    x[npts-1] = end;
    x
}

pub fn logspace(start: f64, end: f64, npts: usize) -> Vec<f64> {
    //! evenly-spaced npts between start and end, evaluated 
    //! such that there are equal intervals between log(start)
    //! and log(end)
    let start_log = start.log10();
    let end_log = end.log10();
    let x_log = linspace(start_log, end_log, npts);
    let mut x = vec![start; npts];
    x[npts-1] = end;
    for i in 1..(npts-1) {
        x[i] = 10.0_f64.powf(x_log[i]);
    }
    x
}

pub fn get_num_threads() -> usize {
    match std::thread::available_parallelism() {
        Ok(result) => std::cmp::max(result.get() / 2, 2),
        Err(_) => 2
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use crate::tools::{logspace, linspace};

    #[test]
    fn test_linspace() {
        let x = linspace(1.0, 10.0, 3);
        assert_eq!(x.len(), 3);
        assert_eq!(x, vec![1.0, 5.5, 10.0])
    }

    #[test]
    fn test_linspace_2() {
        let x = linspace(1.0, 10.0, 2);
        assert_eq!(x.len(), 2);
        assert_eq!(x, vec![1.0, 10.0])
    }

    #[test]
    fn test_logspace() {
        let x = logspace(1.0, 10.0, 3);
        assert_eq!(x.len(), 3);
        assert_relative_eq!(x[1], 3.1622776601683795, epsilon=1e-8)
    }

    #[test]
    fn test_logspace_2() {
        let x = logspace(1.0, 10.0, 2);
        assert_eq!(x.len(), 2);
        assert_eq!(x, vec![1.0, 10.0])
    }
}