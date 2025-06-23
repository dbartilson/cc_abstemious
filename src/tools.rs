/*!
Methods with no better place
*/

pub fn linspace(start: f64, end: f64, npts: usize) -> Vec<f64> {
    //! evenly-spaced npts between start and end
    let dx = (end - start) / ((npts - 1) as f64);
    let mut x = vec![start; npts];
    for (i, xi) in x.iter_mut().enumerate().take(npts - 1).skip(1) {
        *xi = start + dx * (i as f64);
    }
    x[npts - 1] = end;
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
    x[npts - 1] = end;
    for i in 1..(npts - 1) {
        x[i] = 10.0_f64.powf(x_log[i]);
    }
    x
}

pub fn get_num_threads() -> usize {
    //! get the number of threads to use. Default to nthreads / 2, unless that fails, then default to 2
    match std::thread::available_parallelism() {
        Ok(result) => std::cmp::max(result.get() / 2, 2),
        Err(_) => 2,
    }
}

pub fn get_threadpool() -> rayon::ThreadPool {
    //! return Rayon threadpool
    let nthreads = get_num_threads();
    rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap()
}

struct MemRange {
    limit: usize,
    divisor: usize,
    unit: &'static str,
}

static MEM_RANGES: [MemRange; 4] = [
    MemRange {
        limit: 1e3 as usize,
        divisor: 1,
        unit: "B",
    },
    MemRange {
        limit: 1e6 as usize,
        divisor: 1e3 as usize,
        unit: "KB",
    },
    MemRange {
        limit: 1e9 as usize,
        divisor: 1e6 as usize,
        unit: "MB",
    },
    MemRange {
        limit: 1e12 as usize,
        divisor: 1e9 as usize,
        unit: "GB",
    },
];

/// Convert from a usize number of bytes to a pair of f64 and string, something like
/// 1024 Bytes => (1.024, "KB")
pub fn report_memory(nbytes: usize) -> (f64, String) {
    for mr in &MEM_RANGES {
        if nbytes < mr.limit {
            return (nbytes as f64 / mr.divisor as f64, mr.unit.to_string());
        }
    }
    (nbytes as f64, MEM_RANGES[0].unit.to_string())
}

#[cfg(test)]
mod tests {
    use crate::tools::{linspace, logspace};
    use approx::assert_relative_eq;

    #[test]
    fn linspace_3pts() {
        let x = linspace(1.0, 10.0, 3);
        assert_eq!(x.len(), 3);
        assert_eq!(x, vec![1.0, 5.5, 10.0])
    }

    #[test]
    fn linspace_2pts() {
        let x = linspace(1.0, 10.0, 2);
        assert_eq!(x.len(), 2);
        assert_eq!(x, vec![1.0, 10.0])
    }

    #[test]
    fn logspace_3pts() {
        let x = logspace(1.0, 10.0, 3);
        assert_eq!(x.len(), 3);
        assert_relative_eq!(x[1], 3.1622776601683795, epsilon = 1e-8)
    }

    #[test]
    fn logspace_2pts() {
        let x = logspace(1.0, 10.0, 2);
        assert_eq!(x.len(), 2);
        assert_eq!(x, vec![1.0, 10.0])
    }
}
