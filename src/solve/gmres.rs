use na::{DMatrix, DVector, Normed};
use crate::Cplx;

enum ExitFlag {
    Error,
    Tolerance,
    Iterations
}

pub fn solve_gmresk(a: &DMatrix::<Cplx>, x: &mut DVector<Cplx>, 
                    mut max_it: usize, mut thresh: f64) {

    info!(" Using iterative (GMRES) solver...");
    let b = x.clone();
    // zero as initial guess
    x.fill(Cplx::new(0.0, 0.0));
    if max_it <= 0 {max_it = 1000;}
    if thresh <= 0.0 {thresh = 1.0e-5;}
    let k = (max_it as f64).sqrt().ceil() as usize;
    let max_restarts = max_it / k;
    for i in 0..max_restarts {
        info!("  GMRES restart: {}", i);
        let flag = gmres(a, x, &b, k, thresh);
        match flag {
            ExitFlag::Tolerance => {break;},
            _ => {}
        }
    }
}

fn gmres(a: &DMatrix::<Cplx>, x: &mut DVector<Cplx>, b: &DVector<Cplx>, 
             max_it: usize, thresh: f64) -> ExitFlag {
    
    let mut flag = ExitFlag::Error;
    let n = x.len();
    let m = max_it;

    let c_one = Cplx::new(1.0, 0.0);
    let c_zero = Cplx::new(0.0, 0.0);
    // use x as initial guess
    let b_norm = b.norm();
    let mut r = b.clone();
    r.gemv(-c_one, a, x, c_one);

    let alpha = a.norm();
    let mut error = backward_error(r.norm(), alpha, x, b_norm);

    let mut sn = DVector::<Cplx>::from_element(m, c_zero);
    let mut cs = sn.clone();
    let mut e = vec![error];
    let mut q = DMatrix::<Cplx>::from_element(n, m+1, c_zero);
    q.set_column(0, &r.normalize());
    let mut beta = DVector::<Cplx>::from_element(m+1, c_zero);
    beta[0] = Cplx::new(r.norm(), 0.0);
    let mut h = DMatrix::<Cplx>::from_element(n, m, c_zero);
    for k in 0..m {
        let kk = k+1;
        info!("   Iteration: {}, Error: {}", k, error);
        let (mut hk1, qk1) = arnoldi(a, &q, k);
        apply_givens_rotation(&mut hk1, &mut cs, &mut sn, k);
        h.set_column(k, &hk1);
        q.set_column(k+1, &qk1);
        beta[k+1] = -sn[k].conj() * beta[k];
        beta[k] *= cs[k];

        error = backward_error(beta[k+1].norm(), alpha, x, b_norm);
        e.push(error);

        if error < thresh {
            // recompute backward error using exact arithmetic
            r = b.clone();
            get_x(&h, &q, &beta, kk, x);
            r.gemv(-c_one, a, x, c_one);
            error = backward_error(r.norm(), alpha, x, b_norm);
            *e.last_mut().unwrap() = error;
            if error < thresh {     
                flag = ExitFlag::Tolerance;
                info!("  GMRES tolerance acheived! ({} < {})", error, thresh);
                break;
            }
        }
        if k == m-1 {
            flag = ExitFlag::Iterations;
        }
    }
    return flag;
}

fn backward_error(r: f64, alpha: f64, x: &DVector::<Cplx>, beta: f64) -> f64 {
    r / (alpha * x.norm() + beta)
}

fn get_x(h: &DMatrix::<Cplx>, q: &DMatrix::<Cplx>, beta: &DVector::<Cplx>, k: usize, 
         x: &mut DVector::<Cplx>) {
    let c_one = Cplx::new(1.0, 0.0);
    let hkk = h.view((0, 0), (k, k));
    let hkklu = hkk.lu();
    let betak = beta.rows(0, k);
    let y = hkklu.solve(&betak);
    x.gemv(c_one, &q.columns(0,k), y.as_ref().unwrap(), c_one);
}

fn arnoldi(a: &DMatrix::<Cplx>, q: &DMatrix<Cplx>, k: usize) 
        -> (DVector::<Cplx>, DVector::<Cplx>) {
    let n = a.shape().0;
    let c_zero = Cplx::new(0.0, 0.0);
    let c_one = Cplx::new(1.0, 0.0);
    let mut qk1 = DVector::<Cplx>::from_element(n, c_zero);
    qk1.gemv(c_one, a, &q.column(k), c_zero);
    let mut hk1 = DVector::<Cplx>::from_element(n, c_zero);
    for i in 0..k+1 {
        hk1[i] = qk1.dotc(&q.column(i));
        qk1.axpy(-hk1[i], &q.column(i), c_one);
    }
    hk1[k+1] = Cplx::new(qk1.norm(), 0.0);
    let scale = c_one / hk1[k+1];
    for qi in &mut qk1 {
        *qi *= scale;
    }
    return (hk1, qk1)
}

fn apply_givens_rotation(hk1: &mut DVector<Cplx>, cs: &mut DVector<Cplx>, 
                         sn: &mut DVector<Cplx>, k: usize) {

    for i in 0..k {
        let temp = cs[i] * hk1[i] + sn[i] * hk1[i+1];
        hk1[i+1] = -sn[i].conj() * hk1[i] + cs[i] * hk1[i+1];
        hk1[i] = temp;
    }
    (cs[k], sn[k]) = givens_rotation(&hk1[k], &hk1[k+1]);
    hk1[k] = cs[k] * hk1[k] + sn[k] * hk1[k+1];
    hk1[k+1] = Cplx::new(0.0, 0.0);
}

fn givens_rotation(v1: &Cplx, v2: &Cplx) -> (Cplx, Cplx) {
    let c_zero = Cplx::new(0.0, 0.0);
    let c_one = Cplx::new(1.0, 0.0);
    let signf = if *v1 == c_zero {c_one} else {v1 / v1.norm()};
    let r = (v1.norm_sqr() + v2.norm_sqr()).sqrt();
    let cs = Cplx::new(v1.norm() / r, 0.0);
    let sn = signf * v2.conj() / r;
    return (cs, sn)
}