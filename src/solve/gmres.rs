use na::{DMatrix, DVector, Normed};
use crate::Cplx;

pub fn gmres(a: &DMatrix::<Cplx>, x: &mut DVector<Cplx>, b: &DVector<Cplx>, 
             max_it: usize, thresh: f64) {
    let n = x.len();
    let m = max_it;

    let c_one = Cplx::new(1.0, 0.0);
    let c_zero = Cplx::new(0.0, 0.0);
    // use x as initial guess
    let mut r = b.clone();
    r.gemv(c_one, a, x, c_one);
    let b_norm = b.norm();
    let mut error = r.norm() / b_norm;

    let mut sn = DVector::<Cplx>::from_element(m, c_zero);
    let mut cs = sn.clone();
    let mut e1 = DVector::<Cplx>::from_element(m+1, c_zero);
    e1[0] = c_one;
    let mut e = vec![error];
    let r_norm = r.norm();
    let mut q = DMatrix::<Cplx>::from_element(n, 1, c_zero);
    q.set_column(0, &r);
    q.scale_mut(1.0 / r_norm);
    let mut beta = e1.clone().scale(r_norm);
    let mut h = DMatrix::<Cplx>::from_element(0, 0, c_zero);
    let mut kk: usize = 0;
    for k in 0..m {
        kk += 1;
        let (mut hk1, qk1) = arnoldi(a, &q, k);
        h.set_column(k, &hk1);
        q.set_column(k+1, &qk1);
        apply_givens_rotation(&mut hk1, &mut cs, &mut sn, k);
        beta[k+1] = -sn[k] * beta[k];
        beta[k] = cs[k] * beta[k];
        error = beta[k+1].norm() / b_norm;
        e.push(error);

        if error < thresh {break;}
    }

    let hkk = h.view((0, 0), (kk, kk));
    let hkklu = hkk.lu();
    let y = hkklu.solve(&beta);
    x.gemv(c_one, &q.view((0,0), (n, kk)), y.as_ref().unwrap(), c_one)
}

fn arnoldi(a: &DMatrix::<Cplx>, q: &DMatrix<Cplx>, k: usize) 
        -> (DVector::<Cplx>, DVector::<Cplx>) {
    let n = a.shape().0;
    let c_zero = Cplx::new(0.0, 0.0);
    let c_one = Cplx::new(1.0, 0.0);
    let mut qk1 = DVector::<Cplx>::from_element(n, c_zero);
    qk1.gemv(c_one, a, &q.column(k), c_zero);
    let mut hk1 = DVector::<Cplx>::from_element(n, c_zero);
    for i in 0..k {
        hk1[i] = qk1.tr_dot(&q.column(i));
        qk1.axpy(-hk1[i], &q.column(i), c_one);
    }
    hk1[k] = Cplx::new(qk1.norm(), 0.0);
    qk1.scale_mut(1.0 / qk1.norm());
    return (hk1, qk1)
}

fn apply_givens_rotation(hk1: &mut DVector<Cplx>, cs: &mut DVector<Cplx>, 
                         sn: &mut DVector<Cplx>, k: usize) {

    for i in 0..k-1 {
        let temp = cs[i] * hk1[i] + sn[i] * hk1[i+1];
        hk1[i+1] = -sn[i] * hk1[i] + cs[i] * hk1[i+1];
        hk1[i] = temp;
    }
    (cs[k], sn[k]) = givens_rotation(&hk1[k], &hk1[k+1]);
    hk1[k] = cs[k] * hk1[k] + sn[k] * hk1[k+1];
    hk1[k+1] = Cplx::new(0.0, 0.0);
}

fn givens_rotation(v1: &Cplx, v2: &Cplx) -> (Cplx, Cplx) {
    let t = (v1.norm_sqr() + v2.norm()).sqrt();
    let cs = v1 / t;
    let sn = v2 / t;
    return (cs, sn)
}