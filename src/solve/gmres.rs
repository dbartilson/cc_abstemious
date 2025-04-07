use na::{DMatrix, DVector, Normed};
use crate::Cplx;

use super::h_matrix;

enum ExitFlag {
    Error,
    Tolerance,
    Iterations
}

pub struct GMRES{
    max_it: usize,
    max_it_per_restart: usize,
    thresh: f64,
    num_mv: usize,
    pub a: Option<DMatrix::<Cplx>>,
    pub hmatrix: Option<h_matrix::HMatrix>
}

impl GMRES{
    pub fn new(max_it: usize, thresh: f64) -> GMRES {
        GMRES {
            max_it: max_it,
            max_it_per_restart: 0,
            thresh: thresh,
            num_mv: 0,
            a: None,
            hmatrix: None
        }
    }
    /// Solve the system of equations in-place for a given RHS 'x'
    /// This uses the GMRES(k) method, restarting after k iterations
    /// where k ~ sqrt(max_iterations)
    pub fn solve(&mut self, x: &mut DVector<Cplx>) {
        info!(" Using iterative (GMRES) solver...");
        let b = x.clone();
        // zero as initial guess
        x.fill(Cplx::new(0.0, 0.0));
        // if max_it is zero, set to lesser of {default} or num_eqn
        if self.max_it <= 0 {self.max_it = 1000;}
        self.max_it = std::cmp::min(self.max_it, self.get_num_eqn());
        // if thresh is zero, set to {default}
        if self.thresh <= 0.0 {self.thresh = 1.0e-5;}
        // set max_it_per_restart to ~ sqrt(max_iterations)
        let k = (self.max_it as f64).sqrt().ceil() as usize;
        let max_restarts = self.max_it / k;
        self.max_it_per_restart = k;
        for i in 0..max_restarts {
            info!("  GMRES restart: {}", i);
            // run GMRES algorithm, return updated solution
            let flag = self.gmres(x, &b);
            match flag {
                ExitFlag::Tolerance => {break;},
                _ => {}
            }
        }
        info!("  Number of matrix-vector products: {}", self.num_mv);
    }
    /*
    Set up gemv, get_num_eqn, and get_norm as methods which can access either the 
    full (dense) representation or the approximate (ACA) representation
    */
    fn gemv(&mut self, alpha: Cplx, x: &DVector::<Cplx>, beta: Cplx, b: &mut DVector::<Cplx>) {
        // Computes b = alpha * self * x + beta * b, where a is a matrix, x a vector, and alpha, beta two scalars
        self.num_mv += 1; // keep track of number of matrix-vector products
        if self.a.is_some() {b.gemv(alpha, &self.a.as_ref().unwrap(), x, beta); return}
        if self.hmatrix.is_some() {self.hmatrix.as_ref().unwrap().gemv(alpha, x, beta, b); return}
    }
    fn get_num_eqn(&self) -> usize {
        if self.a.is_some() { return self.a.as_ref().unwrap().shape().0} else {}
        if self.hmatrix.is_some() { return self.hmatrix.as_ref().unwrap().get_num_eqn()} else {}
        return 0
    }
    fn get_norm(&self) -> f64 {
        if self.a.is_some() { return self.a.as_ref().unwrap().norm()} else {}
        if self.hmatrix.is_some() { return self.hmatrix.as_ref().unwrap().get_norm()} else {}
        return 0.0
    }
    fn gmres(&mut self, x: &mut DVector<Cplx>, b: &DVector<Cplx>) -> ExitFlag {
        // see the example code at 
        // https://en.wikipedia.org/wiki/Generalized_minimal_residual_method#Regular_GMRES_(MATLAB_/_GNU_Octave)
        let mut flag = ExitFlag::Error;
        let n = x.len();
        let m = self.max_it_per_restart;
        let c_one = Cplx::new(1.0, 0.0);
        let c_zero = Cplx::new(0.0, 0.0);
        // use x as initial guess
        let b_norm = b.norm();
        // calculate initial residual r = b - Ax
        let mut r = b.clone();
        self.gemv(-c_one, x, c_one, &mut r);

        let alpha = self.get_norm();
        let mut error = Self::backward_error(r.norm(), alpha, x, b_norm);

        let mut sn = DVector::<Cplx>::from_element(m, c_zero);
        let mut cs = sn.clone();
        let mut e = vec![error];
        let mut q = DMatrix::<Cplx>::from_element(n, m+1, c_zero);
        q.set_column(0, &r.normalize());
        let mut qk1 = DVector::<Cplx>::from_element(n, c_zero);
        let mut beta = DVector::<Cplx>::from_element(m+1, c_zero);
        beta[0] = Cplx::new(r.norm(), 0.0);
        let mut h = DMatrix::<Cplx>::from_element(m+1, m, c_zero);
        let mut hk1 = DVector::<Cplx>::from_element(m+1, c_zero);
        for k in 0..m {
            if error.is_nan() { error!("GMRES error is NaN!"); }
            info!("   Iteration: {}, Error: {:10.3e}", k, error);
            self.arnoldi(&q, k, &mut qk1, &mut hk1);
            Self::apply_givens_rotation(&mut hk1, &mut cs, &mut sn, k);
            for i in 0..k+2 {
                h[(i, k)] = hk1[i];
            }
            q.set_column(k+1, &qk1);
            beta[k+1] = -sn[k].conj() * beta[k];
            beta[k] *= cs[k];

            error = Self::backward_error(beta[k+1].norm(), alpha, x, b_norm);
            e.push(error);

            if error < self.thresh || k == m-1 {
                // recompute backward error using exact arithmetic
                r = b.clone();
                Self::get_x(&h, &q, &beta, k+1, x);
                self.gemv(-c_one, x, c_one, &mut r);
                error = Self::backward_error(r.norm(), alpha, x, b_norm);
                *e.last_mut().unwrap() = error;
                if error < self.thresh {     
                    flag = ExitFlag::Tolerance;
                    info!("  GMRES tolerance acheived! ({:10.3e} < {:10.3e})", error, self.thresh);
                    break;
                }
            }
            if k == m-1 {
                flag = ExitFlag::Iterations;
            }
        }
        return flag;
    }
    /// Calculate backward error, see https://doi.org/10.1145/1067967.1067970
    fn backward_error(r: f64, alpha: f64, x: &DVector::<Cplx>, beta: f64) -> f64 {
        r / (alpha * x.norm() + beta)
    }
    /// After GMRES iterations, calculate the result from the Hessenberg matrix
    fn get_x(h: &DMatrix::<Cplx>, q: &DMatrix::<Cplx>, beta: &DVector::<Cplx>, k: usize, 
             x: &mut DVector::<Cplx>) {
        let c_one = Cplx::new(1.0, 0.0);
        let hkk = h.view((0, 0), (k, k));
        let hkklu = hkk.lu();
        let betak = beta.rows(0, k);
        // y = H^-1 * beta
        let y = hkklu.solve(&betak);
        // x = x + Q * y
        x.gemv(c_one, &q.columns(0,k), y.as_ref().unwrap(), c_one);
    }
    ///
    fn arnoldi(&mut self, q: &DMatrix<Cplx>, k: usize, 
            qk1: &mut DVector<Cplx>, hk1: &mut DVector<Cplx>)  {
        let c_zero = Cplx::new(0.0, 0.0);
        let c_one = Cplx::new(1.0, 0.0);
        // zero out and set qk1 = A * Q (Krylov vector)
        let qk = q.column(k).clone_owned();
        self.gemv(c_one, &qk, c_zero, qk1);
        hk1.fill(c_zero);
        // Modified Gram-Schmidt, keeping the Hessenberg matrix
        for i in 0..k+1 {
            hk1[i] = qk1.dotc(&q.column(i));
            qk1.axpy(-hk1[i], &q.column(i), c_one);
        }
        hk1[k+1] = Cplx::new(qk1.norm(), 0.0);
        let scale = c_one / hk1[k+1];
        for qi in qk1 {
            *qi *= scale;
        }
    }
    /// Use the Givens rotation angles to zero out the k+1 element of hk1
    fn apply_givens_rotation(hk1: &mut DVector<Cplx>, cs: &mut DVector<Cplx>, 
                             sn: &mut DVector<Cplx>, k: usize) {
        for i in 0..k {
            let temp = cs[i] * hk1[i] + sn[i] * hk1[i+1];
            hk1[i+1] = -sn[i].conj() * hk1[i] + cs[i] * hk1[i+1];
            hk1[i] = temp;
        }
        // update the next sin cos values for rotation
        (cs[k], sn[k]) = Self::givens_rotation(&hk1[k], &hk1[k+1]);
        // eliminate H(i + 1, i)
        hk1[k] = cs[k] * hk1[k] + sn[k] * hk1[k+1];
        hk1[k+1] = Cplx::new(0.0, 0.0);
    }
    /// Get givens rotation angles (cos, sin) for a specified [v1, v2] RHS to eliminate v2
    fn givens_rotation(v1: &Cplx, v2: &Cplx) -> (Cplx, Cplx) {
        // see https://www.netlib.org/lapack/lawnspdf/lawn148.pdf for more information
        // on reliable calculation of Givens rotation angles
        let c_zero = Cplx::new(0.0, 0.0);
        let c_one = Cplx::new(1.0, 0.0);
        let signf = if *v1 == c_zero {c_one} else {v1 / v1.norm()};
        let r = (v1.norm_sqr() + v2.norm_sqr()).sqrt();
        let cs = Cplx::new(v1.norm() / r, 0.0);
        let sn = signf * v2.conj() / r;
        return (cs, sn)
    }
}
