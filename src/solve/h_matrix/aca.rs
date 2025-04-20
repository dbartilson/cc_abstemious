/*!
Adaptive Cross Approximation, used to compress admissible blocks
*/

use std::collections::HashSet;
use na::{ComplexField, DMatrix, DVector};
use rand::{distr::{Distribution, Uniform}, SeedableRng, rngs::StdRng};
use crate::Cplx;

/// U and V terms for pseudo SVD representation
#[derive(Debug)]
struct UV {
    u: DVector::<Cplx>,
    v: DVector::<Cplx>
}
/// Approximation of a matrix/block
#[derive(Debug)]
pub struct ACA
{
    num_rows: usize,
    num_columns: usize,
    uv: Vec<UV>, // vector of UV, see below for how this corresponds to full matrix
    norm: f64,
}

impl ACA
{
    /// Decompose a matrix using ACA using functional and tolerance
    pub fn new<F,G>(tol: f64, m: usize, n: usize, get_row: &F, get_column: &G) -> ACA
        where F: Fn(usize) -> Vec::<Cplx>,
              G: Fn(usize) -> Vec::<Cplx>  {
        let mut a = ACA {
            num_rows: m,
            num_columns: n,
            uv: Vec::new(), 
            norm: 0.0};
        a.decompose(tol, &get_row, &get_column);
        return a
    }
    fn get_max_rank(&self) -> usize { std::cmp::min(self.num_rows, self.num_columns) }
    /// Get the current residual for the given row
    fn get_residual_row<F>(&self, get_row: &F, i: usize) -> DVector::<Cplx>
    where F: Fn(usize) -> Vec::<Cplx> {
        let rv = get_row(i);
        let mut r = DVector::<Cplx>::from_column_slice(&rv);
        for uv in &self.uv {
            r.axpy(-uv.u[i], &uv.v, Cplx::new(1.0, 0.0));
        }
        return r;
    }
    /// Get the current residual for the given column
    fn get_residual_column<G>(&self, get_column: &G, j: usize) -> DVector::<Cplx>
    where G: Fn(usize) -> Vec::<Cplx> {
        let rv = get_column(j);
        let mut r = DVector::<Cplx>::from_column_slice(&rv);
        for uv in &self.uv {
            r.axpy(-uv.v[j], &uv.u, Cplx::new(1.0, 0.0));
        }
        return r;
    }
    /// Get number of UV vectors in decomposition
    pub fn get_num_uv(&self) -> usize { return self.uv.len();}
    /// Calculate the Adaptive Cross Approximation for the given matrix
    fn decompose<F,G>(&mut self, tol: f64, get_row: &F, get_column: &G) 
    where F: Fn(usize) -> Vec::<Cplx>,
          G: Fn(usize) -> Vec::<Cplx> {
        //! largely adapted from https://tbenthompson.com/book/tdes/low_rank.html
        //! also see https://doi.org/10.1007/s00607-004-0103-1
        debug!("  Assembling Adaptive Cross Approximation matrix ({} x {})...", self.num_rows, self.num_columns);
        let n = self.get_max_rank();
        let mut rng = rand::rngs::StdRng::seed_from_u64(10);
        let drange = Uniform::new_inclusive(0, n-1).expect("failed to create rand distr");
        // get randomized starting row/col
        let mut iref = drange.sample(&mut rng);
        let mut jref = drange.sample(&mut rng);
        // list of previous row/column indices to exclude
        let mut istar_list = HashSet::<usize>::new();
        let mut jstar_list = istar_list.clone();
        // calculate residual for reference row/col
        let mut riref = self.get_residual_row(get_row, iref);
        let mut rjref = self.get_residual_column(get_column, jref);
        #[allow(unused_assignments)]
        let mut ristar = DVector::<Cplx>::from_element(1, Cplx::new(0.0, 0.0));
        #[allow(unused_assignments)]
        let mut rjstar = DVector::<Cplx>::from_element(1, Cplx::new(0.0, 0.0));
        let mut norm_f2_total = 0.0;
        let max_vec = std::cmp::min(self.num_columns, self.num_rows);
        for k in 0..max_vec{
            // get column index of maximum value of row
            let mut jstar = Self::index_max_exclude(&riref, &jstar_list);
            let jstar_val = riref[jstar].abs();
            // get row index of maximum value of column
            let mut istar = Self::index_max_exclude(&rjref, &istar_list);
            let istar_val = rjref[istar].abs();
            // choose whichever index set has larger value
            if istar_val > jstar_val {
                // get the residual row corresponding to istar
                ristar = self.get_residual_row(get_row,istar);
                // get the column index of maximum of ristar
                jstar = Self::index_max_exclude(&ristar, &jstar_list);
                // get the residual for this column
                rjstar = self.get_residual_column(get_column, jstar);
            }
            else {
                rjstar = self.get_residual_column(get_column, jstar);
                istar = Self::index_max_exclude(&rjstar, &istar_list);
                ristar = self.get_residual_row(get_row,istar);
            }
            istar_list.insert(istar);
            jstar_list.insert(jstar);
            // simplest way to scale, u_kj = ristar / ristar[rjstar]
            let factor = Cplx::new(1.0, 0.0) / ristar[jstar];
            for rs in &mut ristar {
                *rs *= factor;
            }
            // v_ki = rjstar
            self.uv.push(UV {u: rjstar, v: ristar});

            let norm_f2_k = self.update_norm_estimate(&mut norm_f2_total);
            let step_size = norm_f2_k.sqrt() / norm_f2_total.sqrt();
            debug!("   Step: {}, Row: {}, Column: {}, Step size: {}", k, istar, jstar, step_size);
            if step_size < tol / 50.0 { break; }
            // If pivoting occurred on reference row (iref == istar), get new random row
            if iref == istar {
                Self::get_random_index_exclude(&mut iref, &n, &istar_list, &mut rng, &drange);
                riref = self.get_residual_row(get_row, iref);
            }
            else {
                // update reference row with new approximation
                riref.axpy(-self.uv[k].u[iref], &self.uv[k].v, Cplx::new(1.0, 0.0));
            }
            if jref == jstar {
                Self::get_random_index_exclude(&mut jref, &n, &jstar_list, &mut rng, &drange);
                rjref = self.get_residual_column(get_column, jref);
            }
            else {
                // update reference column
                rjref.axpy(-self.uv[k].v[jref], &self.uv[k].u, Cplx::new(1.0, 0.0));
            }
        }
        self.norm = norm_f2_total.sqrt();
        self.svd_recompression(tol);
    }
    fn svd_recompression(&mut self, tol: f64) {
        //! Largely adapted from https://tbenthompson.com/book/tdes/low_rank.html
        //! Build U and V as matrix represenations for QR
        let czero = Cplx::new(0.0,0.0);
        let mut u = DMatrix::<Cplx>::from_element(self.num_rows, self.get_num_uv(), czero);
        let mut v = DMatrix::<Cplx>::from_element(self.num_columns, self.get_num_uv(), czero);
        for (i, uv) in self.uv.iter().enumerate() {
            u.set_column(i, &uv.u);
            v.set_column(i, &uv.v);
        }
        // Perform QR decomposition on U and V => U V^T = Q_u * R_u * R_v^T * Q_v^T
        let u_qr = u.qr();
        let q_u = u_qr.q();
        let r_u = u_qr.r();
        let v_qr = v.qr();
        let q_v = v_qr.q();
        let r_v_t = v_qr.r().transpose();
        // Compute R_u * R_v^T 
        let rr = r_u * r_v_t;
        // Compute SVD of R_u * R_v^T  => W Sigma Z^H
        let svd = rr.svd(true, true);
        let w = &svd.u.unwrap();
        let sigma = &svd.singular_values;
        let z_h = svd.v_t.unwrap();
        let z = z_h.transpose();
        // Find index of singular value below eps * cumulative sum
        let mut index: usize = 0;
        let mut sum = 0.0_f64;
        for s in sigma {
            sum += s;
            if *s < tol*sum {break;}
            index += 1;
        }
        // Deal with zero-index for size and iter
        index = std::cmp::max(1, std::cmp::min(sigma.len(), index+1));
        // Reassemble U * V^T = (Q_u * W * Sigma) * (Z^H * Q_v^T)
        self.uv.truncate(index);
        for i in 0..index {
            self.uv[i].u.gemv(Cplx::new(sigma[i], 0.0), &q_u, &w.column(i), czero);
            self.uv[i].v.gemv(Cplx::new(1.0, 0.0), &q_v, &z.column(i), czero);
        }
        self.calculate_norm_estimate();
    }
    /// Computes b = alpha * self * x + beta * b, where a is a matrix, x a vector, and alpha, beta two scalars
    pub fn gemv(&self, alpha: Cplx, x: &DVector::<Cplx>, beta: Cplx, b: &mut DVector::<Cplx>)  {
        if self.num_rows != b.len() || self.num_columns != x.len() {
            error!("Dimension mismatch in ACA gemv");
        }
        // M ~ sum u_k v_k^H
        // do first step separately to use beta scaling of b vector
        let dot_product = alpha * self.uv[0].v.dot(x);
        b.axpy(dot_product, &self.uv[0].u, beta);
        // now do the other vectors
        let c_one = Cplx::new(1.0, 0.0);
        for k in 1..self.uv.len() {
            let dot_product = alpha * self.uv[k].v.dot(x);
            b.axpy(dot_product, &self.uv[k].u, c_one);
        }
    }
    /// For testing/diagnostic, transform ACA representation to full (lossy)
    #[allow(dead_code)]
    pub fn to_full(&self) -> DMatrix::<Cplx> {
        let mut a = DMatrix::<Cplx>::from_element(
            self.num_rows, self.num_columns, Cplx::new(0.0,0.0));
        for uv in &self.uv {
            for i in 0..self.num_rows {
                for j in 0..self.num_columns {
                    a[(i,j)] += uv.u[i] * uv.v[j];
                }
            }
        }
        return a;
    }
    /// Return the Frobenius norm estimate from the ACA approximation
    pub fn get_norm(&self) -> f64 { self.norm }
    /// Calculate norm estimate ab initio
    fn calculate_norm_estimate(&mut self) {
        let mut norm_f2_total = 0.0;
        for (i, uv) in self.uv.iter().enumerate() {
            norm_f2_total += uv.u.norm_squared() * uv.v.norm_squared();
            for j in 0..i {
                let uvj = &self.uv[j];
                norm_f2_total += 2.0 * uv.u.dot(&uvj.u).abs() * uv.v.dot(&uvj.v).abs();
            }
        }
        self.norm = norm_f2_total.sqrt();
    }
    /// Update existing norm estimate using the most recent UV vectors
    fn update_norm_estimate(&self, norm_f2_total: &mut f64) -> f64 {
        // see Eq. 9 of https://doi.org/10.3970/cmes.2009.043.149 for Frobenius norm estimate
        match self.uv.last() {
            Some(uv) => {
                // squared norm of current step
                let norm_f2_k = uv.u.norm_squared() * uv.v.norm_squared();
                // add to total norm estimate
                *norm_f2_total += norm_f2_k;
                let n = self.uv.len();
                for i in 0..n-1 {
                    let uvi = &self.uv[i];
                    *norm_f2_total += 2.0 * uvi.u.dot(&uv.u).abs() * uvi.v.dot(&uv.v).abs();
                }
                return norm_f2_k;
            },
            None => {error!("UV last not found"); return 0.0}
        }
    }
    /// Check if desired index is already in not_allowed HashSet, 
    /// if so, randomly find a new index
    fn get_random_index_exclude(i: &mut usize, n: &usize, not_allowed: &HashSet<usize>, 
        rng: &mut StdRng, drange: &Uniform<usize>) {
        if not_allowed.contains(i) {
            *i = drange.sample(rng);
            for _ in 0..*n-1 {
                *i = (*i + 1) % n;
                if !not_allowed.contains(i) {break;}
            }
        }
    }
    /// Get the index corresponding to the max(abs()) of the array while also
    /// excluding any indices which are in not_allowed HashSet
    fn index_max_exclude(array: &DVector::<Cplx>, not_allowed: &HashSet<usize>) -> usize {
        let mut index = array.len();
        let mut val_abs: f64 = -1.0;
        for (i, &val) in array.iter().enumerate() {
            if not_allowed.contains(&i) { continue; }
            let ival = val.abs();
            if val_abs < ival {
                index = i;
                val_abs = ival;
            }
        }
        return index
    }
}

#[cfg(test)]
mod tests {
    extern crate approx;
    use na::DVector;

    use crate::solve::h_matrix::aca::ACA;
    use crate::Cplx;
    use crate::solve::tests::generate_random_ab;

    #[test]
    fn mult_aca() {
        let m = 20;
        let n = 100;
        let (mut a, b) = generate_random_ab(m, n, 13);
        //a.fill(Cplx::new(0.0,0.0));
        for i in 0..m {
            for j in 0..n {
                a[(i,j)] *= 1.0e-8;
            }
        }
        for k in 0..13 {
            let (_, c) = generate_random_ab(m, m, k);
            let (_, d) = generate_random_ab(n, n, k+2);
            for i in 0..m {
                for j in 0..n {
                    a[(i,j)] += c[i] * d[j];
                }
            }
        }
        let get_row = |i: usize| -> Vec<Cplx> {a.clone().row(i).iter().cloned().collect()};
        let get_col = |i: usize| -> Vec<Cplx> {a.clone().column(i).iter().cloned().collect()};
        // build ACA of matrix and compare norms
        let aca = ACA::new(1.0e-6, m, n, &get_row, &get_col);
        approx::assert_relative_eq!(aca.get_norm(), a.norm(), max_relative = 0.05);
        // compare matrix multiplication against random vector for both
        let mut x1 = DVector::<Cplx>::from_element(m, Cplx::new(0.0,0.0));
        x1.gemv(Cplx::new(1.0, 0.0), &a, &b, Cplx::new(0.0, 0.0));
        let mut x2 = x1.clone();
        aca.gemv(Cplx::new(1.0, 0.0), &b, Cplx::new(0.0, 0.0), &mut x2);
        // calculate norm of difference
        let f = (x2 - x1.clone()).norm() / x1.norm();
        approx::assert_relative_eq!(f, 0.0, epsilon = 1.0e-8);
    }
}