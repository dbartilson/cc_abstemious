/*!
Hierarchical matrix decomposition

This build a matrix decomposition using the following steps
    1. Building the cluster tree from node locations
    2. Building a block tree indicating admissibility of the blocks 
        (i.e., whether the clusters are distant enough apart to have their 
         block of equations/interactions approximated by an ACA approximation)
    3. Flattening the block tree into a block list
    4. Forming the admissible and inadmissible blocks into the HMatrix in parallel
        - Admissible blocks are decomposed by ACA
        - Inadmissible blocks (think diagonal, close together points) are left as full/dense representation

This results in a much smaller memory footprint and faster matrix-vector multiplies
*/

pub mod block;
pub mod aca;

use std::{collections::HashMap, rc::Rc, sync::{Arc, Mutex}};
use na::{DMatrix, DVector};
use scoped_threadpool::Pool;
use crate::{preprocess::mesh::CollocationPoint, tools, Cplx};
use block::cluster::Cluster;
use block::{BlockTree, BlockList};

/// Matrix represented in reduced format: ACA
#[derive(Debug)]
struct AdmissibleBlock {
    rows: Vec<usize>,
    columns: Vec<usize>,
    values: aca::ACA
}
impl AdmissibleBlock {
    fn new<F>(rows: Vec<usize>, columns: Vec<usize>, get_row_or_column: F, tolerance: f64) -> AdmissibleBlock
    where F: Fn(Vec<usize>, Vec<usize>) -> Vec::<Cplx> {
        let aca = aca::ACA::new(tolerance, rows.len(), columns.len(), 
        &|i| get_row_or_column(vec![rows[i]], columns.clone()),
        &|j| get_row_or_column(rows.clone(), vec![columns[j]]));
        AdmissibleBlock {
            rows: rows,
            columns: columns,
            values: aca
        }
    }
    /// Perform gather, multiply, scatter for this contribution to b = alpha * self * x + beta * b
    /// Note that beta is not used, assumed to be done on top level of hierarchical matrix
    fn gemv(&self, alpha: Cplx, x: &DVector::<Cplx>, _beta: Cplx, b: &mut DVector::<Cplx>) {
        // gather x
        let czero = Cplx::new(0.0, 0.0);
        let mut x1 = DVector::<Cplx>::from_element(self.columns.len(), czero);
        for (j, column) in self.columns.iter().enumerate() {
            x1[j] = x[*column];
        }
        let mut b1 = DVector::<Cplx>::from_element(self.rows.len(), czero);
        // multiply
        self.values.gemv(alpha, &x1, czero, &mut b1);
        // scatter
        for (i, row) in self.rows.iter().enumerate() {
            b[*row] += b1[i];
        }
    }
    /// For testing/debugging
    fn to_full(&self, a: &mut DMatrix::<Cplx>) {
        let ai = self.values.to_full();
        for (i, row) in self.rows.iter().enumerate() {
            for (j, column) in self.columns.iter().enumerate() {
                a[(*row, *column)] += ai[(i,j)];
            }
        }
    }
}

/// Must be represented fully (in dense form)
#[derive(Debug)]
struct InadmissibleBlock {
    rows: Vec<usize>,
    columns: Vec<usize>,
    values: DMatrix::<Cplx>
}
impl InadmissibleBlock {
    fn new<F>(rows: Vec<usize>, columns: Vec<usize>, get_row: &F) -> InadmissibleBlock
    where F: Fn(usize) -> Vec::<Cplx> {
        let mut values = DMatrix::<Cplx>::from_element(rows.len(), columns.len(), Cplx::new(0.0,0.0));
        for (i, row_index) in rows.iter().enumerate() {
            let row = get_row(*row_index);
            for (j, val) in row.iter().enumerate() {
                values[(i,j)] = *val;
            }
        }
        InadmissibleBlock {
            rows: rows,
            columns: columns,
            values: values
        }
    }
    /// Perform gather, multiply, scatter for this contribution to b = alpha * self * x + beta * b
    /// Note that beta is not used, assumed to be done on top level of hierarchical matrix
    fn gemv(&self, alpha: Cplx, x: &DVector::<Cplx>, _beta: Cplx, b: &mut DVector::<Cplx>) {
        // gather x
        let czero = Cplx::new(0.0, 0.0);
        let mut x1 = DVector::<Cplx>::from_element(self.columns.len(), czero);
        for (j, column) in self.columns.iter().enumerate() {
            x1[j] = x[*column];
        }
        let mut b1 = DVector::<Cplx>::from_element(self.rows.len(), czero);
        // multiply
        b1.gemv(alpha, &self.values, &x1, czero);
        // scatter
        for (i, row) in self.rows.iter().enumerate() {
            b[*row] += b1[i];
        }
    }
    #[allow(dead_code)]
    fn to_full(&self, a: &mut DMatrix::<Cplx>) {
        for (i, row) in self.rows.iter().enumerate() {
            for (j, column) in self.columns.iter().enumerate() {
                a[(*row, *column)] += self.values[(i,j)];
            }
        }
    }
}

/// Hierarchical matrix built from block tree
pub struct HMatrix {
    num_eqn: usize,
    norm: f64,
    admissible_blocks: Vec<AdmissibleBlock>,
    inadmissible_blocks: Vec<InadmissibleBlock>
}

impl HMatrix {
    /// Default 
    pub fn new() -> HMatrix {
        return HMatrix{
            num_eqn: 0,
            norm: 0.0,
            admissible_blocks: Vec::new(),
            inadmissible_blocks: Vec::new()
        }
    }
    /// Generate a hierarchical matrix from collocation points and functional
    pub fn new_from<F>(n: usize, 
                   get_row_or_column: &F, 
                   cpts: &Vec<CollocationPoint>, 
                   eqn_map: &HashMap::<usize, usize>,
                   leaf_cardinality: usize, 
                   tolerance: f64) -> HMatrix 
    where F: Fn(Vec<usize>, Vec<usize>) -> Vec::<Cplx> + std::marker::Sync {
        info!("  Building hierarchical matrix decomposition...");
        let cluster_tree = Rc::new(Cluster::new_from(&cpts, (0..cpts.len()).collect(), leaf_cardinality, eqn_map));
        let block_tree = BlockTree::new_from(cluster_tree.clone(), cluster_tree.clone(), 4.0);
        let block_list = BlockList::new_from(&block_tree);
        let mut mat = HMatrix {
            num_eqn: n,
            norm: 0.0,
            admissible_blocks: Vec::new(),
            inadmissible_blocks: Vec::new()
        };
        mat.load_from(block_list, get_row_or_column, tolerance);
        mat.update_norm();
        mat.print_stats();
        return mat;
    }
    pub fn get_num_eqn(&self) -> usize { self.num_eqn }
    /// Get matrix norm, used for iterative stopping criterion
    pub fn get_norm(&self) -> f64 { self.norm }
    /// Process block tree into admissible and inadmissible blocks
    fn load_from<F>(&mut self, blocklist: BlockList, get_row_or_column: &F, tolerance: f64)
    where F: Fn(Vec<usize>, Vec<usize>) -> Vec::<Cplx> + std::marker::Sync {
        let num_threads = tools::get_num_threads();
        // use a parallel pool of threads
        info!("  Using {} threads...", num_threads);
        let mut pool = Pool::new(num_threads as u32);
        let ad = Arc::new(Mutex::new(Vec::<AdmissibleBlock>::new()));
        let iad = Arc::new(Mutex::new(Vec::<InadmissibleBlock>::new()));
        pool.scoped(|scope| {
            for block in blocklist.get_list() {
                scope.execute(|| {
                    if block.is_admissible() {
                        let a = AdmissibleBlock::new(
                            block.get_row_indices().clone(),
                            block.get_column_indices().clone(),
                            get_row_or_column,
                            tolerance
                        );
                        let mut ad = ad.lock().unwrap();
                        ad.push(a);
                    }
                    else {
                        let ia = InadmissibleBlock::new(
                            block.get_row_indices().clone(),
                            block.get_column_indices().clone(),
                            &|i: usize| -> Vec<Cplx> {get_row_or_column(vec![i], 
                                block.get_column_indices().clone())}
                        );
                        let mut iad = iad.lock().unwrap();
                        iad.push(ia);
                    }
                });
            }
        });
        self.admissible_blocks = Arc::try_unwrap(ad).unwrap().into_inner().unwrap();
        self.inadmissible_blocks = Arc::try_unwrap(iad).unwrap().into_inner().unwrap();
    }
    fn update_norm(&mut self) {
        for block in &self.inadmissible_blocks {
            self.norm += block.values.norm();
        }
        for block in &self.admissible_blocks {
            self.norm += block.values.get_norm();
        }
    }
    fn print_stats(&self) {
        let adl = self.admissible_blocks.len();
        let iadl = self.inadmissible_blocks.len();
        let admissible_ratio = 100.0 * (adl as f64) / ((iadl + adl) as f64);
        // calculate compression ratio
        // equal to storage size for ACA + full blocks
        // divided by total size of matrix in full form
        let mut numerator: usize = 0;
        let denominator = self.num_eqn * self.num_eqn;
        for a in &self.admissible_blocks {
            numerator += a.values.get_num_uv() * (a.columns.len() + a.rows.len());
        }
        for ia in &self.inadmissible_blocks {
            numerator += ia.columns.len() * ia.rows.len();
        }
        let compression_ratio = 100.0 * (1.0 - (numerator as f64 / denominator as f64));
        info!("  Decomposition info:");
        info!("   Admissible block ratio: {:4.1}%", admissible_ratio);
        info!("   Compression ratio: {:4.1}%", compression_ratio);
    }
    /// Computes b = alpha * self * x + beta * b, where a is a matrix, x a vector, and alpha, beta two scalars
    pub fn gemv(&self, alpha: Cplx, x: &DVector::<Cplx>, beta: Cplx, b: &mut DVector::<Cplx>) {
        if beta != Cplx::new(1.0,0.0) {
            for i in 0..b.len() {
                b[i] *= beta;
            }
        }
        if self.inadmissible_blocks.is_empty() && self.admissible_blocks.is_empty() {return;}
        if self.num_eqn != b.len() || self.num_eqn != x.len() {
            error!("Dimension mismatch in H matrix gemv");
        }
        for block in &self.inadmissible_blocks {
            block.gemv(alpha, &x, Cplx::new(1.0, 0.0), b);
        }
        for block in &self.admissible_blocks {
            block.gemv(alpha, &x, Cplx::new(1.0, 0.0), b);
        }
    }
    #[allow(dead_code)]
    fn to_full(&self) -> DMatrix::<Cplx> {
        let mut h = DMatrix::<Cplx>::from_element(
            self.num_eqn, self.num_eqn, Cplx::new(0.0,0.0));
        for block in &self.inadmissible_blocks {
            block.to_full(&mut h);
        }
        for block in &self.admissible_blocks {
            block.to_full(&mut h);
        }
        return h;
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use na::{DMatrix, DVector, Vector3};
    use crate::{preprocess::mesh::CollocationPoint, solve::{h_matrix, tests::generate_random_ab}, Cplx};

    fn get_hyperbolic_matrix(n: usize) -> (Vec::<CollocationPoint>, HashMap::<usize, usize>,
        DMatrix::<Cplx>, DVector::<Cplx>) {
        let mut hmap = HashMap::<usize, usize>::new();
        let mut nodes = Vec::<CollocationPoint>::new();
        let mut i = 0;
        for j in 0..n {
            for k in 0..n {
                nodes.push(CollocationPoint { 
                    id: i, 
                    coords: Vector3::new(j as f64 / n as f64, k as f64 / n as f64, 0.0),
                    normal: Vector3::from_element(0.0),
                    area: 0.0,
                    wt: 0.0});
                hmap.insert(i, i);
                i += 1;
            }
        }
        let n2 = n*n;
        let (mut a, b) = generate_random_ab(n2, n2, 10);
        for i in 0..n2 {
            for j in 0..n2 {
                let ii = &nodes[i].coords;
                let jj = &nodes[j].coords;
                let dist = (ii-jj).norm_squared();
                a[(i,j)] = Cplx::new(1.0 / f64::max(1e-5, dist), 0.0);
            }
        }
        return (nodes, hmap, a, b)
    }
    fn get_row_or_column(a: &DMatrix::<Cplx>, i: Vec<usize>, j: Vec<usize>) -> Vec<Cplx> {
        if i.len() == 1 {
             let mut b = Vec::<Cplx>::new();
             for ji in j {
                b.push(a[(i[0],ji)]);
             }
             return b;
        }
        else {
            let mut b = Vec::<Cplx>::new();
            for ii in i {
               b.push(a[(ii,j[0])]);
            }
            return b;
        }; 
    }
    #[test]
    fn get_hmatrix() {
        let (cpts, hmap, a, b) = get_hyperbolic_matrix(10);
        let get_row_or_column = |i,j| get_row_or_column(&a, i, j);
        // build ACA of matrix and compare norms
        let hm = h_matrix::HMatrix::new_from(b.len(), &get_row_or_column, &cpts, &hmap, 20, 1e-4);
        let a_hm = hm.to_full();
        let eps = (a_hm - a.clone()).norm() / a.norm();
        approx::assert_relative_eq!(eps, 0.0, epsilon = 1.0e-6);
    }
    #[test]
    fn mult_hmatrix() {
        let (cpts, hmap, a, b) = get_hyperbolic_matrix(20);
        let get_row_or_column = |i,j| get_row_or_column(&a, i, j);
        // build ACA of matrix and compare norms
        let hm = h_matrix::HMatrix::new_from(b.len(), &get_row_or_column, &cpts, &hmap, 20, 1e-4);
        // compare matrix multiplication against random vector for both
        let mut x1 = b.clone();
        x1.gemv(Cplx::new(1.0, 0.0), &a, &b, Cplx::new(0.0, 0.0));
        let mut x2 = b.clone();
        hm.gemv(Cplx::new(1.0, 0.0), &b, Cplx::new(0.0, 0.0), &mut x2);
        // calculate norm of difference
        let f = (x2 - x1.clone()).norm() / x1.norm();
        approx::assert_relative_eq!(f, 0.0, epsilon = 1.0e-6);
    }
}
