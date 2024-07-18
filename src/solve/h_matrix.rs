use std::{collections::HashMap, rc::Rc};
use na::{DMatrix, DVector};
use crate::{preprocess::mesh_data::Node, Cplx};

use super::aca::ACA;

#[derive(Clone)]
pub struct Cluster {
    u_bound: [f64;3],
    l_bound: [f64;3],
    diameter: f64,
    indices_contained: Vec<usize>, // node indices, not eqn indices
    sons: Vec<Rc<Cluster>>
}

impl Cluster {
    pub fn new_from(nodes: &Vec<Node>, 
                    indices_contained: Vec<usize>, 
                    leaf_cardinality: usize, 
                    eqn_map: &HashMap::<usize, usize>) -> Cluster {
        let mut cluster = Cluster {
            u_bound: [f64::NEG_INFINITY;3],
            l_bound: [f64::INFINITY;3],
            diameter: 0.0,
            indices_contained: indices_contained,
            sons: Vec::new()
        };
        cluster.process_cluster(nodes, leaf_cardinality, eqn_map);
        return cluster;
    }
    fn process_cluster(&mut self, 
                       nodes: &Vec<Node>, 
                       leaf_cardinality: usize, 
                       eqn_map: &HashMap::<usize, usize>) {
        self.update_bounds(nodes);
        self.update_diameter();
        self.map_nodes_to_eqns(eqn_map);
        if self.indices_contained.len() <= leaf_cardinality {
            // Once the cluster is finalized, map from node indices to eqn indices
            self.map_nodes_to_eqns(eqn_map);
            return;
        }
        // Determine which direction to split the cluster
        let mut jmax: usize = 0;
        let mut diff = 0.0;
        let alpha = &mut self.l_bound;
        let beta = &mut self.u_bound;
        for j in 0_usize..3 {
            let c = beta[j] - alpha[j];
            if c > diff {jmax = j; diff = c;}
        }
        // Split cluster at gamma in jmax direction
        let gamma = 0.5 * (alpha[jmax] + beta[jmax]);
        let mut indx1 = Vec::<usize>::new();
        let mut indx2 = indx1.clone();
        for index in &self.indices_contained {
            if nodes[*index].coords[jmax] <= gamma {
                indx1.push(*index);
            }
            else {
                indx2.push(*index);
            }
        }
        //self.indices_contained = Vec::new();
        self.sons.push(Rc::new(Cluster::new_from(nodes, indx1, leaf_cardinality, eqn_map)));
        self.sons.push(Rc::new(Cluster::new_from(nodes, indx2, leaf_cardinality, eqn_map)));
    }
    fn is_leaf(&self) -> bool { return self.sons.is_empty()}
    fn get_diameter(&self) -> f64 { return self.diameter;}
    fn update_bounds(&mut self, nodes: &Vec<Node>) {
        let alpha = &mut self.l_bound;
        let beta = &mut self.u_bound;
        for index in &self.indices_contained {
            for j in 0_usize..3 {
                alpha[j] = f64::min(nodes[*index].coords[j], alpha[j]);
                beta[j] = f64::max(nodes[*index].coords[j], beta[j]);
            }
        }
    }
    fn update_diameter(&mut self) {
        let diam = &mut self.diameter;
        *diam = 0.0;
        for j in 0..3 {
            *diam += f64::powi(self.u_bound[j] - self.l_bound[j],2);
        }
        *diam = diam.sqrt();
    }
    fn get_distance(c1: &Cluster, c2: &Cluster) -> f64 {
        let mut dist = 0.0;
        for j in 0..3 {
            dist += f64::powi(f64::max(0.0, c1.l_bound[j] - c2.u_bound[j]),2);
            dist += f64::powi(f64::max(0.0, c2.l_bound[j] - c1.u_bound[j]),2);
        }
        return dist.sqrt();
    }
    fn map_nodes_to_eqns(&mut self, eqn_map: &HashMap::<usize, usize>) {
        for idx in &mut self.indices_contained {
            if let Some(new_idx) = eqn_map.get(idx) {*idx = *new_idx;}
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, rc::Rc};
    use na::{DMatrix, Vector3};
    use crate::{preprocess::mesh_data::Node, solve::{h_matrix, tests::generate_random_ab}, Cplx};
    use super::{Block, Cluster};

    #[test]
    fn build_cluster_tree() {
        let mut hmap = HashMap::<usize, usize>::new();
        let mut nodes = Vec::<Node>::new();
        let mut i = 0;
        for j in 0..25 {
            for k in 0..25 {
                nodes.push(Node { id: i, coords: Vector3::new(j as f64, k as f64, 0.0)});
                hmap.insert(i, i);
                i += 1;
            }
        }
        let tree = Cluster::new_from(&nodes, (0..nodes.len()).collect(), 32, &hmap);
        approx::assert_relative_eq!(tree.get_diameter(), 33.94, epsilon = 1e-2);
    }
    #[test]
    fn build_block_tree() {
        let mut hmap = HashMap::<usize, usize>::new();
        let mut nodes = Vec::<Node>::new();
        let mut i = 0;
        for j in 0..25 {
            for k in 0..25 {
                nodes.push(Node { id: i, coords: Vector3::new(j as f64, k as f64, 0.0)});
                hmap.insert(i, i);
                i += 1;
            }
        }
        let cluster_tree = Rc::new(Cluster::new_from(&nodes, (0..nodes.len()).collect(), 32, &hmap));
        let block_tree = Block::new_from(cluster_tree.clone(), cluster_tree.clone(), 4.0);
        assert_eq!(block_tree.children[1].children[1].children[1].admissible, true)
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
    fn mult_hmatrix() {
        let mut hmap = HashMap::<usize, usize>::new();
        let mut nodes = Vec::<Node>::new();
        let n = 20;
        let mut i = 0;
        for j in 0..n {
            for k in 0..n {
                nodes.push(Node { id: i, coords: Vector3::new(j as f64, k as f64, 0.0)});
                hmap.insert(i, i);
                i += 1;
            }
        }
        let n2 = n*n;
        let (mut a, b) = generate_random_ab(n2, 10);
        for i in 0..n2 {
            for j in 0..n2 {
                let ii = i as f64;
                let jj = j as f64;
                a[(i,j)] += 1.0 / f64::max(0.5, (ii-jj)*(ii-jj));
            }
        }
        let get_row_or_column = |i,j| get_row_or_column(&a, i, j);
        // build ACA of matrix and compare norms
        let hm = h_matrix::HMatrix::new_from(n, get_row_or_column, &nodes, 32, &hmap);
        // compare matrix multiplication against random vector for both
        let mut x1 = b.clone();
        x1.gemv(Cplx::new(1.0, 0.0), &a, &b, Cplx::new(0.0, 0.0));
        let mut x2 = b.clone();
        hm.gemv(Cplx::new(1.0, 0.0), &b, Cplx::new(0.0, 0.0), &mut x2);
        // calculate norm of difference
        let f = (x2 - x1.clone()).norm() / x1.norm();
        approx::assert_relative_eq!(f, 0.0, epsilon = 1.0e-10);
    }
}

pub struct Block {
    rows: Rc<Cluster>,
    columns: Rc<Cluster>,
    admissible: bool,
    children: Vec<Block>
}

impl Block {
    pub fn new_from(rows: Rc<Cluster>, columns: Rc<Cluster>, eta: f64) -> Block {
        let mut tree = Block {
            rows: rows,
            columns: columns,
            admissible: false,
            children: Vec::new()
        };
        tree.process_block(eta);
        return tree;
    }
    /// Check if block can be divided (is not admissible, neither cluster is a leaf)
    #[inline]
    fn is_divisible(&self) -> bool {
        !(self.admissible || self.rows.is_leaf() || self.columns.is_leaf())
    }
    fn process_block(&mut self, eta: f64) {
        // update admissibility 
        let diam1 = self.rows.get_diameter();
        let diam2 = self.columns.get_diameter();
        let dist12 = Cluster::get_distance(&self.rows, &self.columns);
        self.admissible = f64::min(diam1, diam2) <= eta * dist12;
        // check divisibility, if divisible, do so
        if self.is_divisible() {
            // divide into 4 blocks according to sons of row/son clusters
            for rson in &self.rows.sons {
                for cson in &self.columns.sons {
                    self.children.push(Block::new_from(rson.clone(),cson.clone(), eta));
                }
            }
        }
    }
}


/// Can be represented in reduced format (e.g., ACA)
struct AdmissibleBlock {
    rows: Vec<usize>,
    columns: Vec<usize>,
    values: ACA
}
impl AdmissibleBlock {
    fn new<F>(rows: Vec<usize>, columns: Vec<usize>, get_row_or_column: F) -> AdmissibleBlock
    where F: Fn(Vec<usize>, Vec<usize>) -> Vec::<Cplx> {
        let aca = ACA::new(1.0e-5, 10, rows.len(), columns.len(), 
        |i| get_row_or_column(vec![i], columns.clone()),
        |j| get_row_or_column(rows.clone(), vec![j]));
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
}

/// Must be represented fully (in dense form)
struct InadmissibleBlock {
    rows: Vec<usize>,
    columns: Vec<usize>,
    values: DMatrix::<Cplx>
}
impl InadmissibleBlock {
    fn new<F>(rows: Vec<usize>, columns: Vec<usize>, get_row: F) -> InadmissibleBlock
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
}

/// Hierarchical matrix built from block tree
pub struct HMatrix {
    num_eqn: usize,
    norm: f64,
    admissible_blocks: Vec<AdmissibleBlock>,
    inadmissible_blocks: Vec<InadmissibleBlock>
}

impl HMatrix {
    pub fn new_from<F>(n: usize, 
                   get_row_or_column: F, 
                   nodes: &Vec<Node>, 
                   leaf_cardinality: usize, 
                   eqn_map: &HashMap::<usize, usize>) -> HMatrix 
    where F: Fn(Vec<usize>, Vec<usize>) -> Vec::<Cplx> {
        let cluster_tree = Rc::new(Cluster::new_from(&nodes, (0..nodes.len()-1).collect(), leaf_cardinality, eqn_map));
        let block_tree = Block::new_from(cluster_tree.clone(), cluster_tree.clone(), 4.0);
        let mut mat = HMatrix {
            num_eqn: n,
            norm: 0.0,
            admissible_blocks: Vec::new(),
            inadmissible_blocks: Vec::new()
        };
        mat.load_from(block_tree, &get_row_or_column);
        mat.update_norm();
        return mat;
    }
    pub fn get_num_eqn(&self) -> usize { self.num_eqn }
    pub fn get_norm(&self) -> f64 { self.norm }
    /// Process block tree into flat vectors of admissible and inadmissible blocks
    fn load_from<F>(&mut self, block: Block, get_row_or_column: &F)
    where F: Fn(Vec<usize>, Vec<usize>) -> Vec::<Cplx> {
        if block.admissible {
            self.admissible_blocks.push(
                AdmissibleBlock::new(
                    block.rows.indices_contained.clone(),
                    block.columns.indices_contained.clone(),
                    get_row_or_column
                )
            );
        }
        else if block.children.is_empty() {
            self.inadmissible_blocks.push(
                InadmissibleBlock::new(
                    block.rows.indices_contained.clone(),
                    block.columns.indices_contained.clone(),
                    |i: usize| -> Vec<Cplx> {get_row_or_column(vec![i], block.columns.indices_contained.clone())}
                )
            );
        }
        else {
            for child in block.children {
                self.load_from(child, get_row_or_column);
            }
        }
    }
    fn update_norm(&mut self) {
        for block in &self.inadmissible_blocks {
            self.norm += block.values.norm();
        }
        for block in &self.admissible_blocks {
            self.norm += block.values.get_norm();
        }
    }
    /// Computes b = alpha * self * x + beta * b, where a is a matrix, x a vector, and alpha, beta two scalars
    pub fn gemv(&self, alpha: Cplx, x: &DVector::<Cplx>, beta: Cplx, b: &mut DVector::<Cplx>) {
        for i in 0..b.len() {
            b[i] *= beta;
        }
        for block in &self.inadmissible_blocks {
            block.gemv(alpha, &x, Cplx::new(1.0, 0.0), b);
        }
        for block in &self.admissible_blocks {
            block.gemv(alpha, &x, Cplx::new(1.0, 0.0), b);
        }
    }
}