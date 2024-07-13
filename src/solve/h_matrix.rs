use std::rc::Rc;

use crate::preprocess::mesh_data::Node;

const LEAF_CARDINALITY: usize = 32;

#[derive(Clone)]
pub struct Cluster {
    u_bound: [f64;3],
    l_bound: [f64;3],
    diameter: f64,
    indices_contained: Vec<usize>,
    sons: Vec<Rc<Cluster>>
}

impl Cluster {
    fn new() -> Cluster {
        Cluster {
            u_bound: [f64::NEG_INFINITY;3],
            l_bound: [f64::INFINITY;3],
            diameter: 0.0,
            indices_contained: Vec::new(),
            sons: Vec::new()
        }
    }
    pub fn build_tree(nodes: &Vec<Node>) -> Cluster {
        let mut tree = Cluster::new();
        tree.indices_contained = (0..nodes.len()).collect();
        tree.process_cluster(nodes);
        return tree
    }
    fn process_cluster(&mut self, nodes: &Vec<Node>) {
        self.update_bounds(nodes);
        self.update_diameter();
        if self.indices_contained.len() <= LEAF_CARDINALITY {return;}
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
        let mut son1 = Cluster::new();
        let mut son2 = Cluster::new();
        for index in &self.indices_contained {
            if nodes[*index].coords[jmax] <= gamma {
                son1.indices_contained.push(*index);
            }
            else {
                son2.indices_contained.push(*index);
            }
        }
        //self.indices_contained = Vec::new();
        son1.process_cluster(nodes);
        son2.process_cluster(nodes);
        self.sons.push(Rc::new(son1));
        self.sons.push(Rc::new(son2));
    }
    pub fn is_leaf(&self) -> bool { return self.sons.is_empty()}
    pub fn get_sons(&self) -> &Vec<Rc<Cluster>> { return &self.sons;}
    pub fn get_bounds(&self) -> (&[f64;3], &[f64;3]) {return (&self.l_bound, &self.u_bound)}
    pub fn get_indices(&self) -> &Vec<usize> { return &self.indices_contained;}
    pub fn get_diameter(&self) -> f64 { return self.diameter;}
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
    }
    pub fn get_distance(c1: &Cluster, c2: &Cluster) -> f64 {
        let mut dist = 0.0;
        for j in 0..3 {
            dist += f64::powi(f64::max(0.0, c1.l_bound[j] - c2.u_bound[j]),2);
            dist += f64::powi(f64::max(0.0, c2.l_bound[j] - c1.u_bound[j]),2);
        }
        return dist.sqrt();
    }
}

#[cfg(test)]
mod tests {
    use std::rc::Rc;
    use na::Vector3;
    use crate::preprocess::mesh_data::Node;
    use super::{Block, Cluster};

    #[test]
    fn build_cluster_tree() {
        let mut nodes = Vec::<Node>::new();
        let mut i = 0;
        for j in 0..25 {
            for k in 0..25 {
                nodes.push(Node { id: i, coords: Vector3::new(j as f64, k as f64, 0.0)});
                i += 1;
            }
        }
        let tree = Cluster::build_tree(&nodes);
        assert_eq!(tree.get_bounds().1[0], 24.0)
    }
    #[test]
    fn build_block_tree() {
        let mut nodes = Vec::<Node>::new();
        let mut i = 0;
        for j in 0..25 {
            for k in 0..25 {
                nodes.push(Node { id: i, coords: Vector3::new(j as f64, k as f64, 0.0)});
                i += 1;
            }
        }
        let cluster_tree = Rc::new(Cluster::build_tree(&nodes));
        let block_tree = Block::new(cluster_tree.clone(), cluster_tree.clone(), 4.0);
    }
}

pub struct Block {
    rows: Rc<Cluster>,
    columns: Rc<Cluster>,
    admissible: bool,
    children: Vec<Block>
}

impl Block {
    pub fn new(rows: Rc<Cluster>, columns: Rc<Cluster>, eta: f64) -> Block {
        let mut tree = Block {
            rows: rows,
            columns: columns,
            admissible: false,
            children: Vec::new()
        };
        tree.process_block(eta);
        return tree;
    }
    fn update_admissibility(&mut self, eta: f64) {
        let diam1 = self.rows.get_diameter();
        let diam2 = self.columns.get_diameter();
        let dist12 = Cluster::get_distance(&self.rows, &self.columns);
        self.admissible = f64::min(diam1, diam2) <= eta * dist12 
    }
    /// Check if block can be divided (is not admissible, neither cluster is a leaf)
    #[inline]
    fn is_divisible(&self) -> bool {
        !(self.admissible || self.rows.is_leaf() || self.columns.is_leaf())
    }
    fn process_block(&mut self, eta: f64) {
        self.update_admissibility(eta);
        if self.is_divisible() {
            for rson in &self.rows.sons {
                for cson in &self.columns.sons {
                    let cblock = Block::new(
                        rson.clone(),cson.clone(), eta);
                    self.children.push(cblock);
                }
            }
        }
    }
}


struct DenseBlock {
    rows: Vec<usize>,
    columns: Vec<usize>
}

struct CompressedBlock {
    rows: Vec<usize>,
    columns: Vec<usize>
}
pub struct HMatrix {
    num_eqn: usize,
    dense_blocks: Vec<DenseBlock>,
    compressed_blocks: Vec<CompressedBlock>
}

impl HMatrix {
    pub fn new(n: usize) -> HMatrix {
        HMatrix {
            num_eqn: n,
            dense_blocks: Vec::new(),
            compressed_blocks: Vec::new()
        }
    }
}