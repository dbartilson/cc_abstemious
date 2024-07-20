use std::rc::Rc;

use super::cluster::Cluster;

/// Block Tree, which is used to compare distances between clusters and determine
/// admissibility, i.e., are they separated enough to approximate their interaction
pub struct Block {
    rows: Rc<Cluster>,
    columns: Rc<Cluster>,
    admissible: bool,
    children: Vec<Block>
}

impl Block {
    /// Build 'Block' of interaction between rows and columns, each represented by a cluster
    /// eta controls the required separation for admissibility, though it is not particularly
    /// sensitive. Recommend eta [4-10] based on http://dx.doi.org/10.3970/cmes.2009.043.149
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
    #[inline]
    pub fn is_admissible(&self) -> bool {return self.admissible;}
    pub fn get_children(&self) -> &Vec<Block> {return &self.children;}
    pub fn get_row_indices(&self) -> &Vec<usize> { return self.rows.get_indices();}
    pub fn get_column_indices(&self) -> &Vec<usize> { return &self.columns.get_indices();}
    fn process_block(&mut self, eta: f64) {
        // update admissibility 
        let diam1 = self.rows.get_diameter();
        let diam2 = self.columns.get_diameter();
        let dist12 = Cluster::get_distance(&self.rows, &self.columns);
        // admissibiliy criterion 
        self.admissible = f64::min(diam1, diam2) <= eta * dist12;
        // check divisibility, if divisible, do so
        if self.is_divisible() {
            // divide into 4 blocks according to sons of row/son clusters
            for rson in self.rows.get_sons() {
                for cson in self.columns.get_sons() {
                    self.children.push(Block::new_from(rson.clone(),cson.clone(), eta));
                }
            }
        }
    }
}#[cfg(test)]
mod tests {
    use std::{collections::HashMap, rc::Rc};
    use na::Vector3;
    use crate::{preprocess::mesh_data::Node, solve::h_matrix::{block::Block, cluster::Cluster}};

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
}