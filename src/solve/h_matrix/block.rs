use std::rc::Rc;

use super::cluster::Cluster;

/// Block Tree, which is used to compare distances between clusters and determine
/// admissibility, i.e., are they separated enough to approximate their interaction
pub struct BlockTree {
    rows: Rc<Cluster>,
    columns: Rc<Cluster>,
    admissible: bool,
    children: Vec<BlockTree>
}

impl BlockTree {
    /// Build 'Block' of interaction between rows and columns, each represented by a cluster
    /// eta controls the required separation for admissibility, though it is not particularly
    /// sensitive. Recommend eta [4-10] based on http://dx.doi.org/10.3970/cmes.2009.043.149
    pub fn new_from(rows: Rc<Cluster>, columns: Rc<Cluster>, eta: f64) -> BlockTree {
        let mut tree = BlockTree {
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
    fn is_admissible(&self) -> bool {return self.admissible;}
    fn get_children(&self) -> &Vec<BlockTree> {return &self.children;}
    fn get_row_indices(&self) -> &Vec<usize> { return self.rows.get_indices();}
    fn get_column_indices(&self) -> &Vec<usize> { return &self.columns.get_indices();}
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
                    self.children.push(BlockTree::new_from(rson.clone(), cson.clone(), eta));
                }
            }
        }
    }
}

pub struct Block {
    rows: Vec<usize>,
    columns: Vec<usize>,
    admissible: bool
}

impl Block {
    pub fn get_row_indices(&self) -> &Vec<usize> { return &self.rows; }
    pub fn get_column_indices(&self) -> &Vec<usize> { return &self.columns;}
    pub fn is_admissible(&self) -> bool { return self.admissible; }
}

// Flattened version of block tree, all held in one vector
pub struct BlockList {
    list: Vec<Block>
}

impl BlockList {
    /// Create a block list from a block tree (flatten)
    pub fn new_from(block: &BlockTree) -> BlockList {
        let mut result = BlockList { list: Vec::new() };
        result.load_from(block);
        return result;
    }
    /// Recursive call to load children from the block tree into vector
    fn load_from(&mut self, block: &BlockTree) {
        if block.get_children().is_empty() {
            self.list.push(Block { 
                rows: block.get_row_indices().clone(), 
                columns: block.get_column_indices().clone(), 
                admissible: block.is_admissible() });
        }
        else {
            for child in block.get_children() {
                self.load_from(child);
            }
        }
    }
    pub fn get_list(&self) -> &Vec<Block> { return &self.list }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, rc::Rc};
    use na::Vector3;
    use crate::{preprocess::mesh_data::CollocationPoint, solve::h_matrix::{block::{BlockList, BlockTree}, cluster::Cluster}};

    #[test]
    fn build_block_tree() {
        let mut hmap = HashMap::<usize, usize>::new();
        let mut cpts = Vec::<CollocationPoint>::new();
        let mut i = 0;
        for j in 0..25 {
            for k in 0..25 {
                cpts.push(CollocationPoint { 
                    id: i, 
                    coords: Vector3::new(j as f64, k as f64, 0.0),
                    normal: Vector3::from_element(0.0),
                    area: 0.0,
                    wt: 0.0});
                hmap.insert(i, i);
                i += 1;
            }
        }
        let cluster_tree = Rc::new(Cluster::new_from(&cpts, (0..cpts.len()).collect(), 32, &hmap));
        let block_tree = BlockTree::new_from(cluster_tree.clone(), cluster_tree.clone(), 4.0);
        assert_eq!(block_tree.children[1].children[1].children[1].admissible, true)
    }
    #[test]
    fn flatten_block_tree() {
        let mut hmap = HashMap::<usize, usize>::new();
        let mut cpts = Vec::<CollocationPoint>::new();
        let mut i = 0;
        for j in 0..25 {
            for k in 0..25 {
                cpts.push(CollocationPoint { 
                    id: i, 
                    coords: Vector3::new(j as f64, k as f64, 0.0),
                    normal: Vector3::from_element(0.0),
                    area: 0.0,
                    wt: 0.0});
                hmap.insert(i, i);
                i += 1;
            }
        }
        let cluster_tree = Rc::new(Cluster::new_from(&cpts, (0..cpts.len()).collect(), 32, &hmap));
        let block_tree = BlockTree::new_from(cluster_tree.clone(), cluster_tree.clone(), 4.0);
        let block_list = BlockList::new_from(&block_tree);
        assert_eq!(block_list.list[5].admissible, true);
    }
}