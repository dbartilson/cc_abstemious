use std::{collections::HashMap, rc::Rc};

use crate::preprocess::mesh_data::Node;

/// Cluster Tree, used to partition the surface in space until
/// they contain fewer nodes than the 'cardinality' parameter of the tree
#[derive(Clone)]
pub struct Cluster {
    u_bound: [f64;3],
    l_bound: [f64;3],
    diameter: f64,
    indices_contained: Vec<usize>, // node indices, not eqn indices
    sons: Vec<Rc<Cluster>>
}

impl Cluster {
    /// Build cluster tree from nodal data
    /// Leaf cardinality is the approximate minimum size of clusters
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
    /// Can be called recursively to process current cluster, split if applicable,
    /// then process those clusters
    fn process_cluster(&mut self, 
                       nodes: &Vec<Node>, 
                       leaf_cardinality: usize, 
                       eqn_map: &HashMap::<usize, usize>) {
        // Largely adapted from https://doi.org/10.1016/S0955-7997(02)00152-2
        // Example 2.1
        self.update_bounds(nodes);
        self.update_diameter();
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
    pub fn is_leaf(&self) -> bool { return self.sons.is_empty()}
    pub fn get_indices(&self) -> &Vec<usize> {return &self.indices_contained;}
    pub fn get_diameter(&self) -> f64 { return self.diameter;}
    pub fn get_sons(&self) -> &Vec<Rc<Cluster>> { return &self.sons;}
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
    /// From example 2.2
    fn update_diameter(&mut self) {
        let diam = &mut self.diameter;
        *diam = 0.0;
        for j in 0..3 {
            *diam += f64::powi(self.u_bound[j] - self.l_bound[j],2);
        }
        *diam = diam.sqrt();
    }
    /// Get distance between two clusters, used for building block tree
    /// From Example 2.2
    pub fn get_distance(c1: &Cluster, c2: &Cluster) -> f64 {
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
    use std::collections::HashMap;
    use na::Vector3;

    use crate::{preprocess::mesh_data::Node, solve::h_matrix::cluster::Cluster};

    #[test]
    fn build_cluster_tree() {
        let mut hmap = HashMap::<usize, usize>::new();
        let mut nodes = Vec::<Node>::new();
        let mut i = 0;
        for j in 0..25 {
            for k in 0..25 {
                nodes.push(Node { 
                    id: i, 
                    coords: Vector3::new(j as f64, k as f64, 0.0),
                    normal: Vector3::from_element(0.0)});
                hmap.insert(i, i);
                i += 1;
            }
        }
        let tree = Cluster::new_from(&nodes, (0..nodes.len()).collect(), 32, &hmap);
        approx::assert_relative_eq!(tree.get_diameter(), 33.94, epsilon = 1e-2);
    }
}