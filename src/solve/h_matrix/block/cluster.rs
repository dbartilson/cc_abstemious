/*!
Cluster tree, representing a hierarchical clustering of physical points

The cluster tree is used to partition the surface in space until
the leaves contain fewer nodes than the 'cardinality' parameter of the tree
*/

use std::{collections::HashMap, rc::Rc};
use crate::preprocess::mesh::CollocationPoint;

/// Cluster Tree, used to partition the surface in space until
/// they contain fewer nodes than the 'cardinality' parameter of the tree
#[derive(Clone)]
pub struct Cluster {
    /// Physical upper bound of the coordinates in the cluster
    u_bound: [f64;3],
    /// Physical lower bound of the coordinates in the cluster
    l_bound: [f64;3],
    /// See update_diameter for definition
    diameter: f64,
    /// node indices, not eqn indices
    indices_contained: Vec<usize>, 
    /// References to 'sons' (i.e., subclusters this was split into)
    sons: Vec<Rc<Cluster>>
}

impl Cluster {
    /// Build cluster tree from nodal data.
    /// Leaf cardinality is the approximate minimum size of clusters
    pub fn new_from(cpts: &Vec<CollocationPoint>, 
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
        cluster.process_cluster(cpts, leaf_cardinality, eqn_map);
        return cluster;
    }
    /// Can be called recursively to process current cluster, split if applicable,
    /// then process those clusters
    fn process_cluster(&mut self, 
                       nodes: &Vec<CollocationPoint>, 
                       leaf_cardinality: usize, 
                       eqn_map: &HashMap::<usize, usize>) {
        // Largely adapted from https://doi.org/10.1016/S0955-7997(02)00152-2
        // Example 2.1
        self.update_bounds(nodes);
        self.update_diameter();
        if self.indices_contained.len() <= leaf_cardinality {
            // Once the cluster is finalized, map from node indices to eqn indices
            self.map_cpts_to_eqns(eqn_map);
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
    /// Check if a cluster is a leaf (has no sons, is not further split)
    pub fn is_leaf(&self) -> bool { return self.sons.is_empty()}
    /// Return reference to indices contained in this cluster (includes sons)
    pub fn get_indices(&self) -> &Vec<usize> {return &self.indices_contained;}
    /// Return diameter of this cluster
    pub fn get_diameter(&self) -> f64 { return self.diameter;}
    /// Return reference to the sons of this cluster
    pub fn get_sons(&self) -> &Vec<Rc<Cluster>> { return &self.sons;}
    /// Update the bounds (maxima/minima in each coordinate direction)
    fn update_bounds(&mut self, cpts: &Vec<CollocationPoint>) {
        let alpha = &mut self.l_bound;
        let beta = &mut self.u_bound;
        for index in &self.indices_contained {
            for j in 0_usize..3 {
                alpha[j] = f64::min(cpts[*index].coords[j], alpha[j]);
                beta[j] = f64::max(cpts[*index].coords[j], beta[j]);
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
    /// Use equation map to go from cpt indices to equations
    fn map_cpts_to_eqns(&mut self, eqn_map: &HashMap::<usize, usize>) {
        for idx in &mut self.indices_contained {
            if let Some(new_idx) = eqn_map.get(idx) {*idx = *new_idx;}
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use na::Vector3;

    use crate::{preprocess::mesh::CollocationPoint, solve::h_matrix::block::cluster::Cluster};

    #[test]
    fn build_cluster_tree() {
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
                    wt: 1.0 } );
                hmap.insert(i, i);
                i += 1;
            }
        }
        let tree = Cluster::new_from(&cpts, (0..cpts.len()).collect(), 32, &hmap);
        approx::assert_relative_eq!(tree.get_diameter(), 33.94, epsilon = 1e-2);
    }
}