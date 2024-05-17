pub mod preprocess {
    use std::collections::HashMap;
    use crate::model_data::model_data::MeshData;

    pub fn get_eqn_map(meshdata: &MeshData, body_id: usize) -> (usize, HashMap::<usize, usize>) {
        let mut eqn_map = HashMap::<usize, usize>::new();

        let nnode = meshdata.nodes.len();
        let ibody = &meshdata.bodies[body_id-1];
        let mut eqn_vector = vec![false; nnode];
        // Mark each possible node (eqn) as used or not
        for element_id in &ibody.element_ids {
            let element = &meshdata.elements[*element_id-1];
            for node_id in &element.node_ids {
                eqn_vector[*node_id] = true;
            }
        }
        // Scroll through all used eqns in order and put in map
        let mut eqn_number: usize = 0;
        for i in 0..eqn_vector.len() {
            if eqn_vector[i] {
                eqn_map.insert(i, eqn_number);
                eqn_number += 1;
            }
        }
        return (eqn_number, eqn_map);
    }
}