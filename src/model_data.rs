pub mod text_reader;

pub mod model_data {

    use std::path::Path;
    use std::io::Write;
    use super::text_reader;
    use na::Vector3;

    pub type Coords = Vector3<f64>;

    pub enum ElementType {
        Null,
        Point,
        Line,
        Tri,
        Quad
    }

    pub struct Element {
        pub id: usize,
        pub body_id: usize,
        pub etype: ElementType,
        pub node_ids: Vec<usize>
    }

    pub enum BodyType {
        Null,
        Point,
        Line,
        Surface,
        Volume
    }
    pub struct Body {
        pub id: usize,
        pub element_ids: Vec<usize>
    }

    pub struct Node {
        pub id: usize,
        pub coords: Coords,
    }

    pub struct MeshData {
        pub nodes: Vec<Node>,
        pub elements: Vec<Element>,
        pub bodies: Vec<Body>,
    }
    impl Default for MeshData {
        fn default() -> MeshData {
            MeshData {
                nodes: Vec::new(),
                elements: Vec::new(),
                bodies: Vec::new()
            }
        }
    }
    impl MeshData {
        pub fn read_from_vtk(&mut self, path: &Path) -> std::io::Result<()> {
            print!(" Reading VTK file '{}' ...", path.display().to_string());
            std::io::stdout().flush().unwrap();
            let mut reader = text_reader::text_reader::BufReader::open(path)?;
            let mut buffer = String::new();

            let nodes = &mut self.nodes;
            let elements = &mut self.elements;
            while let Some(_line) = reader.read_line(&mut buffer) {
                //println!("{}", line?.trim());
                let mut sline = buffer.trim().split_whitespace();
                let first_word = sline.next();
                if first_word == Some("POINTS") {
                    let npts: usize = sline.next().as_ref().unwrap().parse().unwrap();
                    for i in 1..npts+1 {
                        reader.read_line(&mut buffer);
                        sline = buffer.split_whitespace();
                        let mut node_temp: Node = Node{id: i, coords:Coords::from_element(0.0)};
                        for j in 0..3 {
                            node_temp.coords[j] = sline.next().as_ref().unwrap().parse().unwrap();
                        }
                        nodes.push(node_temp);
                    }
                }
                else if first_word == Some("CELLS") {
                    let nelem: usize = sline.next().as_ref().unwrap().parse().unwrap();
                    for i in 1..nelem+1 {
                        reader.read_line(&mut buffer);
                        sline = buffer.split_whitespace();
                        let body_id: usize = sline.next().as_ref().unwrap().parse().unwrap();
                        let mut elem_temp: Element = Element{id: i, body_id: body_id as usize,
                            etype: ElementType::Null, node_ids: Vec::new()};
                        for slinej in sline {
                            elem_temp.node_ids.push(slinej.parse().unwrap());
                        }
                        elements.push(elem_temp)
                    }
                }
                else if first_word == Some("CELL_TYPES") {
                    let nelem: usize = sline.next().as_ref().unwrap().parse().unwrap();
                    for i in 1..nelem+1 {
                        reader.read_line(&mut buffer);
                        sline = buffer.split_whitespace();
                        let etype: usize = sline.next().as_ref().unwrap().parse().unwrap();
                        elements[i-1].etype = match etype {
                            1 => ElementType::Point,
                            3 => ElementType::Line,
                            5 => ElementType::Tri,
                            _ => ElementType::Null
                        }
                    }
                }
            }

            // Now assign elements to bodies
            let bodies = &mut self.bodies;
            for el in &mut *elements {
                let body_id = &el.body_id;
                if (*body_id) > bodies.len() {
                    bodies.push(Body {id: *body_id, element_ids: Vec::new()});
                }
                bodies[(body_id-1) as usize].element_ids.push(el.id);
            }
            println!(" Done!");
            println!(" Read {} bodies, {} elements, {} nodes", bodies.len(), elements.len(), nodes.len());
            Ok(())
        }
    }
    
}