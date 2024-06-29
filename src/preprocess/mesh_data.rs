mod text_reader {
    // module for reading text files line by line to a buffer
    use std::{
        fs::File,
        io::{self, prelude::*},
    };

    pub struct BufReader {
        reader: io::BufReader<File>,
    }

    impl BufReader {
        pub fn open(path: impl AsRef<std::path::Path>) -> io::Result<Self> {
            let file = File::open(path)?;
            let reader = io::BufReader::new(file);
            Ok(Self { reader })
        }

        pub fn read_line<'buf>(&mut self,
            buffer: &'buf mut String,
        ) -> Option<io::Result<&'buf mut String>> {
            buffer.clear();

            self.reader
                .read_line(buffer)
                .map(|u| if u == 0 { None } else { Some(buffer) })
                .transpose()
        }
    }
}

use std::path::Path;
use std::io::Write;
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

//pub enum BodyType {
//    Null,
//    Point,
//    Line,
//    Surface,
//    Volume
//}
pub struct Body {
    pub id: usize,
    pub element_ids: Vec<usize>
}

pub struct Node {
    pub id: usize,
    pub coords: Coords,
}

pub struct Mesh {
    pub nodes: Vec<Node>,
    pub elements: Vec<Element>,
    pub bodies: Vec<Body>,
}
impl Default for Mesh {
    fn default() -> Mesh {
        Mesh {
            nodes: Vec::new(),
            elements: Vec::new(),
            bodies: Vec::new()
        }
    }
}
impl Mesh {
    pub fn read_from_vtk(&mut self, path: &Path) -> std::io::Result<()> {
        // read mesh from VTK (ASCII) format
        info!(" Reading VTK (ASCII) file '{}' ...", path.display().to_string());
        std::io::stdout().flush().unwrap();
        let mut reader = text_reader::BufReader::open(path)?;
        let mut buffer = String::new();

        let nodes = &mut self.nodes;
        let elements = &mut self.elements;
        while let Some(_line) = reader.read_line(&mut buffer) {
            let mut sline = buffer.trim().split_whitespace();
            let first_word = sline.next();
            // Read POINTS data block
            if first_word == Some("POINTS") {
                let npts: usize = sline.next().as_ref().unwrap().parse().unwrap();
                for i in 0..npts {
                    reader.read_line(&mut buffer);
                    sline = buffer.split_whitespace();
                    let mut node_temp: Node = Node{id: i, coords:Coords::from_element(0.0)};
                    for j in 0..3 {
                        node_temp.coords[j] = sline.next().as_ref().unwrap().parse().unwrap();
                    }
                    nodes.push(node_temp);
                }
            }
            // Read CELL data block
            else if first_word == Some("CELLS") {
                let nelem: usize = sline.next().as_ref().unwrap().parse().unwrap();
                for i in 0..nelem {
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
            // Read CELL_TYPES data block
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
                        9 => ElementType::Quad,
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
                bodies.push(Body {id: *body_id-1, element_ids: Vec::new()});
            }
            bodies[body_id-1].element_ids.push(el.id);
        }
        info!(" Read {} bodies, {} elements, {} nodes", bodies.len(), elements.len(), nodes.len());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn vtk_reader() {
        // test mesh reader (VTK ASCII) capability
        use std::path::Path;

        let mut mesh: crate::preprocess::mesh_data::Mesh = Default::default();
        let _result = mesh.read_from_vtk(Path::new("./src/tests/sphere.vtk"));

        assert_eq!(mesh.elements.len(), 336);
    }
}