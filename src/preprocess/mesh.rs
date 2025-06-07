/*!
Mesh file processing, currently only VTK ASCII
*/

mod text_reader {
    /// For reading text files line by line to a buffer
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

        pub fn get_reader(&mut self) -> &mut std::io::BufReader<File> {
            &mut self.reader
        }
    }
}

use std::{io::Read, path::{self, Path}};
use na::Vector3;

pub type Coords = Vector3<f64>;

/// VTK element types
#[derive(Clone, PartialEq)]
pub enum ElementType {
    Null,
    Point,
    Line,
    Tri,
    Quad
}
/// Element struct
pub struct Element {
    /// Element ID #
    pub id: usize,
    /// Body ID associated with element
    pub body_id: usize,
    /// Type
    pub etype: ElementType,
    /// Global node IDs
    pub node_ids: Vec<usize>,
    /// Integration data
    pub intpts: Vec<CollocationPoint>
}

//pub enum BodyType {
//    Null,
//    Point,
//    Line,
//    Surface,
//    Volume
//}
/// Body, having an ID and vector ofelement IDs
pub struct Body {
    pub id: usize,
    pub element_ids: Vec<usize>
}
/// Node, having ID and coordinates
#[derive(Clone)]
pub struct Node {
    pub id: usize,
    pub coords: Coords
}

/// Collocation points, used for integration
pub struct CollocationPoint {
    /// Collocation point ID #
    pub id: usize,
    /// Global coordinates
    pub coords: Coords,
    /// Normal vector 
    pub normal: Vector3<f64>,
    /// Area associated with this (detj)
    pub area: f64,
    /// Gauss integration weight
    pub wt: f64 
}

/// Mesh, from VTK
#[derive(Default)]
pub struct Mesh {
    /// All nodes
    pub nodes: Vec<Node>,
    /// All elements
    pub elements: Vec<Element>,
    /// All bodies 
    pub bodies: Vec<Body>,
    /// Collocation points (added in preprocessing)
    pub cpts: Vec<CollocationPoint>
}

#[allow(clippy::upper_case_acronyms)]
#[derive(PartialEq)]
enum VTKFormat {
    Binary,
    ASCII, 
    Null
}

impl Mesh {
    /// Set up mesh from VTK ASCII
    pub fn read_from_vtk(&mut self, path: &Path) -> std::io::Result<()> {
        if !path.exists() { 
            error!("Mesh file not found at specified path: {}", path::absolute(path).unwrap().display())
        }
        // read mesh from VTK format
        info!(" Reading VTK file '{}' ...", path.display());
        let mut reader = text_reader::BufReader::open(path)?;
        let mut str_buf = String::new();

        // Header info
        reader.read_line(&mut str_buf); 
        if !( str_buf.contains("# vtk DataFile Version 1.0") || str_buf.contains("# vtk DataFile Version 2.0") || str_buf.contains("# vtk DataFile Version 3.0")) {
            error!("Invalid VTK mesh file version")
        }
        reader.read_line(&mut str_buf); // Discard second line, header
        reader.read_line(&mut str_buf);
        
        let mut format = VTKFormat::Null;
        if str_buf.contains("ASCII") {
            format = VTKFormat::ASCII
        }
        else if str_buf.contains("BINARY") {
            format = VTKFormat::Binary
        }
        else {
            error!("Only BINARY or ASCII allowed for VTK mesh file (line 3)")
        }
        reader.read_line(&mut str_buf);
        if !str_buf.contains("DATASET UNSTRUCTURED_GRID") {
            error!("Only UNSTRUCTURED_GRID dataset allowed for VTK mesh file")
        }

        // Loop over lines
        let nodes = &mut self.nodes;
        let elements = &mut self.elements;
        while let Some(_line) = reader.read_line(&mut str_buf) {
            let mut sline = str_buf.split_whitespace();
            let first_word = sline.next();
            
            // Read POINTS data block
            if first_word == Some("POINTS") {
                let npts: usize = sline.next().as_ref().unwrap().parse().unwrap();
                if format == VTKFormat::ASCII {
                    for i in 0..npts {
                        reader.read_line(&mut str_buf);
                        sline = str_buf.split_whitespace();
                        let mut node_temp: Node = Node{
                            id: i, 
                            coords:Coords::from_element(0.0)
                        };
                        for j in 0..3 {
                            node_temp.coords[j] = sline.next().as_ref().unwrap().parse().unwrap();
                        }
                        nodes.push(node_temp);
                    }
                }
                else if format == VTKFormat::Binary {
                    let mut bin_buf = [0u8; std::mem::size_of::<f64>()];
                    for i in 0..npts {
                        let mut node_temp: Node = Node{
                            id: i, 
                            coords:Coords::from_element(0.0)
                        };
                        for j in 0..3 {
                            reader.get_reader().read_exact(&mut bin_buf).expect("Unable to process binary VTK: POINTS");
                            node_temp.coords[j] = f64::from_be_bytes(bin_buf);
                        }
                        nodes.push(node_temp);
                    }
                }
            }
            // Read CELL data block
            else if first_word == Some("CELLS") {
                let nelem: usize = sline.next().as_ref().unwrap().parse().unwrap();
                if format == VTKFormat::ASCII {
                    for i in 0..nelem {
                        reader.read_line(&mut str_buf);
                        sline = str_buf.split_whitespace();
                        let npts: usize = sline.next().as_ref().unwrap().parse().unwrap();
                        let mut elem_temp: Element = 
                            Element{id: i, 
                                    body_id: 0,
                                    etype: ElementType::Null, 
                                    node_ids: Vec::new(),
                                    intpts: Vec::new()};
                        for slinej in sline {
                            elem_temp.node_ids.push(slinej.parse().unwrap());
                        }
                        if elem_temp.node_ids.len() != npts {
                            error!("Incorrect number of nodes read for element: {}", i)
                        }
                        elements.push(elem_temp)
                    }
                }
                else if format == VTKFormat::Binary {
                    let mut bin_buf = [0u8; std::mem::size_of::<u32>()];
                    for i in 0..nelem {
                        reader.get_reader().read_exact(&mut bin_buf).expect("Unable to process binary VTK: CELLS NPTS");
                        let npts: u32 = u32::from_be_bytes(bin_buf);
                        let mut elem_temp: Element = 
                            Element{id: i, 
                                    body_id: 0,
                                    etype: ElementType::Null, 
                                    node_ids: Vec::new(),
                                    intpts: Vec::new()};
                        for _j in 0..npts {
                            reader.get_reader().read_exact(&mut bin_buf).expect("Unable to process binary VTK: CELLS NODE_IDS");
                            let uresult = u32::from_be_bytes(bin_buf);
                            elem_temp.node_ids.push(usize::try_from(uresult).expect("Convert usize from u32 failed"))
                        }
                        elements.push(elem_temp)
                    }
                    
                }
            }
            // Read CELL_TYPES data block
            else if first_word == Some("CELL_TYPES") {
                let nelem: usize = sline.next().as_ref().unwrap().parse().unwrap();
                if format == VTKFormat::ASCII {
                    for elem in elements.iter_mut().take(nelem) {
                        reader.read_line(&mut str_buf);
                        sline = str_buf.split_whitespace();
                        let etype: usize = sline.next().as_ref().unwrap().parse().unwrap();
                        elem.etype = match etype {
                            1 => ElementType::Point,
                            3 => ElementType::Line,
                            5 => ElementType::Tri,
                            9 => ElementType::Quad,
                            _ => ElementType::Null
                        }
                    }
                } 
                else if format == VTKFormat::Binary {
                    let mut bin_buf = [0u8; std::mem::size_of::<u32>()];
                    for elem in elements.iter_mut().take(nelem) {
                        reader.get_reader().read_exact(&mut bin_buf).expect("Unable to process binary VTK: CELL_TYPES");
                        let uresult = usize::try_from(u32::from_be_bytes(bin_buf)).expect("Covert usize from u32 failed");
                        elem.etype = match uresult {
                            1 => ElementType::Point,
                            3 => ElementType::Line,
                            5 => ElementType::Tri,
                            9 => ElementType::Quad,
                            _ => ElementType::Null
                        }
                    }
                }
            }
        }

        // Now assign elements to bodies
        let bodies = &mut self.bodies;
        // Create 3 bodies: one each for 0d, 1d, and 2d elements; one for all others
        for i in 0..4 {
            bodies.push(Body {id: i, element_ids: Vec::new()});
        }
        for el in &mut *elements {
            let body_id = match el.etype {
                ElementType::Point => 0,
                ElementType::Line => 1,
                ElementType::Tri | ElementType::Quad => 2,
                ElementType::Null => 3
            };
            el.body_id = body_id;
            bodies[body_id].element_ids.push(el.id);
        }
        info!(" Read {} bodies, {} elements, {} nodes", bodies.len(), elements.len(), nodes.len());
        if !bodies[3].element_ids.is_empty() {
            warn!(" Unable to process {} elements", bodies[3].element_ids.len())
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::preprocess::mesh::ElementType;

    #[test]
    fn vtk_reader_ascii() {
        // test mesh reader (VTK ASCII) capability
        use std::path::Path;

        let mut mesh: crate::preprocess::mesh::Mesh = Default::default();
        let _result = mesh.read_from_vtk(Path::new("./src/tests/sphere_ascii.vtk"));

        assert_eq!(mesh.bodies.len(), 4);
        assert_eq!(mesh.nodes.len(), 164);
        assert_eq!(mesh.elements.len(), 336);
        assert_eq!(mesh.nodes[0].coords[2], 0.5);
        assert_eq!(mesh.elements[12].node_ids[1], 122);
        assert!(mesh.elements[12].etype == ElementType::Tri);
        assert_eq!(mesh.elements[12].body_id, 2);
        assert_eq!(mesh.bodies[2].element_ids[0], 12)
    }

    #[test]
    fn vtk_reader_binary() {
        // test mesh reader (VTK Binary) capability
        use std::path::Path;

        let mut mesh: crate::preprocess::mesh::Mesh = Default::default();
        let _result = mesh.read_from_vtk(Path::new("./src/tests/sphere_1.vtk"));

        assert_eq!(mesh.bodies.len(), 4);
        assert_eq!(mesh.nodes.len(), 164);
        assert_eq!(mesh.elements.len(), 336);
        assert_eq!(mesh.nodes[0].coords[2], 0.5);
        assert_eq!(mesh.elements[12].node_ids[1], 122);
        assert!(mesh.elements[12].etype == ElementType::Tri);
        assert_eq!(mesh.elements[12].body_id, 2);
        assert_eq!(mesh.bodies[2].element_ids[0], 12)
    }
}