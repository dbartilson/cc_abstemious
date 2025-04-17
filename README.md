# cc_abstemious <=> BEM-ACOUSTICS

`cc_abstemious` is a numerical acoustics simulation software based on the boundary element method. The main features are:

* Direct boundary integral formulation
* Collocation representation with analytical integration for singular integrals
* Support for VTK mesh files
* Multiple solver options:
   * Dense matrix, direct (LU)
   * Dense matrix, iteratative (GMRES)
   * Hierarchical matrix decomposition (ACA) with iterative solver
* Parallel processing of surface influence matrix (dense) and hierarchical decompositions

**License:** MIT

**Author:** Daniel T. Bartilson

## Build and test

``cargo`` is used for the build and test phases. If you install Rust via `rustup`, you should have `cargo` as well.

### Automated tests

``cargo test`` will build ``cc_abstemious`` in debug mode, then run the unit and integration tests.

``cargo test --release`` will do the same, but in release mode.

### Build

``cargo build --release`` will build ``cc_abstemious`` in release mode.

## Usage

`cc_abstemious` supports JSON input files. This describes an analysis that is run in the typical command line way: 
```
cc_abstemious(.exe) my_input.json
```

### Input file format

An example input file is available in the `src/tests` directory. 

The input file includes:
* The path to the mesh file, in `.vtk` format (ASCII) [^1][^2][^3]
* Body index in the mesh file [^3]
* Vector of analysis (drive) frequencies [^4]
* Sound speed and mass density of the acoustic fluid
* Problem type (interior or exterior domain)
* Method type (classical or Burton-Miller)
* Solver specification (direct/iterative, tolerance and max iterations for iterative solvers)
* Incident acoustic fields (plane waves and spherical waves)
* Surface boundary conditions (velocity potential, normal velocity, or impedance) [^5]
* Output information, including file name, pressure/velocity potential, scattered/total field, and field point locations [^6]

A file (`input_schema.json`) describing the input file JSON schema is automatically generated during `cargo test` (using the `schemars` crate).

## Limitations & Future Work

[^1]: Only linear, 3-noded triangles and linear, 4-noded quadrilateral surface elements are supported
[^2]: Only 3D bodies and analysis are supported
[^3]: The user must specify the acoustic body index from the VTK mesh file
[^4]: Only one body supported per analysis 
[^5]: Only one surface boundary condition may be utilized, i.e., the pressure, normal velocity, or impedance B.C. must be applied to the whole surface
[^6]: There is no explicit check of whether a field point is interior to the body for interior problems, and similar for exterior problems



