# cc_abstemious <=> BEM-ACOUSTICS

`cc_abstemious` is a numerical acoustics simulation software based on the boundary element method. The main features are:

* Direct boundary integral formulation
* Galerkin representation (i.e., integration over surface elements), avoiding singular integration and allowing easier interpolation of surface results
* Support for VTK mesh files
* Multiple solver options:
   * Dense matrix, direct (LU)
   * Dense matrix, iteratative (GMRES)
   * Hierarchical matrix decomposition with iterative (GMRES)
* Parallel processing of surface influence matrix (dense) and hierarchical decompositions

**License:** MIT

**Author:** Daniel T. Bartilson

## Usage

`cc_abstemious` supports `.json` input files which specify the following analysis fields:

* Mesh file in `.vtk` format (ASCII) [^1][^2][^3]
* Body index in the mesh file [^3]
* Vector of analysis (drive) frequencies [^4]
* Acoustic fluid sound speed and mass density
* Problem type (interior or exterior domain)
* Solver specification (direct/iterative, tolerance and max iterations for iterative solver)
* Incident acoustic fields (plane waves and spherical waves)
* Surface boundary conditions (velocity potential, normal velocity, or impedance) [^5]
* Output information, including file name, pressure/velocity potential, scattered/total field, and field point locations [^6]

## Limitations & Future Work

[^1]: Only linear, 3-noded triangles and linear, 4-noded quadrilateral surface elements are supported
[^2]: Only 3D bodies and analysis are supported
[^3]: The user must specify the acoustic body index from the VTK mesh file
[^4]: Only one body supported per analysis 
[^5]: Only one surface boundary condition may be utilized, i.e., the pressure, normal velocity, or impedance B.C. must be applied to the whole surface
[^6]: There is no explicit check of whether a field point is interior to the body for interior problems, and similar for exterior problems



