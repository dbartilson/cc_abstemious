# CC-ABSTEMIOUS <=> BEM-ACOUSTICS

`CC-ABSTEMIOUS` is a numerical acoustics simulation software based on the boundary element method. The main features are:

* Direct boundary integral formulation
* Formulated using velocity potential
* Support for interior and exterior problems
* Collocation points are located at element nodes, avoiding singular integration and allowing easier interpolation
* Support for VTK mesh files

**License:** MIT

**Author:** Daniel T. Bartilson

## Usage

`CC-ABSTEMIOUS` supports `.json` input files which specify the following analysis fields:

* Mesh file in `.vtk` format [^1][^2][^3]
* Body index in the mesh file [^3]
* Analysis (drive) frequency [^4]
* Acoustic fluid sound speed and mass density
* Problem type (interior or exterior domain)
* Field point locations [^5]
* Incident acoustic fields (plane waves and spherical waves)
* Surface boundary conditions (velocity potential, normal velocity, or impedance) [^6]

## Limitations & Future Work

[^1] Only linear, 3-noded triangles are supported
[^2] Only 3D bodies and analysis is supported
[^3] The user must specify the acoustic body index from the VTK mesh file
[^4] Only one body and one frequency is supported per analysis 
[^5] There is no explicit check of whether a field point is interior to the body for interior problems, and similar for exterior problems
[^6] Only one surface boundary condition may be utilized, i.e., the pressure, normal velocity, or impedance B.C. must be applied to the whole surface



