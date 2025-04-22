# Theory

## The acoustic wave and Helmholtz equations

The linear wave equation is satisfied over a field, at a given location $`\mathbf{x}`$ and at a time $`t`$ [[1]](#1)[[2]](#2)

```math
\nabla ^2 \psi(\mathbf{x},t) = \frac{1}{c^2} \frac{\partial^2}{\partial t^2} \psi(\mathbf{x},t)
```

where $`\psi(\mathbf{x},t)`$ is the velocity potential field and $`c`$ is the propagation velocity. The vector velocity is given $`\mathbf{v}(\mathbf{x},t) = \nabla \psi(\mathbf{x},t)`$ and the acoustic pressure is $`p = -\rho \frac{\partial}{\partial t} \psi(\mathbf{x},t)`$ with $`\rho`$ as the fluid mass density.

Using a harmonic, steady-state assumption, the velocity potential field can be rewritten as [[3]](#3)

```math
\psi(\mathbf{x},t) = \text{Re}[\phi(\mathbf{x}) e^{-i\omega t}]
```

given that the angular drive frequency is $`\omega = 2 \pi f`$, where $`f`$ is the frequency in Hertz. 

The harmonic form of the Helmholtz equation for the fluid is then [[1]](#1)[[2]](#2)

```math
\nabla ^2 \phi + k^2 \phi= \phi_{I}
```

where $`k=\omega / c`$ is the wavenumber. The directional derivative, along a specified normal direction, will also be important. This is expressed as [[2]](#2)

```math
v_n(\mathbf{x}) = \nabla \phi (\mathbf{x}) \cdot \mathbf{e}_n(\mathbf{x}) = \frac{\partial \phi(\mathbf{x})}{\partial n(\mathbf{x})}
```

### Integral solution

The integral solution to the Helmholtz equation is written as [[0]](#0)[[4]](#4):

```math
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) - v_n(\mathbf{y}) g(\mathbf{x}, \mathbf{y}) d\mathbf{y} =
\left\lbrace 
\begin{array}{l l}
\phi(\mathbf{x})-\phi_I(\mathbf{x}) & \text{at field point, ${x}'$}\\
\phi(\mathbf{x})/2-\phi_I(\mathbf{x}) & \text{on surface for exterior problem}  \\
-\phi_I(\mathbf{x}) & \text{on surface for interior problem}
\end{array}
\right.
```

where $`\phi_I`$ is the incident wave (free-field) velocity potential. $`g`$ and its normal derivative are given as [[0]](#0)[[5]](#5):

```math
g(\mathbf{x}, \mathbf{y}) = \frac{e^{ikr}}{4 \pi r} \\
h(\mathbf{x}, \mathbf{y}) = \frac{\partial g(\mathbf{x}, \mathbf{y})}{\partial n(\mathbf{y})} = \left(ik - \frac{1}{r}\right) g(\mathbf{x}, \mathbf{y}) (-\mathbf{e}_r \cdot \mathbf{e}_n(\mathbf{y}))
```

with $`\mathbf{r} = \mathbf{x} - \mathbf{y}`$ as the vector pointing from $`\mathbf{y}`$ to $`\mathbf{x}`$ and $`r = \| \mathbf{r} \|`$ as the magnitude (distance) and $`\mathbf{e}_r = \mathbf{r}/r`$ as the unit vector form.

# Implementation

## Boundary conditions and source (incident field) terms

Three classes of boundary condition may be used on the surface, and only one type may be used on any partition of the surface [[0]](#0):

1. **Velocity potential or pressure B.C.**: $`\mathbf{\phi} = \bar{\mathbf{\phi}}`$. This can represent a known non-zero surface pressure, or if set to zero, it represents a sound-soft (infinitely absorptive) boundary. ( $`\beta = 0, \gamma/\alpha = \bar{\mathbf{\phi}}`$ )
2. **Surface normal velocity B.C.**: $`\mathbf{v}_n = \bar{\mathbf{v}}_n`$. This can represent a known non-zero surface motion, or if set to zero, it represents a sound-hard (reflective) boundary. ( $`\alpha = 0, \gamma/\beta = \bar{\mathbf{v}}_n`$ )
3. **Impedance B.C.**: A known, linear relation between the surface normal velocity and the fluid pressure: $`(i \omega \rho) \phi = p = Z v_n`$ where $`p`$ is the pressure and $`Z`$ is the impedance [[0]](#0)[[3]](#3). ( $`\gamma = 0, \alpha / \beta = -Z / (i \omega \rho)`$ )

These may be stated generally as variants of a Robin condition condition $`\alpha \phi + \beta v_n = \gamma`$, as shown above.

The incident sound sources (incident acoustic fields) can be defined similarly:

1. **Plane wave**: Has a defined strength (amplitude: $`A`$) and propagates along a unit vector $`\mathbf{e}_d`$. The phase is assumed to be zero at the origin, i.e. $`\mathbf{x} \cdot \mathbf{e}_d = 0`$ [[2]](#2): 

```math
\phi_I(\mathbf{x}) = A e^{ik(\mathbf{x} \cdot \mathbf{e}_d)}
```

```math
v_{nI}(\mathbf{x}) = \frac{\partial \phi_I(\mathbf{x})}{\partial n(\mathbf{x})} = ik \left(\frac{\partial \mathbf{x}}{\partial n(\mathbf{x})}\cdot \mathbf{e}_d+\mathbf{x} \cdot \frac{\partial \mathbf{e}_d}{\partial n(\mathbf{x})} \right) \phi_I(\mathbf{x}) = ik (\mathbf{e}_n(\mathbf{x}) \cdot \mathbf{e}_d) \phi_I(\mathbf{x})
```

2. **Point source**: Also known as a volume source or spherical source. Has a defined amplitude $`A`$ and position vector $`\mathbf{x}_I`$ [[2]](#2):

```math
\phi_I(\mathbf{x}) = \frac{A}{4 \pi r}e^{ikr}; \qquad r = \mathbf{x} - \mathbf{x}_I
```

```math
v_{nI}(\mathbf{x}) = \frac{\partial \phi_I(\mathbf{x})}{\partial n(\mathbf{x})} = \left(ik - \frac{1}{r}\right) (\mathbf{e}_r \cdot \mathbf{e}_n(\mathbf{x})) \phi_I(\mathbf{x})
```

## Boundary elements and discretization

The collocation method is used, and thus, the primary variables are assumed to be constant within the element: 

```math
\phi(\mathbf{y}_e) \approx \phi_e \qquad v_n(\mathbf{y}_e) \approx v_{ne}
```

where $`\mathbf{y}_e`$ is the ordinate within element $`e`$, and $`\phi_e, v_{ne}`$ are the assumed-constant primary variables for that element, usually assigned to a collocation node at the center of the element.

The integral form can be discretized over the surface this way as

```math
\mathbf{H} \mathbf{\phi} - \mathbf{G} \mathbf{v}_n = 
\left\lbrace 
\begin{array}{l l}
\mathbf{\phi}/2-\mathbf{\phi}_I & \text{on surface for exterior problem}  \\
-\mathbf{\phi}_I & \text{on surface for interior problem}
\end{array}
\right.
```

Examining the $`\mathbf{H}\mathbf{\phi}`$ term, each row $`i`$ of $`\mathbf{H}`$ corresponds to the contribution (influence) of the velocity potentials at all surface points to the velocity potential at node $`i`$, which can be split into an integral over the elements [[2]](#2): 

```math
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) d\mathbf{y} \approx
\sum_e \int_{S_e} \phi (\mathbf{y_e}) h(\mathbf{x}, \mathbf{y_e}) d\mathbf{y}_e = 
\left(\sum _e \int_{S_e} h(\mathbf{x}, \mathbf{y}_e)d\mathbf{y}_e \right)\mathbf{\phi}
```

Note that the collocation values of velocity potential are assembled into the vector $`\mathbf{\phi}`$ and the in-element integral can be calculated using quadrature rules.

The velocity potentials must first be solved on the surface, giving the matrix equation (in an exterior problem) [[2]](#2):

```math
\left[\mathbf{H} - \frac{1}{2} \mathbf{I} \right] \mathbf{\phi} = \mathbf{G} \mathbf{v}_n - \mathbf{\phi}_I
```

For the three possible boundary conditions, this is solved as:

```math
\begin{array}{l l}
\mathbf{G} \mathbf{v}_n = \mathbf{\phi}_I -\left[\mathbf{H} - \frac{1}{2} \mathbf{I} \right] \bar{\mathbf{\phi}}  & \text{pressure B.C.} \\
\left[\mathbf{H} - \frac{1}{2} \mathbf{I} \right] \mathbf{\phi} = \mathbf{G} \bar{\mathbf{v}}_n - \mathbf{\phi}_I & \text{normal velocity B.C.} \\
\left[\mathbf{H} - \frac{1}{2} \mathbf{I} - \frac{i \omega \rho}{Z} \mathbf{G}\right] \mathbf{\phi} = -\mathbf{\phi}_I  & \text{impedance B.C.}
\end{array}
```

## Burton-Miller Method

In its classical form, the boundary integral formulation suffers from problems where interior mode frequencies can be present in exterior analyses. The Burton-Miller method is one technique for avoiding this issue. First, consider the 'hypersingular' formulation of the boundary integral equation which is equal to the derivative of the classical formulation with respect to the normal vector at $`\mathbf{x}`$ [[0]](#0)[[2]](#2):

```math
\int_S \phi (\mathbf{y}) dh(\mathbf{x}, \mathbf{y}) - v_n(\mathbf{y}) dg(\mathbf{x}, \mathbf{y}) d\mathbf{y} = \frac{v_{n}(\mathbf{x})}{2} - v_{nI}(\mathbf{x})
```

where [[0]](#0)[[5]](#5)
```math
dg(\mathbf{x}, \mathbf{y}) = \frac{\partial g(\mathbf{x}, \mathbf{y})}{\partial n(\mathbf{x})} = \left(ik - \frac{1}{r}\right) g(\mathbf{x}, \mathbf{y}) (\mathbf{e}_r \cdot \mathbf{e}_n(\mathbf{x})) 
```
```math
\begin{align*}
dh(\mathbf{x}, \mathbf{y}) = \frac{\partial h(\mathbf{x}, \mathbf{y}) }{\partial n(\mathbf{x})} = \frac{g(\mathbf{x}, \mathbf{y})}{r} &\left[  -\left(ik - \frac{1}{r}\right)\left( \mathbf{e}_n(\mathbf{x}) \cdot \mathbf{e}_n(\mathbf{y}) \right) \right. + \\
&\left.  \left( k^2 r + 3 \left(ik - \frac{1}{r}\right) \right) \left( \mathbf{e}_r \cdot \mathbf{e}_n(\mathbf{x}) \right) \left( \mathbf{e}_r \cdot \mathbf{e}_n(\mathbf{y}) \right) \right]
\end{align*}
```

Discretized over the elements and assembled into matrix form, similar to the classical formulation, yields:
```math
[d\mathbf{H}] \mathbf{\phi} = \left[[d\mathbf{G}] + \frac{1}{2} \mathbf{I} \right]\mathbf{v}_n - \mathbf{v}_{nI}
```
Summing together the classical and hypersingular equations using a coupling factor $`\beta`$ yields the following form [[0]](#0)[[2]](#2):

```math
\left[\mathbf{H} + \beta [d\mathbf{H}] - \frac{1}{2} \mathbf{I}\right] \mathbf{\phi} = \left[\mathbf{G} + \beta[d\mathbf{G}] + \frac{\beta}{2} \mathbf{I} \right]\mathbf{v}_n - \mathbf{\phi}_I - \beta\mathbf{v}_{nI}
```
where $`\beta`$ is traditionally set to $`i/k`$ [[0]](#0)[[6]](#6).

## Field solution

Once the velocity potential and normal velocity fields are known on the surface, the velocity potential for an arbitrary point $`\mathbf{x}`$ in the interior or exterior field can be found from the general solution via: 

```math
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) - v_n(\mathbf{y}) g(\mathbf{x}, \mathbf{y}) d\mathbf{y} =
\phi(\mathbf{x})-\phi_I(\mathbf{x}) \\
\mathbf{\phi}_{fp} = \mathbf{M} \mathbf{\phi} - \mathbf{L} \mathbf{v}_n + \mathbf{\phi}_I
```

where $`\mathbf{M}`$ and $`\mathbf{G}`$ are analagously constructed to $`\mathbf{H}`$ and $`\mathbf{G}`$, but are (generally) rectangular matrices of dimension $`n_{fp} \times n_s`$ where $`n_{fp}`$ is the number of field points and $`n_s`$ is the number of surface points. $`\mathbf{M}`$ and $`\mathbf{G}`$ represent the influence of the surface fields on the velocity potential at each field point. Note that $`\mathbf{\phi}_I`$ contains the vector of $`n_{fp}`$ incident wave velocity potentials at each field point.

## Numerical integration and singular integrals

For generating the influence matrices, both for the surface and field problems, requires integrating over surface elements, with a common form. For example, the $`(i,j)`$ component of the $`\mathbf{G}`$ matrix corresponds to the influence of collocation point (equivalently, element) $`j`$ on point $`i`$

```math
G_{ij} = \int_{S_j} g(\mathbf{x}_i, \mathbf{y}_j)d\mathbf{y}_j
```
with $`g(\mathbf{x}, \mathbf{y}) = \frac{e^{ikr}}{4 \pi r}`$ and $`r = \|\mathbf{x} - \mathbf{y}\|`$. This is straight-forward to compute numerically using Gaussian quadrature:
```math
G_{ij} \approx \sum_k g(\mathbf{x}_i, \mathbf{y}(\mathbf{\xi}_k))| \frac{\partial \mathbf{y}}{\partial \mathbf{\xi}} |_k w_k
```
where $`\xi_k`$ is the $`k^{th}`$ Gauss point, expressed in natural coordinate space $`\xi`$ with weight $`w_k`$. $`| \frac{\partial \mathbf{y}}{\partial \mathbf{\xi}} |`$ is the Jacobian of the natural/physical coordinate transformation, analogous to the area associated with the Gauss point.

One can note that $`r`$ goes to 0 as $`\mathbf{y}`$ approaches $`\mathbf{x}_i`$. This is problematic and leads to singular integrals. There are a variety of strategies to ameliorate this problem, but the chosen approach is based on Chertock's analytical approximations [[5]](#5). Thus, when $`i=j`$ the following approximations are used:

```math
G_{jj} \approx \frac{1}{2 \pi b_j} \\
(dH)_{jj} \approx \frac{-1}{2 \pi b_j^3}
```
where $`\pi b_j^2`$ is the area of the element, i.e., $`b_j = \sqrt \frac{S_j}{\pi}`$.

For the other matrices, the approximant is inversely proportional to the mean radius of curvature of element $`j`$, $`c_j`$

```math
H_{jj} \approx (dG)_{jj} \approx \frac{ikb_j - 1}{4 \pi b_j c_j}
```
For flat elements, the average curvature is infinite, and thus the approximant terms are zero [[5]](#5).

## Radiated and incident power

The temporal average power for a point force aligned with a normal velocity is [[3, eq. 2.16]](#3)
```math
W = \frac{1}{2} \text{Re}(f^*v)
```
where $`f^*`$ is the conjugate of the complex force amplitude, and $`v`$ is the (aligned) complex velocity amplitude.

Similarly, the acoustic intensity is [[4, eq. 28]](#4):
```math
I = \frac{1}{2} \text{Re}(p^*v_n)
```
Thus, the total power over the surface is
```math
W = \frac{1}{2} \int _S \text{Re}(p^*v_n)dS
```

# Numerical solution

Naïve solution to the surface system of equations has the form of an $`\mathbf{A} \mathbf{x} = \mathbf{b}`$ linear system of equations, with the note that $`\mathbf{A}`$ is, in general:
- non-symmetric
- complex-valued
- dense

Three solvers are available for this problem:

## LU (Direct, Dense)

The LU solver is practical for small problems and can be more robust than iterative methods. Note that the $`\mathbf{A}`$ matrix storage is $`n^2`$ and the computational complexity of the LU factorization phase is $`\sim n^3`$ and the solution phase is $`\sim n^2`$, where $`n`$ is the number of equations.

## GMRES(k) (Iterative, Dense)

The restarted GMRES iterative iterative solver [[wiki]](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method) removes some of the problems inherent to the LU solver. The matrix is no longer factored, but repeated matrix-vector multiplies are required. The computational complexity is reduced to an approximate $`\sim n^2`$ (for the matrix-vector multiplies), while the memory demand remains the same since the  dense $`\mathbf{A}`$ matrix is still required.

## Hierarchical matrix decomposition (Iterative, Decomposed)

For practical problems, holding the full $`\mathbf{A}`$ matrix becomes impossible. Hierarchical matrix decomposition, often referred to as H-matrix decomposition (but this will be avoided due to confusion with the $`\mathbf{H}`$ matrix defined previously) reduces the memory requirements and compute complexity by a significant margin. 

The theory and implementation is complex [[7-9]](#7), but the main idea is to compress the off-diagonal 'blocks' of the matrix, which represent more distant points, using an SVD-like compression method. The blocks along or close to the diagonal, representing spatially-close pairs of points, are not good candidates for approximation, since they contain strongly-interacting sets. However, much of the matrix can be compressed in this way, leading to large improvements for computation.

For a broad overview, the surface is first processed into a hierarchical clustering of points based on coordinates. This clustering tree then forms the basis for forming another tree: this one representing the pair-wise interactions of clusters. This directly translates to the row and columns of the $`\mathbf{A}`$ matrix, and determining whether the points represented by the rows of the block are sufficiently distant from the points represented by the columns of the block. If they are, the block is deemed 'admissible', in which case it is decomposed by the Adaptive Cross Approximation (ACA) algorithm, which is similar to SVD. If not, the block is inadmissible and it is computed in full, dense format. 

The block tree is then flattened into a list of blocks for easier parallel computing. The matrix-vector multiply stage, essential to the iterative solution, is then computed via a sum of chunked multiplies against the blocks.

## References

<a id="0">[0]</a> Liu, Y. (2019). On the BEM for acoustic wave problems. *Engineering Analysis with Boundary Elements*, 107, 53-62.

<a id="1">[1]</a> Marburg, S., & Nolte, B. (2008). *Computational acoustics of noise propagation in fluids: finite and boundary element methods* (Vol. 578). Berlin: Springer.

<a id="2">[2]</a> Kreuzer, W., Pollack, K., Brinkmann, F., & Majdak, P. (2024). NumCalc: An open-source BEM code for solving acoustic scattering problems. *Engineering Analysis with Boundary Elements*, 161, 157-178.

<a id="3">[3]</a> Cremer, L., & Heckl, M. (2005). *Structure-borne sound: structural vibrations and sound radiation at audio frequencies*. Springer Science & Business Media.

<a id="4">[4]</a> Everstine, G. C., & Henderson, F. M. (1990). Coupled finite element/boundary element approach for fluid–structure interaction. *The Journal of the Acoustical Society of America*, 87(5), 1938-1947.

<a id="5">[5]</a> Chertock, G. (1972). Integral equation methods in sound radiation and scattering from arbitrary surfaces. *The Journal of the Acoustical Society of America*, 52(6A), 1588-1588.

<a id="6">[6]</a> Marburg, S. (2016). The Burton and Miller method: Unlocking another mystery of its coupling parameter. *Journal of Computational Acoustics*, 24(01), 1550016.

<a id="7">[7]</a> Börm, S., Grasedyck, L., & Hackbusch, W. (2003). Introduction to hierarchical matrices with applications. Engineering analysis with boundary elements, 27(5), 405-422.

<a id="8">[8]</a> Thompson, T. Ben. (2021). *Integral equation tutorials: 4. Low rank approximation of BEM matrices with adaptive cross approximation (ACA)*. Available at https://tbenthompson.com/book/tdes/low_rank.html, accessed 19 April 2025.

<a id="9">[9]</a> Grasedyck, L. (2005). Adaptive recompression of H-matrices for BEM. *Computing*, 74, 205-223.

<!---

 ```
 \left[ 
    \begin{matrix}
    \mathbf{H}_{11} & \mathbf{H}_{12} \\
    \mathbf{H}_{21} & \mathbf{H}_{22}
    \end{matrix}
 \right]
 \left\lbrace
    \begin{matrix}
    \mathbf{\phi}_1 \\
    \mathbf{\phi}_2
    \end{matrix}
 \right\rbrace
 +
  \left[ 
    \begin{matrix}
    \mathbf{G}_{11} & \mathbf{G}_{12} \\
    \mathbf{G}_{21} & \mathbf{G}_{22}
    \end{matrix}
 \right]
 \left\lbrace
    \begin{matrix}
    \mathbf{v}_1 \\
    \mathbf{v}_2
    \end{matrix}
 \right\rbrace
 =
  \left\lbrace
    \begin{matrix}
    \mathbf{\phi}^I_1 \\
    \mathbf{\phi}^I_2
    \end{matrix}
 \right\rbrace
 ```

 where $`\mathbf{v}_1 = \bar{\mathbf{v}}`$ (velocity B.C.s), $`\mathbf{\phi_2} = \bar{\mathbf{\phi}} + Z\mathbf{v}_2`$ (pressure and impedance B.C.s)

 ```
 \mathbf{H}_{21} \mathbf{\phi}_1 + \mathbf{H}_{22} (\bar{\mathbf{\phi}} + Z\mathbf{v}_2) + \mathbf{G}_{21} \bar{\mathbf{v}} + \mathbf{G}_{22} \mathbf{v}_2 = \mathbf{\phi}^I_2 \\
\mathbf{v}_2 = [\mathbf{G}_{22} + \mathbf{H}_{22} Z]^{-1}(\mathbf{\phi}^I_2 - \mathbf{H}_{21} \mathbf{\phi}_1 - \mathbf{H}_{22} \bar{\mathbf{\phi}} - \mathbf{G}_{21} \bar{\mathbf{v}})
 ```

 ```
 \mathbf{H}_{11} \mathbf{\phi}_1 + \mathbf{H}_{12} \bar{\mathbf{\phi}} + [\mathbf{H}_{12} Z + \mathbf{G}_{12}] [\mathbf{G}_{22} + \mathbf{H}_{22} Z]^{-1}(\mathbf{\phi}^I_2 - \mathbf{H}_{21} \mathbf{\phi}_1 - \mathbf{H}_{22} \bar{\mathbf{\phi}} - \mathbf{G}_{21} \bar{\mathbf{v}}) + \mathbf{G}_{11} \bar{\mathbf{v}}  = \mathbf{\phi}^I_1 \\
 [\mathbf{H}_{11} - \mathbf{\beta} \mathbf{H}_{21}] \mathbf{\phi}_1 = (\mathbf{\phi}^I_1 - \mathbf{\beta} \mathbf{\phi}^I_2) + [\mathbf{H}_{12} - \mathbf{\beta}\mathbf{H}_{22}] \bar{\mathbf{\phi}} + [\mathbf{G}_{11} - \mathbf{\beta} \mathbf{G}_{21}] \bar{\mathbf{v}}
 ```

 $`\mathbf{\beta} = [\mathbf{H}_{12} Z + \mathbf{G}_{12}] [\mathbf{G}_{22} + \mathbf{H}_{22} Z]^{-1}$

 -->
