The linear wave equation is satisfied over a field, at a given location $`\mathbf{x}`$ and at a time $`t$

```math
\nabla ^2 \psi(\mathbf{x},t) = \frac{1}{c^2} \frac{\partial^2}{\partial t^2} \psi(\mathbf{x},t)
```

where $`\psi(\mathbf{x},t)`$ is the velocity potential field and $`c`$ is the propagation velocity. The vector velocity is given $`\mathbf{v}(\mathbf{x},t) = \nabla \psi(\mathbf{x},t)`$ and the acoustic pressure is $`p = -\rho \frac{\partial}{\partial t} \psi(\mathbf{x},t)`$ with $`\rho`$ as the fluid mass density.

Using a harmonic, steady-state assumption, the velocity potential field can be rewritten as 

```math
\psi(\mathbf{x},t) = Re[\phi(\mathbf{x}) e^{-i\omega t}]
```

given that the angular drive frequency is $`\omega = 2 \pi f`$, where $`f`$ is the frequency in Hertz. 

The harmonic form of the Helmholtz equation for the fluid is then

```math
\nabla ^2 \psi + k^2 \psi= \psi_{inc}
```

where $`k=\omega / c`$ is the wavenumber. Three classes of boundary condition may be used on the surface, and only one type may be used on any partition of the surface:

1. **Velocity potential or pressure B.C.**: $`\mathbf{\phi} = \bar{\mathbf{\phi}}`$. This can represent a known non-zero surface pressure, or if set to zero, it represents a sound-soft (infinitely absorptive) boundary. ( $`\beta = 0, \gamma/\alpha = \bar{\mathbf{\phi}}`$ )
2. **Surface normal velocity B.C.**: $`\mathbf{v}_n = \bar{\mathbf{v}}_n$. This can represent a known non-zero surface motion, or if set to zero, it represents a sound-hard (reflective) boundary. ( $`\alpha = 0, \gamma/\beta = \bar{\mathbf{v}}_n`$ )
3. **Impedance B.C.**: A known, linear relation between the surface normal velocity and the fluid pressure: $`(i \omega \rho) \phi = p = Z v_n`$ where $`p`$ is the pressure and $`Z`$ is the impedance. ( $`\gamma = 0, \alpha / \beta = -Z / (i \omega \rho)`$ )

These may be stated generally as variants of a Robin condition condition $`\alpha \phi + \beta v_n = \gamma`$, as shown above.

The incident sound sources (incident acoustic fields) can be defined similarly:

1. **Plane wave**: Has a defined strength (amplitude: $`A$) and propagates along a unit vector $`\mathbf{e}_d`$. The phase is assumed to be zero at the origin, i.e. $`\mathbf{x} \cdot \mathbf{e}_d = 0$: 

```math
\phi_I(\mathbf{x}) = A e^{ik(\mathbf{x} \cdot \mathbf{e}_d)}
```

2. **Point source**: Also known as a volume source or spherical source. Has a defined amplitude $`A`$ and position vector $`\mathbf{x}_I$:

```math
\phi_I(\mathbf{x}) = \frac{A}{4 \pi r}e^{ikr}; \qquad r = \mathbf{x} - \mathbf{x}_I
```

The directional derivative, along a specified normal direction, will also be important. This is expressed as

```math
v_n(\mathbf{x}) = \nabla \phi (\mathbf{x}) \cdot \mathbf{n}(\mathbf{x}) = \frac{\partial \phi(\mathbf{x})}{\partial n(\mathbf{x})}
```

The integral solution to the Helmholtz equation is written as

```math
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) - v_n(\mathbf{y}) g(\mathbf{x}, \mathbf{y}) d\mathbf{y} =
\left\lbrace 
\begin{array}{l l}
\phi(\mathbf{x})-\phi_I(\mathbf{x}) & \text{at field point, $`{\mathbf{x}}$}\\
\phi(\mathbf{x})/2-\phi_I(\mathbf{x}) & \text{on surface for exterior problem}  \\
-\phi_I(\mathbf{x}) & \text{on surface for interior problem}
\end{array}
\right.
```

where $`\phi_I`$ is the incident wave (free-field) velocity potential. $`g`$ and its normal derivative are given as:

```math
g(\mathbf{x}, \mathbf{y}) = \frac{e^{ikr}}{4 \pi r} \\
h(\mathbf{x}, \mathbf{y}) = \frac{\partial g(\mathbf{x}, \mathbf{y})}{\partial n(\mathbf{y})} = \left(ik - \frac{1}{r}\right) g(\mathbf{x}, \mathbf{y}) (\mathbf{e}_r \cdot \mathbf{n}(\mathbf{y}))
```

with $`\mathbf{r} = \mathbf{x} - \mathbf{y}`$ as the vector pointing from $`\mathbf{y}`$ to $`\mathbf{x}`$ and $`r = || \mathbf{r} ||`$ as the magnitude (distance) and $`\mathbf{e}_r = \mathbf{r}/r`$ as the unit vector form.

Within an element, the shape functions are used to represent the variation in value in the domain, such that

```math
\phi(\mathbf{\xi}) = \sum_i N_i^e(\mathbf{\xi}) \phi^e_i
```

with $`\mathbf{\xi}`$ as the natural domain coordinate within the element, $`\phi^e_i`$ as the value at element node $`i`$, and $`N_i^e(\mathbf{\xi})`$ as the shape function corresponding to node $`i`$ evaluated at $`\mathbf{\xi}`$.

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


Examining the $`\mathbf{H}\mathbf{\phi}`$ term, each row $`i`$ of $`\mathbf{H}`$ corresponds to the contribution (influence) of the velocity potentials at all surface points to the velocity potential at node $`i$, which can be split into an integral over the elements: 

```math
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) d\mathbf{y} \approx
\sum_e \int_{S_e} \phi (\mathbf{y_e}) h(\mathbf{x}, \mathbf{y_e}) d\mathbf{y}_e = 
\left(\sum _e \int_{S_e} \sum _i N_i(\mathbf{\xi})h(\mathbf{x}, \mathbf{y}_e(\mathbf{\xi}))d\mathbf{y}_e \right)\mathbf{\phi}
```

Note that the nodal values of velocity potential are assembled into the vector $`\mathbf{\phi}`$ and the in-element integral can be calculated using quadrature rules.

The velocity potentials must first be solved on the surface, giving the matrix equation (in an exterior problem):

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

Once the velocity potential and normal velocity fields are known on the surface, the velocity potential for an arbitrary point $`\mathbf{x}`$ in the interior or exterior field can be found from the general solution via: 

```math
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) - v_n(\mathbf{y}) g(\mathbf{x}, \mathbf{y}) d\mathbf{y} =
\phi(\mathbf{x})-\phi_I(\mathbf{x}) \\
\mathbf{\phi}_{fp} = -\mathbf{M} \mathbf{\phi} + \mathbf{L} \mathbf{v}_n + \mathbf{\phi}_I
```

where $`\mathbf{M}`$ and $`\mathbf{G}`$ are analagously constructed to $`\mathbf{H}`$ and $`\mathbf{G}$, but are (generally) rectangular matrices of dimension $`n_{fp} \times n_s`$ where $`n_{fp}`$ is the number of field points and $`n_s`$ is the number of surface points. $`\mathbf{M}`$ and $`\mathbf{G}`$ represent the influence of the surface fields on the velocity potential at each field point. Note that $`\mathbf{\phi}_I`$ contains the vector of $`n_{fp}`$ incident wave velocity potentials at each field point.

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
