The linear wave equation is satisfied over a field, at a given location $\bm{x}$ and at a time $t$
$$
\nabla ^2 \psi(\bm{x},t) = \frac{1}{c^2} \frac{\partial^2}{\partial t^2} \psi(\bm{x},t)
$$
where $\psi(\bm{x},t)$ is the velocity potential field and $c$ is the propagation velocity. The vector velocity is given $\bm{v}(\bm{x},t) = \nabla \psi(\bm{x},t)$ and the acoustic pressure is $p = -\rho \frac{\partial}{\partial t} \psi(\bm{x},t)$ with $\rho$ as the fluid mass density.

Using a harmonic, steady-state assumption, the velocity potential field can be rewritten as 
$$
\psi(\bm{x},t) = Re[\phi(\bm{x}) e^{-i\omega t}]
$$
given that the angular drive frequency is $\omega = 2 \pi f$, where $f$ is the frequency in Hertz. 

The harmonic form of the Helmholtz equation for the fluid is then
$$
\nabla ^2 \psi + k^2 \psi= \psi_{inc}
$$
where $k=\omega / c$ is the wavenumber. Three classes of boundary condition may be used on the surface, and only one type may be used on any partition of the surface:

1. **Velocity potential or pressure B.C.**: $\bm{\phi} = \bar{\bm{\phi}}$. This can represent a known non-zero surface pressure, or if set to zero, it represents a sound-soft (infinitely absorptive) boundary. ($\beta = 0, \gamma/\alpha = \bar{\bm{\phi}}$)
2. **Surface normal velocity B.C.**: $\bm{v}_n = \bar{\bm{v}}_n$. This can represent a known non-zero surface motion, or if set to zero, it represents a sound-hard (reflective) boundary. ($\alpha = 0, \gamma/\beta = \bar{\bm{v}}_n$)
3. **Impedance B.C.**: A known, linear relation between the surface normal velocity and the fluid pressure: $(i \omega \rho) \phi = p = Z v_n$ where $p$ is the pressure and $Z$ is the impedance. ($\gamma = 0, \alpha / \beta = -Z / (i \omega \rho)$)

These may be stated generally as variants of a Robin condition condition $\alpha \phi + \beta v_n = \gamma$, as shown above.

The incident sound sources (incident acoustic fields) can be defined similarly:

1. **Plane wave**: Has a defined strength (amplitude: $A$) and propagates along a unit vector $\bm{e}_d$. The phase is assumed to be zero at the origin, i.e. $\bm{x} \cdot \bm{e}_d = 0$: 
$$
\phi_I(\bm{x}) = A e^{ik(\bm{x} \cdot \bm{e}_d)}
$$
2. **Point source**: Also known as a volume source or spherical source. Has a defined amplitude $A$ and position vector $\bm{x}_I$:
$$
\phi_I(\bm{x}) = \frac{A}{4 \pi r}e^{ikr}; \qquad r = \bm{x} - \bm{x}_I
$$

The directional derivative, along a specified normal direction, will also be important. This is expressed as
$$
v_n(\bm{x}) = \nabla \phi (\bm{x}) \cdot \bm{n}(\bm{x}) = \frac{\partial \phi(\bm{x})}{\partial n(\bm{x})}
$$

The integral solution to the Helmholtz equation is written as
$$
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) - v_n(\bm{y}) g(\mathbf{x}, \mathbf{y}) d\bm{y} =
\left\lbrace 
\begin{array}{l l}
\phi(\mathbf{x})-\phi_I(\mathbf{x}) & \text{at field point, ${\bm{x}}$}\\
\phi(\mathbf{x})/2-\phi_I(\mathbf{x}) & \text{on surface for exterior problem}  \\
-\phi_I(\mathbf{x}) & \text{on surface for interior problem}
\end{array}
\right.
$$
where $\phi_I$ is the incident wave (free-field) velocity potential. $g$ and its normal derivative are given as:
$$
g(\mathbf{x}, \mathbf{y}) = \frac{e^{ikr}}{4 \pi r} \\
h(\mathbf{x}, \mathbf{y}) = \frac{\partial g(\mathbf{x}, \mathbf{y})}{\partial n(\mathbf{y})} = \left(ik - \frac{1}{r}\right) g(\bm{x}, \bm{y}) (\bm{e}_r \cdot \bm{n}(\bm{y}))
$$
with $\bm{r} = \bm{x} - \bm{y}$ as the vector pointing from $\bm{y}$ to $\bm{x}$ and $r = || \bm{r} ||$ as the magnitude (distance) and $\bm{e}_r = \bm{r}/r$ as the unit vector form.

Within an element, the shape functions are used to represent the variation in value in the domain, such that
$$
\phi(\bm{\xi}) = \sum_i N_i^e(\bm{\xi}) \phi^e_i
$$
with $\bm{\xi}$ as the natural domain coordinate within the element, $\phi^e_i$ as the value at element node $i$, and $N_i^e(\bm{\xi})$ as the shape function corresponding to node $i$ evaluated at $\bm{\xi}$.

The integral form can be discretized over the surface this way as
$$
\bm{H} \bm{\phi} - \bm{G} \bm{v}_n = 
\left\lbrace 
\begin{array}{l l}
\bm{\phi}/2-\bm{\phi}_I & \text{on surface for exterior problem}  \\

-\bm{\phi}_I & \text{on surface for interior problem}
\end{array}
\right.
$$

Examining the $\bm{H}\bm{\phi}$ term, each row $i$ of $\bm{H}$ corresponds to the contribution (influence) of the velocity potentials at all surface points to the velocity potential at node $i$, which can be split into an integral over the elements: 
$$
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) d\bm{y} \approx
\sum _e \int_{S_e} \phi (\mathbf{y_e}) h(\mathbf{x}, \mathbf{y_e}) d\bm{y}_e =
\left(\sum _e \int_{S_e} \sum _i N_i(\bm{\xi})h(\bm{x}, \bm{y}_e(\bm{\xi}))d\bm{y}_e \right)\bm{\phi}
$$
Note that the nodal values of velocity potential are assembled into the vector $\bm{\phi}$ and the in-element integral can be calculated using quadrature rules.

The velocity potentials must first be solved on the surface, giving the matrix equation (in an exterior problem):
$$
\left[\bm{H} - \frac{1}{2} \bm{I} \right] \bm{\phi} = \bm{G} \bm{v}_n - \bm{\phi}_I
$$
For the three possible boundary conditions, this is solved as:
$$
\begin{array}{l l}
\bm{G} \bm{v}_n = \bm{\phi}_I -\left[\bm{H} - \frac{1}{2} \bm{I} \right] \bar{\bm{\phi}}  & \text{pressure B.C.} \\
\left[\bm{H} - \frac{1}{2} \bm{I} \right] \bm{\phi} = \bm{G} \bar{\bm{v}}_n - \bm{\phi}_I & \text{normal velocity B.C.} \\
\left[\bm{H} - \frac{1}{2} \bm{I} - \frac{i \omega \rho}{Z} \bm{G}\right] \bm{\phi} = -\bm{\phi}_I  & \text{impedance B.C.}
\end{array}
$$

Once the velocity potential and normal velocity fields are known on the surface, the velocity potential for an arbitrary point $\mathbf{x}$ in the interior or exterior field can be found from the general solution via: 

$$
\int_S \phi (\mathbf{y}) h(\mathbf{x}, \mathbf{y}) - v_n(\bm{y}) g(\mathbf{x}, \mathbf{y}) d\bm{y} =
\phi(\mathbf{x})-\phi_I(\mathbf{x}) \\
\bm{\phi}_{fp} = -\bm{M} \bm{\phi} + \bm{L} \bm{v}_n + \bm{\phi}_I
$$

where $\bm{M}$ and $\bm{G}$ are analagously constructed to $\bm{H}$ and $\bm{G}$, but are (generally) rectangular matrices of dimension $n_{fp} \times n_s$ where $n_{fp}$ is the number of field points and $n_s$ is the number of surface points. $\bm{M}$ and $\bm{G}$ represent the influence of the surface fields on the velocity potential at each field point. Note that $\bm{\phi}_I$ contains the vector of $n_{fp}$ incident wave velocity potentials at each field point.

<!---

 $$
 \left[ 
    \begin{matrix}
    \bm{H}_{11} & \bm{H}_{12} \\
    \bm{H}_{21} & \bm{H}_{22}
    \end{matrix}
 \right]
 \left\lbrace
    \begin{matrix}
    \bm{\phi}_1 \\
    \bm{\phi}_2
    \end{matrix}
 \right\rbrace
 +
  \left[ 
    \begin{matrix}
    \bm{G}_{11} & \bm{G}_{12} \\
    \bm{G}_{21} & \bm{G}_{22}
    \end{matrix}
 \right]
 \left\lbrace
    \begin{matrix}
    \bm{v}_1 \\
    \bm{v}_2
    \end{matrix}
 \right\rbrace
 =
  \left\lbrace
    \begin{matrix}
    \bm{\phi}^I_1 \\
    \bm{\phi}^I_2
    \end{matrix}
 \right\rbrace
 $$

 where $\bm{v}_1 = \bar{\bm{v}}$ (velocity B.C.s), $\bm{\phi_2} = \bar{\bm{\phi}} + Z\bm{v}_2$ (pressure and impedance B.C.s)

 $$
 \bm{H}_{21} \bm{\phi}_1 + \bm{H}_{22} (\bar{\bm{\phi}} + Z\bm{v}_2) + \bm{G}_{21} \bar{\bm{v}} + \bm{G}_{22} \bm{v}_2 = \bm{\phi}^I_2 \\
\bm{v}_2 = [\bm{G}_{22} + \bm{H}_{22} Z]^{-1}(\bm{\phi}^I_2 - \bm{H}_{21} \bm{\phi}_1 - \bm{H}_{22} \bar{\bm{\phi}} - \bm{G}_{21} \bar{\bm{v}})
 $$

 $$
 \bm{H}_{11} \bm{\phi}_1 + \bm{H}_{12} \bar{\bm{\phi}} + [\bm{H}_{12} Z + \bm{G}_{12}] [\bm{G}_{22} + \bm{H}_{22} Z]^{-1}(\bm{\phi}^I_2 - \bm{H}_{21} \bm{\phi}_1 - \bm{H}_{22} \bar{\bm{\phi}} - \bm{G}_{21} \bar{\bm{v}}) + \bm{G}_{11} \bar{\bm{v}}  = \bm{\phi}^I_1 \\
 [\bm{H}_{11} - \bm{\beta} \bm{H}_{21}] \bm{\phi}_1 = (\bm{\phi}^I_1 - \bm{\beta} \bm{\phi}^I_2) + [\bm{H}_{12} - \bm{\beta}\bm{H}_{22}] \bar{\bm{\phi}} + [\bm{G}_{11} - \bm{\beta} \bm{G}_{21}] \bar{\bm{v}}
 $$

 $\bm{\beta} = [\bm{H}_{12} Z + \bm{G}_{12}] [\bm{G}_{22} + \bm{H}_{22} Z]^{-1}$

 -->