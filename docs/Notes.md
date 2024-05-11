Starting with typical structural dynamics equation for an *n* DoF system:
$$
\mathbf{M} \ddot{\mathbf{u}}+\mathbf{C} \dot{\mathbf{u}}+\mathbf{K} \mathbf{u}=\mathbf{f}
$$
Let's make the harmonic, steady-state assumption:
$$
\mathbf{f}=\hat{\mathbf{f}}e^{i\omega t} \\
\mathbf{u}=\hat{\mathbf{u}}e^{i\omega t}
$$
which gives
$$
[-\omega^2 \mathbf{M} + i\omega \mathbf{C} +\mathbf{K}] \hat{\mathbf{u}}=\hat{\mathbf{f}}
$$
Neglect the hats and introduce another forcing term from the effect of the fluid on the structural surface:
$$
\mathbf{Z}(\omega)\dot{\mathbf{u}}=\mathbf{f}-\mathbf{N}\mathbf{A}\mathbf{p}
$$
where $\mathbf{Z}(\omega)$ is the *n* by *n* structural impedance matrix ($\mathbf{Z}(\omega)=\frac{1}{i\omega}[-\omega^2 \mathbf{M} + i\omega \mathbf{C} +\mathbf{K}] $), $\mathbf{p}$ is the vector of *m* pressures on the boundary, $\mathbf{N}$ is the *n* by *m* matrix which selects the structural boundary DoF and orients in the outward normal direction, and $\mathbf{A}$ is the *m* by *m* boundary area matrix. A positive pressure leads to a negative force, since the convention uses an outward pointing normal vector.

The structural solution can be separated between the 'dry' and 'wet' (or total) solution:
$$
\dot{\mathbf{u}}=\dot{\mathbf{u}}_{dry}-\mathbf{Z}^{-1}\mathbf{N}\mathbf{A}\mathbf{p}
$$
where $\dot{\mathbf{u}}_{dry}=\mathbf{Z}^{-1} \mathbf{f}$ can be solved independently and prior to solving the fluid portion. Now the fluid portion will be expressed before forming the coupled equation.

The harmonic form of the Helmholtz equation for the exterior fluid is
$$
\nabla ^2 p+k^2 p=0
$$
where *p* is the pressure field and $k=\omega/ c$ is the wavenumber. The solution is
$$
\int_S p(\mathbf{x}) \frac{\partial g(r)}{\partial n} - \frac{\partial p(\mathbf{x})}{\partial n}g(r)dS =
\left\lbrace 
\begin{array}{l l}
p(\mathbf{x}')/2-p_I & \mathbf{x}' \text{ on boundary} \\
p(\mathbf{x}')-p_I & \mathbf{x}' \text{ outside boundary}\\
-p_I & \mathbf{x}' \text{ inside boundary}
\end{array}
\right.
$$
where $p_I$ is the incident pressure on the surface and  $r=|\mathbf{x}' - \mathbf{x}|$ is the distance from a point on the surface $\mathbf{x}$ to a point of interest $\mathbf{x}'$. $g(r)$ is the Green's function
$$
g(r)=e^{-ikr}/4\pi r \\
\frac{\partial g(r)}{\partial n} = -(ik+1/r) g(r) \frac{\partial r}{\partial n}
$$
with $\frac{\partial r}{\partial n}$ as the cosine of the angle between the surface normal vector at $\mathbf{x}$ and the vector $\mathbf{r}=\mathbf{x}' - \mathbf{x}$. 

Momentum considerations give
$$
\frac{\partial p(\mathbf{x})}{\partial n} = -i\omega \rho \dot{u}_n(\mathbf{x})
$$
where $\dot{u}_n$ is the normal velocity on the boundary and $\rho$ is the fluid mass density. Thus, the pressure at every point on the surface can be related to the pressure and surface normal velocity on every other point on the surface. Thus, the integral on the surface can be written:
$$
p(\mathbf{x}')/2 +\int_S p(\mathbf{x}) (ik+1/r)\frac{e^{-ikr}}{4\pi r} \frac{\partial r}{\partial n} dS = i\omega \rho \int_S  \dot{u}_n(\mathbf{x})  \frac{e^{-ikr}}{4\pi r}dS + p_I
$$
These integrals can be discretized over the surface and precomputed
$$
\mathbf{\Gamma}\mathbf{p}=\mathbf{\Lambda}\dot{\mathbf{u}}_n+\mathbf{p}_I
$$
which is the typical discretized form of the BEM equations, with $\mathbf{\Gamma}$ and $\mathbf{\Lambda}$ as the dense, complex, unsymmetric BEM influence matrices.

Now, going back to the structural equations, the normal velocities on the surface can be related to the full set of structural velocity DoF:
$$
\dot{\mathbf{u}}_n = \mathbf{N}^T \dot{\mathbf{u}}
$$
Combining this with equations \ref{eq:BEM} and \ref{eq:partitioned_response}:
$$
[\mathbf{\Gamma} + \mathbf{N}^T \mathbf{Z}^{-1}\mathbf{N}\mathbf{A}]\mathbf{p}=\mathbf{\Lambda}\mathbf{N}^T\dot{\mathbf{u}}_{dry}+\mathbf{p}_I
$$
Which can be solved to find the submerged surface pressures from the precomputed structural impedance $\mathbf{Z}$, dry velocity responses $\dot{\mathbf{u}}_{dry}$, and the BEM influence matrices $\mathbf{\Gamma}$ and $\mathbf{\Lambda}$. Computationally, it is advantageous to compute $\mathbf{Z}^{-1}_n =\mathbf{N}^T \mathbf{Z}^{-1}\mathbf{N}$ and $\mathbf{N}^T\dot{\mathbf{u}}_{dry}$ directly by only saving the surface-normal portions, rather than computing the full responses and multiplying by $\mathbf{N}$.

The submerged surface velocities can be recovered either by modifying the dry responses (Eq. $\ref{eq:partitioned_response}$) or by re-solving the structural response with the known surface pressures (Eq. $\ref{eq:structure_with_fluid}$). Computationally, the second method is preferred since $\mathbf{Z}$ can be stored (in a factored state) on disk and this avoids needing to hold $\dot{\mathbf{u}}_{dry}$ after computing the submerged response.

With the known surface pressures and submerged surface normal velocities, the pressure or velocity at any point in the field can be computed from Eq. $\ref{eq:Helmholtz_integral}$ with analogous intermediate steps to Eqs. $\ref{eq:BEM_integral}$  and $\ref{eq:BEM}$. Note that this stage doesn't require matrix factorization/solution, only multiplications.

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