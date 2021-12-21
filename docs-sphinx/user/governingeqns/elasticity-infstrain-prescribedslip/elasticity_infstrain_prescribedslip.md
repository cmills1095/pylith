## Elasticity with Infinitesimal Strain and Prescribed Slip on Faults

For each fault, which is an internal interface, we add a boundary condition to
the elasticity equation prescribing the jump in the displacement field across
the fault, $$\begin{gathered}
  \label{eqn:bc:prescribed_slip}
  \vec{u}^+(\vec{x},t) - \vec{u}^-(\vec{x},t) - \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,\end{gathered}$$
where $\vec{u}^+$ is the displacement vector on the &ldquo;positive&rdquo;
side of the fault, $\vec{u}^-$ is the displacement vector on the
&ldquo;negative&rdquo; side of the fault, $\vec{d}$ is the slip vector on the
fault, and $\vec{n}$ is the fault normal which points from the negative side
of the fault to the positive side of the fault. Note that as an alternative to
prescribing the jump in displacement across the fault, we can also prescribe
the jump in velocity or acceleration across the fault in terms of slip rate or
slip acceleration, respectively.

Using a domain decomposition approach for constraining the fault slip yields
additional Neumann-like boundary conditions on the fault surface,
$$\begin{gathered}
  \bm{\sigma} \cdot \vec{n} = -\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^+}, \\
  \bm{\sigma} \cdot \vec{n} = +\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^-},\end{gathered}$$
where $\vec{\lambda}$ is the vector of Lagrange multipliers corresponding to
the tractions applied to the fault surface to generate the prescribed slip.

<div id="tab:notation:elasticity:prescribed:slip">

| **Category**                   |   **Symbol**    | **Description**                                                                                   |
|:-------------------------------|:---------------:|:--------------------------------------------------------------------------------------------------|
| Unknowns                       |    $\vec{u}$    | Displacement field                                                                                |
|                                |    $\vec{v}$    | Velocity field                                                                                    |
|                                | $\vec{\lambda}$ | Lagrange multiplier field                                                                         |
| Derived quantities             |  $\bm{\sigma}$  | Cauchy stress tensor                                                                              |
|                                | $\bm{\epsilon}$ | Cauchy strain tensor                                                                              |
| Common constitutive parameters |     $\rho$      | Density                                                                                           |
|                                |      $\mu$      | Shear modulus                                                                                     |
|                                |       $K$       | Bulk modulus                                                                                      |
| Source terms                   |    $\vec{f}$    | Body force per unit volume, for example $\rho \vec{g}$                                            |
|                                |    $\vec{d}$    | Slip vector field on the fault corresponding to a jump in the displacement field across the fault |

Mathematical notation for elasticity equation with infinitesimal strain and
prescribed slip on faults.

</div>

### Quasistatic

As in the case of elasticity without faults, we first consider the quasistatic
case in which we neglect the inertial term
($\rho \frac{\partial \vec{v}}{\partial t} \approx \vec{0}$). We place all of
the terms in the elasticity equation on the LHS, consistent with implicit time
stepping. Our solution vector consists of both displacements and Lagrange
multipliers, and the strong form for the system of equations is
$$\begin{gathered}
  % Solution
  \vec{s}^T = \left( \vec{u} \quad \vec{\lambda} \right)^T \\
  % Elasticity
  \vec{f}(\vec{x},t) + \bm{\nabla} \cdot \bm{\sigma}(\vec{u}) = \vec{0} \text{ in }\Omega, \\
  % Neumann
  \bm{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on }\Gamma_\tau, \\
  % Dirichlet
  \vec{u} = \vec{u}_0(\vec{x},t) \text{ on }\Gamma_u, \\
  % Prescribed slip
  \vec{u}^+(\vec{x},t) - \vec{u}^-(\vec{x},t) - \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,  \\
  \bm{\sigma} \cdot \vec{n} = -\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^+}, \text{ and}\\
  \bm{\sigma} \cdot \vec{n} = +\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^-}.\end{gathered}$$
We create the weak form by taking the dot product with the trial function
${\vec{\psi}_\mathit{trial}^{u}}$ or ${\vec{\psi}_\mathit{trial}^{\lambda}}$
and integrating over the domain. After using the divergence theorem and
incorporating the Neumann boundary and fault interface conditions, we have
$$\begin{gathered}
  % Elasticity
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\bm{\sigma}(\vec{u}) \, d\Omega
  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma,
  + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{u^+}} \cdot \left(-\vec{\lambda}(\vec{x},t)\right)
  + {\vec{\psi}_\mathit{trial}^{u^-}} \cdot \left(+\vec{\lambda}(\vec{x},t)\right)\, d\Gamma = 0\\
  % Prescribed slip
  \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \left(
    -\vec{u}^+(\vec{x},t) + \vec{u}^-(\vec{x},t) + \vec{d}(\vec{x},t) \right) \, d\Gamma = 0.\end{gathered}$$
We have chosen the signs in the last equation so that the Jacobian will be
symmetric with respect to the Lagrange multiplier. We solve the system of
equations using implicit time stepping, making use of residuals functions and
Jacobians for the LHS.

#### Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have $$\begin{aligned}
 % Fu
F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}\mathclap{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} : %
  {\color{blue}%
  \underbrace{\color{black}-\bm{\sigma}(\vec{u})}_{\color{blue}\mathclap{\bm{f^u_1}}}} \, d\Omega
  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}\mathclap{\vec{f}^u_0}}} \, d\Gamma
  + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{u^+}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\left(-\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}\mathclap{\vec{f}^u_0}}}
  + {\vec{\psi}_\mathit{trial}^{u^-}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\left(+\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}\mathclap{\vec{f}^u_0}}}\, d\Gamma \\
  % Fl
  F^\lambda(t,s,\dot{s}) &= \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\left(
    -\vec{u}^+(\vec{x},t) + \vec{u}^-(\vec{x},t) + \vec{d}(\vec{x},t) \right)}_{\color{blue}\mathclap{\vec{f}^\lambda_0}}} \, d\Gamma, \\
  % Gu
  G^u(t,s) &= 0 \\
  % Gl
  G^\lambda(t,s) &= 0\end{aligned}$$ Compared to the quasistatic elasticity
case without a fault, we have simply added additional pointwise functions
associated with the fault. Our fault implementation does not change the
formulation for the materials or external Dirichlet or Neumann boundary
conditions.

#### Jacobian Pointwise Functions

The LHS Jacobians are: $$\begin{aligned}
  % J_F uu
  J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{u}}
      = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\bm{C} : \frac{1}{2}(\nabla + \nabla^T){\vec{\psi}_\mathit{basis}^{u}}
\, d\Omega
      = \int_\Omega {\psi_\mathit{trial}^{u}}_{i,k} \, %
  {\color{blue}%
  \underbrace{\color{black}\left( -C_{ikjl} \right)}_{\color{blue}\mathclap{J_{f3}^{uu}}}} \, {\psi_\mathit{basis}^{u}}_{j,l}\, d\Omega \\
  % J_F ul
  J_F^{u\lambda} &= \frac{\partial F^u}{\partial \lambda} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{\lambda}}
      = \int_{\Gamma_{f}} {\psi_\mathit{trial}^{u^+}}_i %
  {\color{blue}%
  \underbrace{\color{black}\left(-\delta_{ij}\right)}_{\color{blue}\mathclap{J^{u\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j
                   + {\psi_\mathit{trial}^{u^-}}_i %
  {\color{blue}%
  \underbrace{\color{black}\left(+\delta_{ij}\right)}_{\color{blue}\mathclap{J^{u\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j\, d\Gamma, \\
  % J_F lu
  J_F^{\lambda u} &= \frac{\partial F^\lambda}{\partial u} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{u}}
      = \int_{\Gamma_{f}} {\psi_\mathit{trial}^{\lambda}}_i
                    %
  {\color{blue}%
  \underbrace{\color{black}\left(-\delta_{ij}\right)}_{\color{blue}\mathclap{J^{\lambda u}_{f0}}}} {\psi_\mathit{basis}^{u^+}}_j
                    + {\psi_\mathit{trial}^{\lambda}}_i %
  {\color{blue}%
  \underbrace{\color{black}\left(+\delta_{ij}\right)}_{\color{blue}\mathclap{J^{\lambda u}_{f0}}}} {\psi_\mathit{basis}^{u^-}}_j \, d\Gamma, \\
  % J_F ll
  J_F^{\lambda \lambda} &= \bm{0}\end{aligned}$$ This LHS Jacobian has the
structure
$$J_F = \left( \begin{array} {cc} J_F^{uu} & J_F^{u\lambda} \\ J_F^{\lambda u} & 0 \end{array} \right)
      = \left( \begin{array} {cc} J_F^{uu} & C^T \\ C & 0 \end{array} \right),$$
where $C$ contains entries of $\pm 1$ for degrees of freedom on the two sides
of the fault. The Schur complement of $J$ with respect to $J_F^{uu}$ is
$-C\left(J_F^{uu}\right)^{-1}C^T$.

### Dynamic

The equation prescribing fault slip is independent of the Lagrange multiplier,
so we do not have a system of equations that we can put in the form
$\dot{s} = G^*(t,s)$. Instead, we have a differential-algebraic set of
equations (DAEs), which we solve using an implicit-explicit (IMEX) time
integration scheme. As in the case of dynamic elasticity without faults, we
introduce the velocity ($\vec{v}$) as an unknown to turn the elasticity
equation into two first order equations. The strong form for our system of
equations is: $$\begin{gathered}
  % Solution
  \vec{s}^T = \left( \vec{u} \quad \vec{v} \quad \vec{\lambda} \right)^T \\
  % Displacement-velocity
  \frac{\partial \vec{u}}{\partial t} = \vec{v}, \\
  % Elasticity
  \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} =
  \vec{f}(\vec{x},t) + \nabla \cdot \bm{\sigma}(\vec{u}), \\
  % Neumann BC
  \bm{\sigma} \cdot \vec{n} = \vec{\tau} \text{ on } \Gamma_\tau. \\
  % Dirichlet BC
  \vec{u} = \vec{u}_0 \text{ on } \Gamma_u, \\
  % Presribed slip
  \vec{u}^+(\vec{x},t) - \vec{u}^-(\vec{x},t) - \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f,  \\
  \bm{\sigma} \cdot \vec{n} = -\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^+}, \text{ and}\\
  \bm{\sigma} \cdot \vec{n} = +\vec{\lambda}(\vec{x},t) \text{ on }\Gamma_{f^-}.\end{gathered}$$
The differentiation index is 2 because we must take the second time derivative
of the prescribed slip equation to match the order of the time derivative in
the elasticity equation, $$\begin{gathered}
  \frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} -
  \frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} = \vec{0}.\end{gathered}$$
We generate the weak form in the usual way, $$\begin{gathered}
  % Displacement-velocity
  \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \frac{\partial \vec{u}}{\partial t} \, d\Omega =
  \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{v} \, d\Omega, \\
  % Elasticity
  \begin{multlined}
  \int_{\Omega} {\vec{\psi}_\mathit{trial}^{v}} \cdot \rho(\vec{x}) \frac{\partial \vec{v}}{\partial t} \, d\Omega
  = \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\bm{\sigma}(\vec{u}) \, d\Omega
  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma \\
  + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\vec{\lambda}(\vec{x},t)\right)
  + {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(+\vec{\lambda}(\vec{x},t)\right)\, d\Gamma,
  \end{multlined}\\
  % Prescribed slip
  \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \left(
    \frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} -
    \frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} \right) \, d\Gamma = 0.\end{gathered}$$

For compatibility with PETSc TS IMEX implementations, we need $\dot{\vec{s}}$
on the LHS for the explicit part (displacement-velocity and elasticity
equations) and we need $\vec{\lambda}$ in the equation for the implicit part
(prescribed slip equation). We first focus on the explicit part and create a
lumped LHS Jacobian matrix, $M$, so that we have $$\begin{gathered}
  % Displacement-velocity
  \label{eqn:displacement:velocity:prescribed:slip:weak:form}
  \frac{\partial \vec{u}}{\partial t} = M_u^{-1} \int_{\Omega} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{v} \, d\Omega, \\
  % Elasticity
  \label{eqn:elasticity:prescribed:slip:dynamic:weak:form}
  \begin{multlined}
  \frac{\partial \vec{v}}{\partial t}
  = M_v^{-1} \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\bm{\sigma}(\vec{u}) \, d\Omega
  + M_v^{-1} \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma \\
  + M_{v^+}^{-1} \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\vec{\lambda}(\vec{x},t)\right) \, d\Gamma
  + M_{v^-}^{-1} \int_{\Gamma_{f}}{\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(+\vec{\lambda}(\vec{x},t)\right) \, d\Gamma,
\end{multlined}\\
% Mu
M_u = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{u}}_i \delta_{ij} {\psi_\mathit{basis}^{u}}_j \, d\Omega \right), \\
% Mv
M_v = \mathit{Lump}\left( \int_\Omega {\psi_\mathit{trial}^{v}}_i \rho(\vec{x}) \delta_{ij} {\psi_\mathit{basis}^{v}}_j \, d\Omega \right).\end{gathered}$$
Now, focusing on the implicit part we want to introduce $\vec{\lambda}$ into
the prescribed slip equation. We solve the elasticity equation for
$\frac{\partial \vec{v}}{\partial t}$, create the weak form, and substitute
into the prescribed slip equation. Solving the elasticity equation for
$\frac{\partial \vec{v}}{\partial t}$, we have
$$\frac{\partial \vec{v}}{\partial t} = \frac{1}{\rho(x)} \vec{f}(\vec{x},t) + \frac{1}{\rho(x)} \left(\nabla \cdot \bm{\sigma}(\vec{u}) \right),$$
and the corresponding weak form is $$\label{eqn:prescribed:slip:DAE:weak:form}
  \int_{\Omega} {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{\partial \vec{v}}{\partial t} \, d\Omega
  = \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \vec{f}(\vec{x},t) + {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \bm{\sigma}(\vec{u}) \right) \, d\Omega,$$
We apply the divergence theorem,
$$\int_{\Omega} \nabla \cdot \vec{F} \, d\Omega = \int_\Gamma \vec{F} \cdot \vec{n} \, d\Gamma,$$
with
$\vec{F} = {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \bm{\sigma}(\vec{u})\right)$
to get
$$\int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \bm{\sigma}(\vec{u}) \right) \, d\Omega
  = \int_\Gamma {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( \frac{1}{\rho(x)} \bm{\sigma}(\vec{u}) \cdot \vec{n} \right) \, d\Gamma
  + \int_\Omega \nabla{\vec{\psi}_\mathit{trial}^{v}} : \left(-\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right)
  + {\vec{\psi}_\mathit{trial}^{v}} \cdot \left(-\frac{\nabla \rho(\vec{x})}{\rho^2} \cdot \bm{\sigma}(\vec{u}) \right) \, d\Omega.$$
Restricting the trial function to $v^+$ and $v^-$ while recognizing that there
is no overlap between the external Neumann boundary conditions $\Gamma_\tau$
and the fault interfaces $\Gamma_f$, yields $$\begin{gathered}
  \int_\Omega {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \bm{\sigma}(\vec{u}) \right) \, d\Omega
  = \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\frac{1}{\rho(x)} \vec{\lambda} \right) \, d\Gamma
  + \int_\Omega \nabla{\vec{\psi}_\mathit{trial}^{v^+}} : \left(-\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right)
  + {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \, \left(-\frac{\nabla \rho(\vec{x})}{\rho^2} \cdot \bm{\sigma}(\vec{u}) \right) \, d\Omega, \\
  \int_\Omega {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \frac{1}{\rho(x)} \left(\nabla \cdot \bm{\sigma}(\vec{u}) \right) \, d\Omega
  = \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(+\frac{1}{\rho(x)} \vec{\lambda} \right) \, d\Gamma
  + \int_\Omega \nabla{\vec{\psi}_\mathit{trial}^{v^-}} : \left(-\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right)
  + {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \, \left(-\frac{\nabla \rho(\vec{x})}{\rho^2} \cdot \bm{\sigma}(\vec{u}) \right) \, d\Omega. \end{gathered}$$
Picking
${\vec{\psi}_\mathit{trial}^{v}}={\vec{\psi}_\mathit{trial}^{\lambda}}$ and
substituting into equation&nbsp;[\[eqn:prescribed:slip:DAE:weak:form\]][1]
gives $$\begin{gathered}
  \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \left(
    \frac{\partial \vec{v}^+}{\partial t} - \frac{\partial \vec{v}^-}{\partial t} -
    \frac{\partial^2 \vec{d}(\vec{x},t)}{\partial t^2} \right) \, d\Gamma = \\
  \int_\Omega {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left( \frac{1}{\rho(\vec{x})} \vec{f}(\vec{x}, t)
    -\frac{\nabla \rho(\vec{x})}{\rho^2(\vec{x})} \cdot \bm{\sigma}(\vec{u}) \right)
  + \nabla {\vec{\psi}_\mathit{trial}^{v^+}} : \left(-\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right) \, d\Omega
  + \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot \left(-\frac{1}{\rho(\vec{x})} \vec{\lambda} \right) \, d\Gamma \\
  - \int_\Omega {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left( \frac{1}{\rho(\vec{x})} \vec{f}(\vec{x}, t)
  -\frac{\nabla \rho(\vec{x})}{\rho^2(\vec{x})} \cdot \bm{\sigma}(\vec{u}) \right)
  + \nabla {\vec{\psi}_\mathit{trial}^{v^-}} : \left(-\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right) \, d\Omega
  + \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{v^-}} \cdot \left(-\frac{1}{\rho(\vec{x})} \vec{\lambda} \right) \, d\Gamma \\
  - \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{\partial^2 \vec{d}(\vec{x}, t)}{\partial t^2} \, d\Gamma.\end{gathered}$$
We rewrite the integrals over the domain involving the degrees of freedom
adjacent to the fault as integrals over the positive and negative sides of the
fault. These are implemented as integrals over the faces of cells adjacent to
the fault; they involve quantities, such as density, that are defined only
within the domain cells. After collecting and rearranging terms, we have
$$\begin{gathered}
  \label{eqn:elasticity:prescribed:slip:dynamic:DAE:weak:form}
  \int_{\Gamma_{f^+}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{1}{\rho(\vec{x})} \left(
    \vec{\lambda} - \vec{f}(\vec{x},t) + \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \bm{\sigma}(\vec{u}) \right)
  + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : \left(+\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right) \, d\Gamma \\
  + \int_{\Gamma_{f^-}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{1}{\rho(\vec{x})} \left(
    \vec{\lambda} + \vec{f}(\vec{x},t) - \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \bm{\sigma}(\vec{u}) \right)
  + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : \left(-\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right)  \, d\Gamma \\
  + \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot \frac{\partial^2 \vec{d}(\vec{x}, t)}{\partial t^2} \, d\Gamma
  = 0.\end{gathered}$$

### Residual Pointwise Functions

Combining the explicit parts of the weak form in
equations&nbsp;[\[eqn:displacement:velocity:prescribed:slip:weak:form\]][2]
and [\[eqn:elasticity:prescribed:slip:dynamic:weak:form\]][3] with the
implicit part of the weak form in
equation&nbsp;[\[eqn:elasticity:prescribed:slip:dynamic:DAE:weak:form\]][4]
and identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have $$\begin{gathered}
  % Fu
  F^u(t,s,\dot{s}) = \frac{\partial \vec{u}}{\partial t} \\
  % Fv
  F^v(t,s,\dot{s}) = \frac{\partial \vec{v}}{\partial t} \\
  % Fl
  \begin{multlined}
    F^\lambda(t,s,\dot{s}) =   \int_{\Gamma_{f^+}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\frac{1}{\rho(\vec{x})} \left(
    \vec{\lambda} - \vec{f}(\vec{x},t) + \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \bm{\sigma}(\vec{u}) \right)}_{\color{blue}\mathclap{f^\lambda_0}}}
  + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : %
  {\color{blue}%
  \underbrace{\color{black}\left(+\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u})\right)}_{\color{blue}\mathclap{f^\lambda_1}}} \, d\Gamma \\
  + \int_{\Gamma_{f^-}} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\frac{1}{\rho(\vec{x})} \left(
    \vec{\lambda} + \vec{f}(\vec{x},t) - \frac{\nabla\rho(\vec{x})}{\rho(\vec{x})} \cdot \bm{\sigma}(\vec{u}) \right)}_{\color{blue}\mathclap{f^\lambda_0}}}
  + \nabla {\vec{\psi}_\mathit{trial}^{\lambda}} : %
  {\color{blue}%
  \underbrace{\color{black}\left(-\frac{1}{\rho(\vec{x})} \bm{\sigma}(\vec{u}) \right)}_{\color{blue}\mathclap{f^\lambda_1}}}  \, d\Gamma \\
  + \int_{\Gamma_f} {\vec{\psi}_\mathit{trial}^{\lambda}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\frac{\partial^2 \vec{d}(\vec{x}, t)}{\partial t^2}}_{\color{blue}\mathclap{f^\lambda_0}}} \, d\Gamma
  \end{multlined}\\
  % Gu
  G^u(t,s) = \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\vec{v}}_{\color{blue}\mathclap{\vec{g}^u_0}}} \, d\Omega, \\
  % Gv
  G^v(t,s) =  \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}\mathclap{\vec{g}^v_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{v}} : %
  {\color{blue}%
  \underbrace{\color{black}-\bm{\sigma}(\vec{u})}_{\color{blue}\mathclap{\bm{g^v_1}}}} \, d\Omega
  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}\mathclap{\vec{g}^v_0}}} \, d\Gamma,
  + \int_{\Gamma_{f}} {\vec{\psi}_\mathit{trial}^{v^+}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\left(-\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}\mathclap{\vec{g}^v_0}}}
             + {\vec{\psi}_\mathit{trial}^{v^-}} \cdot %
  {\color{blue}%
  \underbrace{\color{black}\left(+\vec{\lambda}(\vec{x},t)\right)}_{\color{blue}\mathclap{\vec{g}^v_0}}} \, d\Gamma, \\
  % Gl
  G^l(t,s) = 0\end{gathered}$$

### Jacobian Pointwise Functions

For the explicit part we have pointwise functions for computing the lumped LHS
Jacobian. These are exactly the same pointwise functions as in the dynamic
case without a fault, $$\begin{aligned}
  % J_F uu
  J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_\mathit{tshift} \frac{\partial F^u}{\partial \dot{u}} =
             \int_\Omega {\psi_\mathit{trial}^{u}}_i %
  {\color{blue}%
  \underbrace{\color{black}s_\mathit{tshift} \delta_{ij}}_{\color{blue}\mathclap{J^{uu}_{f0}}}} {\psi_\mathit{basis}^{u}}_j  \, d\Omega, \\
  % J_F vv
  J_F^{vv} &= \frac{\partial F^v}{\partial v} + s_\mathit{tshift} \frac{\partial F^v}{\partial \dot{v}} =
             \int_\Omega {\psi_\mathit{trial}^{v}}_i %
  {\color{blue}%
  \underbrace{\color{black}\rho(\vec{x}) s_\mathit{tshift} \delta_{ij}}_{\color{blue}\mathclap{J ^{vv}_{f0}}}} {\psi_\mathit{basis}^{v}}_j \, d\Omega\end{aligned}$$
For the implicit part, we have pointwise functions for the LHS Jacobians
associated with the prescribed slip, $$\begin{gathered}
  \begin{multlined}
  % J_F lu
  J_F^{\lambda u} = \frac{\partial F^\lambda}{\partial u} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{u}} = \\
                    \int_{\Gamma_{f^+}} {\psi_\mathit{trial}^{\lambda}}_i %
  {\color{blue}%
  \underbrace{\color{black}\frac{\rho_{,j}(\vec{x})}{\rho^2(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{j,l}}_{\color{blue}\mathclap{J^{\lambda u}_{fX}}}}
                    + {\psi_\mathit{trial}^{\lambda}}_{i,k} %
  {\color{blue}%
  \underbrace{\color{black}+\frac{1}{\rho(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{k,l}}_{\color{blue}\mathclap{J^{\lambda u}_{f3}}}} \, d\Gamma \\
                    +\int_{\Gamma_{f^-}} {\psi_\mathit{trial}^{\lambda}}_i %
  {\color{blue}%
  \underbrace{\color{black}-\frac{\rho_{,j}(\vec{x})}{\rho^2(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{j,l}}_{\color{blue}\mathclap{J^{\lambda u}_{fX}}}}
                    + {\psi_\mathit{trial}^{\lambda}}_{i,k} %
  {\color{blue}%
  \underbrace{\color{black}-\frac{1}{\rho(\vec{x})} C_{ikjl} {\psi_\mathit{basis}^{u}}_{k,l}}_{\color{blue}\mathclap{J^{\lambda u}_{f3}}}} \, d\Gamma
                  \end{multlined} \\
  % J_F ll
  J_F^{\lambda \lambda} = \frac{\partial F^\lambda}{\partial \lambda} + s_\mathit{tshift} \frac{\partial F^\lambda}{\partial \dot{\lambda}} =
             \int_{\Gamma_{f^+}} {\psi_\mathit{trial}^{\lambda}}_i %
  {\color{blue}%
  \underbrace{\color{black}\frac{1}{\rho(\vec{x})} \delta_{ij}}_{\color{blue}\mathclap{J^{\lambda\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j \, d\Gamma
            + \int_{\Gamma_{f^-}} {\psi_\mathit{trial}^{\lambda}}_i %
  {\color{blue}%
  \underbrace{\color{black}\frac{1}{\rho(\vec{x})} \delta_{ij}}_{\color{blue}\mathclap{J^{\lambda\lambda}_{f0}}}} {\psi_\mathit{basis}^{\lambda}}_j \, d\Gamma\end{gathered}$$

  [1]: #eqn:prescribed:slip:DAE:weak:form
  [2]: #eqn:displacement:velocity:prescribed:slip:weak:form
  [3]: #eqn:elasticity:prescribed:slip:dynamic:weak:form
  [4]: #eqn:elasticity:prescribed:slip:dynamic:DAE:weak:form
