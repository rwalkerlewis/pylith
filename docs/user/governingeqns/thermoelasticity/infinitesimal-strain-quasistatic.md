# Quasistatic

If we neglect the inertial term ($\rho \frac{\partial \vec{v}}{\partial t} \approx \vec{0}$), then time dependence only arises from the temperature evolution and boundary conditions.
Our solution vector is the displacement vector and temperature, and the thermoelasticity equations reduce to
%
```{math}
:label: eqn:thermoelasticity:strong:form:quasistatic
\begin{gathered}
\vec{f}(\vec{x},t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u},T) = \vec{0} \text{ in }\Omega, \\
%
\rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) = 0 \text{ in }\Omega, \\
%
\boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on }\Gamma_\tau, \\
%
\vec{u} = \vec{u}_0(\vec{x},t) \text{ on }\Gamma_u, \\
%
T = T_0(\vec{x},t) \text{ on }\Gamma_T, \\
%
-k \nabla T \cdot \vec{n} = q_0(\vec{x},t) \text{ on }\Gamma_q.
\end{gathered}
```
%
Because we will use implicit time stepping, we place all of the terms in the thermoelasticity equations on the LHS.
We create the weak form by taking the dot product with the trial functions ${\vec{\psi}_\mathit{trial}^{u}}$ and ${\psi_\mathit{trial}^{T}}$ and integrating over the domain:
%
\begin{gather}
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \vec{f}(t) + \boldsymbol{\nabla}\cdot \boldsymbol{\sigma} (\vec{u},T) \right) \, d\Omega = 0, \\
\int_\Omega {\psi_\mathit{trial}^{T}} \left( \rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) \right) \, d\Omega = 0.
\end{gather}
%
Using the divergence theorem and incorporating the Neumann boundary conditions, we have
%
\begin{gather}
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{\sigma}(\vec{u},T) \, d\Omega  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = 0, \\
\int_\Omega {\psi_\mathit{trial}^{T}} \left( \rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) \right) + \nabla {\psi_\mathit{trial}^{T}} \cdot (-k \nabla T) \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{T}} q_0(\vec{x},t) \, d\Gamma = 0.
\end{gather}

## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have
%
\begin{align}
% Fu
F^u(t,s,\dot{s}) &=  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot{\color{blue}\underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} :{\color{blue} \underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},T)}_{\color{blue}{\boldsymbol{f^u_1}}}} \, d\Omega  + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} \, d\Gamma, \\
% Gu
G^u(t,s) &= 0, \\
% FT
F^T(t,s,\dot{s}) &= \int_\Omega {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}\left( \rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) \right)}_{\color{blue}{f^T_0}}} + \nabla {\psi_\mathit{trial}^{T}} \cdot {\color{blue}\underbrace{\color{black}(-k \nabla T)}_{\color{blue}{\vec{f}^T_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}q_0(\vec{x},t)}_{\color{blue}{f^T_0}}} \, d\Gamma, \\
% GT
G^T(t,s) &= 0.
\end{align}
%
Note that we have multiple $\vec{f}_0$ functions, each associated with a trial function and an integral over a different domain or boundary.
Each material and boundary condition (except Dirichlet) contribute pointwise functions.
With $G=0$ it is clear that we have a formulation that will use implicit time stepping algorithms.

## Jacobian Pointwise Functions

We have Jacobians for the LHS for both displacement and temperature:
%
```{math}
:label: eqn:thermoelasticity:quasistatic:jacobian:pointwise
\begin{aligned}
J_F^{uu} &= \frac{\partial F^u}{\partial u} = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial u}(-\boldsymbol{\sigma}) \, d\Omega  = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{C} : \frac{1}{2}(\nabla + \nabla^T){\vec{\psi}_\mathit{basis}^{u}}\, d\Omega, \\
J_F^{uT} &= \frac{\partial F^u}{\partial T} = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : \frac{\partial}{\partial T}(-\boldsymbol{\sigma}) \, d\Omega = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{u}} : \boldsymbol{C} : \alpha \mathbf{I} {\psi_\mathit{basis}^{T}} \, d\Omega, \\
J_F^{Tu} &= \frac{\partial F^T}{\partial u} = 0, \\
J_F^{TT} &= \frac{\partial F^T}{\partial T} + s_{tshift} \frac{\partial F^T}{\partial \dot{T}} = \int_\Omega \nabla {\psi_\mathit{trial}^{T}} \cdot \frac{\partial}{\partial T}(-k \nabla T) \, d\Omega + s_{tshift} \int_\Omega {\psi_\mathit{trial}^{T}} \rho c {\psi_\mathit{basis}^{T}} \, d\Omega.
\end{aligned}
```