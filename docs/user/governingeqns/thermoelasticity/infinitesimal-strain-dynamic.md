# Dynamic

For the dynamic case we retain the inertial term and solve
%
```{math}
:label: eqn:thermoelasticity:strong:form:dynamic
\begin{gathered}
\rho \frac{\partial^2\vec{u}}{\partial t^2} - \vec{f}(\vec{x},t) - \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u},T) = \vec{0} \text{ in }\Omega, \\
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
We rewrite the equation as a system of first order equations in time by introducing the velocity field, $\vec{v}$.
%
\begin{gather}
\rho \frac{\partial \vec{v}}{\partial t} - \vec{f}(\vec{x},t) - \boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u},T) = \vec{0} \text{ in }\Omega, \\
\frac{\partial \vec{u}}{\partial t} - \vec{v} = \vec{0} \text{ in }\Omega.
\end{gather}
%
We form the weak form in the same way as we did for the quasistatic case.
%
\begin{gather}
\int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( \rho \frac{\partial \vec{v}}{\partial t} - \vec{f}(\vec{x},t) - \boldsymbol{\nabla} \cdot \boldsymbol{\sigma}(\vec{u},T) \right) \, d\Omega = 0, \\
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \frac{\partial \vec{u}}{\partial t} - \vec{v} \right) \, d\Omega = 0, \\
\int_\Omega {\psi_\mathit{trial}^{T}} \left( \rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) \right) \, d\Omega = 0.
\end{gather}
%
Applying the divergence theorem to the first and third equations and incorporating the Neumann boundary conditions, we have
%
\begin{gather}
\int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \left( \rho \frac{\partial \vec{v}}{\partial t} - \vec{f}(\vec{x},t) \right) + \nabla{\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{\sigma}(\vec{u},T) \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = 0, \\
\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \frac{\partial \vec{u}}{\partial t} - \vec{v} \right) \, d\Omega = 0, \\
\int_\Omega {\psi_\mathit{trial}^{T}} \left( \rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) \right) + \nabla {\psi_\mathit{trial}^{T}} \cdot (-k \nabla T) \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{T}} q_0(\vec{x},t) \, d\Gamma = 0.
\end{gather}
%

## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$, we have
%
\begin{align}
% Fv
F^v(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}\underbrace{\color{black}\left( \rho \frac{\partial \vec{v}}{\partial t} - \vec{f}(\vec{x},t) \right)}_{\color{blue}{\vec{f}^v_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{v}} : {\color{blue}\underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},T)}_{\color{blue}{\boldsymbol{f}^v_1}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}\underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^v_0}}} \, d\Gamma, \\
% Gv
G^v(t,s) &= 0, \\
% Fu
F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}\underbrace{\color{black}\left( \frac{\partial \vec{u}}{\partial t} - \vec{v} \right)}_{\color{blue}{\vec{f}^u_0}}} \, d\Omega, \\
% Gu
G^u(t,s) &= 0, \\
% FT
F^T(t,s,\dot{s}) &= \int_\Omega {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}\left( \rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) \right)}_{\color{blue}{f^T_0}}} + \nabla {\psi_\mathit{trial}^{T}} \cdot {\color{blue}\underbrace{\color{black}(-k \nabla T)}_{\color{blue}{\vec{f}^T_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}q_0(\vec{x},t)}_{\color{blue}{f^T_0}}} \, d\Gamma, \\
% GT
G^T(t,s) &= 0.
\end{align}

## Jacobian Pointwise Functions

For the dynamic case, Jacobians couple velocity, displacement, and temperature:
%
\begin{align}
J_F^{vv} &= \frac{\partial F^v}{\partial v} + s_{tshift} \frac{\partial F^v}{\partial \dot{v}} = s_{tshift} \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \rho {\vec{\psi}_\mathit{basis}^{v}} \, d\Omega, \\
J_F^{vu} &= \frac{\partial F^v}{\partial u} = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{C} : \frac{1}{2}(\nabla + \nabla^T) {\vec{\psi}_\mathit{basis}^{u}} \, d\Omega, \\
J_F^{vT} &= \frac{\partial F^v}{\partial T} = \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{v}} : \boldsymbol{C} : \alpha \mathbf{I} {\psi_\mathit{basis}^{T}} \, d\Omega, \\
J_F^{uu} &= \frac{\partial F^u}{\partial u} + s_{tshift} \frac{\partial F^u}{\partial \dot{u}} = s_{tshift} \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\vec{\psi}_\mathit{basis}^{u}} \, d\Omega, \\
J_F^{uv} &= \frac{\partial F^u}{\partial v} = -\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\vec{\psi}_\mathit{basis}^{v}} \, d\Omega, \\
J_F^{TT} &= \frac{\partial F^T}{\partial T} + s_{tshift} \frac{\partial F^T}{\partial \dot{T}} = \int_\Omega \nabla {\psi_\mathit{trial}^{T}} \cdot \frac{\partial}{\partial T}(-k \nabla T) \, d\Omega + s_{tshift} \int_\Omega {\psi_\mathit{trial}^{T}} \rho c {\psi_\mathit{basis}^{T}} \, d\Omega.
\end{align}