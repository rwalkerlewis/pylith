# Quasistatic

For ease of solution in the quasistatic case, we introduce volumetric strain ($\epsilon_v$) as a third variable in addition to displacement and pressure, and we include temperature as the fourth field.
The strong form of the problem may be expressed as
%
\begin{gather}
% Solution
  \vec{s}^{T} = \left(\vec{u} \quad p \quad T \quad \epsilon_v\right), \\
% Elasticity
  \vec{f}(t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u},p,T) = \vec{0} \text{ in } \Omega, \\
% Pressure
  \frac{\partial \zeta(\vec{u},p,T)}{\partial t} - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p,T) = 0 \text{ in } \Omega, \\
% Temperature
  \rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T = 0 \text{ in } \Omega, \\
% Vol. Strain
  \nabla \cdot \vec{u} - \epsilon_{v} = 0 \text{ in } \Omega, \\
% Neumann traction
  \boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on } \Gamma_{\tau}, \\
% Neumann flow
  \vec{q} \cdot \vec{n} = q_0(\vec{x}, t) \text{ on } \Gamma_{q}, \\
% Neumann heat flux
  -k \nabla T \cdot \vec{n} = h_0(\vec{x},t) \text{ on } \Gamma_h, \\
% Dirichlet displacement
  \vec{u} = \vec{u}_0(\vec{x}, t) \text{ on } \Gamma_{u}, \\
% Dirichlet pressure
  p = p_0(\vec{x},t) \text{ on } \Gamma_{p}, \text{ and } \\
% Dirichlet temperature
  T = T_0(\vec{x},t) \text{ on } \Gamma_T.
\end{gather}
%
We place all terms for the elasticity, pressure, temperature, and volumetric strain equations on the left-hand-side, consistent with PETSc TS implicit time stepping.

We create the weak form by taking the dot product with the trial functions ${\vec{\psi}_\mathit{trial}^{u}}$, ${\psi_\mathit{trial}^{p}}$, ${\psi_\mathit{trial}^{T}}$, and ${\psi_\mathit{trial}^{\epsilon_{v}}}$ and integrating over the domain:
%
\begin{gather}
% Weak conservation of momentum
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left( \vec{f}(\vec{x},t) + \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u},p,T) \right) \, d\Omega = 0, \\
% Weak conservation of mass
  \int_\Omega  {\psi_\mathit{trial}^{p}} \left( \frac{\partial \zeta(\vec{u},p,T)}{\partial t} - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p,T)\right) \, d\Omega = 0,\\
% Weak energy balance
  \int_\Omega {\psi_\mathit{trial}^{T}} \left( \rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T \right) \, d\Omega = 0, \\
% Weak vol. strain
  \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}}\cdot \left( \nabla \cdot \vec{u} - \epsilon_v \right) \, d\Omega = 0.
\end{gather}
%
Applying the divergence theorem to the first three equations and incorporating the Neumann boundary conditions yields
%
\begin{gather}
% Weak conservation of momentum
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{f}(\vec{x},t) + \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{\sigma}(\vec{u},p,T) \,
  d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = 0, \\
% Weak conservation of mass
  \int_\Omega  {\psi_\mathit{trial}^{p}} \left( \frac{\partial \zeta(\vec{u},p,T)}{\partial t} - \gamma(\vec{x},t)\right)
  + \nabla {\psi_\mathit{trial}^{p}} \cdot \left(-\vec{q}(p,T)\right) \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} q_0(\vec{x},t) \, d\Gamma = 0, \\
% Weak energy balance
  \int_\Omega {\psi_\mathit{trial}^{T}} \left( \rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T \right) + \nabla {\psi_\mathit{trial}^{T}} \cdot (-k \nabla T) \, d\Omega + \int_{\Gamma_h} {\psi_\mathit{trial}^{T}} h_0(\vec{x},t) \, d\Gamma = 0, \\
% Weak vol. strain
  \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot \left(\nabla \cdot \vec{u} - \epsilon_{v} \right) d\Omega = 0.
\end{gather}
%

## Residual Pointwise Functions

Identifying $F(t,s,\dot{s})$ and $G(t,s)$ we have
%
\begin{align}
  % LHS displacement
  F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{f}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{u}} : {\color{blue}  \underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},p,T)}_{\color{blue}{\boldsymbol{f}^u_1}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}  \underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^u_0}}} \, d\Gamma, \\
% RHS displacement
  G^u(t,s) &= 0, \\
% LHS fluid pressure
  F^p(t,s,\dot{s}) &= \int_\Omega  {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\left( \color{black}\frac{\partial \zeta(\vec{u},p,T)}{\partial t} - \gamma(\vec{x},t)\right)}_{\color{blue}{f^p_0}}} + \nabla {\psi_\mathit{trial}^{p}} \cdot {\color{blue}  \underbrace{\color{black}-\vec{q}(p,T)}_{\color{blue}{\vec{f}^p_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} {\color{blue} \underbrace{\color{black}q_0(\vec{x},t)}_{\color{blue}{f^p_0}}} \, d\Gamma, \\
% RHS fluid pressure
  G^p(t,s) &= 0, \\
% LHS temperature
  F^T(t,s,\dot{s}) &= \int_\Omega {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}\left( \rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T \right)}_{\color{blue}{f^T_0}}} + \nabla {\psi_\mathit{trial}^{T}} \cdot {\color{blue}\underbrace{\color{black}(-k \nabla T)}_{\color{blue}{\vec{f}^T_1}}} \, d\Omega + \int_{\Gamma_h} {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}h_0(\vec{x},t)}_{\color{blue}{f^T_0}}} \, d\Gamma, \\
% RHS temperature
  G^T(t,s) &= 0, \\
% LHS trace strain
  F^{\epsilon_{v}}(t,s,\dot{s}) &= \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \cdot {\color{blue}
  \underbrace{\color{black}\left(\nabla \cdot \vec{u} - \epsilon_{v} \right)}_{\color{blue}{f^{\epsilon_{v}}_{0}}}} \, d\Omega, \\
% RHS trace strain
  G^{\epsilon_v}(t,s) &= 0.
\end{align}
%

## Jacobian Pointwise Functions

Four fields yields potentially 16 Jacobian pointwise functions for the LHS. The key thermal coupling terms are $J_F^{uT}$, $J_F^{pT}$, and $J_F^{Tp}$:
%
\begin{align}
  J_F^{uu} &= \frac{\partial F^u}{\partial u} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\boldsymbol{C}: \frac{1}{2} (\nabla + \nabla^T) {\vec{\psi}_\mathit{basis}^{u}} \ d\Omega, \\
  J_F^{up} &= \frac{\partial F^u}{\partial p} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \alpha \mathbf{I} {\psi_\mathit{basis}^{p}} \ d\Omega, \\
  J_F^{uT} &= \frac{\partial F^u}{\partial T} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : \boldsymbol{C} : \alpha_T \mathbf{I} {\psi_\mathit{basis}^{T}} \ d\Omega, \\
  J_F^{u \epsilon_{v}} &= \frac{\partial F^u}{\partial \epsilon_{v}} = \int_{\Omega} \nabla {\vec{\psi}_\mathit{trial}^{u}} : -\lambda \mathbf{I} {\psi_\mathit{basis}^{\epsilon_{v}}} d\Omega, \\
  J_F^{pu} &= \frac{\partial F^p}{\partial u} = 0, \\
  J_F^{pp} &= \frac{\partial F^p}{\partial p} + s_{tshift} \frac{\partial F^p}{\partial \dot{p}} = \int_{\Omega} \nabla {\psi_\mathit{trial}^{p}} \cdot \frac{\boldsymbol{k}}{\mu_{f}} \nabla {\psi_\mathit{basis}^{p}} \ d\Omega + \int_{\Omega} {\psi_\mathit{trial}^{p}} \left(s_{tshift} \frac{1}{M}\right) {\psi_\mathit{basis}^{p}} \ d\Omega, \\
  J_F^{pT} &= \frac{\partial F^p}{\partial T} + s_{tshift} \frac{\partial F^p}{\partial \dot{T}} = s_{tshift} \int_{\Omega} {\psi_\mathit{trial}^{p}} \beta {\psi_\mathit{basis}^{T}} \ d\Omega, \\
  J_F^{p\epsilon_{v}} &= \frac{\partial F^p}{\partial \epsilon_{v}} + s_{tshift} \frac{\partial F^p}{\partial \dot{\epsilon_{v}}} = \int_{\Omega} {\psi_\mathit{trial}^{p}} \left(s_{tshift} \alpha \right) {\psi_\mathit{basis}^{\epsilon_{v}}} \ d\Omega, \\
  J_F^{Tu} &= \frac{\partial F^T}{\partial u} = 0, \\
  J_F^{Tp} &= \frac{\partial F^T}{\partial p} = \int_{\Omega} {\psi_\mathit{trial}^{T}} \rho_f c_f \left(-\frac{\boldsymbol{k}}{\mu_f} \nabla {\psi_\mathit{basis}^{p}}\right) \cdot \nabla T \ d\Omega, \\
  J_F^{TT} &= \frac{\partial F^T}{\partial T} + s_{tshift} \frac{\partial F^T}{\partial \dot{T}} = \int_{\Omega} \nabla {\psi_\mathit{trial}^{T}} \cdot k \nabla {\psi_\mathit{basis}^{T}} \ d\Omega + s_{tshift} \int_{\Omega} {\psi_\mathit{trial}^{T}} \rho_b c_b {\psi_\mathit{basis}^{T}} \ d\Omega, \\
  J_F^{T\epsilon_v} &= \frac{\partial F^T}{\partial \epsilon_v} = 0, \\
  J_F^{\epsilon_{v}u} &= \frac{\partial F^{\epsilon_{v}}}{\partial u} = \int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} \nabla \cdot {\vec{\psi}_\mathit{basis}^{u}} \ d\Omega, \\
  J_F^{\epsilon_{v}p} &= \frac{\partial F^{\epsilon_{v}}}{\partial p} = 0, \\
  J_F^{\epsilon_{v}T} &= \frac{\partial F^{\epsilon_{v}}}{\partial T} = 0, \\
  J_F^{\epsilon_{v}\epsilon_{v}} &= \frac{\partial F^{\epsilon_{v}}}{\partial \epsilon_{v}} = -\int_{\Omega} {\psi_\mathit{trial}^{\epsilon_{v}}} {\psi_\mathit{basis}^{\epsilon_{v}}} \ d\Omega.
\end{align}
%

````