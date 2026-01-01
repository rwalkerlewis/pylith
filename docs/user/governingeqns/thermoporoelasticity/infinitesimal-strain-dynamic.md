# Dynamic

For the dynamic case we retain the inertial term and introduce velocity as an additional field. The solution vector becomes $(\vec{v}, \vec{u}, p, T, \epsilon_v)$.

The strong form is:
%
\begin{gather}
  \rho_b \frac{\partial \vec{v}}{\partial t} - \vec{f}(t) - \nabla \cdot \boldsymbol{\sigma}(\vec{u},p,T) = \vec{0} \text{ in } \Omega, \\
  \frac{\partial \vec{u}}{\partial t} - \vec{v} = \vec{0} \text{ in } \Omega, \\
  \frac{\partial \zeta(\vec{u},p,T)}{\partial t} - \gamma(\vec{x},t) + \nabla \cdot \vec{q}(p,T) = 0 \text{ in } \Omega, \\
  \rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T = 0 \text{ in } \Omega, \\
  \nabla \cdot \vec{u} - \epsilon_{v} = 0 \text{ in } \Omega.
\end{gather}
%

Applying the weak form and divergence theorem:
%
\begin{gather}
  \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \left(\rho_b \frac{\partial \vec{v}}{\partial t} - \vec{f}(\vec{x},t)\right) + \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{\sigma}(\vec{u},p,T) \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot \vec{\tau}(\vec{x},t) \, d\Gamma = 0, \\
  \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot \left(\frac{\partial \vec{u}}{\partial t} - \vec{v}\right) \, d\Omega = 0, \\
  \int_\Omega {\psi_\mathit{trial}^{p}} \left(\frac{\partial \zeta(\vec{u},p,T)}{\partial t} - \gamma(\vec{x},t)\right) + \nabla {\psi_\mathit{trial}^{p}} \cdot (-\vec{q}(p,T)) \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} q_0(\vec{x},t) \, d\Gamma = 0, \\
  \int_\Omega {\psi_\mathit{trial}^{T}} \left(\rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T\right) + \nabla {\psi_\mathit{trial}^{T}} \cdot (-k \nabla T) \, d\Omega + \int_{\Gamma_h} {\psi_\mathit{trial}^{T}} h_0(\vec{x},t) \, d\Gamma = 0, \\
  \int_\Omega {\psi_\mathit{trial}^{\epsilon_v}} (\nabla \cdot \vec{u} - \epsilon_v) \, d\Omega = 0.
\end{gather}
%

## Residual Pointwise Functions

\begin{align}
  F^v(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}\underbrace{\color{black}\left(\rho_b \frac{\partial \vec{v}}{\partial t} - \vec{f}(\vec{x},t)\right)}_{\color{blue}{\vec{f}^v_0}}} + \nabla {\vec{\psi}_\mathit{trial}^{v}} : {\color{blue}\underbrace{\color{black}-\boldsymbol{\sigma}(\vec{u},p,T)}_{\color{blue}{\boldsymbol{f}^v_1}}} \, d\Omega + \int_{\Gamma_\tau} {\vec{\psi}_\mathit{trial}^{v}} \cdot {\color{blue}\underbrace{\color{black}\vec{\tau}(\vec{x},t)}_{\color{blue}{\vec{f}^v_0}}} \, d\Gamma, \\
  G^v(t,s) &= 0, \\
  F^u(t,s,\dot{s}) &= \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\color{blue}\underbrace{\color{black}\left(\frac{\partial \vec{u}}{\partial t} - \vec{v}\right)}_{\color{blue}{\vec{f}^u_0}}} \, d\Omega, \\
  G^u(t,s) &= 0, \\
  F^p(t,s,\dot{s}) &= \int_\Omega {\psi_\mathit{trial}^{p}} {\color{blue}\underbrace{\color{black}\left(\frac{\partial \zeta}{\partial t} - \gamma\right)}_{\color{blue}{f^p_0}}} + \nabla {\psi_\mathit{trial}^{p}} \cdot {\color{blue}\underbrace{\color{black}(-\vec{q})}_{\color{blue}{\vec{f}^p_1}}} \, d\Omega + \int_{\Gamma_q} {\psi_\mathit{trial}^{p}} {\color{blue}\underbrace{\color{black}q_0}_{\color{blue}{f^p_0}}} \, d\Gamma, \\
  G^p(t,s) &= 0, \\
  F^T(t,s,\dot{s}) &= \int_\Omega {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}\left(\rho_b c_b \frac{\partial T}{\partial t} - Q + \rho_f c_f \vec{q} \cdot \nabla T\right)}_{\color{blue}{f^T_0}}} + \nabla {\psi_\mathit{trial}^{T}} \cdot {\color{blue}\underbrace{\color{black}(-k \nabla T)}_{\color{blue}{\vec{f}^T_1}}} \, d\Omega + \int_{\Gamma_h} {\psi_\mathit{trial}^{T}} {\color{blue}\underbrace{\color{black}h_0}_{\color{blue}{f^T_0}}} \, d\Gamma, \\
  G^T(t,s) &= 0, \\
  F^{\epsilon_v}(t,s,\dot{s}) &= \int_\Omega {\psi_\mathit{trial}^{\epsilon_v}} {\color{blue}\underbrace{\color{black}(\nabla \cdot \vec{u} - \epsilon_v)}_{\color{blue}{f^{\epsilon_v}_0}}} \, d\Omega, \\
  G^{\epsilon_v}(t,s) &= 0.
\end{align}

## Jacobian Pointwise Functions

The dynamic formulation introduces additional coupling through velocity:
%
\begin{align}
  J_F^{vv} &= s_{tshift} \int_\Omega {\vec{\psi}_\mathit{trial}^{v}} \cdot \rho_b {\vec{\psi}_\mathit{basis}^{v}} \, d\Omega, \\
  J_F^{vu} &= \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\boldsymbol{C} : \frac{1}{2}(\nabla + \nabla^T) {\vec{\psi}_\mathit{basis}^{u}} \, d\Omega, \\
  J_F^{vp} &= \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{v}} : \alpha \mathbf{I} {\psi_\mathit{basis}^{p}} \, d\Omega, \\
  J_F^{vT} &= \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{v}} : \boldsymbol{C} : \alpha_T \mathbf{I} {\psi_\mathit{basis}^{T}} \, d\Omega, \\
  J_F^{v\epsilon_v} &= \int_\Omega \nabla {\vec{\psi}_\mathit{trial}^{v}} : -\lambda \mathbf{I} {\psi_\mathit{basis}^{\epsilon_v}} \, d\Omega, \\
  J_F^{uu} &= s_{tshift} \int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\vec{\psi}_\mathit{basis}^{u}} \, d\Omega, \\
  J_F^{uv} &= -\int_\Omega {\vec{\psi}_\mathit{trial}^{u}} \cdot {\vec{\psi}_\mathit{basis}^{v}} \, d\Omega.
\end{align}
%
The Jacobian terms for $p$, $T$, and $\epsilon_v$ follow the same pattern as in the quasistatic case.

````