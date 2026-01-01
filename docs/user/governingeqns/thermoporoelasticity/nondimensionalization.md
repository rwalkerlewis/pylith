# Nondimensionalization

Starting with the thermoporoelasticity equations:
%
\begin{gather}
\rho_b(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u},p,T)=\vec{0}\text{ in }\Omega,\\
\frac{\partial \zeta(\vec{u},p,T)}{\partial t} + \nabla \cdot \vec{q}(p,T) - \gamma(\vec{x},t) = 0 \text{ in }\Omega,\\
\rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T = 0 \text{ in }\Omega,
\end{gather}
%
we define nondimensional values:
%
\begin{align}
\vec{x}^* &= \frac{\vec{x}}{x_o}, \\
\vec{u}^* &= \frac{\vec{u}}{u_o}, \\
p^* &= \frac{p}{p_o}, \\
T^* &= \frac{T - T_0}{\Delta T}, \\
\rho^* &= \frac{\rho}{\rho_o}, \\
\vec{f}^* &= \frac{\vec{f}}{f_o}, \\
\gamma^* &= \frac{\gamma}{\gamma_o}, \\
Q^* &= \frac{Q}{Q_o}, \\
\boldsymbol{\sigma}^* &= \frac{\boldsymbol{\sigma}}{\sigma_o}.
\end{align}
%
We also recognize that $\boldsymbol{\nabla}^* = x_o \boldsymbol{\nabla}$.

Substituting into the equations and grouping terms yields nondimensional groups. For the momentum equation:
%
\begin{equation}
\left( \rho_o \frac{u_o}{t_o^2}\frac{x_o}{\sigma_o}\right) \rho^*(\vec{x}^*) \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - \left(f_o \frac{x_o}{\sigma_o}\right) \vec{f}^*(\vec{x}^*,t^*) -  \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^*(p^*,T^*) = \vec{0}.
\end{equation}
%

For the mass balance equation:
%
\begin{equation}
\left(\frac{x_o^2 \mu_f}{\boldsymbol{k} p_o t_o}\right) \frac{\partial \zeta^*}{\partial t^*} + \nabla^* \cdot \vec{q}^* - \left(\gamma_o \frac{x_o^2 \mu_f}{\boldsymbol{k} p_o}\right) \gamma^* = 0.
\end{equation}
%

For the energy balance equation:
%
\begin{equation}
\left( \rho_o c_o \frac{x_o^2}{k_o t_o} \right) \rho^* c^* \frac{\partial T^*}{\partial t^*} - \left( Q_o \frac{x_o^2}{k_o \Delta T} \right) Q^* - \nabla^* \cdot (k^* \nabla^* T^*) + \left(\rho_f c_f \frac{\boldsymbol{k} p_o x_o}{k_o \mu_f \Delta T}\right) \vec{q}^* \cdot \nabla^* T^* = 0.
\end{equation}

## Nondimensional Groups

### Poroelastic Biot Number

The ratio of fluid diffusion time to problem time scale:
\begin{equation}
\Pi_\mathit{poro} = \frac{x_o^2 \mu_f}{\boldsymbol{k} M t_o}.
\end{equation}
If $\Pi_\mathit{poro} \ll 1$, fluid pressure equilibrates quickly (drained limit). If $\Pi_\mathit{poro} \gg 1$, the process is undrained on the loading time scale.

### Thermal Biot Number

The ratio of thermal diffusion time to problem time scale:
\begin{equation}
\Pi_\mathit{thermal} = \frac{\rho_o c_o x_o^2}{k_o t_o}.
\end{equation}
If $\Pi_\mathit{thermal} \ll 1$, thermal diffusion is fast and temperature equilibrates quickly. If $\Pi_\mathit{thermal} \gg 1$, the process is nearly adiabatic.

### Thermal-Poroelastic Coupling Number

The ratio of advective heat transport to conductive heat transport:
\begin{equation}
\Pi_\mathit{advection} = \frac{\rho_f c_f \boldsymbol{k} p_o x_o}{k_o \mu_f \Delta T}.
\end{equation}
If $\Pi_\mathit{advection} \gg 1$, heat transport by fluid flow dominates over conduction. This is important in highly permeable media with significant fluid flow.

### Thermal-Mechanical Coupling

From the constitutive relation, thermal effects on stress scale as:
\begin{equation}
\frac{\mu_o \alpha_T \Delta T}{p_o \alpha}.
\end{equation}
If this ratio is O(1), thermal stresses are comparable to poroelastic stresses.

### Inertia

The scale of the inertial term is $\Pi_\mathit{inertia} = \frac{\rho_o x_o^2}{\mu_o t_o^2}$.
If this term is O(1), include inertia and solve the dynamic formulation.
If very small, solve the quasistatic formulation.

````