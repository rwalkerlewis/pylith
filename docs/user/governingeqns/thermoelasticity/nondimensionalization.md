````markdown
# Nondimensionalization

Starting with the thermoelasticity equation,
%
\begin{gather}
\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u},T)=\vec{0}\text{ in }\Omega,\\
\rho(\vec{x}) c \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) = 0 \text{ in }\Omega,\\
\boldsymbol{\sigma}(\vec{u},T)\cdot\vec{n}=\vec{\tau}(\vec{x},t)\text{ on }\Gamma_{\tau}\text{,}\\
-k \nabla T \cdot \vec{n} = q_0(\vec{x},t) \text{ on }\Gamma_q,
\end{gather}
%
we define nondimensional values:
%
\begin{align}
\vec{x}^* &= \frac{\vec{x}}{x_o}, \\
\vec{u}^* &= \frac{\vec{u}}{u_o}, \\
T^* &= \frac{T - T_0}{\Delta T}, \\
\rho^* &= \frac{\rho}{\rho_o}, \\
\vec{f}^* &= \frac{\vec{f}}{f_o}, \\
Q^* &= \frac{Q}{Q_o}, \\
\boldsymbol{\sigma}^* &= \frac{\boldsymbol{\sigma}}{\sigma_o}, \\
\vec{\tau}^* &= \frac{\vec{\tau}}{\sigma_o}, \\
q^* &= \frac{q}{q_o}.
\end{align}
%
We also recognize that
%
\begin{equation}
\boldsymbol{\nabla}^* = x_o \boldsymbol{\nabla}.
\end{equation}

Substituting into the equations, we have
%
\begin{gather}
\rho_o  \rho^*(\vec{x}^*) \frac{u_o}{t_o^2} \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - f_o \vec{f}^*(\vec{x}^*,t^*) - \frac{\sigma_o}{x_o} \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^*(T^*) = \vec{0}\text{ in }\Omega,\\
\rho_o c_o \rho^* c^* \frac{\Delta T}{t_o} \frac{\partial T^*}{\partial t^*} - Q_o Q^*(\vec{x}^*,t^*) - \frac{k_o \Delta T}{x_o^2} \nabla^* \cdot (k^* \nabla^* T^*) = 0 \text{ in }\Omega,\\
\sigma_o \boldsymbol{\sigma}^*(\vec{u^*}, T^*) \cdot \vec{n} = \sigma_o \vec{\tau}^*(\vec{x}^*,t^*)\text{ on }\Gamma_{\tau}\text{,}\\
-\frac{k_o \Delta T}{x_o} k^* \nabla^* T^* \cdot \vec{n} = q_o q_0^*(\vec{x}^*,t^*) \text{ on }\Gamma_q.
\end{gather}

For the boundary condition equations, the nondimensional scales for the terms in each equation are consistent, so we will limit our discussion to the first two equations.
Grouping terms and multiplying the first equation by $\frac{x_o}{\sigma_o}$ and the second by $\frac{x_o^2}{k_o \Delta T}$, we have
%
\begin{gather}
\left( \rho_o \frac{u_o}{t_o^2}\frac{x_o}{\sigma_o}\right) \rho^*(\vec{x}^*) \frac{\partial^{2}\vec{u}^*}{\partial {t^*}^{2}} - \left(f_o \frac{x_o}{\sigma_o}\right) \vec{f}^*(\vec{x}^*,t^*) -  \boldsymbol{\nabla}^*\cdot\boldsymbol{\sigma}^*(T^*) = \vec{0}\text{ in }\Omega, \\
\left( \rho_o c_o \frac{x_o^2}{k_o t_o} \right) \rho^* c^* \frac{\partial T^*}{\partial t^*} - \left( Q_o \frac{x_o^2}{k_o \Delta T} \right) Q^*(\vec{x}^*,t^*) - \nabla^* \cdot (k^* \nabla^* T^*) = 0 \text{ in }\Omega.
\end{gather}
%
All terms should be nondimensional, which implies
%
\begin{align}
f_o &= \frac{\sigma_o}{x_o}, \\
\rho_o &= \frac{t_o^2 \sigma_o}{u_o x_o}, \\
Q_o &= \frac{k_o \Delta T}{x_o^2}.
\end{align}

We want to determine the stress scale, $\sigma_o$.
Considering isotropic, linear thermoelasticity we have
%
\begin{equation}
\boldsymbol{\sigma} = \boldsymbol{C} : \boldsymbol{\epsilon} - \boldsymbol{C} : \alpha (T - T_0) \mathbf{I} = \boldsymbol{C} : \frac{1}{2}\left(\boldsymbol{\nabla} + \boldsymbol{\nabla}^T \right) \vec{u} - \boldsymbol{C} : \alpha (T - T_0) \mathbf{I}.
\end{equation}
%
Substituting in our nondimensional values yields
%
\begin{equation}
\sigma_o \boldsymbol{\sigma}^* = \mu_o \boldsymbol{C}^* : \frac{u_o}{x_o} \frac{1}{2}\left(\boldsymbol{\nabla}^* + \boldsymbol{\nabla}^{*^T}\right) \vec{u}^* - \mu_o \boldsymbol{C}^* : \alpha_o \alpha^* \Delta T T^* \mathbf{I}.
\end{equation}
%
We recognize that for the equation to be nondimensional, both mechanical and thermal terms contribute:
\begin{equation}
\sigma_o = \mu_o \frac{u_o}{x_o} = \mu_o \alpha_o \Delta T.
\end{equation}
%
This gives us
\begin{equation}
\frac{u_o}{x_o} = \alpha_o \Delta T,
\end{equation}
which relates the displacement scale to the temperature scale.

Returning to the expression for $\rho_o$ and substituting in the expression for $\sigma_o$, we have
\begin{align}
\rho_o &= \frac{t_o^2 \sigma_o}{u_o x_o}, \\
\rho_o &= \mu_o \frac{t_o^2}{x_o^2}.
\end{align}

## Nondimensional Groups

### Thermal Biot Number

The second equation contains a nondimensional group analogous to the Biot number in poroelasticity:
\begin{equation}
\Pi_\mathit{thermal} = \frac{\rho_o c_o x_o^2}{k_o t_o}.
\end{equation}
This represents the ratio of thermal diffusion time to the problem time scale.
If $\Pi_\mathit{thermal} \ll 1$, thermal diffusion is fast compared to mechanical loading and the temperature field equilibrates quickly.
If $\Pi_\mathit{thermal} \gg 1$, thermal diffusion is slow and the process is nearly adiabatic on the loading time scale.

### Inertia

The scale of the inertial term is $\Pi_\mathit{inertia} = \frac{\rho_o x_o^2}{\mu_o t_o^2}$.
If this term is O(1), then we should include inertia and solve the dynamic form of the thermoelasticity equation.
If this term is very small, then we can neglect inertia and solve the quasistatic form of the thermoelasticity equation.

**If the time scale equals the time it takes the shear wave ($v_o^2 = \frac{\mu_o}{\rho_o}$) to propagate over the length scale**, then $\Pi_\mathit{inertia} = 1$ and inertia must be included.

````