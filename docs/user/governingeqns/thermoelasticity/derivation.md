````markdown
# Derivation of Thermoelasticity Equations

For completeness we start our discussion of the governing equations with a derivation of the coupled thermoelasticity equations.
Consider domain $\Omega$ bounded by boundary $\Gamma$.

## Momentum Balance

Applying a Lagrangian description of the conservation of momentum gives
%
```{math}
:label: eqn:thermoelasticity:momentum:vec
\frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\, d\Omega=\int_{\Omega}\vec{f}(\vec{x},t)\, d\Omega + \int_{\Gamma}\vec{\tau}(\vec{x},t)\, d\Gamma.
```
%
The traction vector field is related to the stress tensor through
%
\begin{equation}
\vec{\tau}(\vec{x},t) = \boldsymbol{\sigma}(\vec{u},T) \cdot \vec{n},
\end{equation}
%
where $\vec{n}$ is the outward normal vector to $\Gamma$.
Substituting into equation {math:numref}`eqn:thermoelasticity:momentum:vec` yields
%
\begin{equation}
\frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\, d\Omega = \int_{\Omega}\vec{f}(\vec{x},t)\, d\Omega+\int_{\Gamma}\boldsymbol{\sigma}(\vec{u},T)\cdot\vec{n}\, d\Gamma.
\end{equation}
%
Applying the divergence theorem,
%
\begin{equation}
\int_{\Omega}\boldsymbol{\nabla}\cdot\vec{a}\: d\Omega=\int_{\Gamma}\vec{a}\cdot\vec{n}\: d\Gamma,
\end{equation}
%
to the boundary integral results in
%
\begin{equation}
\frac{\partial}{\partial t}\int_{\Omega}\rho(\vec{x})\frac{\partial\vec{u}}{\partial t}\, d\Omega=\int_{\Omega}\vec{f}(\vec{x},t)\, d\Omega+\int_{\Omega}\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u},T)\, d\Omega,
\end{equation}
%
which we can rewrite as
%
\begin{equation}
\int_{\Omega}\left(\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u},T)\right)\, d\Omega=\vec{0}.
\end{equation}
%
Because the domain $\Omega$ is arbitrary, the integrand must be the zero vector at every location in the domain.

## Energy Balance

The energy balance with heat conduction is given by
%
\begin{equation}
\int_{\Omega} \rho c \frac{\partial T}{\partial t} \, d\Omega = \int_{\Omega} Q(\vec{x},t) \, d\Omega + \int_{\Gamma} q_n(\vec{x},t) \, d\Gamma,
\end{equation}
%
where $c$ is the specific heat capacity, $Q$ is the volumetric heat source, and $q_n$ is the heat flux normal to the boundary.
The heat flux is related to the temperature gradient by Fourier's law:
%
\begin{equation}
q_n = -k \nabla T \cdot \vec{n},
\end{equation}
%
where $k$ is the thermal conductivity.
Applying the divergence theorem yields
%
\begin{equation}
\int_{\Omega} \left( \rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) \right) \, d\Omega = 0.
\end{equation}

## Strong Form

Combining the momentum and energy equations, we have
%
\begin{gather}
\rho(\vec{x})\frac{\partial^{2}\vec{u}}{\partial t^{2}}-\vec{f}(\vec{x},t)-\boldsymbol{\nabla}\cdot\boldsymbol{\sigma}(\vec{u},T)=\vec{0}\text{ in }\Omega,\\
\rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) = 0 \text{ in }\Omega,\\
\boldsymbol{\sigma}(\vec{u},T)\cdot\vec{n}=\vec{\tau}(\vec{x},t)\text{ on }\Gamma_{\tau},\\
\vec{u}=\vec{u}_0(\vec{x},t)\text{ on }\Gamma_{u},\\
T = T_0(\vec{x},t) \text{ on } \Gamma_T,\\
-k \nabla T \cdot \vec{n} = q_0(\vec{x},t) \text{ on } \Gamma_q.
\end{gather}
%
We specify tractions, $\vec{\tau}$, on boundary $\Gamma_{\tau}$, displacements, $\vec{u}_0$, on boundary $\Gamma_{u}$, temperatures, $T_0$, on boundary $\Gamma_T$, and heat flux, $q_0$, on boundary $\Gamma_q$.

## Constitutive Relation

For linear isotropic thermoelasticity, the stress tensor is
%
\begin{equation}
\boldsymbol{\sigma}(\vec{u},T) = \boldsymbol{C} : (\boldsymbol{\epsilon} - \alpha (T - T_{\mathrm{ref}}) \mathbf{I}),
\end{equation}
%
where $\boldsymbol{C}$ is the elasticity tensor, $\boldsymbol{\epsilon}$ is the strain tensor, $\alpha$ is the coefficient of thermal expansion, $T_{\mathrm{ref}}$ is a reference temperature, and $\mathbf{I}$ is the identity tensor.

````