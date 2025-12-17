# Derivation

We derive the coupled thermoporoelasticity equations from conservation of momentum and energy, along with mass balance for the fluid phase.

## Momentum Balance

From conservation of momentum applied to the bulk porous medium, we have
%
\begin{equation}
\rho_b \frac{\partial \vec{v}}{\partial t} = \vec{f}(\vec{x},t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u},p,T),
\end{equation}
%
where $\rho_b = (1-\phi)\rho_s + \phi \rho_f$ is the bulk density, $\phi$ is the porosity, $\rho_s$ is the solid density, $\rho_f$ is the fluid density, $\vec{v} = \frac{\partial \vec{u}}{\partial t}$ is the velocity, $\vec{f}$ is the body force per unit volume, and $\boldsymbol{\sigma}$ is the Cauchy stress tensor.

## Fluid Mass Balance

Conservation of fluid mass gives
%
\begin{equation}
\frac{\partial \zeta}{\partial t} + \nabla \cdot \vec{q} = \gamma(\vec{x},t),
\end{equation}
%
where $\zeta$ is the variation in fluid content, $\vec{q}$ is the Darcy flux, and $\gamma$ is the rate of injected fluid per unit volume of the porous solid.

## Energy Balance

The energy balance equation for the coupled thermal-fluid-solid system is
%
\begin{equation}
\rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) + \rho_f c_f \vec{q} \cdot \nabla T = 0,
\end{equation}
%
where $c_b$ is the bulk specific heat, $c_f$ is the fluid specific heat, $T$ is temperature, $k$ is the bulk thermal conductivity, $Q$ is a volumetric heat source, and the term $\rho_f c_f \vec{q} \cdot \nabla T$ represents heat advection by the fluid flow.

## Strong Form

Applying the divergence theorem and incorporating boundary conditions, the strong form of the thermoporoelasticity equations becomes:
%
\begin{gather}
\rho_b \frac{\partial \vec{v}}{\partial t} - \vec{f}(\vec{x},t) - \nabla \cdot \boldsymbol{\sigma}(\vec{u},p,T) = \vec{0} \text{ in }\Omega, \\
\frac{\partial \zeta(\vec{u},p,T)}{\partial t} + \nabla \cdot \vec{q}(p,T) - \gamma(\vec{x},t) = 0 \text{ in }\Omega, \\
\rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) + \rho_f c_f \vec{q}(p,T) \cdot \nabla T = 0 \text{ in }\Omega, \\
\boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on }\Gamma_\tau, \\
\vec{u} = \vec{u}_0(\vec{x},t) \text{ on }\Gamma_u, \\
\vec{q} \cdot \vec{n} = q_0(\vec{x},t) \text{ on }\Gamma_q, \\
p = p_0(\vec{x},t) \text{ on }\Gamma_p, \\
T = T_0(\vec{x},t) \text{ on }\Gamma_T, \\
-k \nabla T \cdot \vec{n} = h_0(\vec{x},t) \text{ on }\Gamma_h.
\end{gather}

## Constitutive Relations

### Solid Stress-Strain Relation

For linear isotropic thermoporoelasticity, the effective stress is
%
\begin{equation}
\boldsymbol{\sigma}(\vec{u},p,T) = \boldsymbol{C} : \boldsymbol{\epsilon} - \alpha p \mathbf{I} - \boldsymbol{C} : \alpha_T (T - T_0) \mathbf{I},
\end{equation}
%
where $\boldsymbol{C}$ is the drained elasticity tensor, $\boldsymbol{\epsilon} = \frac{1}{2}(\nabla \vec{u} + \nabla^T \vec{u})$ is the strain tensor, $\alpha$ is the Biot coefficient, $\alpha_T$ is the linear thermal expansion coefficient of the bulk medium, and $T_0$ is a reference temperature.

### Fluid Content

The variation in fluid content couples volumetric strain, pressure, and temperature:
%
\begin{equation}
\zeta(\vec{u},p,T) = \alpha \epsilon_v + \frac{p}{M} + \beta (T - T_0),
\end{equation}
%
where $\epsilon_v = \nabla \cdot \vec{u}$ is the volumetric strain, $M$ is the Biot modulus, and $\beta$ is the thermal expansion coefficient of the fluid-solid system.

### Darcy Flow

Fluid flow follows Darcy's law:
%
\begin{equation}
\vec{q}(p,T) = -\frac{\boldsymbol{k}(T)}{\mu_f(T)}(\nabla p - \vec{f}_f),
\end{equation}
%
where $\boldsymbol{k}$ is the intrinsic permeability (which may be temperature-dependent), $\mu_f$ is the fluid viscosity (which is typically temperature-dependent), and $\vec{f}_f$ is the fluid body force (e.g., $\rho_f \vec{g}$).