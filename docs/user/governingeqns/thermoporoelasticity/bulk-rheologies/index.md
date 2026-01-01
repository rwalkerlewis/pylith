(sec-user-governingeqns-thermoporoelasticity-rheologies)=
# Bulk Rheologies

In this section we describe the mathematical formulations of bulk rheologies for thermoporoelasticity.
The coupled system requires constitutive relations for the solid skeleton, fluid content, and fluid flow, all of which may be temperature-dependent.

## Linear Isotropic Thermoporoelasticity

The fundamental constitutive relations are:

### Effective Stress

\begin{equation}
\boldsymbol{\sigma} = \boldsymbol{C} : \boldsymbol{\epsilon} - \alpha p \mathbf{I} - \boldsymbol{C} : \alpha_T (T - T_0) \mathbf{I},
\end{equation}
%
where $\alpha$ is the Biot coefficient, $\alpha_T$ is the thermal expansion coefficient, and $\boldsymbol{C}$ is the drained elasticity tensor.

For isotropic materials:
%
\begin{equation}
\boldsymbol{\sigma} = 2\mu \boldsymbol{\epsilon} + \lambda \mathop{\mathrm{Tr}}(\boldsymbol{\epsilon}) \mathbf{I} - \alpha p \mathbf{I} - (3\lambda + 2\mu) \alpha_T (T - T_0) \mathbf{I}.
\end{equation}

### Fluid Content

\begin{equation}
\zeta = \alpha \epsilon_v + \frac{p}{M} + \beta (T - T_0),
\end{equation}
%
where $M$ is the Biot modulus and $\beta$ is the thermal expansion coefficient of the fluid content:
%
\begin{equation}
\beta = \phi \beta_f + (\alpha - \phi) \beta_s,
\end{equation}
%
with $\beta_f$ and $\beta_s$ the volumetric thermal expansion coefficients of the fluid and solid, respectively.

### Darcy Flow

\begin{equation}
\vec{q} = -\frac{\boldsymbol{k}(T)}{\mu_f(T)} (\nabla p - \rho_f \vec{g}),
\end{equation}
%
where both permeability and fluid viscosity may depend on temperature.

## Temperature-Dependent Properties

### Fluid Viscosity

Fluid viscosity typically decreases with temperature. For water, an empirical relation is:
%
\begin{equation}
\mu_f(T) = \mu_{f,0} \exp\left(-\beta_\mu (T - T_0)\right),
\end{equation}
%
or more accurately, the Vogel-Fulcher-Tammann equation.

### Permeability

Permeability may change with temperature due to:
- Thermal expansion of the solid matrix (typically small effect)
- Changes in pore structure from thermal stresses
- Phase changes or mineral dissolution/precipitation at higher temperatures

### Bulk Properties

The drained and undrained bulk moduli, Biot coefficient, and Biot modulus all depend on the constituent properties:
%
\begin{align}
K_d &= K_d(T), \\
\alpha &= 1 - \frac{K_d(T)}{K_s(T)}, \\
\frac{1}{M} &= \frac{\alpha - \phi}{K_s(T)} + \frac{\phi}{K_f(T)}.
\end{align}

## Nonlinear Extensions

For large temperature or pressure variations:
- Nonlinear thermal expansion: $\boldsymbol{\epsilon}^{\mathrm{th}} = \alpha(T) (T - T_0) \mathbf{I}$
- Pressure-dependent permeability: $\boldsymbol{k} = \boldsymbol{k}(p, T)$
- Stress-dependent permeability: $\boldsymbol{k} = \boldsymbol{k}(\boldsymbol{\sigma}, T)$
- Temperature-dependent storage: $M = M(T)$
```