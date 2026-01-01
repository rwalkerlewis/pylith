# Infinitesimal Strain and No Faults

We base this formulation for thermoporoelasticity on extensions of {cite:t}`Detournay:Cheng:Poroelasticity:1993` and {cite:t}`Cheng:2016:Poroelasticity` to include thermal coupling.
We assume a slightly compressible fluid that completely saturates a porous solid, undergoing infinitesimal strain with coupled thermal effects.

We begin with the conservation of linear momentum, including inertia:
%
```{math}
:label: eqn:thermoporoelasticity:momentum
\rho_b\frac{\partial^2 \vec{u}}{\partial t^2} = \vec{f}(t) + \nabla \cdot \boldsymbol{\sigma}(\vec{u},p,T).
```
%
Enforcing mass balance of the fluid gives
%
```{math}
:label: eqn:thermoporoelasticity:mass
\begin{gather}
  \frac{\partial \zeta(\vec{u},p,T)}{\partial t} + \nabla \cdot \vec{q}(p,T) =
  \gamma(\vec{x},t) \text{ in } \Omega, \\
%
  \vec{q} \cdot \vec{n} = q_0(\vec{x},t) \text{ on }\Gamma_q, \\
%
  p = p_0(\vec{x},t) \text{ on }\Gamma_p,
\end{gather}
```
%
where $\zeta$ is the variation in fluid content, $\vec{q}$ is the rate of fluid volume crossing a unit area of the porous solid, $\gamma$ is the rate of injected fluid per unit volume of the porous solid, $q_0$ is the outward fluid velocity normal to the boundary $\Gamma_q$, and $p_0$ is the fluid pressure on boundary $\Gamma_p$.

The energy balance equation is
%
```{math}
:label: eqn:thermoporoelasticity:energy
\begin{gather}
\rho_b c_b \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) + \rho_f c_f \vec{q} \cdot \nabla T = 0 \text{ in }\Omega, \\
%
T = T_0(\vec{x},t) \text{ on }\Gamma_T, \\
%
-k \nabla T \cdot \vec{n} = h_0(\vec{x},t) \text{ on }\Gamma_h,
\end{gather}
```
%
where $c_b$ is the bulk specific heat, $c_f$ is the fluid specific heat, $k$ is the thermal conductivity, $Q$ is a volumetric heat source, and the advection term $\rho_f c_f \vec{q} \cdot \nabla T$ represents heat transport by fluid flow.

We require the fluid flow to follow Darcy's law (Navier-Stokes equation neglecting inertial effects),
%
\begin{equation}
  \vec{q}(p,T) = -\frac{\boldsymbol{k}(T)}{\mu_{f}(T)}(\nabla p - \vec{f}_f),
\end{equation}
%
where $\boldsymbol{k}$ is the intrinsic permeability (potentially temperature-dependent), $\mu_f$ is the viscosity of the fluid (typically temperature-dependent), $p$ is the fluid pressure, and $\vec{f}_f$ is the body force in the fluid.
If gravity is included in a problem, then usually $\vec{f}_f = \rho_f \vec{g}$, where $\rho_f$ is the density of the fluid and $\vec{g}$ is the gravitational acceleration vector.

## Constitutive Behavior

We assume linear elasticity for the solid phase with thermal and poroelastic coupling, so the constitutive behavior can be expressed as
%
\begin{equation}
  \boldsymbol{\sigma}(\vec{u},p,T) = \boldsymbol{C} : \boldsymbol{\epsilon} - \alpha p \boldsymbol{I} - \boldsymbol{C} : \alpha_T (T - T_0) \boldsymbol{I},
\end{equation}
%
where $\boldsymbol{\sigma}$ is the stress tensor, $\boldsymbol{C}$ is the drained tensor of elasticity constants, $\alpha$ is the Biot coefficient (effective stress coefficient), $\alpha_T$ is the linear thermal expansion coefficient, $\boldsymbol{\epsilon}$ is the strain tensor, $T_0$ is a reference temperature, and $\boldsymbol{I}$ is the identity tensor.

For the constitutive behavior of the fluid, we use the volumetric strain, pressure, and temperature to couple the fluid-solid-thermal behavior:
%
\begin{gather}
  \zeta(\vec{u},p,T) = \alpha \mathop{\mathrm{Tr}}({\boldsymbol{\epsilon}}) + \frac{p}{M} + \beta (T - T_0), \\
%
  \frac{1}{M} = \frac{\alpha-\phi}{K_s} + \frac{\phi}{K_f},
\end{gather}
%
where $1/M$ is the specific storage coefficient at constant strain and temperature, $K_s$ is the bulk modulus of the solid, $K_f$ is the bulk modulus of the fluid, and $\beta$ is the thermal expansion coefficient of the fluid content.
We can write the trace of the strain tensor as the dot product of the gradient and displacement field, so we have
%
\begin{equation}
  \zeta(\vec{u},p,T) = \alpha (\nabla \cdot \vec{u}) + \frac{p}{M} + \beta (T - T_0).
\end{equation}

```{table} Mathematical notation for thermoporoelasticity with infinitesimal strain.
:name: tab:notation:thermoporoelasticity

| **Category**                   |       **Symbol**        | **Description**                                                                                               |
| :----------------------------- | :---------------------: | :------------------------------------------------------------------------------------------------------------ |
| Unknowns                       |        $\vec{u}$        | Displacement field                                                                                            |
|                                |        $\vec{v}$        | Velocity field                                                                                                |
|                                |           $p$           | Pressure field (corresponds to pore fluid pressure)                                                           |
|                                |           $T$           | Temperature field                                                                                             |
|                                |     $\epsilon_{v}$      | Volumetric (trace) strain                                                                                     |
| Derived quantities             |  $\boldsymbol{\sigma}$  | Cauchy stress tensor                                                                                          |
|                                | $\boldsymbol{\epsilon}$ | Cauchy strain tensor                                                                                          |
|                                |         $\zeta$         | Variation of fluid content, $\alpha \epsilon_{v} + \frac{p}{M} + \beta (T - T_0)$                             |
|                                |       $\rho_{b}$        | Bulk density, $\left(1 - \phi\right) \rho_{s} + \phi \rho_{f}$                                                |
|                                |        $\vec{q}$        | Darcy flux, $-\frac{\boldsymbol{k}}{\mu_{f}} \cdot \left(\nabla p - \vec{f}_{f}\right)$                       |
| Common constitutive parameters |       $\rho_{s}$        | Solid (matrix) density                                                                                        |
|                                |       $\rho_{f}$        | Fluid density                                                                                                 |
|                                |        $\mu_{f}$        | Fluid viscosity                                                                                               |
|                                |         $\phi$          | Porosity                                                                                                      |
|                                |          $\mu$          | Shear modulus                                                                                                 |
|                                |         $K_{d}$         | Drained bulk modulus                                                                                          |
|                                |        $\alpha$         | Biot coefficient, $1 - \frac{K_{d}}{K_{s}}$                                                                   |
|                                |      $\alpha_T$         | Linear thermal expansion coefficient                                                                          |
|                                |           $M$           | Biot modulus                                                                                                  |
|                                |        $\beta$          | Thermal expansion coefficient of fluid content                                                                |
|                                |    $\boldsymbol{k}$     | Permeability                                                                                                  |
|                                |           $k$           | Thermal conductivity                                                                                          |
|                                |         $c_b$           | Bulk specific heat capacity                                                                                   |
|                                |         $c_f$           | Fluid specific heat capacity                                                                                  |
| Source terms                   |        $\vec{f}$        | Body force per unit volume, for example: $\rho_{b} \vec{g}$                                                   |
|                                |      $\vec{f}_{f}$      | Fluid body force, for example: $\rho_{f} \vec{g}$                                                             |
|                                |        $\gamma$         | Source density; rate of injected fluid per unit volume of the porous solid                                    |
|                                |           $Q$           | Volumetric heat source                                                                                        |
```

:::{toctree}
infinitesimal-strain-quasistatic.md
infinitesimal-strain-dynamic.md
:::

````