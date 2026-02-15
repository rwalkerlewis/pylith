````markdown
# Infinitesimal Strain and No Faults

We begin with the thermoelasticity system including the inertial term,

```{math}
:label: eqn:thermoelasticity:strong:form
\rho \frac{\partial^2\vec{u}}{\partial t^2} - \vec{f}(\vec{x},t) - \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} (\vec{u},T) = \vec{0} \text{ in }\Omega,
```

```{math}
:label: eqn:thermoelasticity:heat:strong:form
\rho c \frac{\partial T}{\partial t} - Q(\vec{x},t) - \nabla \cdot (k \nabla T) = 0 \text{ in }\Omega,
```

```{math}
:label: eqn:thermoelasticity:bc:Neumann
\boldsymbol{\sigma} \cdot \vec{n} = \vec{\tau}(\vec{x},t) \text{ on }\Gamma_\tau,
```

```{math}
:label: eqn:thermoelasticity:bc:Dirichlet
\vec{u} = \vec{u}_0(\vec{x},t) \text{ on }\Gamma_u,
```

```{math}
:label: eqn:thermoelasticity:bc:Temp:Dirichlet
T = T_0(\vec{x},t) \text{ on }\Gamma_T,
```

```{math}
:label: eqn:thermoelasticity:bc:Heat:Neumann
-k \nabla T \cdot \vec{n} = q_0(\vec{x},t) \text{ on }\Gamma_q,
```

where $\vec{u}$ is the displacement vector, $T$ is the temperature, $\rho$ is the mass density, $c$ is the specific heat capacity, $\vec{f}$ is the body force vector, $Q$ is the volumetric heat source, $k$ is the thermal conductivity, $\boldsymbol{\sigma}$ is the Cauchy stress tensor, $\vec{x}$ is the spatial coordinate, and $t$ is time. We specify tractions $\vec{\tau}$ on boundary $\Gamma_\tau$, displacements $\vec{u}_0$ on boundary $\Gamma_u$, temperatures $T_0$ on boundary $\Gamma_T$, and heat flux $q_0$ on boundary $\Gamma_q$.
Because both $\vec{\tau}$ and $\vec{u}$ are vector quantities, there can be some spatial overlap of boundaries $\Gamma_\tau$ and $\Gamma_u$; however, a degree of freedom at any location cannot be associated with both prescribed displacements (Dirichlet) and traction (Neumann) boundary conditions simultaneously. Similar constraints apply to thermal boundaries $\Gamma_T$ and $\Gamma_q$.

## Kinematics

The small-strain tensor is
%
\begin{equation}
\boldsymbol{\epsilon} = \tfrac{1}{2}(\nabla\vec{u} + (\nabla\vec{u})^T).
\end{equation}

## Constitutive Behavior

For linear isotropic thermoelasticity, the stress tensor is
%
\begin{equation}
\boldsymbol{\sigma}(\vec{u},T) = \boldsymbol{C} : (\boldsymbol{\epsilon} - \alpha (T - T_{\mathrm{ref}}) \mathbf{I}),
\end{equation}
%
where $\boldsymbol{C}$ is the elasticity tensor, $\alpha$ is the coefficient of thermal expansion, $T_{\mathrm{ref}}$ is a reference temperature, and $\mathbf{I}$ is the identity tensor.

```{table} Mathematical notation for thermoelasticity equation with infinitesimal strain.
:name: tab:notation:thermoelasticity

| **Category**                   |       **Symbol**        | **Description**                                        |
| :----------------------------- | :---------------------: | :----------------------------------------------------- |
| Unknowns                       |        $\vec{u}$        | Displacement field                                     |
|                                |        $\vec{v}$        | Velocity field                                         |
|                                |           $T$           | Temperature field                                      |
| Derived quantities             |  $\boldsymbol{\sigma}$  | Cauchy stress tensor                                   |
|                                | $\boldsymbol{\epsilon}$ | Cauchy strain tensor                                   |
| Common constitutive parameters |         $\rho$          | Density                                                |
|                                |          $\mu$          | Shear modulus                                          |
|                                |           $K$           | Bulk modulus                                           |
|                                |        $\alpha$         | Coefficient of thermal expansion                       |
|                                |           $c$           | Specific heat capacity                                 |
|                                |           $k$           | Thermal conductivity                                   |
|                                |      $T_{\mathrm{ref}}$ | Reference temperature                                  |
| Source terms                   |        $\vec{f}$        | Body force per unit volume, for example $\rho \vec{g}$ |
|                                |           $Q$           | Volumetric heat source                                 |
```

:::{toctree}
infinitesimal-strain-quasistatic.md
infinitesimal-strain-dynamic.md
:::

````