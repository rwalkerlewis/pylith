````markdown
# Infinitesimal Strain and Prescribed Fault Slip

For each fault, we add a boundary condition prescribing the jump in the displacement field across the fault, following the same approach as in thermoelasticity and poroelasticity:
%
```{math}
:label: eqn:bc:prescribed:slip:thermoporoelasticity
\begin{gathered}
  -\vec{u}^+ + \vec{u}^- + \vec{d}(\vec{x},t) = \vec{0} \text{ on }\Gamma_f.
\end{gathered}
```
%
We enforce the jump in displacements using a Lagrange multiplier corresponding to equal and opposite tractions on the two sides of the fault.

## Frictional Heating

As in thermoelasticity, fault slip generates heat through frictional dissipation:
\begin{equation}
  \Phi_f = \vec{\lambda} \cdot \frac{\partial \vec{d}}{\partial t},
\end{equation}
where $\vec{\lambda}$ is the fault traction and $\frac{\partial \vec{d}}{\partial t}$ is the slip rate.
This heating term can be incorporated as a volumetric heat source in the fault zone or as a boundary condition on the fault surface.

## Pore Pressure Effects on Faults

In thermoporoelasticity, pore pressure directly affects the effective stress on the fault:
\begin{equation}
  \boldsymbol{\sigma}_{\mathrm{eff}} = \boldsymbol{\sigma} + \alpha p \mathbf{I}.
\end{equation}
Elevated pore pressure reduces the effective normal stress, potentially promoting fault slip. The coupling between slip, pore pressure evolution, and temperature is critical in earthquake mechanics.

## Permeability Changes

Fault slip can alter the permeability structure of the fault zone:
- Compaction during slip may reduce permeability
- Dilatancy may increase permeability
- Fracturing during dynamic rupture can create high-permeability damage zones

These permeability changes affect subsequent fluid flow and pressure evolution, creating complex feedback between mechanical, hydraulic, and thermal processes.

## Thermal Pressurization

Frictional heating in low-permeability fault zones can cause thermal pressurization: as temperature rises, thermal expansion of the pore fluid increases pore pressure if drainage is restricted. This positive feedback can lead to significant weakening during earthquake rupture.

```{table} Mathematical notation for thermoporoelasticity with prescribed slip on faults.
:name: tab:notation:thermoporoelasticity:prescribed:slip
| Category                       |         Symbol          | Description                                                                                       |
| :----------------------------- | :---------------------: | :------------------------------------------------------------------------------------------------ |
| Unknowns                       |        $\vec{u}$        | Displacement field                                                                                |
|                                |        $\vec{v}$        | Velocity field                                                                                    |
|                                |           $p$           | Pore fluid pressure field                                                                         |
|                                |           $T$           | Temperature field                                                                                 |
|                                |     $\vec{\lambda}$     | Lagrange multiplier field                                                                         |
| Derived quantities             |  $\boldsymbol{\sigma}$  | Cauchy stress tensor                                                                              |
|                                | $\boldsymbol{\epsilon}$ | Cauchy strain tensor                                                                              |
|                                |        $\vec{q}$        | Darcy flux                                                                                        |
| Common constitutive parameters |       $\rho_b$          | Bulk density                                                                                      |
|                                |       $\rho_f$          | Fluid density                                                                                     |
|                                |          $\mu$          | Shear modulus                                                                                     |
|                                |        $\alpha$         | Biot coefficient                                                                                  |
|                                |      $\alpha_T$         | Linear thermal expansion coefficient                                                              |
|                                |           $M$           | Biot modulus                                                                                      |
|                                |    $\boldsymbol{k}$     | Permeability                                                                                      |
|                                |        $\mu_f$          | Fluid viscosity                                                                                   |
|                                |           $k$           | Thermal conductivity                                                                              |
|                                |         $c_b$           | Bulk specific heat capacity                                                                       |
|                                |         $c_f$           | Fluid specific heat capacity                                                                      |
| Source terms                   |        $\vec{f}$        | Body force per unit volume                                                                        |
|                                |        $\gamma$         | Fluid source density                                                                              |
|                                |           $Q$           | Volumetric heat source                                                                            |
|                                |        $\vec{d}$        | Slip vector field on the fault                                                                    |
|                                |      $\Phi_f$           | Frictional heating rate per unit area on fault                                                    |
```

````