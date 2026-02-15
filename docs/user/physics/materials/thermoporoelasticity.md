````markdown
(sec-user-physics-thermoporoelasticity)=
# Thermoporoelasticity

You can use the `Thermoporoelasticity` component to solve coupled thermo-poroelastic problems.
The component couples the poroelastic (displacement + pore pressure) equations with the heat (temperature) equation.
Whether inertia or body forces are included is determined by the `Thermoporoelasticity` property settings.
Gravitational body forces are included if the `gravity_field` is set in the `Problem`.
{numref}`tab:thermoporoelasticity:rheologies` lists the thermoporoelastic bulk rheologies implemented for the thermoporoelasticity equation.

```{table} Thermoporoelasticity bulk rheologies
:name: tab:thermoporoelasticity:rheologies
| Bulk Rheology                          | Description                                |
|:--------------------------------------:|:-------------------------------------------|
| `IsotropicLinearThermoporoelasticity` | Isotropic, linear thermo-poroelasticity    |
```

```{table} Properties defining thermoporoelasticity auxiliary subfields
:name: tab:thermoporoelasticity:auxiliary:subfields
| Subfield                         | Required | Components / Notes                            |
|:--------------------------------:|:--------:|:----------------------------------------------|
| `solid_density`                  |    X     |                                               |
| `fluid_density`                  |    X     |                                               |
| `fluid_viscosity`                |    X     |                                               |
| `porosity`                       |    X     |                                               |
| `biot_coefficient`               |    X     |                                               |
| `biot_modulus`                   |    X     |                                               |
| `drained_bulk_modulus`           |    X     | provided by rheology                           |
| `shear_modulus`                  |    X     | provided by rheology                           |
| `isotropic_permeability`         |    X     | or tensor permeability (xx, yy, zz, xy, yz, xz)|
| `reference_temperature`          |    X     |                                               |
| `thermal_expansion_coefficient`  |    X     |                                               |
| `fluid_thermal_expansion`        |    X     |                                               |
| `thermal_conductivity`           |    X     | scalar or tensor                               |
| `specific_heat`                  |    X     | heat capacity per unit mass                    |
| `reference_stress`               |    O     | xx, yy, zz, xy, yz, xz                         |
| `reference_strain`               |    O     | xx, yy, zz, xy, yz, xz                         |
| `body_force`                     |    O     | x, y, z                                        |
| `gravitational_acceleration`     |    O     | x, y, z                                        |
| `source_density`                 |    O     | (optional volumetric mass source)              |
| `heat_source`                    |    O     | volumetric heat source                          |
```

X: required value in auxiliary field spatial database  
O: optional value in auxiliary field spatial database

```{table} Derived subfields that are available for output for thermoporoelasticity bulk rheologies.
:name: tab:thermoporoelasticity:derived:subfields
| Subfield        |  Available  | Components               |
|:---------------:|:-----------:|:-------------------------|
| `cauchy_stress` |     ✓       | xx, yy, zz, xy, yz, xz   |
| `cauchy_strain` |     ✓       | xx, yy, zz, xy, yz, xz   |
| `bulk_density`  |     ✓       |                         |
| `water_content` |     ✓       |                         |
```

When porosity is enabled as a state variable, it will be included in the output along with the derived subfields.

:::{seealso}
See [`Thermoporoelasticity` Component](../../components/materials/Thermoporoelasticity.md) for the Pyre properties and facilities and configuration examples.
:::

````