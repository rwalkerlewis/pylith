````markdown
(sec-user-physics-thermoelasticity)=
# Thermoelasticity

You can use the `Thermoelasticity` component to solve coupled thermo-mechanical problems.
The component couples the elasticity (displacement) equation with the heat (temperature) equation.
Whether inertia or body forces are included is determined by the `Thermoelasticity` property settings.
Gravitational body forces are included if the `gravity_field` is set in the `Problem`.
{numref}`tab:thermoelasticity:rheologies` lists the thermoelastic bulk rheologies implemented for the thermoelasticity equation.

```{table} Thermoelasticity bulk rheologies
:name: tab:thermoelasticity:rheologies
| Bulk Rheology                    | Description                            |
|:---------------------------------|:---------------------------------------|
| `IsotropicLinearThermoelasticity`| Isotropic, linear thermoelasticity     |
```

```{table} Properties defining thermoelasticity auxiliary subfields
:name: tab:thermoelasticity:auxiliary:subfields
| Subfield                      |  Required  | Components / Notes                         |
|:-----------------------------:|:----------:|:-------------------------------------------|
| `density`                     |     X      |                                            |
| `specific_heat`               |     X      | heat capacity per unit mass                 |
| `thermal_conductivity`        |     X      | scalar or tensor (xx, yy, zz, xy, yz, xz)   |
| `reference_temperature`       |     X      |                                            |
| `thermal_expansion_coefficient`|    X      | scalar                                      |
| `shear_modulus`               |     X      | provided by rheology                        |
| `bulk_modulus`                |     X      | provided by rheology                        |
| `body_force`                  |     O      | x, y, z                                     |
| `heat_source`                 |     O      | volumetric heat source                      |
| `gravitational_acceleration`  |     O      | x, y, z (used with `gravity_field`)         |
```

X: required value in auxiliary field spatial database  
O: optional value in auxiliary field spatial database

```{table} Derived subfields that are available for output for thermoelasticity bulk rheologies.
:name: tab:thermoelasticity:derived:subfields
|      Subfield    |  Available  | Components               |
|:----------------:|:-----------:|:-------------------------|
| `cauchy_stress`  |     ✓       | xx, yy, zz, xy, yz, xz   |
| `cauchy_strain`  |     ✓       | xx, yy, zz, xy, yz, xz   |
```

:::{seealso}
See [`Thermoelasticity` Component](../../components/materials/Thermoelasticity.md) for the Pyre properties and facilities and configuration examples.
:::

````