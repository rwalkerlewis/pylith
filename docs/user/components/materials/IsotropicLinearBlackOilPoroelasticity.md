# IsotropicLinearBlackOilPoroelasticity

## Overview

The `IsotropicLinearBlackOilPoroelasticity` rheology extends the standard isotropic linear poroelasticity with pressure-dependent fluid properties commonly used in petroleum reservoir simulation. This is known as the "black oil" model.

## Black Oil Model Physics

The black oil formulation introduces pressure-dependent fluid properties:

### Pressure-Dependent Fluid Viscosity

$$\mu_{eff}(p) = \mu_{ref} \cdot \exp(\alpha_\mu (p - p_{ref}))$$

where:
- $\mu_{ref}$ is the reference viscosity (`fluid_viscosity`)
- $\alpha_\mu$ is the viscosity coefficient (`viscosity_coefficient`)
- $p_{ref}$ is the reference pressure (`reference_pressure`)

### Pressure-Dependent Fluid Compressibility

$$c_f(p) = c_{f0} + c_{f1} (p - p_{ref})$$

where:
- $c_{f0}$ is the reference compressibility (`fluid_compressibility`)
- $c_{f1}$ is the compressibility coefficient (`fluid_compressibility_coefficient`)

### Effective Biot Modulus

$$\frac{1}{M_{eff}} = \frac{1}{M_{solid}} + \phi \cdot c_f(p)$$

where:
- $M_{solid}$ is the solid matrix Biot modulus (`biot_modulus`)
- $\phi$ is the porosity (`porosity`)

## Auxiliary Subfields

The black oil rheology requires the following auxiliary subfields:

| Subfield | Description | Units |
|----------|-------------|-------|
| `shear_modulus` | Shear modulus G | Pa |
| `drained_bulk_modulus` | Drained bulk modulus K_d | Pa |
| `biot_coefficient` | Biot coefficient α | - |
| `biot_modulus` | Biot modulus M | Pa |
| `reference_pressure` | Reference pressure for fluid property correlations | Pa |
| `fluid_compressibility` | Fluid compressibility at reference pressure | 1/Pa |
| `fluid_compressibility_coefficient` | Rate of change of compressibility with pressure | 1/Pa² |
| `viscosity_coefficient` | Rate of change of log viscosity with pressure | 1/Pa |
| `isotropic_permeability` | Isotropic permeability k | m² |

## Configuration Example

```cfg
[pylithapp.problem.materials]
poroelastic_mat.bulk_rheology = pylith.materials.IsotropicLinearBlackOilPoroelasticity

[pylithapp.problem.materials.poroelastic_mat]
description = Black oil poroelastic material
label_value = 1

db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Black oil properties
db_auxiliary_field.iohandler.filename = blackoil_properties.spatialdb

[pylithapp.problem.materials.poroelastic_mat.bulk_rheology]
use_reference_state = False
use_tensor_permeability = False

auxiliary_subfields.shear_modulus.basis_order = 0
auxiliary_subfields.drained_bulk_modulus.basis_order = 0
auxiliary_subfields.biot_coefficient.basis_order = 0
auxiliary_subfields.biot_modulus.basis_order = 0
auxiliary_subfields.reference_pressure.basis_order = 0
auxiliary_subfields.fluid_compressibility.basis_order = 0
auxiliary_subfields.fluid_compressibility_coefficient.basis_order = 0
auxiliary_subfields.viscosity_coefficient.basis_order = 0
auxiliary_subfields.isotropic_permeability.basis_order = 0
```

## Spatial Database Example

```spatialdb
#SPATIAL.ascii 1
SimpleDB {
  num-values = 13

  value-names =  solid_density  fluid_density  fluid_viscosity  porosity  shear_modulus  drained_bulk_modulus  biot_coefficient  fluid_bulk_modulus  reference_pressure  fluid_compressibility  fluid_compressibility_coefficient  viscosity_coefficient  isotropic_permeability
  value-units =  kg/m**3  kg/m**3  Pa*s  none  Pa  Pa  none  Pa  Pa  1/Pa  1/Pa**2  1/Pa  m**2

  num-locs = 1
  data-dim = 0
  space-dim = 2

  cs-data = cartesian {
    to-meters = 1.0
    space-dim = 2
  }
}
  0.0  0.0  2500.0  1000.0  1.0e-3  0.2  3.0e+10  8.0e+10  0.8  2.0e+9  1.0e+7  5.0e-10  1.0e-18  1.0e-8  1.0e-14
```

## Python API

```python
from pylith.materials.IsotropicLinearBlackOilPoroelasticity import IsotropicLinearBlackOilPoroelasticity

rheology = IsotropicLinearBlackOilPoroelasticity()
rheology.useReferenceState = False
rheology.useTensorPermeability = False
```

## Notes

1. The black oil model introduces nonlinearity through the pressure-dependent properties. For significant pressure variations, use the nonlinear solver.

2. For small pressure changes near the reference pressure, the model behaves similarly to standard linear poroelasticity.

3. The viscosity coefficient typically has values around 1e-8 to 1e-7 1/Pa for petroleum applications.

4. The compressibility coefficient can be zero for constant compressibility (linear case).
