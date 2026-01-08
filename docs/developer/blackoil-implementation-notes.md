# Black Oil Poroelasticity Implementation Notes

## Overview

This document describes the implementation of the black oil formulation for poroelastic materials in PyLith. The black oil model extends the standard isotropic linear poroelasticity model by incorporating pressure-dependent fluid properties.

## Physics

The black oil formulation models:

1. **Pressure-dependent viscosity**: 
   ```
   μ_eff(p) = μ_ref * exp(c_μ * (p - p_ref))
   ```

2. **Pressure-dependent compressibility**:
   ```
   c_eff(p) = c_ref * exp(c_c * (p - p_ref))
   ```

3. **Effective Biot modulus**:
   ```
   M_eff = 1 / (φ * c_eff + (α - φ) / K_s)
   ```

where:
- `μ_ref` = reference fluid viscosity at `p_ref`
- `c_ref` = reference fluid compressibility at `p_ref`
- `p_ref` = reference pressure
- `c_μ` = viscosity coefficient (pressure sensitivity)
- `c_c` = compressibility coefficient (pressure sensitivity)
- `φ` = porosity
- `α` = Biot coefficient
- `K_s` = solid grain bulk modulus

## Implementation Files

### C++ Core (libsrc)

1. **`libsrc/pylith/materials/IsotropicLinearBlackOilPoroelasticity.hh/.cc`**
   - Main rheology class extending `RheologyPoroelasticity`
   - Manages black oil specific auxiliary fields
   - Returns appropriate kernel functions for residuals and Jacobians

2. **`libsrc/pylith/materials/AuxiliaryFactoryBlackOilPoroelastic.hh/.cc`**
   - Factory for creating black oil specific auxiliary subfields:
     - `reference_pressure`
     - `fluid_compressibility`
     - `fluid_compressibility_coefficient`
     - `viscosity_coefficient`

3. **`libsrc/pylith/fekernels/IsotropicLinearBlackOilPoroelasticity.hh`**
   - Finite element kernels implementing the black oil physics
   - Includes residual and Jacobian calculations for 2D (plane strain) and 3D

### SWIG Interface (modulesrc)

1. **`modulesrc/materials/IsotropicLinearBlackOilPoroelasticity.i`**
   - SWIG interface exposing C++ class to Python

### Python Wrappers (pylith)

1. **`pylith/materials/IsotropicLinearBlackOilPoroelasticity.py`**
   - Python wrapper class with Pyre configuration

2. **`pylith/materials/AuxSubfieldsIsotropicLinearBlackOilPoroelasticity.py`**
   - Pyre facility for auxiliary subfield configuration

### Tests

1. **Fullscale Tests** (`tests/fullscale/poroelasticity/blackoil/`)
   - `test_pylith.py` - Test runner
   - `TestBlackOil.py` - Test case module
   - `blackoil_soln.py` - Analytical solution
   - `pylithapp.cfg`, `blackoil.cfg`, etc. - Configuration files

### Documentation

1. **`docs/user/components/materials/IsotropicLinearBlackOilPoroelasticity.md`**
   - User documentation for the black oil component

## Building

The implementation requires PyLith to be built with:
- PETSc
- Spatialdata
- Pythia

Build files modified:
- `libsrc/pylith/Makefile.am` - Added C++ source files
- `libsrc/pylith/materials/Makefile.am` - Added header files
- `libsrc/pylith/fekernels/Makefile.am` - Added kernel header
- `modulesrc/materials/Makefile.am` - Added SWIG interface
- `modulesrc/materials/materials.i` - Added includes

## Configuration Example

```ini
[pylithapp.problem.materials.mat_poroelastic]
bulk_rheology = pylith.materials.IsotropicLinearBlackOilPoroelasticity

[pylithapp.problem.materials.mat_poroelastic.bulk_rheology]
auxiliary_subfields.reference_pressure.basis_order = 0
auxiliary_subfields.fluid_compressibility.basis_order = 0
auxiliary_subfields.fluid_compressibility_coefficient.basis_order = 0
auxiliary_subfields.viscosity_coefficient.basis_order = 0
```

## Spatial Database Example

```
#SPATIAL_GRID.ascii 1
SimpleDB {
  num-values = 4
  value-names = reference_pressure fluid_compressibility fluid_compressibility_coefficient viscosity_coefficient
  value-units = Pa 1/Pa 1/Pa 1/Pa
  num-locs = 1
  data-dim = 0
  space-dim = 2
  cs-data = cartesian {
    space-dim = 2
  }
}
0.0  0.0   1.0e7  1.0e-10  1.0e-8  1.0e-8
```

## Testing Strategy

### Unit Tests (libtests)
- Test auxiliary factory
- Test kernel computations

### MMS Tests (mmstests)
- Verify convergence rates
- Test with manufactured solutions

### Fullscale Tests
- Integration tests with realistic configurations
- Compare with known solutions

## Known Limitations

1. Only implicit time-stepping is currently supported
2. Dynamic (inertia) terms not yet implemented
3. State variables not currently tracked for black oil specific fields

## Future Work

1. Add explicit time-stepping support
2. Implement state variable tracking
3. Add MMS test cases specific to pressure-dependent properties
4. Support for additional black oil parameters (e.g., gas-oil ratio)
