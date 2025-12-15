# Thermoelasticity 2D Examples

This example suite demonstrates thermoelastic coupling with and without faults in a 2D domain. Thermoelasticity couples the mechanical and thermal fields:

- **Temperature → Stress**: Temperature changes cause thermal strain via thermal expansion
- **Heat Conduction**: Fourier's law governs heat diffusion

## Domain

The domain is a 40 km wide × 20 km deep vertical cross-section through the crust:

- -20 km ≤ x ≤ 20 km
- -20 km ≤ y ≤ 0 km (y=0 is the surface)

For fault examples, a vertical strike-slip fault runs along x=0 from the surface (y=0) to 15 km depth (y=-15 km).

```
         surface (y=0, T=300 K)
     --------------------
     |         |        |
     |         |fault   |
     |         |        |
     |         *        |  <- fault tip at 15 km depth (fault examples only)
     |                  |
     --------------------
         bottom (y=-20 km, T=600 K)
```

## Physics

### Thermoelasticity

The material is governed by coupled thermoelasticity:

- **Thermal Strain**: σ = C:(ε - α(T-T_ref)I)
- **Heat Conduction**: Governed by Fourier's law with isotropic thermal conductivity

### Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Density | 2700 | kg/m³ |
| S-wave velocity | 3500 | m/s |
| P-wave velocity | 6000 | m/s |
| Specific heat | 900 | J/(kg·K) |
| Thermal conductivity | 3.0 | W/(m·K) |
| Reference temperature | 300 | K |
| Thermal expansion coefficient | 2.4×10⁻⁵ | 1/K |

### Boundary Conditions

- **Displacement**: 
  - Left/Right boundaries: Fixed in x and y
  - Bottom boundary: Fixed in y
  - Top boundary: Free (traction-free)
- **Temperature**: 
  - Surface (y=0): 300 K (27°C)
  - Bottom (y=-20 km): 600 K (327°C)
  - Geothermal gradient: 15 K/km

### Initial Conditions

- Linear temperature gradient from surface to bottom (geothermal gradient)

## Simulations

### Step 1: Thermoelasticity Without Fault (`step01_no_fault.cfg`)

This simulation models thermal equilibration in a domain without any faults:
- Initial temperature perturbation (hot region at depth)
- Heat conduction toward equilibrium
- Thermal stresses from temperature changes

**To run:**
```bash
pylith step01_no_fault.cfg
```

### Step 2: Thermoelasticity With Fault (`step02_fault.cfg`)

This simulation adds a kinematic fault with prescribed slip:
- Vertical strike-slip fault at x=0
- Prescribed coseismic slip of 2 m (left-lateral)
- Slip distribution: maximum at 7.5 km depth, tapering to zero at surface and fault tip
- Demonstrates coupled mechanical-thermal response

**To run:**
```bash
pylith step02_fault.cfg
```

### Step 3: Thermal Diffusion Post-Earthquake (`step03_postseismic.cfg`)

This simulation models post-seismic thermal evolution:
- Follows Step 2 with continued heat diffusion
- No additional fault slip (locked fault)
- Long-term thermal equilibration

**To run:**
```bash
pylith step03_postseismic.cfg
```

## Mesh Generation

The mesh is generated using Gmsh. To create the mesh:

```bash
./generate_gmsh.py --write
```

This creates `mesh_tri.msh` with triangular elements and refined mesh near the fault.

## Files

| File | Description |
|------|-------------|
| `generate_gmsh.py` | Gmsh mesh generation script |
| `mesh_tri.msh` | Generated mesh file |
| `pylithapp.cfg` | Common simulation parameters |
| `step01_no_fault.cfg` | Thermoelasticity without fault |
| `step02_fault.cfg` | Thermoelasticity with fault |
| `step03_postseismic.cfg` | Post-seismic thermal diffusion |
| `fault_slip.spatialdb` | Fault slip distribution |
| `initial_temperature.spatialdb` | Initial geothermal gradient |
| `initial_temperature_perturbed.spatialdb` | Initial temperature with perturbation |

## Visualization

Output is written in HDF5 format to the `output/` directory. You can visualize results using:

- **ParaView**: Open the `.xmf` files
- **PyLith viz tools**: Use the `pylith_viz` module

Key fields to examine:
- `temperature`: Shows thermal evolution
- `displacement`: Shows thermally-induced deformation
- `cauchy_stress`: Shows thermal stresses

## Scientific Background

### Thermal Stresses

When a material experiences non-uniform temperature, thermal expansion causes internal stresses. For a constrained material:

σ_thermal = -α E ΔT / (1 - ν)

where:
- α = thermal expansion coefficient
- E = Young's modulus
- ΔT = temperature change
- ν = Poisson's ratio

### Thermal Diffusion

Heat diffuses according to the heat equation:

∂T/∂t = κ ∇²T

where κ = k/(ρc) is the thermal diffusivity. For crustal rocks with κ ≈ 1 mm²/s:

- 1 m diffusion: ~30 years
- 10 m diffusion: ~3000 years
- 1 km diffusion: ~30 million years
