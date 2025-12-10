# Earthquake Thermoelasticity 2D Example

This example demonstrates thermoelastic coupling during and after an earthquake, including:

1. **Coseismic slip** on a strike-slip fault
2. **Shear heating** from fault slip
3. **Heat conduction** away from the fault zone

## Domain

The domain is a 40 km wide × 20 km deep vertical cross-section through the crust:

- -20 km ≤ x ≤ 20 km
- -20 km ≤ y ≤ 0 km (y=0 is the surface)

A vertical strike-slip fault runs along x=0 from the surface (y=0) to 15 km depth (y=-15 km).

```
         surface (y=0, T=300 K)
     --------------------
     |         |        |
     |         |fault   |
     |  slip<--|        |
     |         *        |  <- fault tip at 15 km depth
     |                  |
     --------------------
         bottom (y=-20 km, T=600 K)
```

## Physics

### Thermoelasticity

The material is governed by coupled thermoelasticity:

- **Thermal → Mechanical**: Temperature changes cause thermal strain: σ = C:(ε - α(T-T_ref)I)
- **Heat conduction**: Fourier's law with isotropic thermal conductivity

### Shear Heating

During fault slip, frictional work generates heat proportional to:

Q = τ × slip_rate

where τ is the shear stress and slip_rate is the velocity of slip. This heat source is localized on the fault.

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

- **Displacement**: Fixed on all boundaries
- **Temperature**: 
  - Surface (y=0): 300 K (27°C)
  - Bottom (y=-20 km): 600 K (327°C)
  - This gives a geothermal gradient of 15 K/km

### Initial Conditions

- Linear temperature gradient from surface to bottom (geothermal gradient)

## Simulations

### Step 1: Coseismic Slip with Shear Heating (`step01_coseismic.cfg`)

This simulation models:
- Instantaneous coseismic slip of up to 2 m (left-lateral)
- Slip distribution: maximum at 7.5 km depth, tapering to zero at surface and fault tip
- Shear heating from the slip
- Initial post-seismic heat diffusion (100 years)

**To run:**
```bash
pylith step01_coseismic.cfg
```

### Step 2: Post-seismic Heat Diffusion (`step02_postseismic.cfg`)

This simulation models:
- Long-term (1000 years) heat diffusion after the earthquake
- No additional fault slip (locked fault)
- Thermal relaxation of the shear heating anomaly

**To run:**
```bash
pylith step02_postseismic.cfg
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
| `step01_coseismic.cfg` | Coseismic slip simulation |
| `step02_postseismic.cfg` | Post-seismic heat diffusion |
| `fault_slip.spatialdb` | Fault slip distribution |
| `initial_temperature.spatialdb` | Initial geothermal gradient |
| `initial_temperature_heated.spatialdb` | Initial temperature with shear heating anomaly |

## Visualization

Output is written in HDF5 format to the `output/` directory. You can visualize results using:

- **ParaView**: Open the `.xmf` files
- **PyLith viz tools**: Use the `pylith_viz` module

Key fields to examine:
- `temperature`: Shows heating near fault and subsequent diffusion
- `displacement`: Shows coseismic deformation
- `cauchy_stress`: Shows stress changes from thermal and mechanical loading

## Scientific Background

### Shear Heating in Earthquakes

During an earthquake, the rapid slip on a fault generates heat through friction. The temperature rise ΔT can be estimated as:

ΔT = (τ × d) / (ρ × c × w)

where:
- τ = shear stress (~10-100 MPa)
- d = slip (~1-10 m)
- ρ = density (~2700 kg/m³)
- c = specific heat (~900 J/kg/K)
- w = fault zone width (~mm to cm)

For narrow fault zones, this can produce temperature rises of 100s to 1000s of degrees, potentially causing:
- Thermal pressurization of pore fluids
- Flash melting
- Mineral transformations

### Thermal Diffusion Time Scale

Heat diffuses away from the fault over time. The characteristic diffusion time scale is:

t ~ L² / κ

where L is the length scale and κ = k/(ρc) is the thermal diffusivity. For crustal rocks with κ ≈ 1 mm²/s:

- 1 m diffusion: ~30 years
- 10 m diffusion: ~3000 years
- 1 km diffusion: ~30 million years

This example demonstrates these concepts in a simplified 2D geometry.

## References

1. Rice, J.R. (2006). Heating and weakening of faults during earthquake slip. Journal of Geophysical Research, 111, B05311.
2. Lachenbruch, A.H. (1980). Frictional heating, fluid pressure, and the resistance to fault motion. Journal of Geophysical Research, 85, 6097-6112.
