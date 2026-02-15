# Thermoporoelasticity 2D Examples

This example suite demonstrates fully coupled thermal-hydraulic-mechanical (THM) simulations using thermoporoelasticity with and without faults in a 2D domain.

Thermoporoelasticity couples three physical processes:

1. **Mechanics**: Elasticity with pore pressure and thermal effects
2. **Hydraulics**: Fluid flow (Darcy's law) with thermal effects
3. **Thermal**: Heat conduction

## Domain

The domain is a 20 km wide × 10 km deep vertical cross-section:

- -10 km ≤ x ≤ 10 km
- -10 km ≤ y ≤ 0 km (y=0 is the surface)

For fault examples, a vertical fault runs along x=0 from the surface to 7.5 km depth.

```
         surface (y=0, T=300 K, p=0)
     --------------------
     |         |        |
     |         |fault   |
     |         |        |
     |         *        |  <- fault tip at 7.5 km depth (fault examples only)
     |                  |
     --------------------
         bottom (y=-10 km, T=450 K)
```

## Physics

### Thermoporoelasticity Governing Equations

1. **Momentum Balance**:
   ∇·σ + ρ_b g = 0
   
   where σ = σ' - α p I (effective stress with Biot coefficient α)

2. **Fluid Mass Balance**:
   ∂ζ/∂t + ∇·q = Q_f
   
   where q = -k/μ (∇p - ρ_f g) (Darcy flow)

3. **Energy Balance**:
   ρc ∂T/∂t + ∇·q_T = Q_T
   
   where q_T = -κ ∇T (Fourier's law)

### Material Properties

| Property | Value | Unit |
|----------|-------|------|
| Solid density | 2700 | kg/m³ |
| Fluid density | 1000 | kg/m³ |
| Fluid viscosity | 1.0×10⁻³ | Pa·s |
| Porosity | 0.1 | - |
| Permeability | 1.0×10⁻¹⁵ | m² |
| Drained bulk modulus | 10 | GPa |
| Shear modulus | 6 | GPa |
| Biot coefficient | 0.8 | - |
| Biot modulus | 10 | GPa |
| Specific heat | 900 | J/(kg·K) |
| Thermal conductivity | 3.0 | W/(m·K) |
| Solid thermal expansion | 2.4×10⁻⁵ | 1/K |
| Fluid thermal expansion | 2.1×10⁻⁴ | 1/K |
| Reference temperature | 300 | K |

### Boundary Conditions

- **Displacement**: 
  - Left/Right boundaries: Fixed in x
  - Bottom boundary: Fixed in y
  - Top boundary: Free (traction-free)
- **Pressure**: 
  - Top boundary: Zero pressure (drained)
  - Other boundaries: No flow
- **Temperature**: 
  - Surface (y=0): 300 K
  - Bottom (y=-10 km): 450 K
  - Geothermal gradient: 15 K/km

## Simulations

### Step 1: THM Without Fault (`step01_no_fault.cfg`)

This simulation models coupled THM processes without any fault:
- Initial thermal and pressure perturbation
- Coupled thermal-hydraulic-mechanical response
- Heat and fluid pressure diffusion

**To run:**
```bash
pylith step01_no_fault.cfg
```

### Step 2: THM With Fault (`step02_fault.cfg`)

This simulation adds a kinematic fault:
- Vertical strike-slip fault at x=0
- Prescribed slip induces mechanical deformation
- Coupled pressure and temperature response to slip
- Demonstrates fault-THM interaction

**To run:**
```bash
pylith step02_fault.cfg
```

### Step 3: Injection-Induced Response (`step03_injection.cfg`)

This simulation models fluid injection:
- Localized fluid source at depth
- Pressure diffusion and induced deformation
- Thermal effects from injection
- No fault

**To run:**
```bash
pylith step03_injection.cfg
```

## Mesh Generation

The mesh is generated using Gmsh:

```bash
./generate_gmsh.py --write
```

This creates `mesh_tri.msh` with triangular elements.

## Files

| File | Description |
|------|-------------|
| `generate_gmsh.py` | Gmsh mesh generation script |
| `mesh_tri.msh` | Generated mesh file |
| `pylithapp.cfg` | Common simulation parameters |
| `step01_no_fault.cfg` | THM without fault |
| `step02_fault.cfg` | THM with fault |
| `step03_injection.cfg` | Injection-induced response |
| `fault_slip.spatialdb` | Fault slip distribution |
| `mat_thermoporoelastic.spatialdb` | Material properties |
| `initial_temperature.spatialdb` | Initial temperature |
| `initial_pressure.spatialdb` | Initial pressure |

## Visualization

Output is written in HDF5 format to the `output/` directory. Visualize using:

- **ParaView**: Open the `.xmf` files
- **PyLith viz tools**: Use the `pylith_viz` module

Key fields to examine:
- `temperature`: Thermal evolution
- `pressure`: Pore pressure changes
- `displacement`: Deformation
- `cauchy_stress`: Stress state

## Scientific Background

### Coupled THM Processes

In thermoporoelasticity, the three fields (displacement, pressure, temperature) are coupled:

1. **T → σ**: Temperature changes cause thermal stress
2. **T → p**: Fluid thermal expansion affects pore pressure
3. **p → σ**: Pore pressure reduces effective stress
4. **σ → p**: Volumetric strain affects pore pressure (Biot coupling)
5. **Conduction**: Heat diffuses through the medium
6. **Darcy flow**: Fluid flows down pressure gradients

### Diffusion Time Scales

The system has two characteristic diffusion times:

- **Thermal diffusivity**: κ_T = k_T/(ρc) ≈ 1×10⁻⁶ m²/s
- **Hydraulic diffusivity**: κ_H = k/(μS) where S is storage

For crustal rocks, hydraulic diffusion is typically faster than thermal diffusion.

### Applications

THM coupling is important for:
- Geothermal energy extraction
- CO₂ sequestration
- Nuclear waste disposal
- Induced seismicity from injection
- Earthquake-generated thermal and pressure anomalies
