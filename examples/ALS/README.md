# ALS Example - Minimum Working Example for PyLith

This example demonstrates a basic 2D elasticity problem using coordinates
based on the ALS mesh data. The mesh covers approximately 19 km x 21 km
in a projected coordinate system (e.g., UTM).

## Directory Contents

* `generate_gmsh.py` - Python script to generate the mesh using Gmsh
* `pylithapp.cfg` - Common PyLith parameters
* `step01_axialdisp.cfg` - Axial displacement boundary conditions
* `mat_elastic.spatialdb` - Material properties database

## Mesh Generation

Generate the mesh using the `generate_gmsh.py` script:

```bash
# Generate triangle mesh (default)
./generate_gmsh.py --write

# Generate quadrilateral mesh
./generate_gmsh.py --write --cell quad
```

This creates `mesh_tri.msh` or `mesh_quad.msh` depending on the cell type.

## Running the Simulation

```bash
# Run with default (triangle) mesh
pylith step01_axialdisp.cfg

# Run with quad mesh
pylith step01_axialdisp.cfg mesh_quad.cfg
```

## Domain Description

The domain is a 2D rectangular region:
- X: 370,000 m to 389,000 m (19 km)
- Y: 6,527,000 m to 6,548,000 m (21 km)

The coordinates are in a projected coordinate system typical of UTM zones.

### Boundary Conditions

* `boundary_xneg` (-x face): Displacement Ux = 0
* `boundary_xpos` (+x face): Displacement Ux = 1.0 m
* `boundary_yneg` (-y face): Displacement Uy = 0

### Material Properties

Linear elastic material with:
- Density: 2500 kg/m³
- Vp: 5.0 km/s
- Vs: 3.0 km/s
