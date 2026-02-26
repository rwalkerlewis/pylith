# ALS Paleoseismic Example

This example demonstrates a 2D elasticity problem with a fault, using
coordinates and material properties derived from ALS geospatial survey
data (EPSG 28351 – GDA94 / MGA zone 51).

The domain covers approximately 19 km × 21 km.  A NE–SW-trending fault
trace, digitised from the ALS fault layers, splits the domain into a
western and an eastern material block.

## Directory Contents

* `generate_gmsh.py` – Python script to generate the Gmsh mesh
* `pylithapp.cfg` – Common PyLith parameters (mesh, materials, fault, BCs)
* `step01_axialdisp.cfg` – Axial extension with zero fault slip
* `step02_coseismic.cfg` – Coseismic left-lateral fault slip (1 m)
* `mesh_quad.cfg` – Override to use quad mesh
* `20260206_190857_pylith_inputs/` – Original ALS-generated spatial-database
  files and configuration template
* `20260206_190857_pylith_model.geo` – Original GIS-generated Gmsh `.geo` file
* `run_metadata.json` – Metadata from the ALS processing run

## Mesh Generation

```bash
# Generate triangle mesh (default)
python generate_gmsh.py --write --cell tri

# Generate quadrilateral mesh
python generate_gmsh.py --write --cell quad
```

This creates `mesh_tri.msh` or `mesh_quad.msh`.

## Running the Simulations

```bash
# Step 01 – Axial extension (triangle mesh)
pylith step01_axialdisp.cfg

# Step 02 – Coseismic fault slip (triangle mesh)
pylith step02_coseismic.cfg

# Using the quad mesh
pylith step02_coseismic.cfg mesh_quad.cfg
```

## Domain Description

The domain is a 2D rectangular region:

* **X**: 370 000 m to 389 000 m (19 km)
* **Y**: 6 527 000 m to 6 548 000 m (21 km)

Coordinates are in EPSG 28351 (GDA94 / MGA zone 51).

### Fault

A simplified fault trace runs from (388 000, 6 527 000) on the south
boundary to (375 500, 6 548 000) on the north boundary, broadly
following the main NE–SW lineament mapped in the ALS survey data.

### Boundary Conditions

| Boundary | Label | Tag | Step 01 (axial disp) | Step 02 (coseismic) |
|----------|-------|-----|----------------------|---------------------|
| West  (`-x`) | `boundary_xneg` | 10 | Ux = 0, Uy = 0 | Ux = 0, Uy = 0 |
| East  (`+x`) | `boundary_xpos` | 11 | Ux = +1 m       | Ux = 0, Uy = 0 |
| South (`-y`) | `boundary_yneg` | 12 | Uy = 0          | Ux = 0, Uy = 0 |
| North (`+y`) | `boundary_ypos` | 13 | (free)          | Ux = 0, Uy = 0 |

### Material Properties

Both material blocks use uniform elastic properties based on
Christensen & Mooney (1995) average continental-crust values:

| Property | Value |
|----------|-------|
| Density  | 2 700 kg/m³ |
| Vp       | 6 000 m/s |
| Vs       | 3 464 m/s |

### Geospatial Data

The `20260206_190857_pylith_inputs/` directory contains the original
ALS-generated spatial-database files (materials, faults, boundaries,
gravity) and a manifest describing each file.  These files informed
the material properties, fault geometry, and boundary conditions used
in this example.
