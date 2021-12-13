# Examples: 2-D black oil

This suite of examples serves as a testbed for the black oil multiphase
poroelasticity functionality in PyLith. The mesh from the 2D box example is used,
with some modifications.

## Step01: Axial extension with Dirichlet boundary conditions

Axial extension with Dirichlet boundary conditions on the +x, -x, and
-y boundaries. Features used in this simulation include:

* Static simulation
* UniformDB spatial database for specifying values for properties and
  boundary conditions

The simulation parameters are in the `pylithapp.cfg` and
`step01_axialdisp.cfg` files.

To run the example:
```
pylith step01_axialdisp.cfg
```
