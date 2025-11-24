# PyLith Quick Reference

A concise reference for common development tasks and concepts.

## Directory Structure Quick Reference

```
pylith/
├── libsrc/pylith/         # C++ core (performance-critical code)
│   ├── bc/                # Boundary conditions
│   ├── faults/            # Fault implementations
│   ├── feassemble/        # Finite element assembly
│   ├── fekernels/         # Pointwise kernel functions
│   ├── materials/         # Material models
│   ├── meshio/            # Mesh I/O and output
│   ├── problems/          # Problem formulations
│   ├── topology/          # Mesh topology (DMPlex wrapper)
│   └── utils/             # Utilities
│
├── modulesrc/             # SWIG interfaces (C++ ↔ Python bridge)
│   └── */                 # Mirrors libsrc structure
│
├── pylith/                # Python layer (configuration, user interface)
│   ├── apps/              # Application entry points
│   └── */                 # Mirrors libsrc structure
│
├── tests/                 # Test suite
│   ├── libtests/          # C++ unit tests (Catch2)
│   ├── pytests/           # Python unit tests
│   ├── mmstests/          # Method of manufactured solutions
│   └── fullscale/         # Full simulation tests
│
└── examples/              # Example problems
```

## Common Commands

### Building

```bash
# Configure
./configure [OPTIONS]

# Common options
./configure \
    --enable-swig \            # Regenerate Python bindings
    --enable-testing \         # Enable tests
    --enable-hdf5 \           # HDF5 output
    --with-petsc-dir=PATH \   # PETSc location
    PYTHON=python3            # Python interpreter

# Build
make -j$(nproc)               # Parallel build

# Install (optional)
make install

# Run tests
make check
```

### Running PyLith

```bash
# Basic run
pylith simulation.cfg

# Multiple config files (later overrides earlier)
pylith base.cfg problem.cfg output.cfg

# Command-line overrides
pylith sim.cfg --problem.solver=nonlinear

# Get help
pylith --help
pylith --problem.help              # Help for problem component
pylith --problem.help-components   # List sub-components

# Debugging
pylith sim.cfg --petsc.log_view                    # Performance log
pylith sim.cfg --petsc.snes_monitor                # Nonlinear solver
pylith sim.cfg --petsc.ksp_monitor                 # Linear solver
pylith sim.cfg --petsc.snes_test_jacobian          # Check Jacobian
```

## Key Classes

### Problem Hierarchy

```
Problem (abstract)
└── TimeDependent (quasi-static, dynamic)
└── GreensFns (Green's functions)
```

### Material Hierarchy

```
Physics (abstract)
└── Material (abstract)
    ├── Elasticity
    │   ├── IsotropicLinearElasticity
    │   ├── IsotropicLinearMaxwell (viscoelastic)
    │   └── IsotropicPowerLaw
    ├── IncompressibleElasticity
    │   └── IsotropicLinearIncompElasticity
    └── Poroelasticity
        └── IsotropicLinearPoroelasticity
```

### Boundary Conditions

```
BoundaryCondition (abstract)
├── DirichletTimeDependent (prescribed displacement/velocity)
├── NeumannTimeDependent (prescribed traction)
└── AbsorbingDampers (absorbing BC)
```

### Faults

```
FaultCohesive (abstract)
├── FaultCohesiveKin (prescribed slip)
└── FaultCohesiveImpulses (Green's functions)

KinSrc (slip time function)
├── KinSrcStep
├── KinSrcRamp
├── KinSrcBrune
└── KinSrcLiuCos
```

## Component Lifecycle

```python
# 1. Construction
obj = Component()

# 2. Configuration (_configure)
# Pyre reads .cfg files and sets properties

# 3. Pre-initialization (preinitialize)
obj.preinitialize(problem)

# 4. Verification (verifyConfiguration)
obj.verifyConfiguration()

# 5. Initialization (initialize)
obj.initialize()

# 6. Run (run, poststep)
obj.run()

# 7. Finalization (finalize)
obj.finalize()

# 8. Destruction (deallocate)
obj.deallocate()
```

## Finite Element Assembly

### Weak Form

```
Find u such that:
∫_Ω (f0·v + f1·∇v) dV = 0  for all test functions v
```

### Kernels

```cpp
// Residual kernels
void f0(/* ... */, PylithScalar f0[]) {
    // f0 = term integrated with test function
}

void f1(/* ... */, PylithScalar f1[]) {
    // f1 = term integrated with test function gradient
    // Example: f1 = stress for elasticity
}

// Jacobian kernels (J = ∂F/∂u)
void Jf0(/* ... */, PylithScalar J[]);  // ∂f0/∂u
void Jf1(/* ... */, PylithScalar J[]);  // ∂f0/∂(∇u)
void Jf2(/* ... */, PylithScalar J[]);  // ∂f1/∂u
void Jf3(/* ... */, PylithScalar J[]);  // ∂f1/∂(∇u)
```

### Kernel Signature

```cpp
typedef void (*PetscPointFn)(
    PetscInt dim,              // Spatial dimension
    PetscInt Nf,               // Number of solution fields
    PetscInt NfAux,            // Number of auxiliary fields
    const PetscInt uOff[],     // Solution offsets
    const PetscInt uOff_x[],   // Solution gradient offsets
    const PetscScalar u[],     // Solution values
    const PetscScalar u_t[],   // Time derivatives
    const PetscScalar u_x[],   // Gradients
    const PetscInt aOff[],     // Auxiliary offsets
    const PetscInt aOff_x[],   // Auxiliary gradient offsets
    const PetscScalar a[],     // Auxiliary values
    const PetscScalar a_t[],   // Auxiliary time derivatives
    const PetscScalar a_x[],   // Auxiliary gradients
    PetscReal t,               // Time
    const PetscScalar x[],     // Coordinates
    PetscInt numConstants,     // Number of constants
    const PetscScalar constants[], // Constants
    PetscScalar f[]            // Output
);
```

## PETSc Integration

### Data Structures

```cpp
DM        // Distributed mesh (DMPlex)
Vec       // Vector (global, local)
Mat       // Matrix (sparse)
PetscSection  // DOF layout

TS        // Time stepper
SNES      // Nonlinear solver
KSP       // Linear solver
PC        // Preconditioner
```

### Common PETSc Options

```cfg
[pylithapp.petsc]
# Time stepping
ts_type = beuler           # Backward Euler
ts_monitor = true          # Monitor time steps

# Nonlinear solver
snes_monitor = true        # Monitor nonlinear iterations
snes_converged_reason = true
snes_rtol = 1.0e-10        # Relative tolerance

# Linear solver
ksp_type = gmres           # GMRES
ksp_monitor = true
ksp_converged_reason = true
ksp_rtol = 1.0e-10

# Preconditioner
pc_type = ilu              # ILU for serial
# pc_type = asm            # Additive Schwarz for parallel
# pc_type = gamg           # Algebraic multigrid

# Field split (for multi-physics)
pc_type = fieldsplit
pc_fieldsplit_type = schur
```

## Configuration File Syntax

### Basic Syntax

```cfg
# Comments
[section]
property = value
```

### Component Hierarchy

```cfg
[pylithapp]
# Top level

[pylithapp.problem]
# Problem component

[pylithapp.problem.materials]
# Array of materials
slab = pylith.materials.Elasticity

[pylithapp.problem.materials.slab]
# Individual material

[pylithapp.problem.materials.slab.rheology]
# Sub-component (facility)
```

### Units

```cfg
# Pyre handles unit conversion
length = 5.0*km           # Converted to meters
time = 100.0*year         # Converted to seconds
stress = 10.0*MPa         # Converted to Pascals

# Available units: m, km, cm, s, year, Pa, MPa, GPa, kg, etc.
```

## Spatial Databases

### SimpleDB (ASCII)

```
#SPATIAL.ascii 1
SimpleDB {
  num-values = 3
  value-names = density vs vp
  value-units = kg/m**3 m/s m/s
  num-locs = 4
  data-dim = 1
  space-dim = 3
  cs-data = cartesian {
    to-meters = 1.0
    space-dim = 3
  }
}
// Columns: x y z density vs vp
0.0  0.0  0.0   2500.0  3000.0  5200.0
1.0  0.0  0.0   2500.0  3000.0  5200.0
0.0  1.0  0.0   2500.0  3000.0  5200.0
0.0  0.0  1.0   2500.0  3000.0  5200.0
```

### UniformDB (Python config)

```cfg
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.description = Uniform properties
db_auxiliary_field.values = [density, vs, vp]
db_auxiliary_field.data = [2500.0*kg/m**3, 3000.0*m/s, 5200.0*m/s]
```

## Output Configuration

```cfg
[pylithapp.problem]
solution_observers = [domain, boundary, points]

[pylithapp.problem.solution_observers.domain]
writer = pylith.meshio.DataWriterHDF5
writer.filename = output/solution-domain.h5
data_fields = [displacement, velocity]

trigger = pylith.meshio.OutputTriggerStep
trigger.num_skip = 1  # Output every step

[pylithapp.problem.solution_observers.points]
writer.filename = output/solution-points.h5
data_fields = [displacement]
reader = pylith.meshio.PointsList
reader.filename = output_points.txt
```

## Debugging Tips

### Check Residual/Jacobian

```bash
pylith sim.cfg --petsc.snes_test_jacobian --petsc.snes_test_jacobian_view
# Should show Jacobian errors ~ 1e-10 or smaller
```

### Increase Verbosity

```cfg
[pylithapp]
# Global settings
timedependent.progress_monitor = pylith.problems.ProgressMonitorTime

[pylithapp.timedependent.progress_monitor]
t_units = year
update_percent = 5.0

[pylithapp.journal.info]
timedependent = 1
solution = 1
problem = 1
```

### Use Python Debugger

```cfg
[pylithapp]
start_python_debugger = True  # Drops into pdb
```

### GDB (C++ debugging)

```bash
gdb --args python $(which pylith) sim.cfg
(gdb) break pylith::problems::Problem::initialize
(gdb) run
```

### Valgrind (Memory leaks)

```bash
valgrind --leak-check=full --track-origins=yes \
    python $(which pylith) sim.cfg
```

## Common Errors and Solutions

### Error: "PETSc not found"

```bash
# Solution: Set PETSC_DIR and PETSC_ARCH
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-name
./configure
```

### Error: "Jacobian test failed"

```bash
# Solution: Check kernel implementations
# Use simpler problem to isolate issue
pylith sim.cfg --petsc.snes_test_jacobian --petsc.snes_test_jacobian_view
```

### Error: "Linear solve diverged"

```cfg
# Solution: Try different preconditioner
[pylithapp.petsc]
pc_type = gamg
# Or use direct solver (small problems only)
# ksp_type = preonly
# pc_type = lu
```

### Error: "Memory allocation failed"

```bash
# Solution: Increase memory or use fewer processors
# Check matrix preallocation
pylith sim.cfg --petsc.log_view  # Look for "malloc" warnings
```

## Performance Optimization

### Profiling

```bash
pylith sim.cfg --petsc.log_view > performance.log
# Look for:
# - Time in assembly vs solve
# - Memory usage
# - Communication time
```

### Solver Tuning

```cfg
[pylithapp.petsc]
# For linear problems
solver = linear

# Field split for multi-physics
pc_type = fieldsplit
pc_fieldsplit_type = schur
pc_fieldsplit_schur_fact_type = full
pc_fieldsplit_schur_precondition = selfp

# Increase Krylov subspace
ksp_gmres_restart = 200
```

### Parallel Scaling

```bash
# Test weak/strong scaling
mpirun -n 1 pylith sim.cfg
mpirun -n 2 pylith sim.cfg
mpirun -n 4 pylith sim.cfg
mpirun -n 8 pylith sim.cfg
```

## Useful Links

- **Documentation**: https://pylith.readthedocs.io
- **Source Code**: https://github.com/geodynamics/pylith
- **Forum**: https://community.geodynamics.org/c/pylith/
- **PETSc Documentation**: https://petsc.org/release/docs/
- **Issue Tracker**: https://github.com/geodynamics/pylith/issues

## File Naming Conventions

### C++
- Headers: `ClassName.hh`
- Implementation: `ClassName.cc`
- Inline: `ClassName.icc`
- Forward declarations: `modulefwd.hh`

### Python
- Modules: `ClassName.py`
- Tests: `TestClassName.py`

### SWIG
- Interfaces: `ClassName.i`

### Configuration
- Main: `pylithapp.cfg`
- Step-specific: `step01_problem.cfg`

### Data
- Spatial databases: `*.spatialdb`
- Meshes: `*.msh` (Gmsh), `*.exo` (Cubit)

## Quick Checklist for New Material

- [ ] C++ header (`MyMaterial.hh`)
- [ ] C++ implementation (`MyMaterial.cc`)
- [ ] Kernel functions (`fekernels/MyMaterial.hh/cc`)
- [ ] SWIG interface (`modulesrc/materials/MyMaterial.i`)
- [ ] Python wrapper (`pylith/materials/MyMaterial.py`)
- [ ] Update `Makefile.am` files
- [ ] C++ unit tests (`tests/libtests/materials/TestMyMaterial.cc`)
- [ ] Python unit tests (`tests/pytests/materials/TestMyMaterial.py`)
- [ ] Full-scale test (`tests/fullscale/MyMaterial/`)
- [ ] Documentation (docstrings, comments)
- [ ] Example configuration

---

**Version**: PyLith 5.0.0dev  
**Last Updated**: 2025-11-24
