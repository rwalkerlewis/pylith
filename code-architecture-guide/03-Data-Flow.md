# Data Flow in PyLith

## Table of Contents
1. [Simulation Lifecycle Data Flow](#simulation-lifecycle-data-flow)
2. [Mesh Data Flow](#mesh-data-flow)
3. [Field Data Flow](#field-data-flow)
4. [Assembly Data Flow](#assembly-data-flow)
5. [Solver Data Flow](#solver-data-flow)
6. [Output Data Flow](#output-data-flow)

## Simulation Lifecycle Data Flow

### Complete Flow Diagram

```
┌──────────────────────────────────────────────────────────────────┐
│ 1. STARTUP                                                        │
│    User Config → Pyre → PyLithApp → Component Instantiation      │
└───────────────────────────────┬──────────────────────────────────┘
                                ▼
┌──────────────────────────────────────────────────────────────────┐
│ 2. MESHING                                                        │
│    Mesh File → MeshIO → DMPlex → Distribute → Adjusted Mesh      │
│    (Cubit/ASCII)        (PETSc)   (Parallel)   (+ Cohesive)      │
└───────────────────────────────┬──────────────────────────────────┘
                                ▼
┌──────────────────────────────────────────────────────────────────┐
│ 3. PROBLEM SETUP                                                  │
│    ┌──────────────┐    ┌─────────────┐    ┌─────────────┐       │
│    │ Materials    │ → │ Solution    │ ← │ Spatial DBs │       │
│    │ (Properties) │    │ Field       │    │ (Data)      │       │
│    └──────────────┘    └─────────────┘    └─────────────┘       │
│    ┌──────────────┐    ┌─────────────┐                          │
│    │ BCs          │    │ Auxiliary   │                          │
│    │ (Constraints)│ → │ Fields      │                          │
│    └──────────────┘    └─────────────┘                          │
│    ┌──────────────┐    ┌─────────────┐                          │
│    │ Faults       │    │ Derived     │                          │
│    │ (Slip)       │ → │ Fields      │                          │
│    └──────────────┘    └─────────────┘                          │
└───────────────────────────────┬──────────────────────────────────┘
                                ▼
┌──────────────────────────────────────────────────────────────────┐
│ 4. TIME STEPPING LOOP (for each time step)                       │
│                                                                   │
│    ┌───────────────────────────────────────────┐                │
│    │ Compute Residual F(t, s)                  │                │
│    │  ├─ Material contributions                │                │
│    │  ├─ BC contributions                      │                │
│    │  └─ Fault contributions                   │                │
│    └──────────────┬────────────────────────────┘                │
│                   ▼                                              │
│    ┌───────────────────────────────────────────┐                │
│    │ Compute Jacobian J(t, s) = ∂F/∂s          │                │
│    │  ├─ Material Jacobians                    │                │
│    │  ├─ BC Jacobians                          │                │
│    │  └─ Fault Jacobians                       │                │
│    └──────────────┬────────────────────────────┘                │
│                   ▼                                              │
│    ┌───────────────────────────────────────────┐                │
│    │ Solve J·Δs = -F (PETSc)                   │                │
│    │  ├─ KSP (linear solver)                   │                │
│    │  └─ PC (preconditioner)                   │                │
│    └──────────────┬────────────────────────────┘                │
│                   ▼                                              │
│    ┌───────────────────────────────────────────┐                │
│    │ Update s = s + Δs                         │                │
│    └──────────────┬────────────────────────────┘                │
│                   ▼                                              │
│    ┌───────────────────────────────────────────┐                │
│    │ Post-step Processing                      │                │
│    │  ├─ Update state variables                │                │
│    │  ├─ Compute derived fields                │                │
│    │  └─ Check convergence                     │                │
│    └──────────────┬────────────────────────────┘                │
│                   ▼                                              │
│    ┌───────────────────────────────────────────┐                │
│    │ Output (if triggered)                     │                │
│    │  └─ Write HDF5/VTK files                  │                │
│    └───────────────────────────────────────────┘                │
│                                                                   │
└───────────────────────────────┬──────────────────────────────────┘
                                ▼
┌──────────────────────────────────────────────────────────────────┐
│ 5. FINALIZATION                                                   │
│    Write metadata → Cleanup PETSc → Destroy objects              │
└──────────────────────────────────────────────────────────────────┘
```

## Mesh Data Flow

### Mesh Import and Processing

```
┌─────────────────────────────────────────────────────────────────┐
│ Mesh Source                                                      │
│  • Cubit/Trelis (EXODUS II format)                              │
│  • ASCII format                                                  │
│  • Gmsh                                                          │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ MeshIO::read()                                                   │
│  • Read vertices (coordinates)                                   │
│  • Read cells (connectivity)                                     │
│  • Read groups (labels)                                          │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Create PETSc DMPlex                                              │
│  • DMPlexCreateFromCellListPetsc()                               │
│  • Store vertices and cell connectivity                          │
│  • Build mesh topology (edges, faces)                            │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Set Coordinate System                                            │
│  • Proj library for geographic coordinates                       │
│  • Store in DMPlex aux data                                      │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Adjust for Faults (if present)                                   │
│  • Identify fault vertices/edges                                 │
│  • Insert cohesive cells                                         │
│  • Duplicate vertices on fault                                   │
│  • Create negative/positive side labels                          │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Parallel Distribution                                            │
│  • DMPlexDistribute()                                            │
│  • Partition mesh (ParMetis/Chaco)                               │
│  • Migrate cells to processors                                   │
│  • Create overlap (ghost cells)                                  │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Final Mesh Topology                                              │
│  • Vertices, edges, faces, cells                                 │
│  • Labels for materials, BCs, faults                             │
│  • Coordinate field                                              │
│  • Ready for field creation                                      │
└─────────────────────────────────────────────────────────────────┘
```

### Mesh Data Structures

**DMPlex Structure** (PETSc data structure):
```
Mesh Points (0-based indexing):
├─ Cells: [0, numCells)
├─ Vertices: [numCells, numCells+numVertices)
├─ Edges: [numCells+numVertices, numCells+numVertices+numEdges)
└─ Faces: [numCells+numVertices+numEdges, ...)

Labels (for grouping):
├─ "material-id": Cell groups by material
├─ "boundary_xpos": Vertices/edges on +x boundary
├─ "fault": Cohesive cells for fault
└─ "fault_edge": Edges bounding fault
```

## Field Data Flow

### Field Creation and Population

```
┌─────────────────────────────────────────────────────────────────┐
│ Field Definition (Python/C++)                                    │
│  field.addSubfield("displacement", 3)                            │
│  field.addSubfield("velocity", 3)                                │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Set Discretization                                               │
│  • Basis order (0=constant, 1=linear, 2=quadratic)              │
│  • Quadrature order                                              │
│  • DMSetField() - attach to DMPlex                               │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Create Section (Layout)                                          │
│  • DMCreateDS() - create discrete system                         │
│  • Compute DOF for each mesh point                               │
│  • PetscSection: describes data layout                           │
│  │   Point 0: offset=0, dof=3                                    │
│  │   Point 1: offset=3, dof=3                                    │
│  │   ...                                                          │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Allocate Storage                                                 │
│  • DMCreateGlobalVector() - global (distributed) vector          │
│  • DMCreateLocalVector() - local (with ghost) vector             │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Populate Values                                                  │
│  • From spatial database: query at coordinates                   │
│  • From function: evaluate user function                         │
│  • From file: read from HDF5/VTK                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Field Data Access Pattern

**Global Vector** (distributed):
```python
# Each processor owns part of the global vector
Proc 0: [u_0, u_1, u_2, ...]  # Owned points
Proc 1: [u_k, u_k+1, ...]     # Owned points
```

**Local Vector** (with ghosts):
```python
# Each processor has owned + ghost points
Proc 0: [u_0, u_1, u_2, ..., u_ghost_1, u_ghost_2, ...]
        [----owned-------]  [------ghosts------]
```

**Scatter Operation**:
```cpp
// Global to Local (before computation)
DMGlobalToLocalBegin(dm, globalVec, INSERT_VALUES, localVec);
DMGlobalToLocalEnd(dm, globalVec, INSERT_VALUES, localVec);

// Local to Global (after computation)
DMLocalToGlobalBegin(dm, localVec, ADD_VALUES, globalVec);
DMLocalToGlobalEnd(dm, localVec, ADD_VALUES, globalVec);
```

### Auxiliary Field Population from Spatial Database

```
┌─────────────────────────────────────────────────────────────────┐
│ Spatial Database (spatialdata library)                          │
│  • SimpleDB: ASCII file with coordinates and values              │
│  • UniformDB: Uniform values                                     │
│  • SimpleGridDB: Structured grid                                 │
│  • SCEC CVM: 3D seismic velocity models                          │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Query Loop (for each mesh point)                                │
│  FOR each point p in mesh:                                       │
│    1. Get coordinates: x = mesh.getCoordinates(p)                │
│    2. Query database: values = db.query(x)                       │
│    3. Store in auxiliary field: auxField.set(p, values)          │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Auxiliary Field (Local Vector)                                   │
│  [ρ_0, μ_0, K_0, ρ_1, μ_1, K_1, ...]                            │
│  Ready for use in FE kernels                                     │
└─────────────────────────────────────────────────────────────────┘
```

## Assembly Data Flow

### Residual Assembly

```
┌─────────────────────────────────────────────────────────────────┐
│ Input: Solution vector s, time t                                 │
│ Output: Residual vector F                                        │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Zero Residual Vector                                             │
│  VecZeroEntries(residual)                                        │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ FOR each Integrator (Material, BC, Fault):                       │
│                                                                   │
│  ┌──────────────────────────────────────────────┐               │
│  │ Get label and value for integration domain   │               │
│  └──────────────┬───────────────────────────────┘               │
│                 ▼                                                │
│  ┌──────────────────────────────────────────────┐               │
│  │ FOR each cell c in domain:                   │               │
│  │                                               │               │
│  │  ┌────────────────────────────────────┐      │               │
│  │  │ 1. Get cell geometry (Jacobian)    │      │               │
│  │  └────────────┬───────────────────────┘      │               │
│  │               ▼                               │               │
│  │  ┌────────────────────────────────────┐      │               │
│  │  │ 2. Get solution values at cell     │      │               │
│  │  │    (via PetscSection)              │      │               │
│  │  └────────────┬───────────────────────┘      │               │
│  │               ▼                               │               │
│  │  ┌────────────────────────────────────┐      │               │
│  │  │ 3. Get auxiliary field values      │      │               │
│  │  │    (material properties)           │      │               │
│  │  └────────────┬───────────────────────┘      │               │
│  │               ▼                               │               │
│  │  ┌────────────────────────────────────┐      │               │
│  │  │ 4. FOR each quadrature point q:    │      │               │
│  │  │     a. Evaluate basis functions     │      │               │
│  │  │     b. Compute u, ∇u at q          │      │               │
│  │  │     c. Call kernel: f0, f1 = kernel│      │               │
│  │  │     d. Integrate:                   │      │               │
│  │  │        ∫ f0·v + f1·∇v dV           │      │               │
│  │  └────────────┬───────────────────────┘      │               │
│  │               ▼                               │               │
│  │  ┌────────────────────────────────────┐      │               │
│  │  │ 5. Add to residual vector          │      │               │
│  │  │    (via PetscSection offsets)      │      │               │
│  │  └────────────────────────────────────┘      │               │
│  │                                               │               │
│  └───────────────────────────────────────────────┘               │
│                                                                   │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Assembly Complete                                                │
│  • Local contributions summed                                    │
│  • Parallel reduction (if needed)                                │
│  • Residual vector ready                                         │
└─────────────────────────────────────────────────────────────────┘
```

### Kernel Data Flow (Pointwise)

```
Input to Kernel at Quadrature Point:
┌─────────────────────────────────────────────────────────────────┐
│ • dim: Spatial dimension (2 or 3)                                │
│ • t: Current time                                                │
│ • x[]: Coordinates [x, y, z]                                     │
│ • numS: Number of solution fields                                │
│ • numA: Number of auxiliary fields                               │
│ • sOff[]: Offsets into solution array                            │
│ • aOff[]: Offsets into auxiliary array                           │
│ • s[]: Solution values [u_x, u_y, u_z, v_x, v_y, v_z, ...]      │
│ • s_t[]: Time derivatives [u̇_x, u̇_y, u̇_z, ...]                 │
│ • s_x[]: Gradients [∂u_x/∂x, ∂u_x/∂y, ..., ∂u_y/∂x, ...]        │
│ • a[]: Auxiliary values [ρ, μ, K, ...]                           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼ Kernel Computation
┌─────────────────────────────────────────────────────────────────┐
│ Example: Elasticity f1 Kernel                                    │
│                                                                   │
│ 1. Extract material properties:                                  │
│    μ = a[aOff[0]]                                                │
│    K = a[aOff[1]]                                                │
│                                                                   │
│ 2. Compute strain from gradient:                                 │
│    ε = 0.5(∇u + ∇u^T)                                            │
│                                                                   │
│ 3. Compute stress (constitutive law):                            │
│    σ = λ(tr ε)I + 2μ dev(ε)                                      │
│    where λ = K - 2μ/3                                            │
│                                                                   │
│ 4. Store stress in f1[]:                                         │
│    f1[0..8] = [σ_xx, σ_xy, σ_xz, σ_yx, σ_yy, ...]               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Output: f0[], f1[] arrays                                        │
│ • Used by PETSc to compute ∫ f0·v + f1·∇v dV                     │
└─────────────────────────────────────────────────────────────────┘
```

## Solver Data Flow

### PETSc Integration

```
┌─────────────────────────────────────────────────────────────────┐
│ Problem Setup                                                     │
│  ├─ Create TS (time stepper)                                     │
│  ├─ Set IFunction (residual)                                     │
│  ├─ Set IJacobian (Jacobian)                                     │
│  └─ Create SNES (nonlinear solver)                               │
│      └─ Create KSP (linear solver)                               │
│          └─ Create PC (preconditioner)                           │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ TSSolve() [PETSc Time Stepping]                                 │
│                                                                   │
│  FOR each time step:                                             │
│                                                                   │
│   ┌──────────────────────────────────────────┐                  │
│   │ SNESSolve() [Nonlinear Solver]           │                  │
│   │                                           │                  │
│   │  WHILE not converged:                    │                  │
│   │                                           │                  │
│   │   ┌──────────────────────────────────┐   │                  │
│   │   │ Call IFunction                   │   │                  │
│   │   │  → Problem::computeRHSResidual() │   │                  │
│   │   │  → Problem::computeLHSResidual() │   │                  │
│   │   │  → F = LHS - RHS                 │   │                  │
│   │   └──────────────────────────────────┘   │                  │
│   │                                           │                  │
│   │   ┌──────────────────────────────────┐   │                  │
│   │   │ Call IJacobian                   │   │                  │
│   │   │  → Problem::computeLHSJacobian() │   │                  │
│   │   │  → J = ∂F/∂s                     │   │                  │
│   │   └──────────────────────────────────┘   │                  │
│   │                                           │                  │
│   │   ┌──────────────────────────────────┐   │                  │
│   │   │ KSPSolve() [Linear Solver]       │   │                  │
│   │   │  Solve: J·Δs = -F                │   │                  │
│   │   │                                   │   │                  │
│   │   │  ┌───────────────────────────┐   │   │                  │
│   │   │  │ PCApply() [Preconditioner]│   │   │                  │
│   │   │  │  • ILU, GAMG, FieldSplit  │   │   │                  │
│   │   │  └───────────────────────────┘   │   │                  │
│   │   └──────────────────────────────────┘   │                  │
│   │                                           │                  │
│   │   s = s + Δs                             │                  │
│   │                                           │                  │
│   └───────────────────────────────────────────┘                  │
│                                                                   │
│   ┌──────────────────────────────────────────┐                  │
│   │ TSPostStep()                             │                  │
│   │  → Problem::poststep()                   │                  │
│   │  → Update state variables                │                  │
│   └──────────────────────────────────────────┘                  │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
```

### Matrix Assembly (Jacobian)

Similar to residual, but computes derivatives:

```
Jacobian: J_ij = ∂F_i/∂s_j

FOR each cell:
  FOR each quadrature point:
    Call Jacobian kernel: J0, J1, J2, J3
    Compute: ∫ (J0·v·u + J1·v·∇u + J2·∇v·u + J3·∇v·∇u) dV
    Add to global Jacobian matrix
```

## Output Data Flow

### Solution Output Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│ OutputTrigger Evaluation                                         │
│  • Time-based: every N seconds                                   │
│  • Step-based: every N steps                                     │
│  IF trigger activated → proceed to output                        │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Prepare Data                                                     │
│  • Solution field (displacement, velocity, pressure)             │
│  • Derived fields (stress, strain)                               │
│  • Auxiliary fields (material properties)                        │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Select Subfields to Output                                       │
│  data_fields = [displacement, cauchy_stress, density]            │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ DataWriter (HDF5 or VTK)                                         │
│                                                                   │
│  ┌──────────────────────────────────────────┐                   │
│  │ HDF5 Format (parallel-friendly)          │                   │
│  │  ├─ /topology/cells                      │                   │
│  │  ├─ /topology/vertices                   │                   │
│  │  ├─ /geometry/vertices                   │                   │
│  │  ├─ /vertex_fields/displacement          │                   │
│  │  └─ /cell_fields/cauchy_stress           │                   │
│  │  Companion .xmf (Xdmf) for ParaView       │                   │
│  └──────────────────────────────────────────┘                   │
│                                                                   │
│  ┌──────────────────────────────────────────┐                   │
│  │ VTK Format (serial, human-readable)      │                   │
│  │  • Legacy or XML format                   │                   │
│  │  • Unstructured grid                      │                   │
│  │  • Point data and cell data               │                   │
│  └──────────────────────────────────────────┘                   │
└────────────────┬────────────────────────────────────────────────┘
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│ Write to Disk                                                    │
│  • One file per time step (or time history in single file)      │
│  • Parallel I/O (MPI-IO) for HDF5                                │
│  • Master process for VTK                                        │
└─────────────────────────────────────────────────────────────────┘
```

### Observer Pattern for Output

```python
# Output configuration in Python
[pylithapp.problem]
solution_observers = [domain, boundary, points]

[pylithapp.problem.solution_observers.domain]
data_fields = [displacement, velocity]
trigger = pylith.meshio.OutputTriggerStep
trigger.num_skip = 1  # Output every step

[pylithapp.problem.solution_observers.points]
data_fields = [displacement]
points = [[0, 0, 0], [10, 0, 0]]  # Monitor specific points
```

**Flow**:
```
Problem::poststep()
  → ObserversSoln::notifyObservers()
    → FOR each Observer:
        IF trigger.shouldWrite(t):
          observer.update(t, solution)
            → DataWriter::write()
```

## Summary

Data flow in PyLith follows a clear pipeline:

1. **Configuration → Components**: User config creates component hierarchy
2. **Files → Mesh**: Mesh files loaded into DMPlex
3. **Spatial DBs → Fields**: Material properties queried and stored
4. **Fields → Kernels**: Data passed to pointwise functions
5. **Kernels → Assembly**: Pointwise results integrated
6. **Assembly → Vectors/Matrices**: Residual and Jacobian constructed
7. **PETSc → Solution**: Solvers compute solution update
8. **Solution → Output**: Results written to files

Key principles:
- **Locality**: Most operations work on local data (with ghosts)
- **Parallel**: MPI communication via PETSc abstractions
- **Streaming**: Large data processed incrementally
- **Abstraction**: PETSc handles low-level details

Next: [04-Build-System.md](04-Build-System.md) - Compilation and build process
