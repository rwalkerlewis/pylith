# Architecture Overview

## Table of Contents
1. [High-Level Architecture](#high-level-architecture)
2. [Design Philosophy](#design-philosophy)
3. [Layer Architecture](#layer-architecture)
4. [Key Subsystems](#key-subsystems)
5. [Execution Flow](#execution-flow)

## High-Level Architecture

PyLith follows a **hybrid architecture** combining C++ for computational performance with Python for configuration, scripting, and user interaction. This design provides both efficiency and flexibility.

```
┌─────────────────────────────────────────────────────────────┐
│                    User Interface Layer                      │
│                      (Python + Pyre)                         │
│  - Configuration files (.cfg, .json)                         │
│  - Command-line interface                                    │
│  - Parameter management                                      │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                 Python Application Layer                     │
│                    (pylith/*.py)                             │
│  - PyLithApp: Main application                               │
│  - Problem setup and configuration                           │
│  - Component factories                                       │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼ SWIG Bindings
┌─────────────────────────────────────────────────────────────┐
│                  C++ Core Implementation                     │
│                   (libsrc/pylith/)                           │
│  ┌───────────────┬──────────────┬──────────────┐           │
│  │   Problem     │  Materials   │  Faults      │           │
│  │  Management   │  & Rheology  │              │           │
│  └───────────────┴──────────────┴──────────────┘           │
│  ┌───────────────┬──────────────┬──────────────┐           │
│  │ FE Assembly   │   Topology   │   MeshIO     │           │
│  │               │              │              │           │
│  └───────────────┴──────────────┴──────────────┘           │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              External Libraries Layer                        │
│  ┌──────────┬─────────┬─────────┬──────────┐               │
│  │  PETSc   │   MPI   │  HDF5   │  NetCDF  │               │
│  │ (Solvers)│(Parallel)│(Output) │ (Cubit)  │               │
│  └──────────┴─────────┴─────────┴──────────┘               │
└─────────────────────────────────────────────────────────────┘
```

## Design Philosophy

### 1. Separation of Concerns

PyLith separates different aspects of the simulation:

- **Physics**: Material behavior, fault mechanics (independent of discretization)
- **Discretization**: Finite element assembly, mesh management
- **Solution Strategy**: Time stepping, linear/nonlinear solvers (delegated to PETSc)
- **I/O**: Mesh input, solution output (separate from computation)

### 2. Component-Based Architecture

Uses the **Pyre/Pythia framework** for component management:

```python
# Example: Problem is composed of multiple components
class PyLithApp:
    mesher = MeshImporter()      # Component for mesh generation
    problem = TimeDependent()     # Component for problem type
    
    class Problem:
        materials = [Material1(), Material2()]  # Component array
        bc = [BC1(), BC2()]                     # Component array
        interfaces = [Fault1()]                 # Component array
```

Each component:
- Is independently configurable
- Has a well-defined interface
- Can be replaced with custom implementations
- Manages its own initialization and cleanup

### 3. Template Method Pattern

Core algorithms define the skeleton while allowing customization:

```cpp
// Abstract base class defines algorithm structure
class Problem {
public:
    void run() {
        // Template method
        initialize();
        while (!finished()) {
            computeResidual();    // Virtual - customizable
            computeJacobian();    // Virtual - customizable
            solve();
            updateState();        // Virtual - customizable
        }
        finalize();
    }
protected:
    virtual void computeResidual() = 0;
    virtual void computeJacobian() = 0;
};
```

### 4. Factory Pattern

Components are created through factories, enabling:
- Runtime selection of implementations
- Easy addition of new component types
- Configuration-driven instantiation

```python
# Factory functions in Problem.py
def materialFactory(name):
    return facility(name, family="material", factory=Elasticity)

def bcFactory(name):
    return facility(name, family="boundary_condition", 
                   factory=DirichletTimeDependent)
```

## Layer Architecture

### Layer 1: User Interface (Python + Configuration)

**Purpose**: User-facing configuration and control

**Key Files**:
- `pylith/apps/PyLithApp.py` - Main application class
- Configuration files (`.cfg`) - User parameters
- JSON parameter dumps - Runtime configuration

**Responsibilities**:
- Parse command-line arguments
- Read configuration files
- Initialize Pyre component system
- Setup logging and metadata

### Layer 2: Python Application Layer

**Purpose**: High-level problem setup and orchestration

**Key Directories**:
- `pylith/problems/` - Problem definitions
- `pylith/materials/` - Material wrappers
- `pylith/bc/` - Boundary condition wrappers
- `pylith/faults/` - Fault wrappers

**Responsibilities**:
- Component lifecycle management (preinitialize, initialize, run, finalize)
- Parameter validation
- Factory instantiation
- Python-level data structures

### Layer 3: C++ Core (Computational Engine)

**Purpose**: High-performance numerical computation

**Key Directories**:
- `libsrc/pylith/problems/` - Problem core
- `libsrc/pylith/feassemble/` - FE assembly
- `libsrc/pylith/materials/` - Material implementations
- `libsrc/pylith/topology/` - Mesh and field management

**Responsibilities**:
- Finite element assembly
- Residual and Jacobian computation
- State variable management
- PETSc integration

### Layer 4: External Libraries

**Purpose**: Foundational computational services

**Key Libraries**:
- **PETSc**: Linear/nonlinear solvers, vectors, matrices, time steppers
- **MPI**: Parallel communication
- **HDF5**: Binary output for large datasets
- **NetCDF**: Cubit/EXODUS mesh format
- **Proj**: Geographic coordinate transformations
- **spatialdata**: Spatial database queries

## Key Subsystems

### 1. Problem Management Subsystem

**Central Class**: `pylith::problems::Problem`

Orchestrates the entire simulation:

```cpp
Problem::run() {
    // 1. Setup phase
    preinitialize(mesh);
    verifyConfiguration();
    initialize();
    
    // 2. Time stepping loop
    while (t < t_end) {
        computeResidual();
        computeJacobian();
        solve();
        poststep();
        output();
    }
    
    // 3. Cleanup
    finalize();
}
```

**Problem Types**:
- `TimeDependent` - Most common: quasi-static and dynamic
- `GreensFns` - Specialized for Green's function computation

### 2. Finite Element Assembly Subsystem

**Key Classes**:
- `Integrator` - Base class for FE integration
- `IntegratorDomain` - Volume integration (bulk materials)
- `IntegratorBoundary` - Surface integration (boundary conditions)
- `IntegratorInterface` - Interface integration (faults)

**Kernel Architecture**:
```cpp
// Pointwise functions for FE computations
typedef void (*PetscPointFn)(/* ... */);

struct ResidualKernels {
    PetscPointFn f0;  // f0(x) term
    PetscPointFn f1;  // f1(x) term
};

struct JacobianKernels {
    PetscPointFn J0;  // ∂f0/∂u
    PetscPointFn J1;  // ∂f0/∂∇u
    PetscPointFn J2;  // ∂f1/∂u
    PetscPointFn J3;  // ∂f1/∂∇u
};
```

### 3. Material/Physics Subsystem

**Hierarchy**:
```
Physics (abstract base)
├── Material (abstract base for bulk materials)
│   ├── Elasticity
│   │   ├── IsotropicLinearElasticity
│   │   ├── IsotropicLinearMaxwell (viscoelastic)
│   │   └── IsotropicPowerLaw (power-law rheology)
│   ├── IncompressibleElasticity
│   └── Poroelasticity
└── (Other physics like boundary conditions, faults)
```

Each material provides:
- Auxiliary fields (material properties: density, elastic moduli)
- Derived fields (computed quantities: stress, strain)
- Residual kernels (weak form of governing equations)
- Jacobian kernels (linearization for Newton's method)

### 4. Mesh and Topology Subsystem

**Key Classes**:
- `Mesh` - Wrapper around PETSc DMPlex
- `Field` - Fields defined on mesh (solution, auxiliary, derived)
- `Distributor` - Parallel mesh distribution

**Mesh Representation**:
- Uses PETSc DMPlex (unstructured mesh with DAG)
- Supports 2D (tri, quad) and 3D (tet, hex) cells
- Handles parallel distribution automatically

### 5. I/O Subsystem

**Input**:
- `MeshIOAscii` - ASCII mesh format
- `MeshIOCubit` - EXODUS II format (from Cubit/Trelis)

**Output**:
- `OutputSoln*` - Solution output (domain, boundary, points)
- `DataWriterHDF5` - HDF5/Xdmf output (parallel-friendly)
- `DataWriterVTK` - VTK output (for visualization)

## Execution Flow

### Typical Simulation Lifecycle

```
1. Startup
   ├── Parse command line (pylithapp)
   ├── Load configuration files
   ├── Initialize Pyre components
   └── Setup logging

2. Meshing
   ├── Import or generate mesh
   ├── Adjust for faults (cohesive cells)
   ├── Distribute mesh (parallel)
   └── Setup topology

3. Problem Setup
   ├── Create solution field
   ├── Initialize materials (set properties)
   ├── Initialize boundary conditions
   ├── Initialize faults
   ├── Setup constraints
   └── Verify configuration

4. Problem Initialization
   ├── Setup PETSc TS (time stepper)
   ├── Setup PETSc SNES (nonlinear solver)
   ├── Setup PETSc KSP (linear solver)
   ├── Assemble initial Jacobian
   └── Apply initial conditions

5. Time Stepping Loop
   ├── Compute residual F(t,s)
   ├── Compute Jacobian J(t,s)
   ├── Solve linear system J·Δs = -F
   ├── Update solution s += Δs
   ├── Update state variables
   └── Write output

6. Finalization
   ├── Write final output
   ├── Cleanup PETSc objects
   └── Generate metadata
```

### Problem Formulation

PyLith solves problems of the form:

**F(t, s, ṡ) = G(t, s)**

Where:
- **s** = solution (displacement, pressure, etc.)
- **F** = LHS function (inertia, stiffness)
- **G** = RHS function (forces)
- **t** = time

This is mapped to PETSc TS as:
- **IFunction** = F - G (called "LHS - RHS")
- **RHSFunction** = G (explicit terms)

**Quasi-static**: F = Ku (stiffness only, no time derivatives)
**Dynamic**: F = Mü + Cu̇ + Ku (includes inertia and damping)

## Summary

PyLith's architecture achieves:

✅ **Performance**: C++ core with optimized numerics  
✅ **Flexibility**: Python interface for easy configuration  
✅ **Scalability**: MPI parallelism via PETSc  
✅ **Extensibility**: Component-based design  
✅ **Maintainability**: Clear separation of concerns  

The layered design allows users to:
- Run simulations without touching C++ code
- Extend functionality through Python
- Add new physics through C++ with minimal changes to core

Next: [02-Core-Components.md](02-Core-Components.md) - Detailed component descriptions
