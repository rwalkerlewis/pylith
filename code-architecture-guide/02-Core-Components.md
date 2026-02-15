# Core Components

## Table of Contents
1. [Problem Components](#problem-components)
2. [Material Components](#material-components)
3. [Boundary Condition Components](#boundary-condition-components)
4. [Fault Components](#fault-components)
5. [Assembly Components](#assembly-components)
6. [Topology and Field Components](#topology-and-field-components)

## Problem Components

### Problem (Base Class)

**Location**: `libsrc/pylith/problems/Problem.hh/cc` and `pylith/problems/Problem.py`

**Purpose**: Abstract base class for all boundary value problems. Orchestrates the entire simulation.

**Key Responsibilities**:
- Manage solution field
- Coordinate integrators (materials, BCs, faults)
- Interface with PETSc solvers
- Control time stepping
- Manage observers for output

**Important Methods**:

```cpp
class Problem {
    // Lifecycle
    void preinitialize(const Mesh& mesh);
    void verifyConfiguration();
    void initialize();
    void run();
    void finalize();
    
    // Solver setup
    void setFormulation(FormulationEnum value);
    void setSolverType(SolverTypeEnum value);
    
    // Component management
    void setMaterials(Material* materials[], int numMaterials);
    void setBoundaryConditions(BoundaryCondition* bcs[], int numBCs);
    void setInterfaces(FaultCohesive* faults[], int numFaults);
    
    // PETSc integration
    virtual PetscErrorCode computeRHSResidual(PetscVec residual, PetscVec solution);
    virtual PetscErrorCode computeLHSResidual(PetscVec residual, PetscVec solution);
    virtual PetscErrorCode computeLHSJacobian(PetscMat jacobian, PetscVec solution);
};
```

**Configuration Parameters**:
- `formulation`: "quasistatic", "dynamic", or "dynamic_imex"
- `solver`: "linear" or "nonlinear"
- `scales`: Nondimensionalization scales

### TimeDependent

**Location**: `libsrc/pylith/problems/TimeDependent.hh/cc` and `pylith/problems/TimeDependent.py`

**Purpose**: Concrete implementation for time-dependent problems (both quasi-static and dynamic).

**Features**:
- Time stepping with PETSc TS
- Adaptive time stepping support
- Start time, end time, max time step configuration
- Progress monitoring

**Additional Methods**:

```cpp
class TimeDependent : public Problem {
    void setStartTime(double value);
    void setEndTime(double value);
    void setMaxTimeStep(double value);
    void setInitialTimeStep(double value);
    
    // Time stepping
    PetscErrorCode poststep(PetscReal t, PetscVec solution);
};
```

**Python Configuration**:
```python
[pylithapp.problem]
initial_dt = 1.0*year
start_time = 0.0*year  
end_time = 100.0*year
```

### GreensFns

**Location**: `libsrc/pylith/problems/GreensFns.hh/cc` and `pylith/problems/GreensFns.py`

**Purpose**: Specialized problem for computing Green's functions (impulse responses).

**Use Case**: Compute static deformation due to unit slip on fault patches.

## Material Components

### Material (Base Class)

**Location**: `libsrc/pylith/materials/Material.hh/cc` and `pylith/materials/Material.py`

**Purpose**: Abstract base for all bulk materials. Inherits from `Physics`.

**Key Concepts**:
- **Auxiliary Fields**: Material properties (density, elastic moduli, viscosity)
- **Derived Fields**: Computed quantities (stress, strain)
- **Rheology**: Constitutive model (elastic, viscoelastic, plastic)

**Material Hierarchy**:

```
Material (abstract)
├── Elasticity (governing equation for elasticity)
│   ├── IsotropicLinearElasticity (linear elastic)
│   ├── IsotropicLinearMaxwell (Maxwell viscoelastic)
│   ├── IsotropicLinearGenMaxwell (Generalized Maxwell)
│   └── IsotropicPowerLaw (Power-law rheology)
├── IncompressibleElasticity (incompressible materials)
│   └── IsotropicLinearIncompElasticity
└── Poroelasticity (fluid-saturated porous media)
    └── IsotropicLinearPoroelasticity
```

### Elasticity

**Location**: `libsrc/pylith/materials/Elasticity.hh/cc` and `pylith/materials/Elasticity.py`

**Governing Equation**: ∇·σ + f = ρü (quasi-static: acceleration term dropped)

**Solution Fields**:
- Displacement (u)
- Velocity (v) - for dynamic problems
- Lagrange multipliers - for fault constraints

**Auxiliary Fields** (depends on rheology):
- Density (ρ)
- Shear modulus (μ)
- Bulk modulus (K)
- Reference stress/strain
- Viscosity parameters (for viscoelastic)

**Derived Fields**:
- Cauchy stress (σ)
- Strain (ε)

**Example - IsotropicLinearElasticity**:

```cpp
class IsotropicLinearElasticity : public RheologyElasticity {
    // Constitutive model: σ = λ(tr ε)I + 2μ dev(ε)
    
    // Auxiliary field requirements
    static const char* auxFieldsSetup[];
    // ["density", "shear_modulus", "bulk_modulus"]
    
    // Kernel functions
    static void computeStress(/* ... */);  // Derived field kernel
    static void f1(/* ... */);             // Residual f1 kernel
    static void Jf3uu(/* ... */);          // Jacobian kernel
};
```

### Viscoelastic Materials

**Maxwell Viscoelasticity**:
- **Model**: Spring and dashpot in series
- **Constitutive Law**: σ̇ + (μ/η)σ = 2με̇
- **State Variables**: Viscous strain

**Generalized Maxwell**:
- Multiple Maxwell elements in parallel
- Allows for broader relaxation spectrum

**Implementation**:
```cpp
class IsotropicLinearMaxwell : public RheologyElasticity {
    // Additional auxiliary fields
    // ["maxwell_time", "shear_modulus_ratio", "viscous_strain"]
    
    // State variable update kernel
    static void updateStateVars(/* ... */);
};
```

### Material Configuration (Python)

```python
[pylithapp.problem.materials.slab]
# Material uses Python wrapper + C++ implementation
use_reference_state = False

# Auxiliary field setup
auxiliary_subfields.density.basis_order = 0
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.label = Elastic properties
db_auxiliary_field.iohandler.filename = mat_elastic.spatialdb

# Observers for output
observers.observer.data_fields = [displacement, cauchy_stress]
```

## Boundary Condition Components

### BoundaryCondition (Base Class)

**Location**: `libsrc/pylith/bc/` and `pylith/bc/`

**Types**:
- `DirichletTimeDependent` - Prescribed displacement/velocity
- `NeumannTimeDependent` - Prescribed traction
- `AbsorbingDampers` - Absorbing boundary

### DirichletTimeDependent

**Purpose**: Prescribe displacement or velocity on boundaries.

**Key Features**:
- Time-dependent prescribed values
- Can constrain individual components (x, y, z)
- Uses spatial databases for spatial variation
- Implemented as constraint (not in residual/Jacobian)

**Configuration**:
```python
[pylithapp.problem.bc.bc_xneg]
constrained_dof = [0, 1]  # Constrain x and y components
label = boundary_xneg
label_value = 10

db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.values = [initial_amplitude_x, initial_amplitude_y]
db_auxiliary_field.data = [0.0*m, -2.0*m]
```

**Implementation**:
```cpp
class DirichletTimeDependent : public BoundaryCondition {
    // Uses constraint projector in PETSc
    // No contribution to residual/Jacobian
    // Applied before/after solve
    
    void setConstraintDOF(int* dof, int numDOF);
    void setDB(spatialdata::spatialdb::SpatialDB* db);
};
```

### NeumannTimeDependent

**Purpose**: Apply traction (force per unit area) on boundaries.

**Key Features**:
- Time-dependent tractions
- Integrated into residual (weak form)
- Can be initial (reference) + rate + time history

**Weak Form Contribution**:
```
∫_Γ v · T dS
```

Where T is the prescribed traction, Γ is the boundary, v is the test function.

**Kernels**:
```cpp
// Residual kernel
static void f0(
    const PylithInt dim,
    const PylithScalar t,
    const PylithScalar* x,
    const PylithScalar* n,  // Surface normal
    const PylithScalar* auxiliaryField,
    PylithScalar* f0) {
    // f0 = traction from auxiliary field
}
```

## Fault Components

### FaultCohesive (Base Class)

**Location**: `libsrc/pylith/faults/` and `pylith/faults/`

**Purpose**: Represents fault surfaces with slip using cohesive cells.

**Cohesive Cell Approach**:
- Insert zero-thickness cells along fault surface
- Negative and positive sides of fault
- Lagrange multipliers enforce slip constraint

```
    + side          - side
      |              |
    --|--------------|--
      | cohesive    |
      | cell        |
```

**Types**:
- `FaultCohesiveKin` - Prescribed (kinematic) slip
- `FaultCohesiveImpulses` - For Green's functions
- (Future: spontaneous rupture with friction)

### FaultCohesiveKin

**Purpose**: Prescribe slip history on fault.

**Features**:
- Multiple earthquake ruptures
- Slip time function (step, ramp, Brune, Liu-Cosine, etc.)
- Spatially variable slip amplitude

**Slip Sources (KinSrc)**:
- `KinSrcStep` - Instantaneous slip (quasi-static)
- `KinSrcRamp` - Linear ramp
- `KinSrcBrune` - Brune far-field time function
- `KinSrcLiuCos` - Liu cosine time function
- `KinSrcTimeHistory` - User-defined time history

**Configuration**:
```python
[pylithapp.problem.interfaces.fault]
label = fault
label_value = 20
edge = fault_edge
edge_value = 21

observers.observer.data_fields = [slip, traction_change]

# Kinematic rupture
[pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
db_auxiliary_field = spatialdata.spatialdb.UniformDB
db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]
```

**Implementation**:
```cpp
class FaultCohesiveKin : public FaultCohesive {
    // Solution subfields
    // - Lagrange multiplier (traction on fault)
    
    // Auxiliary subfields
    // - Slip
    // - Slip rate (for dynamic)
    
    // Constraint kernel
    static void constraintLagrange(/* ... */);  // Enforce slip
    static void constraintVelocity(/* ... */);  // For dynamic
};
```

## Assembly Components

### Integrator (Base Class)

**Location**: `libsrc/pylith/feassemble/Integrator.hh/cc`

**Purpose**: Base class for finite element assembly.

**Key Concepts**:
- Operates on a labeled region of mesh
- Manages auxiliary fields
- Provides kernels for residual and Jacobian
- Updates state variables

**Types**:
- `IntegratorDomain` - Volume integration (materials)
- `IntegratorBoundary` - Surface integration (Neumann BCs)
- `IntegratorInterface` - Interface integration (faults)

**Integration Points**:
```cpp
enum EquationPart {
    LHS = 0,        // Left-hand side (stiffness)
    RHS = 1,        // Right-hand side (forces)
    LHS_LUMPED_INV = 2,  // Lumped mass matrix inverse (explicit)
    LHS_WEIGHTED = 3     // Mass matrix (implicit)
};
```

### IntegratorDomain

**Purpose**: Integrate bulk material contributions.

**Residual Computation**:
```cpp
// Weak form: ∫_Ω (f0·v + f1·∇v) dV
struct ResidualKernels {
    PetscPointFn f0;  // f0(x, u, ∇u, auxiliary)
    PetscPointFn f1;  // f1(x, u, ∇u, auxiliary)
};

// Example for elasticity:
// f0 = 0 (no body force term in this kernel)
// f1 = σ (stress, integrated as ∫ σ:∇v dV)
```

**Jacobian Computation**:
```cpp
// Jacobian: ∂F/∂u
struct JacobianKernels {
    PetscPointFn J0;  // ∂f0/∂u
    PetscPointFn J1;  // ∂f0/∂∇u
    PetscPointFn J2;  // ∂f1/∂u
    PetscPointFn J3;  // ∂f1/∂∇u
};

// Example for linear elasticity:
// J3 = C (elasticity tensor)
```

### Kernels

**Location**: `libsrc/pylith/fekernels/`

**Purpose**: Pointwise functions evaluated at quadrature points.

**Structure**:
```cpp
// Example: Elasticity residual f1 kernel
void pylith::fekernels::IsotropicLinearElasticity::f1v(
    const PylithInt dim,
    const PylithScalar t,
    const PylithScalar x[],
    const PylithInt numS,         // Number of solution fields
    const PylithInt numA,         // Number of auxiliary fields  
    const PylithInt sOff[],       // Solution field offsets
    const PylithInt aOff[],       // Auxiliary field offsets
    const PylithScalar s[],       // Solution values
    const PylithScalar s_t[],     // Solution time derivatives
    const PylithScalar s_x[],     // Solution gradients
    const PylithScalar a[],       // Auxiliary values
    PylithScalar f1[]             // Output: f1 value
) {
    // Extract material properties from auxiliary fields
    const PylithScalar shearModulus = a[aOff[0]];
    const PylithScalar bulkModulus = a[aOff[1]];
    
    // Extract strain from solution gradient
    // strain = 0.5(∇u + ∇u^T)
    
    // Compute stress: σ = λ(tr ε)I + 2μ dev(ε)
    // Store in f1
}
```

## Topology and Field Components

### Mesh

**Location**: `libsrc/pylith/topology/Mesh.hh/cc`

**Purpose**: Wrapper around PETSc DMPlex for unstructured mesh.

**Features**:
- Parallel mesh distribution
- Multiple cell types (tri, quad, tet, hex)
- Label management for boundaries, materials, faults
- Coordinate system

**Key Methods**:
```cpp
class Mesh {
    void createFromFile(const char* filename);
    void distribute();  // Parallel distribution
    PetscDM getDM();    // Get underlying PETSc DM
    
    // Dimension queries
    int getDimension();
    void getCells(PetscIS* cells);
    void getVertices(PetscIS* vertices);
    
    // Label access
    void createLabel(const char* name);
    void getStratumIS(const char* label, int value, PetscIS* is);
};
```

### Field

**Location**: `libsrc/pylith/topology/Field.hh/cc`

**Purpose**: Represents fields (solution, auxiliary, derived) on mesh.

**Structure**:
- **Field**: Contains multiple subfields
- **Subfield**: Individual component (e.g., "displacement", "pressure")
- **Component**: Scalar component of a subfield

**Discretization**:
- Basis order (0 = piecewise constant, 1 = piecewise linear, etc.)
- Quadrature order
- Basis type (Lagrange, simplex, etc.)

**Important Methods**:
```cpp
class Field {
    // Subfield management
    void addSubfield(const char* name, int numComponents);
    void setSubfieldBasisOrder(const char* name, int order);
    
    // Allocation
    void allocate();
    void createDiscretization();
    
    // Access
    PetscVec getGlobalVector();
    PetscVec getLocalVector();
    PetscSection getSection();  // Describes layout
    
    // Scatter between global and local
    void scatterLocalToGlobal();
    void scatterGlobalToLocal();
};
```

**Example - Solution Field for Elasticity**:
```cpp
// Create solution field
Field solution(mesh);
solution.addSubfield("displacement", 3);  // 3D vector
solution.addSubfield("velocity", 3);      // For dynamic
solution.setSubfieldBasisOrder("displacement", 1);  // Linear
solution.createDiscretization();
solution.allocate();
```

### Visitors

**Location**: `libsrc/pylith/topology/VisitorMesh.hh`

**Purpose**: Efficient traversal and access to field data.

**Use Case**: Access field values at mesh vertices/cells.

```cpp
VisitorMesh visitor(field);
PetscScalar* array = visitor.getLocalArray();

// Access data at point 'p'
PetscInt off = visitor.getSectionOffset(p);
PetscInt dof = visitor.getSectionDof(p);
// Data is array[off] to array[off+dof-1]
```

## Summary

The core components provide:

1. **Problem**: Orchestration and solver integration
2. **Materials**: Constitutive models and governing equations
3. **Boundary Conditions**: Constraints and surface loads
4. **Faults**: Slip representation
5. **Assembly**: Finite element integration
6. **Topology/Fields**: Mesh and data management

Each component is:
- **Modular**: Independent, swappable
- **Configurable**: Via Python/Pyre
- **Extensible**: Can be subclassed

Next: [03-Data-Flow.md](03-Data-Flow.md) - How data moves through the system
