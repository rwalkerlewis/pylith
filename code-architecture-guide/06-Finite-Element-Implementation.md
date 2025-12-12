# Finite Element Implementation

## Table of Contents
1. [Weak Form Formulation](#weak-form-formulation)
2. [PETSc Integration](#petsc-integration)
3. [Kernel Architecture](#kernel-architecture)
4. [Assembly Process](#assembly-process)
5. [Time Integration](#time-integration)
6. [Example: Elasticity](#example-elasticity)

## Weak Form Formulation

### Strong vs Weak Form

**Strong Form** (PDEs):
```
∇·σ + f = ρü    in Ω
u = g           on Γ_D (Dirichlet)
σ·n = h         on Γ_N (Neumann)
```

**Weak Form** (Variational):
```
Find u ∈ V such that for all v ∈ V:
∫_Ω v·ρü dV + ∫_Ω ∇v:σ dV = ∫_Ω v·f dV + ∫_Γ_N v·h dS
```

Where:
- **u**: Solution (displacement)
- **v**: Test function
- **σ**: Stress (constitutive law)
- **f**: Body force
- **h**: Traction

### Residual Form

PyLith casts problems as:
```
F(t, s, ṡ) = 0
```

Where F is the residual. For elasticity:

```
F = ∫_Ω (f0(u)·v + f1(u)·∇v) dV - ∫_Γ_N g·v dS

f0 = -f            (body force, negative because moved to LHS)
f1 = σ             (stress)
```

### Discretization

**Finite Element Approximation**:
```
u(x) ≈ Σ_i u_i φ_i(x)
v(x) = φ_j(x)
```

Where φ_i are basis functions.

**Discrete Residual**:
```
F_j = ∫_Ω (f0·φ_j + f1·∇φ_j) dV

For each DOF j:
F_j = Σ_cells ∫_cell (f0·φ_j + f1·∇φ_j) dV
    = Σ_cells Σ_quadrature_pts w_q |J_q| (f0·φ_j + f1·∇φ_j)
```

Where:
- **w_q**: Quadrature weight
- **|J_q|**: Jacobian determinant (cell volume scaling)

## PETSc Integration

### PETSc Data Structures

**DMPlex** (Distributed Mesh):
```
DM: Unstructured mesh
├─ Points: Vertices, edges, faces, cells
├─ Sections: DOF layout
├─ Coordinates: Vertex positions
└─ Labels: Grouping (materials, BCs)
```

**Vec** (Vectors):
```
- Global vector: Distributed across processes
- Local vector: Includes ghost points
- Section: Maps mesh points → vector indices
```

**Mat** (Matrices):
```
- Sparse matrix (CSR format)
- Parallel distribution
- Preallocation for efficiency
```

### DMPlex Field Layout

**PetscSection** describes how DOF are distributed:

```
Example: 3D displacement field, 4 vertices

Section:
  Vertex 0: offset=0,  dof=3  → [u_x, u_y, u_z]_0
  Vertex 1: offset=3,  dof=3  → [u_x, u_y, u_z]_1
  Vertex 2: offset=6,  dof=3  → [u_x, u_y, u_z]_2
  Vertex 3: offset=9,  dof=3  → [u_x, u_y, u_z]_3

Global vector: [u_x0, u_y0, u_z0, u_x1, u_y1, u_z1, ...]
```

**Multi-field**: Displacement + pressure

```
Section with 2 fields:
  Field 0 (displacement): 3 components
  Field 1 (pressure): 1 component

  Vertex 0: offset=0, dof=4 → [u_x, u_y, u_z, p]_0
  Vertex 1: offset=4, dof=4 → [u_x, u_y, u_z, p]_1
```

### PETSc Residual and Jacobian Functions

**Function Pointers**:
```cpp
// Residual function
PetscErrorCode (*IFunction)(
    TS ts,           // Time stepper
    PetscReal t,     // Current time
    Vec u,           // Solution
    Vec u_t,         // Time derivative
    Vec F,           // Output: Residual
    void* ctx        // User context
);

// Jacobian function
PetscErrorCode (*IJacobian)(
    TS ts,
    PetscReal t,
    Vec u,
    Vec u_t,
    PetscReal shift,  // For implicit time stepping
    Mat J,            // Output: Jacobian matrix
    Mat Jpre,         // Preconditioner matrix
    void* ctx
);
```

**PyLith Implementation**:
```cpp
class Problem {
    static PetscErrorCode computeIFunction(
        TS ts, PetscReal t, Vec u, Vec u_t, Vec F, void* ctx
    ) {
        Problem* problem = (Problem*)ctx;
        
        // Compute RHS (explicit terms)
        problem->computeRHSResidual(F, t, u);
        
        // Compute LHS (implicit terms)
        Vec Flhs;
        problem->computeLHSResidual(Flhs, t, u, u_t);
        
        // F = LHS - RHS
        VecAXPY(F, -1.0, Flhs);
        
        return 0;
    }
};
```

## Kernel Architecture

### Pointwise Functions

**Purpose**: Evaluate weak form terms at quadrature points.

**Signature**:
```cpp
typedef void (*PetscPointFn)(
    PetscInt dim,              // Spatial dimension (2 or 3)
    PetscInt Nf,               // Number of fields (solution)
    PetscInt NfAux,            // Number of auxiliary fields
    const PetscInt uOff[],     // Offsets for solution fields
    const PetscInt uOff_x[],   // Offsets for solution gradients
    const PetscScalar u[],     // Solution values
    const PetscScalar u_t[],   // Time derivatives
    const PetscScalar u_x[],   // Gradients
    const PetscInt aOff[],     // Offsets for auxiliary fields
    const PetscInt aOff_x[],   // Offsets for auxiliary gradients
    const PetscScalar a[],     // Auxiliary values
    const PetscScalar a_t[],   // Auxiliary time derivatives
    const PetscScalar a_x[],   // Auxiliary gradients
    PetscReal t,               // Time
    const PetscScalar x[],     // Coordinates
    PetscInt numConstants,     // Number of constants
    const PetscScalar constants[], // Constant values
    PetscScalar f[]            // Output: f0 or f1
);
```

### Residual Kernels

**f0 Kernel** (integrates with test function):
```
∫_Ω f0·v dV
```

**f1 Kernel** (integrates with test function gradient):
```
∫_Ω f1·∇v dV
```

**Example - Body Force** (f0):
```cpp
void bodyForce_f0(
    const PylithInt dim,
    const PylithScalar t,
    const PylithScalar x[],
    const PylithInt numS,
    const PylithInt numA,
    const PylithInt sOff[],
    const PylithInt aOff[],
    const PylithScalar s[],
    const PylithScalar s_t[],
    const PylithScalar s_x[],
    const PylithScalar a[],
    PylithScalar f0[]
) {
    // Extract body force from auxiliary field
    const PylithScalar* bodyForce = &a[aOff[1]];  // Assuming offset 1
    
    // Output: negative body force (moved to LHS)
    for (int i = 0; i < dim; ++i) {
        f0[i] = -bodyForce[i];
    }
}
```

**Example - Elasticity Stress** (f1):
```cpp
void elasticity_f1(
    const PylithInt dim,
    const PylithScalar t,
    const PylithScalar x[],
    const PylithInt numS,
    const PylithInt numA,
    const PylithInt sOff[],
    const PylithInt aOff[],
    const PylithScalar s[],
    const PylithScalar s_t[],
    const PylithScalar s_x[],
    const PylithScalar a[],
    PylithScalar f1[]
) {
    // Extract material properties
    const PylithScalar mu = a[aOff[0]];      // Shear modulus
    const PylithScalar lambda = a[aOff[1]];  // Lame parameter
    
    // Compute strain: ε = 0.5(∇u + ∇u^T)
    PylithScalar strain[6];  // [ε_xx, ε_yy, ε_zz, ε_xy, ε_yz, ε_xz]
    computeStrain(strain, s_x, dim);
    
    // Compute stress: σ = λ(tr ε)I + 2μ ε
    PylithScalar stress[6];
    computeStress(stress, strain, lambda, mu, dim);
    
    // Output: stress (flattened tensor)
    for (int i = 0; i < 6; ++i) {
        f1[i] = stress[i];
    }
}
```

### Jacobian Kernels

**Purpose**: Compute derivatives ∂F/∂u for Newton's method.

**Four Jacobian Terms**:
```
J0: ∂f0/∂u       (f0 w.r.t. solution)
J1: ∂f0/∂(∇u)    (f0 w.r.t. gradient)
J2: ∂f1/∂u       (f1 w.r.t. solution)
J3: ∂f1/∂(∇u)    (f1 w.r.t. gradient)
```

**Example - Linear Elasticity** (J3):
```cpp
void elasticity_Jf3(
    const PylithInt dim,
    const PylithScalar t,
    const PylithScalar x[],
    const PylithInt numS,
    const PylithInt numA,
    const PylithInt sOff[],
    const PylithInt aOff[],
    const PylithScalar s[],
    const PylithScalar s_t[],
    const PylithScalar s_x[],
    const PylithScalar a[],
    PylithScalar J[]
) {
    // Extract material properties
    const PylithScalar mu = a[aOff[0]];
    const PylithScalar lambda = a[aOff[1]];
    
    // J3 = elasticity tensor C
    // σ_ij = C_ijkl ε_kl
    // For isotropic: C_ijkl = λδ_ij δ_kl + μ(δ_ik δ_jl + δ_il δ_jk)
    
    // Fill elasticity tensor (9x9 for 3D)
    // J[i*dim*dim + j*dim + k] = ∂σ_ij/∂(∂u_k/∂x_l)
    
    // Diagonal terms
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            int idx = (i*dim + j) * dim*dim + (i*dim + j);
            J[idx] = lambda + 2*mu;  // Normal stress components
        }
    }
    
    // Off-diagonal (shear) terms
    // ... (depends on Voigt notation)
}
```

### Kernel Registration

**In Material Class**:
```cpp
std::vector<ResidualKernels>
IsotropicLinearElasticity::getKernelsResidual(...) const {
    std::vector<ResidualKernels> kernels;
    
    ResidualKernels rKernel;
    rKernel.subfieldName = "displacement";
    rKernel.f0 = NULL;  // No f0 term
    rKernel.f1 = pylith::fekernels::IsotropicLinearElasticity::f1v;
    kernels.push_back(rKernel);
    
    return kernels;
}

std::vector<JacobianKernels>
IsotropicLinearElasticity::getKernelsJacobian(...) const {
    std::vector<JacobianKernels> kernels;
    
    JacobianKernels jKernel;
    jKernel.subfieldName = "displacement";
    jKernel.fieldTrial = "displacement";
    jKernel.J0 = NULL;
    jKernel.J1 = NULL;
    jKernel.J2 = NULL;
    jKernel.J3 = pylith::fekernels::IsotropicLinearElasticity::Jf3uu;
    kernels.push_back(jKernel);
    
    return kernels;
}
```

## Assembly Process

### Residual Assembly

```
┌─────────────────────────────────────────────────────────────┐
│ DMPlexComputeResidual_Internal()  [PETSc]                   │
│                                                              │
│  FOR each integrator (material, BC, fault):                 │
│    label = integrator.getLabel()                            │
│    value = integrator.getLabelValue()                       │
│    kernels = integrator.getKernels()                        │
│                                                              │
│    FOR each cell c in label[value]:                         │
│      ┌──────────────────────────────────────────────┐      │
│      │ 1. Get cell closure (vertices, edges)        │      │
│      │    DMPlexGetClosureIndices()                 │      │
│      └──────────────────────────────────────────────┘      │
│      ┌──────────────────────────────────────────────┐      │
│      │ 2. Get cell geometry                          │      │
│      │    - Compute Jacobian J = ∂x/∂ξ             │      │
│      │    - Compute |det(J)|                         │      │
│      │    - Compute J^{-T} for gradient transform   │      │
│      └──────────────────────────────────────────────┘      │
│      ┌──────────────────────────────────────────────┐      │
│      │ 3. Extract solution at cell vertices          │      │
│      │    u_cell[] = solution[vertex_indices]        │      │
│      └──────────────────────────────────────────────┘      │
│      ┌──────────────────────────────────────────────┐      │
│      │ 4. FOR each quadrature point q:              │      │
│      │                                               │      │
│      │  a. Evaluate basis functions: φ_i(ξ_q)      │      │
│      │  b. Evaluate gradients: ∇φ_i(ξ_q)           │      │
│      │  c. Interpolate u: u_q = Σ u_i φ_i(ξ_q)     │      │
│      │  d. Transform gradient: ∇u_q = J^{-T} ∂u/∂ξ │      │
│      │                                               │      │
│      │  e. Call kernel:                             │      │
│      │     kernel(x_q, u_q, ∇u_q, aux, f)          │      │
│      │                                               │      │
│      │  f. Integrate:                                │      │
│      │     F_i += w_q |J| (f0·φ_i + f1·∇φ_i)       │      │
│      │                                               │      │
│      └──────────────────────────────────────────────┘      │
│      ┌──────────────────────────────────────────────┐      │
│      │ 5. Add cell contribution to global residual  │      │
│      │    VecSetValues(F, indices, values, ADD)     │      │
│      └──────────────────────────────────────────────┘      │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Jacobian Assembly

Similar to residual, but computes:
```
∫_Ω (J0 φ_i φ_j + J1 φ_i ∇φ_j + J2 ∇φ_i φ_j + J3 ∇φ_i ∇φ_j) dV
```

**Result**: Sparse matrix
```
Matrix entry: A[i,j] = ∂F_i/∂u_j
```

### Matrix Preallocation

**Purpose**: Allocate memory for sparse matrix before assembly.

```cpp
// Estimate nonzeros per row
DMPlexPreallocateOperator(
    dm,               // Mesh
    numFields,        // Number of fields
    section,          // DOF layout
    globalSection,    // Global DOF layout
    dnz,              // Diagonal nonzeros per row
    onz,              // Off-diagonal nonzeros per row
    matrix,           // Matrix to preallocate
    fillMatrix        // Actually fill structure
);
```

**Why Important**: Prevents expensive dynamic reallocation during assembly.

## Time Integration

### PETSc TS (Time Stepper)

**Time Stepping Methods**:
- **Explicit**: Forward Euler, Runge-Kutta
- **Implicit**: Backward Euler, BDF, Theta methods
- **IMEX**: Implicit-Explicit (different treatment for different terms)

**PyLith Usage**:
```cpp
// Setup
TSCreate(comm, &ts);
TSSetType(ts, TSBEULER);  // Backward Euler (implicit)
TSSetIFunction(ts, residual, computeIFunction, problem);
TSSetIJacobian(ts, jacobian, jacobianPre, computeIJacobian, problem);

// Time step parameters
TSSetTime(ts, t_start);
TSSetMaxTime(ts, t_end);
TSSetTimeStep(ts, dt);
TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);

// Solve
TSSolve(ts, solution);
```

### Implicit Time Stepping

**Backward Euler**:
```
(u^{n+1} - u^n) / Δt = f(t^{n+1}, u^{n+1})

Rearrange:
u^{n+1} - Δt f(t^{n+1}, u^{n+1}) = u^n

Define residual:
F(u^{n+1}) = u^{n+1} - u^n - Δt f(t^{n+1}, u^{n+1})

Solve: F(u^{n+1}) = 0 (nonlinear system)
```

**Newton's Method**:
```
WHILE not converged:
  1. Compute residual: F(u_k)
  2. Compute Jacobian: J(u_k) = ∂F/∂u
  3. Solve linear system: J Δu = -F
  4. Update: u_{k+1} = u_k + Δu
```

**PETSc Implementation**:
```
TSSolve() → SNESSolve()
  → WHILE not converged:
      → KSPSolve(J, Δu, -F)
```

### Quasi-Static Approximation

**Idea**: Drop time derivatives (ü ≈ 0).

**Equation**:
```
∇·σ + f = 0  (equilibrium)
```

**Time Stepping**: "Pseudo-time" for incremental loading.

```python
[pylithapp.problem]
formulation = quasistatic

# Time steps represent loading increments
initial_dt = 0.1*year
end_time = 100.0*year
```

## Example: Elasticity

### Complete FE Implementation

**Problem**: Static linear elasticity

**Governing Equation**:
```
-∇·σ = f    in Ω
σ = C:ε     (constitutive)
ε = ∇_s u   (kinematic)
u = g       on Γ_D
σ·n = h     on Γ_N
```

**Weak Form**:
```
∫_Ω ∇v:C:∇u dV = ∫_Ω v·f dV + ∫_Γ_N v·h dS
```

### Implementation Steps

#### 1. Define Material Properties (Auxiliary Field)

```cpp
void IsotropicLinearElasticity::addAuxiliarySubfields() {
    _auxiliaryFactory->addDensity();         // ρ
    _auxiliaryFactory->addShearModulus();    // μ
    _auxiliaryFactory->addBulkModulus();     // K
}
```

#### 2. Implement Residual Kernel (f1)

```cpp
void pylith::fekernels::IsotropicLinearElasticity::f1v(
    /* ... parameters ... */
    PylithScalar f1[]
) {
    // Get properties
    const PylithScalar mu = a[aOff[0]];
    const PylithScalar lambda = a[aOff[1]];
    
    // Compute strain from gradient
    // ε_ij = 0.5(∂u_i/∂x_j + ∂u_j/∂x_i)
    
    // Compute stress
    // σ_ij = λ δ_ij tr(ε) + 2μ ε_ij
    
    // Return stress (will be integrated with ∇v)
    f1[...] = stress[...];
}
```

#### 3. Implement Jacobian Kernel (Jf3)

```cpp
void pylith::fekernels::IsotropicLinearElasticity::Jf3uu(
    /* ... parameters ... */
    PylithScalar J[]
) {
    // Get properties
    const PylithScalar mu = a[aOff[0]];
    const PylithScalar lambda = a[aOff[1]];
    
    // Elasticity tensor C
    // C_ijkl = λ δ_ij δ_kl + μ(δ_ik δ_jl + δ_il δ_jk)
    
    // Store in J[] (flattened 4th-order tensor)
    J[...] = C[...];
}
```

#### 4. Assembly (Automatic via PETSc)

```cpp
Problem::computeLHSResidual(Vec residual, PetscReal t, Vec solution) {
    // For each material:
    for (auto material : _materials) {
        auto kernels = material->getKernelsResidual();
        DMPlexComputeResidual(dm, label, value, kernels, residual);
    }
}
```

#### 5. Solve

```cpp
// PETSc solves: J Δu = -F
KSPSolve(ksp, residual, update);
VecAXPY(solution, 1.0, update);  // u += Δu
```

## Summary

PyLith's FE implementation:

✅ **Weak Form**: Variational formulation  
✅ **Kernels**: Pointwise evaluation at quadrature points  
✅ **PETSc Integration**: DMPlex for mesh, TS for time stepping  
✅ **Automatic Assembly**: PETSc handles integration  
✅ **Modular**: Kernels separate from assembly  

Key concepts:
- **f0, f1**: Residual terms
- **J0, J1, J2, J3**: Jacobian terms
- **DMPlex**: Mesh representation
- **PetscSection**: DOF layout
- **TS**: Time integration

Next: [07-Extension-Guide.md](07-Extension-Guide.md) - How to add custom components
