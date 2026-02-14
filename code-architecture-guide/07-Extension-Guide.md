# Extension Guide

## Table of Contents
1. [Adding a Custom Material](#adding-a-custom-material)
2. [Adding a Custom Boundary Condition](#adding-a-custom-boundary-condition)
3. [Adding a Custom Fault Slip Function](#adding-a-custom-fault-slip-function)
4. [Adding Derived Fields](#adding-derived-fields)
5. [Testing Your Extension](#testing-your-extension)

## Adding a Custom Material

### Overview

To add a new material rheology (e.g., Drucker-Prager plasticity, temperature-dependent viscosity):

1. Define C++ rheology class
2. Implement kernel functions
3. Create SWIG interface
4. Add Python wrapper
5. Test and document

### Step-by-Step: Custom Viscoelastic Material

#### Step 1: C++ Header File

Create `libsrc/pylith/materials/MyViscoelastic.hh`:

```cpp
#pragma once

#include "pylith/materials/RheologyElasticity.hh"

namespace pylith {
    namespace materials {
        class MyViscoelastic;
    } // materials
} // pylith

/**
 * My custom viscoelastic rheology.
 *
 * Constitutive equation:
 * σ̇ + (μ/η) σ = 2μ ε̇
 *
 * where:
 *   μ = shear modulus
 *   η = viscosity
 */
class pylith::materials::MyViscoelastic : public RheologyElasticity {
    friend class TestMyViscoelastic; // unit testing

public:
    
    /// Constructor
    MyViscoelastic();
    
    /// Destructor
    virtual ~MyViscoelastic();
    
    /// Deallocate
    void deallocate();
    
    /** Use reference stress and strain in formulation?
     *
     * @param[in] value Flag indicating to include reference stress/strain.
     */
    void useReferenceState(const bool value);
    
    /** Use reference stress and strain in formulation?
     *
     * @returns True if using reference stress/strain, false otherwise.
     */
    bool useReferenceState() const;
    
    /** Add rheology subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields();
    
    /** Get stress kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     * @param[in] gravityField Gravity field (NULL if no gravity).
     * @param[in] isJacobianTerm True if adding terms for Jacobian.
     *
     * @returns Residual kernel for stress.
     */
    std::vector<ResidualKernels> getKernelsResidual(
        const spatialdata::geocoords::CoordSys* coordsys,
        const bool gravityField,
        const bool isJacobianTerm
    ) const;
    
    /** Get elastic constants kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     * @param[in] gravityField Gravity field (NULL if no gravity).
     *
     * @returns Jacobian kernel for elastic constants.
     */
    std::vector<JacobianKernels> getKernelsJacobian(
        const spatialdata::geocoords::CoordSys* coordsys,
        const bool gravityField
    ) const;
    
    /** Get kernels for updating state variables.
     *
     * @returns Project kernels for updating state variables.
     */
    std::vector<ProjectKernels> getKernelsUpdateStateVars() const;
    
private:
    
    bool _useReferenceState; ///< Flag to use reference stress/strain.
    
}; // MyViscoelastic
```

#### Step 2: C++ Implementation

Create `libsrc/pylith/materials/MyViscoelastic.cc`:

```cpp
#include "MyViscoelastic.hh"

#include "pylith/fekernels/MyViscoelastic.hh" // kernels
#include "pylith/materials/AuxiliaryFactoryElastic.hh"

#include <cassert>

// Constructor
pylith::materials::MyViscoelastic::MyViscoelastic() :
    _useReferenceState(false) {
}

// Destructor
pylith::materials::MyViscoelastic::~MyViscoelastic() {
    deallocate();
}

// Deallocate
void pylith::materials::MyViscoelastic::deallocate() {
}

// Use reference state?
void pylith::materials::MyViscoelastic::useReferenceState(const bool value) {
    _useReferenceState = value;
}

bool pylith::materials::MyViscoelastic::useReferenceState() const {
    return _useReferenceState;
}

// Add auxiliary subfields
void pylith::materials::MyViscoelastic::addAuxiliarySubfields() {
    // Add standard elastic fields
    _auxiliaryFactory->addDensity();
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addBulkModulus();
    
    // Add viscoelastic-specific fields
    _auxiliaryFactory->addViscosity();
    _auxiliaryFactory->addMaxwellTime();
    
    // Add state variables
    _auxiliaryFactory->addViscousStrain();  // 6 components
    
    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    }
}

// Get residual kernels
std::vector<pylith::materials::RheologyElasticity::ResidualKernels>
pylith::materials::MyViscoelastic::getKernelsResidual(
    const spatialdata::geocoords::CoordSys* coordsys,
    const bool gravityField,
    const bool isJacobianTerm
) const {
    std::vector<ResidualKernels> kernels;
    
    ResidualKernels rKernel;
    rKernel.subfieldName = "displacement";
    
    if (_useReferenceState) {
        rKernel.f0 = NULL;
        rKernel.f1 = pylith::fekernels::MyViscoelastic::f1v_refstate;
    } else {
        rKernel.f0 = NULL;
        rKernel.f1 = pylith::fekernels::MyViscoelastic::f1v;
    }
    
    kernels.push_back(rKernel);
    return kernels;
}

// Get Jacobian kernels
std::vector<pylith::materials::RheologyElasticity::JacobianKernels>
pylith::materials::MyViscoelastic::getKernelsJacobian(
    const spatialdata::geocoords::CoordSys* coordsys,
    const bool gravityField
) const {
    std::vector<JacobianKernels> kernels;
    
    JacobianKernels jKernel;
    jKernel.subfieldName = "displacement";
    jKernel.fieldTrial = "displacement";
    jKernel.J0 = NULL;
    jKernel.J1 = NULL;
    jKernel.J2 = NULL;
    
    if (_useReferenceState) {
        jKernel.J3 = pylith::fekernels::MyViscoelastic::Jf3uu_refstate;
    } else {
        jKernel.J3 = pylith::fekernels::MyViscoelastic::Jf3uu;
    }
    
    kernels.push_back(jKernel);
    return kernels;
}

// Get state variable update kernels
std::vector<pylith::feassemble::Integrator::ProjectKernels>
pylith::materials::MyViscoelastic::getKernelsUpdateStateVars() const {
    std::vector<ProjectKernels> kernels;
    
    ProjectKernels pKernel;
    pKernel.subfield = "viscous_strain";
    pKernel.f = pylith::fekernels::MyViscoelastic::updateViscousStrain;
    
    kernels.push_back(pKernel);
    return kernels;
}
```

#### Step 3: Kernel Implementation

Create `libsrc/pylith/fekernels/MyViscoelastic.hh`:

```cpp
#pragma once

#include "pylith/fekernels/fekernelsfwd.hh"

class pylith::fekernels::MyViscoelastic {
public:
    
    // Auxiliary field indices
    struct AuxConstants {
        static const int density = 0;
        static const int shearModulus = 1;
        static const int bulkModulus = 2;
        static const int viscosity = 3;
        static const int maxwellTime = 4;
        static const int viscousStrain = 5;
        static const int referenceStress = 11;  // if used
        static const int referenceStrain = 17;  // if used
    };
    
    /** f1 function for elasticity equation: f1 = stress.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), 
     *                    viscosity(1), maxwell_time(1), viscous_strain(6)]
     */
    static void f1v(
        const PylithInt dim,
        const PylithInt numS,
        const PylithInt numA,
        const PylithInt sOff[],
        const PylithInt sOff_x[],
        const PylithScalar s[],
        const PylithScalar s_t[],
        const PylithScalar s_x[],
        const PylithInt aOff[],
        const PylithInt aOff_x[],
        const PylithScalar a[],
        const PylithScalar a_t[],
        const PylithScalar a_x[],
        const PylithReal t,
        const PylithScalar x[],
        const PylithInt numConstants,
        const PylithScalar constants[],
        PylithScalar f1[]
    );
    
    /** Jf3 function for elasticity equation (Jacobian).
     */
    static void Jf3uu(/* ... same signature ... */);
    
    /** Update viscous strain state variable.
     *
     * Exponential integration of viscous strain:
     * ε_v^{n+1} = exp(-Δt/τ) ε_v^n + [1-exp(-Δt/τ)] dev(ε_total)
     */
    static void updateViscousStrain(
        const PylithInt dim,
        const PylithInt numS,
        const PylithInt numA,
        const PylithInt sOff[],
        const PylithInt sOff_x[],
        const PylithScalar s[],
        const PylithScalar s_t[],
        const PylithScalar s_x[],
        const PylithInt aOff[],
        const PylithInt aOff_x[],
        const PylithScalar a[],
        const PylithScalar a_t[],
        const PylithScalar a_x[],
        const PylithReal t,
        const PylithScalar x[],
        const PylithInt numConstants,
        const PylithScalar constants[],
        PylithScalar stateVars[]
    );
};
```

Create `libsrc/pylith/fekernels/MyViscoelastic.cc`:

```cpp
#include "MyViscoelastic.hh"
#include "Tensor.hh"  // Tensor utilities

#include <cassert>

// f1 residual kernel
void pylith::fekernels::MyViscoelastic::f1v(
    const PylithInt dim,
    const PylithInt numS,
    const PylithInt numA,
    const PylithInt sOff[],
    const PylithInt sOff_x[],
    const PylithScalar s[],
    const PylithScalar s_t[],
    const PylithScalar s_x[],
    const PylithInt aOff[],
    const PylithInt aOff_x[],
    const PylithScalar a[],
    const PylithScalar a_t[],
    const PylithScalar a_x[],
    const PylithReal t,
    const PylithScalar x[],
    const PylithInt numConstants,
    const PylithScalar constants[],
    PylithScalar f1[]
) {
    // Extract material properties
    const PylithScalar shearModulus = a[aOff[AuxConstants::shearModulus]];
    const PylithScalar bulkModulus = a[aOff[AuxConstants::bulkModulus]];
    const PylithScalar* viscousStrain = &a[aOff[AuxConstants::viscousStrain]];
    
    const PylithScalar lambda = bulkModulus - 2.0/3.0 * shearModulus;
    
    // Compute total strain from displacement gradient
    PylithScalar totalStrain[6];
    Tensor::symmetrize(totalStrain, s_x, dim);  // ε = 0.5(∇u + ∇u^T)
    
    // Compute deviatoric strain
    const PylithScalar meanStrain = (totalStrain[0] + totalStrain[1] + totalStrain[2]) / 3.0;
    PylithScalar devStrain[6];
    for (int i = 0; i < 3; ++i) {
        devStrain[i] = totalStrain[i] - meanStrain;
    }
    for (int i = 3; i < 6; ++i) {
        devStrain[i] = totalStrain[i];
    }
    
    // Elastic deviatoric strain = total - viscous
    PylithScalar elasticDevStrain[6];
    for (int i = 0; i < 6; ++i) {
        elasticDevStrain[i] = devStrain[i] - viscousStrain[i];
    }
    
    // Compute stress
    // σ = λ tr(ε) I + 2μ (ε_dev - ε_viscous)
    const PylithScalar meanStress = 3.0 * bulkModulus * meanStrain;
    
    f1[0] = meanStress + 2.0 * shearModulus * elasticDevStrain[0]; // σ_xx
    f1[1] = meanStress + 2.0 * shearModulus * elasticDevStrain[1]; // σ_yy
    f1[2] = meanStress + 2.0 * shearModulus * elasticDevStrain[2]; // σ_zz
    f1[3] = 2.0 * shearModulus * elasticDevStrain[3]; // σ_xy
    
    if (dim == 3) {
        f1[4] = 2.0 * shearModulus * elasticDevStrain[4]; // σ_yz
        f1[5] = 2.0 * shearModulus * elasticDevStrain[5]; // σ_xz
    }
}

// Jf3 Jacobian kernel
void pylith::fekernels::MyViscoelastic::Jf3uu(/* ... */) {
    // Extract properties
    const PylithScalar shearModulus = a[aOff[AuxConstants::shearModulus]];
    const PylithScalar bulkModulus = a[aOff[AuxConstants::bulkModulus]];
    
    // Same as linear elasticity for implicit formulation
    // (viscous strain is a state variable, not part of solution)
    
    // Fill elasticity tensor
    Tensor::createElasticityTensor(J, shearModulus, bulkModulus, dim);
}

// Update viscous strain
void pylith::fekernels::MyViscoelastic::updateViscousStrain(/* ... */) {
    // Extract current state
    const PylithScalar shearModulus = a[aOff[AuxConstants::shearModulus]];
    const PylithScalar viscosity = a[aOff[AuxConstants::viscosity]];
    const PylithScalar maxwellTime = viscosity / shearModulus;
    
    const PylithScalar* viscousStrainOld = &a[aOff[AuxConstants::viscousStrain]];
    
    // Get time step (from constants)
    const PylithScalar dt = constants[0];
    
    // Compute total strain
    PylithScalar totalStrain[6];
    Tensor::symmetrize(totalStrain, s_x, dim);
    
    // Compute deviatoric strain
    const PylithScalar meanStrain = (totalStrain[0] + totalStrain[1] + totalStrain[2]) / 3.0;
    PylithScalar devStrain[6];
    for (int i = 0; i < 3; ++i) {
        devStrain[i] = totalStrain[i] - meanStrain;
    }
    for (int i = 3; i < 6; ++i) {
        devStrain[i] = totalStrain[i];
    }
    
    // Exponential integration
    const PylithScalar expFactor = exp(-dt / maxwellTime);
    
    for (int i = 0; i < 6; ++i) {
        stateVars[i] = expFactor * viscousStrainOld[i] + 
                       (1.0 - expFactor) * devStrain[i];
    }
}
```

#### Step 4: SWIG Interface

Create `modulesrc/materials/MyViscoelastic.i`:

```swig
// Interface file for MyViscoelastic

%module myviscoelastic

%include "exception.i"
%exception {
    try {
        $action
    } catch (const std::exception& err) {
        SWIG_exception(SWIG_RuntimeError, err.what());
    }
}

%include "typemaps.i"

%{
#include "pylith/materials/MyViscoelastic.hh"
%}

%include "pylith/materials/materialsfwd.hh"
%include "pylith/materials/MyViscoelastic.hh"
```

Update `modulesrc/materials/Makefile.am` to include the new interface.

#### Step 5: Python Wrapper

Create `pylith/materials/MyViscoelastic.py`:

```python
from pylith.materials.RheologyElasticity import RheologyElasticity
from .materials import MyViscoelastic as ModuleMyViscoelastic


class MyViscoelastic(RheologyElasticity, ModuleMyViscoelastic):
    """
    My custom viscoelastic rheology.
    
    Implements a Maxwell viscoelastic model with exponential
    integration for the viscous strain.
    
    Implements:
        pylith.materials.RheologyElasticity
    """
    
    import pythia.pyre.inventory
    
    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain?"
    
    def __init__(self, name="myviscoelastic"):
        """Constructor."""
        RheologyElasticity.__init__(self, name)
    
    def preinitialize(self, problem):
        """Pre-initialization."""
        RheologyElasticity.preinitialize(self, problem)
        
        # Set flag in C++
        ModuleMyViscoelastic.useReferenceState(self, self.useReferenceState)
    
    def _createModuleObj(self):
        """Create handle to C++ object."""
        ModuleMyViscoelastic.__init__(self)
```

#### Step 6: Configuration Example

User configuration file:

```cfg
[pylithapp.problem.materials.crust]
# Use custom material
rheology = pylith.materials.MyViscoelastic

[pylithapp.problem.materials.crust.rheology]
use_reference_state = False

# Spatial database for properties
db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.description = Viscoelastic properties
db_auxiliary_field.iohandler.filename = mat_viscoelastic.spatialdb
db_auxiliary_field.query_type = linear
```

Spatial database `mat_viscoelastic.spatialdb`:

```
#SPATIAL.ascii 1
SimpleDB {
  num-values = 5
  value-names = density vs vp viscosity maxwell_time
  value-units = kg/m**3 m/s m/s Pa*s s
  num-locs = 4
  data-dim = 1
  space-dim = 3
  cs-data = cartesian {
    to-meters = 1.0
    space-dim = 3
  }
}
// Columns: x y z density vs vp viscosity maxwell_time
0.0  0.0  0.0    2500.0  3000.0  5200.0  1.0e19  1.0e11
1.0  0.0  0.0    2500.0  3000.0  5200.0  1.0e19  1.0e11
0.0  1.0  0.0    2500.0  3000.0  5200.0  1.0e19  1.0e11
0.0  0.0  1.0    2500.0  3000.0  5200.0  1.0e19  1.0e11
```

## Adding a Custom Boundary Condition

### Example: Time-Dependent Velocity BC

#### C++ Implementation

`libsrc/pylith/bc/VelocityTimeDependent.hh/cc`

```cpp
class VelocityTimeDependent : public BoundaryCondition {
public:
    
    // Set constrained DOF
    void setConstraintDOF(const int* dof, const int numDOF);
    
    // Set spatial database for velocity
    void setDB(spatialdata::spatialdb::SpatialDB* db);
    
    // Create constraint
    std::vector<Constraint*> createConstraints();
    
    // Get kernel
    static void kernel(/* ... */);
};
```

#### Python Wrapper

`pylith/bc/VelocityTimeDependent.py`

```python
class VelocityTimeDependent(BoundaryCondition, ModuleVelocityTimeDependent):
    """
    Time-dependent velocity boundary condition.
    """
    
    import pythia.pyre.inventory
    
    constrainedDOF = pythia.pyre.inventory.list("constrained_dof", default=[0])
    constrainedDOF.meta['tip'] = "DOF to constrain (0=x, 1=y, 2=z)"
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pythia.pyre.inventory.facility("db_auxiliary_field", family="spatial_database", factory=SimpleDB)
    
    def preinitialize(self, problem):
        BoundaryCondition.preinitialize(self, problem)
        ModuleVelocityTimeDependent.setConstraintDOF(self, self.constrainedDOF)
        ModuleVelocityTimeDependent.setDB(self, self.db)
```

## Adding a Custom Fault Slip Function

### Example: Exponential Slip Function

#### C++ Kernel

`libsrc/pylith/fekernels/KinSrcExponential.hh/cc`

```cpp
class KinSrcExponential {
public:
    
    /** Compute slip at time t.
     *
     * slip(t) = final_slip * (1 - exp(-(t - t0)/tau))
     *
     * Auxiliary fields:
     *   [0] final_slip
     *   [1] tau (time constant)
     *   [2] t0 (initiation time)
     */
    static void slip(/* PetscPointFunc signature */);
    
    /** Compute slip rate at time t.
     */
    static void slipRate(/* ... */);
};

void KinSrcExponential::slip(/* ... */) {
    const PylithScalar finalSlip = a[aOff[0]];
    const PylithScalar tau = a[aOff[1]];
    const PylithScalar t0 = a[aOff[2]];
    
    if (t < t0) {
        f[0] = 0.0;
    } else {
        f[0] = finalSlip * (1.0 - exp(-(t - t0) / tau));
    }
}
```

#### Python Wrapper

`pylith/faults/KinSrcExponential.py`

```python
class KinSrcExponential(KinSrc, ModuleKinSrcExponential):
    """
    Exponential slip time function.
    """
    
    import pythia.pyre.inventory
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pythia.pyre.inventory.facility(
        "db_auxiliary_field",
        family="spatial_database",
        factory=SimpleDB
    )
    db.meta['tip'] = "Spatial database for slip parameters (final_slip, tau, t0)"
    
    def preinitialize(self, problem):
        KinSrc.preinitialize(self, problem)
```

## Adding Derived Fields

### Example: Compute Von Mises Stress

#### C++ Kernel

`libsrc/pylith/fekernels/DerivedFactoryElasticity.cc`

```cpp
// Add to factory
void DerivedFactoryElasticity::addVonMisesStress() {
    const char* subfieldName = "von_mises_stress";
    const char* components[] = {"von_mises"};
    const PylithInt numComponents = 1;
    const PylithScalar scale = _normalizerStress->getScale();
    const PylithInt basisOrder = 0;  // Constant per cell
    
    _field->subfieldAdd(subfieldName, numComponents, components, scale, basisOrder);
    _factory->setSubfieldQuery(subfieldName);
}

// Kernel implementation
void computeVonMisesStress(/* PetscPointFunc signature */) {
    // Extract stress components from auxiliary field
    const PylithScalar* stress = &a[aOff[stressIndex]];
    
    // σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz
    const PylithScalar sxx = stress[0];
    const PylithScalar syy = stress[1];
    const PylithScalar szz = stress[2];
    const PylithScalar sxy = stress[3];
    const PylithScalar syz = stress[4];
    const PylithScalar sxz = stress[5];
    
    // Von Mises: sqrt(3/2 * dev(σ):dev(σ))
    const PylithScalar meanStress = (sxx + syy + szz) / 3.0;
    const PylithScalar s11 = sxx - meanStress;
    const PylithScalar s22 = syy - meanStress;
    const PylithScalar s33 = szz - meanStress;
    
    const PylithScalar J2 = 0.5 * (s11*s11 + s22*s22 + s33*s33) + 
                            sxy*sxy + syz*syz + sxz*sxz;
    
    f[0] = sqrt(3.0 * J2);
}
```

#### Configuration

```cfg
[pylithapp.problem.materials.crust]
derived_subfields = [cauchy_stress, von_mises_stress]

derived_subfields.von_mises_stress.basis_order = 0
```

## Testing Your Extension

### Unit Tests (C++)

`tests/libtests/materials/TestMyViscoelastic.cc`

```cpp
#include "catch2/catch.hpp"

#include "pylith/materials/MyViscoelastic.hh"

TEST_CASE("MyViscoelastic::constructor", "[MyViscoelastic]") {
    pylith::materials::MyViscoelastic material;
    
    CHECK(!material.useReferenceState());
}

TEST_CASE("MyViscoelastic::useReferenceState", "[MyViscoelastic]") {
    pylith::materials::MyViscoelastic material;
    
    material.useReferenceState(true);
    CHECK(material.useReferenceState());
}

TEST_CASE("MyViscoelastic::addAuxiliarySubfields", "[MyViscoelastic]") {
    // Test that auxiliary fields are added correctly
    // ...
}
```

### Integration Tests (Python)

`tests/pytests/materials/TestMyViscoelastic.py`

```python
import unittest
from pylith.materials.MyViscoelastic import MyViscoelastic


class TestMyViscoelastic(unittest.TestCase):
    
    def test_constructor(self):
        """Test constructor."""
        material = MyViscoelastic()
        self.assertEqual(material.useReferenceState, False)
    
    def test_configure(self):
        """Test configuration."""
        material = MyViscoelastic()
        material.inventory.useReferenceState = True
        material._configure()
        self.assertEqual(material.useReferenceState, True)


if __name__ == "__main__":
    unittest.main()
```

### Full-Scale Test

`tests/fullscale/viscoelasticity/MyViscoelastic/`

```
MyViscoelastic/
├── mesh.msh (Gmsh mesh)
├── pylithapp.cfg (main config)
├── mat_viscoelastic.spatialdb (properties)
├── bc_*.spatialdb (boundary conditions)
├── output/ (expected output for comparison)
└── test_myviscoelastic.py (test driver)
```

Test driver:

```python
import unittest
from pylith.testing.FullTestApp import TestComponent


class TestMyViscoelastic(TestComponent):
    """Full-scale test of MyViscoelastic material."""
    
    def setUp(self):
        self.name = "myviscoelastic"
        self.simDir = __file__.replace("test_myviscoelastic.py", "")
    
    def test_run(self):
        """Run simulation and compare output."""
        self.run_pylith()
        self.check_output()


if __name__ == "__main__":
    unittest.main()
```

## Best Practices

### Code Organization

✅ **Separate concerns**: Kernels in `fekernels/`, classes in respective directories  
✅ **Document thoroughly**: Doxygen for C++, docstrings for Python  
✅ **Follow naming conventions**: CamelCase for classes, lowercase_underscore for methods  
✅ **Add unit tests**: Test each component independently  
✅ **Add integration tests**: Test end-to-end workflows  

### Performance Considerations

✅ **Minimize kernel complexity**: Keep pointwise functions simple and fast  
✅ **Avoid branching in kernels**: Use separate kernels for different cases  
✅ **Preallocate memory**: Use PETSc preallocation for matrices  
✅ **Profile code**: Use PETSc logging to identify bottlenecks  

### Debugging Tips

```bash
# Enable verbose logging
pylith --petsc.log_view

# Debug with gdb
gdb --args python $(which pylith) myconfig.cfg

# Check residual/Jacobian consistency
pylith myconfig.cfg --petsc.snes_test_jacobian

# Monitor solver
pylith myconfig.cfg --petsc.ksp_monitor --petsc.snes_monitor
```

### Contributing Back

To contribute your extension to PyLith:

1. **Fork repository**: https://github.com/geodynamics/pylith
2. **Create branch**: `git checkout -b feature/my-extension`
3. **Add tests**: Ensure full test coverage
4. **Document**: Add user documentation in `docs/`
5. **Submit PR**: Create pull request with description

## Summary

Adding extensions to PyLith:

✅ **Materials**: Inherit from `RheologyElasticity`, implement kernels  
✅ **Boundary Conditions**: Inherit from `BoundaryCondition`  
✅ **Fault Slip**: Inherit from `KinSrc`  
✅ **Derived Fields**: Add factory method + kernel  
✅ **Test**: Unit tests + integration tests  

Key steps:
1. C++ implementation (class + kernels)
2. SWIG interface
3. Python wrapper
4. Configuration example
5. Tests
6. Documentation

**Template repository** available at: https://github.com/geodynamics/pylith_templates

---

This completes the PyLith Code Architecture Guide. For questions, visit the [CIG community forum](https://community.geodynamics.org/c/pylith/).
