# Module Structure and Python-C++ Integration

## Table of Contents
1. [Module Organization](#module-organization)
2. [Python-C++ Bridge](#python-c-bridge)
3. [Component Lifecycle](#component-lifecycle)
4. [Pyre Integration](#pyre-integration)
5. [Example: Creating a Material](#example-creating-a-material)

## Module Organization

### Parallel Hierarchies

PyLith maintains parallel implementations in Python and C++:

```
Python Layer (pylith/)          C++ Layer (libsrc/pylith/)
├── apps/                       
│   └── PyLithApp.py           (Application entry)
├── problems/                   ├── problems/
│   ├── Problem.py             │   ├── Problem.hh/cc
│   └── TimeDependent.py       │   └── TimeDependent.hh/cc
├── materials/                  ├── materials/
│   ├── Material.py            │   ├── Material.hh/cc
│   └── Elasticity.py          │   └── Elasticity.hh/cc
├── bc/                         ├── bc/
│   └── DirichletTimeDependent.py  └── DirichletTimeDependent.hh/cc
└── ...                         └── ...

         Bridge (modulesrc/)
         ├── problems/
         │   ├── Problem.i (SWIG)
         │   └── TimeDependent.i
         ├── materials/
         │   ├── Material.i
         │   └── Elasticity.i
         └── ...
```

### Responsibilities by Layer

**Python Layer** (User-facing):
- Configuration management (Pyre inventory)
- Parameter validation
- Component factories
- High-level logic
- Documentation (docstrings)

**C++ Layer** (Computational):
- Numerical algorithms
- Finite element assembly
- PETSc integration
- Memory management
- Performance-critical code

**SWIG Bridge** (Glue):
- Type conversion (Python ↔ C++)
- Exception handling
- Ownership management

## Python-C++ Bridge

### Dual Inheritance Pattern

Most PyLith components use **dual inheritance**:

```python
# Python wrapper class
from pylith.problems.problems import Problem as ModuleProblem

class Problem(PetscComponent, ModuleProblem):
    """
    Python Problem class.
    Inherits from:
    - PetscComponent: Pyre component base (Python)
    - ModuleProblem: C++ implementation (via SWIG)
    """
    
    def __init__(self, name="problem"):
        PetscComponent.__init__(self, name)
        # ModuleProblem constructor called automatically
    
    def preinitialize(self, mesh):
        """Python-level pre-initialization."""
        self.inventory.scales.preinitialize()  # Python logic
        ModuleProblem.preinitialize(self, mesh)  # Call C++
```

### SWIG Type Mapping

**Basic Types**:
```
C++                  Python
int              →   int
double           →   float
const char*      →   str
bool             →   bool
```

**Arrays**:
```cpp
// C++
void setValues(double* values, int size);

// SWIG typemap (in modulesrc/include/scalartypemaps.i)
%apply (double* INPLACE_ARRAY1, int DIM1) {
    (double* values, int size)
};

# Python
values = numpy.array([1.0, 2.0, 3.0])
obj.setValues(values)  # NumPy array passed directly
```

**Objects**:
```cpp
// C++
class Mesh { ... };
void setMesh(const Mesh& mesh);

// SWIG automatically handles references/pointers

# Python
mesh = Mesh()
obj.setMesh(mesh)  # Python object wraps C++ object
```

### Memory Management

**Ownership Rules**:

1. **Python owns**: Object created in Python
   ```python
   mesh = Mesh()  # Python owns, deleted when no references
   ```

2. **C++ owns**: Object created in C++, returned to Python
   ```python
   mesh = problem.getMesh()  # C++ owns, Python just has reference
   ```

3. **Transfer ownership**: SWIG `%newobject` directive
   ```swig
   %newobject createMesh;
   Mesh* createMesh();  // Python takes ownership
   ```

**SWIG Smart Pointers**:
```cpp
// C++ uses std::shared_ptr
std::shared_ptr<Material> material;

// SWIG handles ref counting automatically
# Python
material = Material()  # Ref count = 1
problem.setMaterial(material)  # Ref count = 2
del material  # Ref count = 1, C++ still has it
```

## Component Lifecycle

### Lifecycle Phases

All PyLith components follow a standard lifecycle:

```
┌─────────────────────────────────────────────────────────┐
│ 1. CONSTRUCTION                                          │
│    Python: __init__()                                    │
│    C++: Constructor                                      │
│    • Create object                                       │
│    • Initialize member variables                         │
└────────────┬────────────────────────────────────────────┘
             ▼
┌─────────────────────────────────────────────────────────┐
│ 2. CONFIGURATION (_configure)                            │
│    Python: _configure()                                  │
│    • Read from Pyre inventory                            │
│    • Set object properties                               │
└────────────┬────────────────────────────────────────────┘
             ▼
┌─────────────────────────────────────────────────────────┐
│ 3. PREINITIALIZE (preinitialize)                         │
│    Python + C++: preinitialize()                         │
│    • Setup dependencies                                  │
│    • Validate configuration                              │
│    • Create auxiliary structures                         │
└────────────┬────────────────────────────────────────────┘
             ▼
┌─────────────────────────────────────────────────────────┐
│ 4. VERIFY (verifyConfiguration)                          │
│    C++: verifyConfiguration()                            │
│    • Check consistency                                   │
│    • Validate inter-component relationships              │
└────────────┬────────────────────────────────────────────┘
             ▼
┌─────────────────────────────────────────────────────────┐
│ 5. INITIALIZE (initialize)                               │
│    C++: initialize()                                     │
│    • Setup fields                                        │
│    • Query databases                                     │
│    • Prepare for computation                             │
└────────────┬────────────────────────────────────────────┘
             ▼
┌─────────────────────────────────────────────────────────┐
│ 6. RUN (run, poststep, etc.)                             │
│    C++: run(), poststep()                                │
│    • Main computation                                    │
│    • State updates                                       │
└────────────┬────────────────────────────────────────────┘
             ▼
┌─────────────────────────────────────────────────────────┐
│ 7. FINALIZE (finalize)                                   │
│    C++: finalize()                                       │
│    • Write final output                                  │
│    • Cleanup                                             │
└────────────┬────────────────────────────────────────────┘
             ▼
┌─────────────────────────────────────────────────────────┐
│ 8. DESTRUCTION (deallocate)                              │
│    C++: deallocate()                                     │
│    Python: __del__()                                     │
│    • Free memory                                         │
│    • Destroy PETSc objects                               │
└─────────────────────────────────────────────────────────┘
```

### Example: Material Lifecycle

```python
# 1. CONSTRUCTION
[pylithapp.problem]
materials = [crust, slab]

[pylithapp.problem.materials.crust]
# Component is instantiated here
# Python: Elasticity.__init__()
# C++: Elasticity constructor

# 2. CONFIGURATION
label_value = 1
# Python: _configure() reads this

# 3. PREINITIALIZE
# Python layer
def preinitialize(self, problem):
    # Setup rheology
    self.rheology.preinitialize(problem)
    # Call C++
    ModuleMaterial.preinitialize(self, problem)

# C++ layer
void Material::preinitialize(const Mesh& mesh) {
    // Register fields
    // Setup integrators
}

# 4-8. Continue through lifecycle...
```

## Pyre Integration

### Pyre Component System

**Pyre/Pythia** provides:
- Configuration management
- Component registry
- Factory pattern
- Property validation

### Component Definition

```python
from pythia.pyre.components import Component

class Material(Component):
    """
    Material component.
    """
    
    # Declare configurable properties
    import pythia.pyre.inventory
    
    # Simple property
    labelValue = pythia.pyre.inventory.int("label_value", default=1)
    labelValue.meta["tip"] = "Material ID"
    
    # Validated property
    formulation = pythia.pyre.inventory.str(
        "formulation",
        default="quasistatic",
        validator=pythia.pyre.inventory.choice(["quasistatic", "dynamic"])
    )
    
    # Sub-component (facility)
    from pylith.materials.IsotropicLinearElasticity import IsotropicLinearElasticity
    rheology = pythia.pyre.inventory.facility(
        "rheology",
        family="rheology_elasticity",
        factory=IsotropicLinearElasticity
    )
    rheology.meta["tip"] = "Bulk rheology for elasticity"
    
    # Lifecycle
    def __init__(self, name="material"):
        Component.__init__(self, name, facility="material")
    
    def _configure(self):
        """Transfer inventory to properties."""
        Component._configure(self)
        # Access configured values
        self.myLabelValue = self.labelValue
```

### Configuration Files

**Format**: `.cfg` files (INI-style)

```cfg
# Global settings
[pylithapp]
# ...

# Component hierarchy
[pylithapp.problem]
formulation = quasistatic

# Sub-component
[pylithapp.problem.materials]
slab = pylith.materials.Elasticity

# Sub-sub-component
[pylithapp.problem.materials.slab]
label_value = 2

# Facility (sub-component with type)
[pylithapp.problem.materials.slab.rheology]
use_reference_state = False

# Spatial database (another component type)
[pylithapp.problem.materials.slab.rheology.db_auxiliary_field]
description = Elastic properties for slab
iohandler.filename = mat_slab.spatialdb
query_type = linear
```

### Component Factories

**Purpose**: Create components by name

```python
def materialFactory(name):
    """Factory for material components."""
    from pythia.pyre.inventory import facility
    from pylith.materials.Elasticity import Elasticity
    
    return facility(name, family="material", factory=Elasticity)

# In Problem.py
import pythia.pyre.inventory

materials = pythia.pyre.inventory.facilityArray(
    "materials",
    itemFactory=materialFactory,
    factory=list
)

# User can then configure:
# [pylithapp.problem]
# materials = [crust, slab]
# [pylithapp.problem.materials.crust]
# # ... crust settings
```

### Property Validation

```python
# Validators
from pythia.pyre.inventory import choice, range, str

# Choice validator
solver = pythia.pyre.inventory.str(
    "solver",
    default="nonlinear",
    validator=choice(["linear", "nonlinear"])
)

# Range validator
maxIterations = pythia.pyre.inventory.int(
    "max_iterations",
    default=100,
    validator=range(1, 1000)
)

# Dimensional validator (with units)
from pythia.pyre.units.length import meter
width = pythia.pyre.inventory.dimensional(
    "width",
    default=10.0*meter
)

# In config file:
# width = 5.0*km  # Pyre converts to meters
```

## Example: Creating a Material

### Step-by-Step Implementation

#### 1. C++ Header (libsrc/pylith/materials/MyMaterial.hh)

```cpp
#pragma once

#include "pylith/materials/RheologyElasticity.hh"

class pylith::materials::MyMaterial : public RheologyElasticity {
public:
    
    MyMaterial();
    virtual ~MyMaterial();
    
    // Provide auxiliary field information
    void addAuxiliarySubfields();
    
    // Provide kernels
    std::vector<ResidualKernels> getKernelsResidual(
        const bool gravityField,
        const bool isJacobianTerm
    ) const;
    
    std::vector<JacobianKernels> getKernelsJacobian(
        const bool gravityField
    ) const;

private:
    
    // Kernel functions
    static void f1(/* ... */);  // Residual
    static void Jf3(/* ... */); // Jacobian
};
```

#### 2. C++ Implementation (libsrc/pylith/materials/MyMaterial.cc)

```cpp
#include "MyMaterial.hh"
#include "pylith/fekernels/MyMaterial.hh"  // Kernels

// Constructor
pylith::materials::MyMaterial::MyMaterial() {
    // Initialize
}

// Destructor
pylith::materials::MyMaterial::~MyMaterial() {
    deallocate();
}

// Specify auxiliary fields needed
void pylith::materials::MyMaterial::addAuxiliarySubfields() {
    RheologyElasticity::addAuxiliarySubfields();
    
    // Add custom fields
    _auxiliaryFactory->addMyProperty();
}

// Provide residual kernels
std::vector<RheologyElasticity::ResidualKernels>
pylith::materials::MyMaterial::getKernelsResidual(...) const {
    std::vector<ResidualKernels> kernels;
    
    ResidualKernels rKernel;
    rKernel.subfieldName = "displacement";
    rKernel.f0 = NULL;  // No f0 term
    rKernel.f1 = pylith::fekernels::MyMaterial::f1;  // Stress
    kernels.push_back(rKernel);
    
    return kernels;
}

// Similar for getKernelsJacobian...
```

#### 3. Kernel Implementation (libsrc/pylith/fekernels/MyMaterial.hh)

```cpp
#pragma once

#include "pylith/fekernels/fekernelsfwd.hh"

class pylith::fekernels::MyMaterial {
public:
    
    // Residual kernel
    static void f1(
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
    );
    
    // Jacobian kernel
    static void Jf3(/* ... */);
};
```

#### 4. SWIG Interface (modulesrc/materials/MyMaterial.i)

```swig
%module mymaterial

%include "pylith/materials/materialsfwd.hh"

%{
#include "pylith/materials/MyMaterial.hh"
%}

%include "pylith/materials/MyMaterial.hh"
```

#### 5. Python Wrapper (pylith/materials/MyMaterial.py)

```python
from pylith.materials.RheologyElasticity import RheologyElasticity
from .materials import MyMaterial as ModuleMyMaterial

class MyMaterial(RheologyElasticity, ModuleMyMaterial):
    """
    My custom material.
    
    Implements: (describe constitutive law)
    """
    
    import pythia.pyre.inventory
    
    # Configuration parameters
    myParameter = pythia.pyre.inventory.float("my_parameter", default=1.0)
    myParameter.meta["tip"] = "Description of parameter"
    
    def __init__(self, name="mymaterial"):
        """Constructor."""
        RheologyElasticity.__init__(self, name)
    
    def _configure(self):
        """Set members from inventory."""
        RheologyElasticity._configure(self)
        # Could call C++ setter if needed
        # ModuleMyMaterial.setParameter(self, self.myParameter)
    
    def preinitialize(self, problem):
        """Pre-initialization."""
        RheologyElasticity.preinitialize(self, problem)
        # C++ preinitialize called automatically via inheritance
```

#### 6. Usage (Configuration)

```cfg
[pylithapp.problem.materials.custom]
# Specify custom material type
rheology = pylith.materials.MyMaterial

[pylithapp.problem.materials.custom.rheology]
my_parameter = 2.5

db_auxiliary_field = spatialdata.spatialdb.SimpleDB
db_auxiliary_field.iohandler.filename = my_properties.spatialdb
```

## Data Transfer Patterns

### Python → C++

**Configuration Parameters**:
```python
# Python
self.labelValue = 2

# Transfer in preinitialize()
ModuleMaterial.setLabelValue(self, self.labelValue)

// C++
void Material::setLabelValue(const int value) {
    _labelValue = value;
}
```

**Complex Objects** (Mesh, Fields):
```python
# Python creates, C++ uses
mesh = Mesh()
# ... setup mesh ...

problem.preinitialize(mesh)
# → C++ receives reference to mesh

// C++
void Problem::preinitialize(const Mesh& mesh) {
    _mesh = &mesh;  // Store reference
    // Use mesh throughout C++ code
}
```

### C++ → Python

**Query Results**:
```cpp
// C++
int Problem::getDimension() const {
    return _mesh->getDimension();
}

# Python
dim = problem.getDimension()  # SWIG handles call
```

**Arrays** (via NumPy):
```cpp
// C++ (with SWIG typemap)
void getCoordinates(double** coords, int* size) const {
    // Return pointer and size
}

# Python
coords = problem.getCoordinates()  # Returns NumPy array
```

## Advanced SWIG Topics

### Exception Handling

```swig
%exception {
    try {
        $action
    } catch (const std::exception& err) {
        PyErr_SetString(PyExc_RuntimeError, err.what());
        SWIG_fail;
    }
}
```

### Directors (Callbacks from C++ to Python)

```swig
%module(directors="1") mymodule

%feature("director") MyBase;

class MyBase {
public:
    virtual void virtualMethod() = 0;
};

# Python can subclass and override
class MyDerived(MyBase):
    def virtualMethod(self):
        print("Python implementation")

# C++ can call Python implementation
obj = MyDerived()
obj.virtualMethod()  # Calls Python
```

### Nested Classes

```swig
%nestedworkaround Outer::Inner;

class Outer {
public:
    class Inner {
        // ...
    };
};

%unnestedworkaround Outer::Inner;
```

## Summary

PyLith's module structure:

✅ **Parallel Hierarchies**: Python + C++ + SWIG bridge  
✅ **Dual Inheritance**: Combines Pyre components with C++ implementation  
✅ **Lifecycle**: Standardized initialization sequence  
✅ **Pyre Integration**: Configuration management and validation  
✅ **SWIG**: Automatic binding generation  

Key patterns:
- **Configuration**: Pyre inventory → Python → C++
- **Computation**: C++ (performance-critical)
- **Extension**: Subclass in Python or C++ (or both)

Next: [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md) - FE assembly details
