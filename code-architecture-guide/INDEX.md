# PyLith Code Architecture Guide - Complete Index

## Overview

This index provides a complete reference to all topics covered in the PyLith Code Architecture Guide.

## Getting Started

### New Users
1. Start with [README.md](README.md) for an overview
2. Read [01-Architecture-Overview.md](01-Architecture-Overview.md) to understand the big picture
3. Use [QUICK-REFERENCE.md](QUICK-REFERENCE.md) for common tasks

### Developers Adding Features
1. Review [02-Core-Components.md](02-Core-Components.md) for component details
2. Study [07-Extension-Guide.md](07-Extension-Guide.md) for implementation patterns
3. Check [QUICK-REFERENCE.md](QUICK-REFERENCE.md) for checklists

### Contributors to Core
1. Understand [03-Data-Flow.md](03-Data-Flow.md) for system execution
2. Deep dive into [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md)
3. Review [04-Build-System.md](04-Build-System.md) and [05-Module-Structure.md](05-Module-Structure.md)

## Complete Topic Index

### A

- **Absorbing Boundary Conditions** → [02-Core-Components.md](02-Core-Components.md#boundary-condition-components)
- **Adaptive Time Stepping** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#time-integration)
- **Assembly Process** → [03-Data-Flow.md](03-Data-Flow.md#assembly-data-flow), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#assembly-process)
- **Autotools** → [04-Build-System.md](04-Build-System.md#autotools-configuration)
- **Auxiliary Fields** → [02-Core-Components.md](02-Core-Components.md#material-components), [03-Data-Flow.md](03-Data-Flow.md#field-data-flow)

### B

- **Backward Euler** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#implicit-time-stepping)
- **Boundary Conditions** → [02-Core-Components.md](02-Core-Components.md#boundary-condition-components), [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-boundary-condition)
- **Build System** → [04-Build-System.md](04-Build-System.md), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#building)
- **Brune Slip Function** → [02-Core-Components.md](02-Core-Components.md#faultcohesivekin)

### C

- **C++ Core** → [01-Architecture-Overview.md](01-Architecture-Overview.md#layer-architecture), [02-Core-Components.md](02-Core-Components.md)
- **Cohesive Cells** → [02-Core-Components.md](02-Core-Components.md#faultcohesive-base-class)
- **Component Lifecycle** → [05-Module-Structure.md](05-Module-Structure.md#component-lifecycle), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#component-lifecycle)
- **Configuration Files** → [05-Module-Structure.md](05-Module-Structure.md#configuration-files), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#configuration-file-syntax)
- **Constitutive Models** → [02-Core-Components.md](02-Core-Components.md#material-components)
- **Constraints** → [02-Core-Components.md](02-Core-Components.md#dirichlettimedependent)

### D

- **Data Flow** → [03-Data-Flow.md](03-Data-Flow.md)
- **Data Structures** → [03-Data-Flow.md](03-Data-Flow.md#mesh-data-flow), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#petsc-integration)
- **Debugging** → [07-Extension-Guide.md](07-Extension-Guide.md#best-practices), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#debugging-tips)
- **Derived Fields** → [02-Core-Components.md](02-Core-Components.md#material-components), [07-Extension-Guide.md](07-Extension-Guide.md#adding-derived-fields)
- **Dirichlet BC** → [02-Core-Components.md](02-Core-Components.md#dirichlettimedependent)
- **Discretization** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#discretization)
- **DMPlex** → [01-Architecture-Overview.md](01-Architecture-Overview.md#mesh-and-topology-subsystem), [03-Data-Flow.md](03-Data-Flow.md#mesh-data-flow), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#petsc-integration)
- **Dual Inheritance** → [05-Module-Structure.md](05-Module-Structure.md#dual-inheritance-pattern)
- **Dynamic Simulation** → [01-Architecture-Overview.md](01-Architecture-overview.md#problem-formulation), [02-Core-Components.md](02-Core-Components.md#timedependent)

### E

- **Elasticity** → [02-Core-Components.md](02-Core-Components.md#elasticity), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#example-elasticity)
- **Exception Handling** → [05-Module-Structure.md](05-Module-Structure.md#advanced-swig-topics)
- **EXODUS Format** → [01-Architecture-Overview.md](01-Architecture-Overview.md#io-subsystem)
- **Extension Guide** → [07-Extension-Guide.md](07-Extension-Guide.md)

### F

- **Factory Pattern** → [01-Architecture-Overview.md](01-Architecture-Overview.md#design-philosophy), [05-Module-Structure.md](05-Module-Structure.md#component-factories)
- **Faults** → [02-Core-Components.md](02-Core-Components.md#fault-components)
- **FE Assembly** → [01-Architecture-Overview.md](01-Architecture-Overview.md#finite-element-assembly-subsystem), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md)
- **Field Layout** → [03-Data-Flow.md](03-Data-Flow.md#field-data-flow), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#dmplex-field-layout)
- **Finite Element** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md)

### G

- **Green's Functions** → [02-Core-Components.md](02-Core-Components.md#greensfns)
- **Gravity** → [02-Core-Components.md](02-Core-Components.md#material-base-class)

### H

- **HDF5 Output** → [01-Architecture-Overview.md](01-Architecture-Overview.md#io-subsystem), [03-Data-Flow.md](03-Data-Flow.md#output-data-flow)

### I

- **Implicit Time Stepping** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#implicit-time-stepping)
- **Incompressible Elasticity** → [02-Core-Components.md](02-Core-Components.md#material-hierarchy)
- **Initial Conditions** → [02-Core-Components.md](02-Core-Components.md#problem-components)
- **Integrator** → [02-Core-Components.md](02-Core-Components.md#integrator-base-class), [03-Data-Flow.md](03-Data-Flow.md#assembly-data-flow)

### J

- **Jacobian** → [03-Data-Flow.md](03-Data-Flow.md#assembly-data-flow), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#jacobian-kernels)

### K

- **Kernels** → [02-Core-Components.md](02-Core-Components.md#kernels), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#kernel-architecture), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#kernels)
- **KinSrc** → [02-Core-Components.md](02-Core-Components.md#faultcohesivekin), [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-fault-slip-function)

### L

- **Lagrange Multipliers** → [02-Core-Components.md](02-Core-Components.md#faultcohesive-base-class)
- **Layer Architecture** → [01-Architecture-Overview.md](01-Architecture-Overview.md#layer-architecture)
- **Lifecycle** → [05-Module-Structure.md](05-Module-Structure.md#component-lifecycle)
- **Linear Elasticity** → [02-Core-Components.md](02-Core-Components.md#elasticity), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#example-elasticity)
- **Libtool** → [04-Build-System.md](04-Build-System.md#libtool-libraries)

### M

- **Makefile** → [04-Build-System.md](04-Build-System.md#directory-structure)
- **Material Properties** → [02-Core-Components.md](02-Core-Components.md#material-components), [03-Data-Flow.md](03-Data-Flow.md#auxiliary-field-population-from-spatial-database)
- **Materials** → [02-Core-Components.md](02-Core-Components.md#material-components), [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-material)
- **Matrix Assembly** → [03-Data-Flow.md](03-Data-Flow.md#matrix-assembly-jacobian)
- **Maxwell Viscoelasticity** → [02-Core-Components.md](02-Core-Components.md#viscoelastic-materials)
- **Memory Management** → [05-Module-Structure.md](05-Module-Structure.md#memory-management)
- **Mesh** → [02-Core-Components.md](02-Core-Components.md#mesh), [03-Data-Flow.md](03-Data-Flow.md#mesh-data-flow)
- **MPI** → [01-Architecture-Overview.md](01-Architecture-Overview.md#external-libraries-layer)
- **Module Structure** → [05-Module-Structure.md](05-Module-Structure.md)

### N

- **Neumann BC** → [02-Core-Components.md](02-Core-Components.md#neumanntimedependent)
- **Newton's Method** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#implicit-time-stepping)
- **Nondimensionalization** → [02-Core-Components.md](02-Core-Components.md#problem-base-class)

### O

- **Observer Pattern** → [03-Data-Flow.md](03-Data-Flow.md#observer-pattern-for-output)
- **Output** → [01-Architecture-Overview.md](01-Architecture-Overview.md#io-subsystem), [03-Data-Flow.md](03-Data-Flow.md#output-data-flow), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#output-configuration)

### P

- **Parallel** → [01-Architecture-Overview.md](01-Architecture-Overview.md#mesh-and-topology-subsystem), [03-Data-Flow.md](03-Data-Flow.md#mesh-data-flow)
- **Performance** → [07-Extension-Guide.md](07-Extension-Guide.md#performance-considerations), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#performance-optimization)
- **PETSc** → [01-Architecture-Overview.md](01-Architecture-Overview.md#external-libraries-layer), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#petsc-integration), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#petsc-integration)
- **Physics** → [02-Core-Components.md](02-Core-Components.md#material-components)
- **Pointwise Functions** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#pointwise-functions)
- **Poroelasticity** → [02-Core-Components.md](02-Core-Components.md#material-hierarchy)
- **Preallocation** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#matrix-preallocation)
- **Problem** → [02-Core-Components.md](02-Core-Components.md#problem-components)
- **Pyre** → [01-Architecture-Overview.md](01-Architecture-Overview.md#component-based-architecture), [05-Module-Structure.md](05-Module-Structure.md#pyre-integration)
- **Python Wrapper** → [05-Module-Structure.md](05-Module-Structure.md#module-organization), [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-material)

### Q

- **Quadrature** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#assembly-process)
- **Quasi-static** → [01-Architecture-Overview.md](01-Architecture-Overview.md#problem-formulation), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#quasi-static-approximation)
- **Quick Reference** → [QUICK-REFERENCE.md](QUICK-REFERENCE.md)

### R

- **Residual** → [03-Data-Flow.md](03-Data-Flow.md#residual-assembly), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#residual-kernels)
- **Rheology** → [02-Core-Components.md](02-Core-Components.md#material-components)

### S

- **Scales** → [02-Core-Components.md](02-Core-Components.md#problem-base-class)
- **Section (PETSc)** → [02-Core-Components.md](02-Core-Components.md#field), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#dmplex-field-layout)
- **Slip Functions** → [02-Core-Components.md](02-Core-Components.md#faultcohesivekin), [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-fault-slip-function)
- **SNES** → [03-Data-Flow.md](03-Data-Flow.md#solver-data-flow), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#petsc-integration)
- **Solution Field** → [02-Core-Components.md](02-Core-Components.md#field), [03-Data-Flow.md](03-Data-Flow.md#field-data-flow)
- **Solver** → [03-Data-Flow.md](03-Data-Flow.md#solver-data-flow), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#common-petsc-options)
- **Spatial Database** → [03-Data-Flow.md](03-Data-Flow.md#auxiliary-field-population-from-spatial-database), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#spatial-databases)
- **State Variables** → [02-Core-Components.md](02-Core-Components.md#viscoelastic-materials), [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-material)
- **Stress** → [02-Core-Components.md](02-Core-Components.md#elasticity), [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#residual-kernels)
- **SWIG** → [04-Build-System.md](04-Build-System.md#swig-integration), [05-Module-Structure.md](05-Module-Structure.md#swig-type-mapping)

### T

- **Template Method Pattern** → [01-Architecture-Overview.md](01-Architecture-Overview.md#design-philosophy)
- **Testing** → [04-Build-System.md](04-Build-System.md#testing-system), [07-Extension-Guide.md](07-Extension-Guide.md#testing-your-extension)
- **Time Stepping** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#time-integration)
- **TimeDependent** → [02-Core-Components.md](02-Core-Components.md#timedependent)
- **Topology** → [02-Core-Components.md](02-Core-Components.md#topology-and-field-components), [03-Data-Flow.md](03-Data-Flow.md#mesh-data-flow)
- **Traction** → [02-Core-Components.md](02-Core-Components.md#neumanntimedependent)

### U

- **Unit Tests** → [04-Build-System.md](04-Build-System.md#testing-system), [07-Extension-Guide.md](07-Extension-Guide.md#testing-your-extension)

### V

- **Validation** → [05-Module-Structure.md](05-Module-Structure.md#property-validation)
- **Vectors** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#petsc-integration)
- **Viscoelasticity** → [02-Core-Components.md](02-Core-Components.md#viscoelastic-materials), [07-Extension-Guide.md](07-Extension-Guide.md#step-by-step-custom-viscoelastic-material)
- **Visitors** → [02-Core-Components.md](02-Core-Components.md#visitors)
- **VTK Output** → [01-Architecture-Overview.md](01-Architecture-Overview.md#io-subsystem), [03-Data-Flow.md](03-Data-Flow.md#output-data-flow)

### W

- **Weak Form** → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#weak-form-formulation), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#finite-element-assembly)

### X

- **Xdmf** → [01-Architecture-Overview.md](01-Architecture-Overview.md#io-subsystem), [03-Data-Flow.md](03-Data-Flow.md#solution-output-pipeline)

## Diagrams and Visualizations

### Architecture Diagrams
- High-level architecture → [01-Architecture-Overview.md](01-Architecture-Overview.md#high-level-architecture)
- Layer architecture → [01-Architecture-Overview.md](01-Architecture-Overview.md#layer-architecture)
- Material hierarchy → [02-Core-Components.md](02-Core-Components.md#material-hierarchy)

### Data Flow Diagrams
- Complete simulation flow → [03-Data-Flow.md](03-Data-Flow.md#complete-flow-diagram)
- Mesh import flow → [03-Data-Flow.md](03-Data-Flow.md#mesh-import-and-processing)
- Field creation flow → [03-Data-Flow.md](03-Data-Flow.md#field-creation-and-population)
- Assembly flow → [03-Data-Flow.md](03-Data-Flow.md#residual-assembly)
- Solver flow → [03-Data-Flow.md](03-Data-Flow.md#petsc-integration)

### Process Diagrams
- Build process → [04-Build-System.md](04-Build-System.md#what-happens-during-build)
- Component lifecycle → [05-Module-Structure.md](05-Module-Structure.md#lifecycle-phases)

## Code Examples

### Complete Implementations
- Custom material → [07-Extension-Guide.md](07-Extension-Guide.md#step-by-step-custom-viscoelastic-material)
- Custom boundary condition → [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-boundary-condition)
- Custom slip function → [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-fault-slip-function)
- Derived fields → [07-Extension-Guide.md](07-Extension-Guide.md#adding-derived-fields)

### Code Snippets
- Kernel functions → [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md#residual-kernels)
- Python wrappers → [05-Module-Structure.md](05-Module-Structure.md#dual-inheritance-pattern)
- Configuration files → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#configuration-file-syntax)

## Configuration Examples

- Material configuration → [02-Core-Components.md](02-Core-Components.md#material-configuration-python)
- Boundary condition configuration → [02-Core-Components.md](02-Core-Components.md#dirichlettimedependent)
- Fault configuration → [02-Core-Components.md](02-Core-Components.md#faultcohesivekin)
- Output configuration → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#output-configuration)
- PETSc options → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#common-petsc-options)

## Best Practices

- Code organization → [07-Extension-Guide.md](07-Extension-Guide.md#best-practices)
- Performance considerations → [07-Extension-Guide.md](07-Extension-Guide.md#performance-considerations)
- Debugging tips → [07-Extension-Guide.md](07-Extension-Guide.md#debugging-tips), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#debugging-tips)

## Reference Materials

- Command reference → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#common-commands)
- Class hierarchy → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#key-classes)
- Configuration syntax → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#configuration-file-syntax)
- PETSc options → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#common-petsc-options)
- File naming conventions → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#file-naming-conventions)

## Common Tasks

- Building PyLith → [04-Build-System.md](04-Build-System.md#compilation-process), [QUICK-REFERENCE.md](QUICK-REFERENCE.md#building)
- Running simulations → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#running-pylith)
- Adding new materials → [07-Extension-Guide.md](07-Extension-Guide.md#adding-a-custom-material)
- Debugging problems → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#debugging-tips)
- Testing code → [07-Extension-Guide.md](07-Extension-Guide.md#testing-your-extension)

## Troubleshooting

- Common errors → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#common-errors-and-solutions)
- Build issues → [04-Build-System.md](04-Build-System.md#common-build-issues)
- Performance problems → [QUICK-REFERENCE.md](QUICK-REFERENCE.md#performance-optimization)

---

**Navigation Tips**:
- Use your browser's search function (Ctrl+F / Cmd+F) to find specific topics
- Follow hyperlinks to jump directly to relevant sections
- Start with the README for an overview of the documentation structure

**Last Updated**: 2025-11-24
