# PyLith Code Architecture Guide

Welcome to the comprehensive PyLith code architecture guide. This documentation provides an in-depth understanding of how PyLith works internally, its design patterns, and how to extend it.

## About PyLith

PyLith is an open-source finite-element code for dynamic and quasi-static simulations of crustal deformation, primarily focused on earthquakes and volcanoes. It is developed by the Computational Infrastructure for Geodynamics (CIG).

**Version**: 5.0.0dev  
**Primary Language**: C++ (core) with Python (interface)  
**Key Dependencies**: PETSc, MPI, HDF5, NetCDF, Pythia/Pyre

## Documentation Structure

This guide is organized into the following sections:

1. **[01-Architecture-Overview.md](01-Architecture-Overview.md)** - High-level architecture and design philosophy
2. **[02-Core-Components.md](02-Core-Components.md)** - Detailed description of core components
3. **[03-Data-Flow.md](03-Data-Flow.md)** - How data flows through the system
4. **[04-Build-System.md](04-Build-System.md)** - Build system and compilation process
5. **[05-Module-Structure.md](05-Module-Structure.md)** - Python-C++ integration via SWIG
6. **[06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md)** - FE assembly and solver integration
7. **[07-Extension-Guide.md](07-Extension-Guide.md)** - How to extend PyLith with custom components

## Quick Start for Developers

### Understanding the Codebase Layout

```
pylith/
├── libsrc/pylith/          # C++ core implementation
│   ├── bc/                 # Boundary conditions
│   ├── faults/             # Fault implementations
│   ├── feassemble/         # Finite element assembly
│   ├── fekernels/          # Finite element kernels (pointwise functions)
│   ├── materials/          # Material models and rheologies
│   ├── meshio/             # Mesh I/O and output
│   ├── problems/           # Problem formulations
│   ├── topology/           # Mesh topology management
│   └── utils/              # Utilities and common code
├── modulesrc/              # SWIG interface files (Python bindings)
├── pylith/                 # Python layer
│   ├── apps/               # Application entry points
│   ├── bc/                 # Python BC wrappers
│   ├── faults/             # Python fault wrappers
│   ├── materials/          # Python material wrappers
│   ├── problems/           # Python problem wrappers
│   └── utils/              # Python utilities
├── tests/                  # Comprehensive test suite
├── examples/               # Example problems
└── docs/                   # User documentation
```

### Key Design Principles

1. **Hybrid Architecture**: C++ core for performance, Python for configuration and scripting
2. **Component-Based Design**: Uses Pyre/Pythia component framework
3. **PETSc Integration**: Leverages PETSc for linear algebra and solvers
4. **Separation of Concerns**: Physics, discretization, and solution strategy are decoupled
5. **Extensibility**: Easy to add new materials, boundary conditions, and fault models

## Getting Started with the Code

### For Users Wanting to Understand Internals
Start with [01-Architecture-Overview.md](01-Architecture-Overview.md) to understand the big picture, then move to [02-Core-Components.md](02-Core-Components.md).

### For Developers Adding Features
Read [07-Extension-Guide.md](07-Extension-Guide.md) for specific instructions on adding new components, then refer to the relevant component documentation.

### For Contributors to Core
Study [03-Data-Flow.md](03-Data-Flow.md) and [06-Finite-Element-Implementation.md](06-Finite-Element-Implementation.md) to understand the execution model and numerical methods.

## Contributing

When contributing to PyLith, please:
1. Follow the existing code style and conventions
2. Add tests for new functionality
3. Update documentation
4. Run the full test suite before submitting changes

## Resources

- **Main Website**: https://geodynamics.org/resources/pylith
- **Documentation**: https://pylith.readthedocs.io
- **Source Code**: https://github.com/geodynamics/pylith
- **Community Forum**: https://community.geodynamics.org/c/pylith/

## Citation

If you use PyLith in your research, please cite:

```bibtex
@Manual{PyLith:software,
  title        = {PyLith Version 5.0.0dev},
  author       = {Aagaard, B. and Knepley, M. and Williams, C.},
  organization = {Computational Infrastructure for Geodynamics (CIG)},
  address      = {University of California, Davis},
  year         = {2023},
  doi          = {10.5281/zenodo.14635926}
}
```

## Authors

- Brad Aagaard (U.S. Geological Survey)
- Matthew Knepley (University at Buffalo)
- Charles Williams (GNS Science, New Zealand)

## License

PyLith is released under the MIT License. See LICENSE.md for details.

---

**Last Updated**: 2025-11-24  
**Document Version**: 1.0
