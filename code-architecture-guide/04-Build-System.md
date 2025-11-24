# Build System

## Table of Contents
1. [Build System Overview](#build-system-overview)
2. [Autotools Configuration](#autotools-configuration)
3. [Directory Structure](#directory-structure)
4. [Compilation Process](#compilation-process)
5. [SWIG Integration](#swig-integration)
6. [Testing System](#testing-system)

## Build System Overview

PyLith uses **GNU Autotools** (Autoconf, Automake, Libtool) for its build system. This provides:
- Cross-platform compatibility
- Dependency checking
- Shared library management
- Parallel build support

### Build Dependencies

**Required**:
- C/C++ compiler (GCC 7+ or Clang 9+, C++14 support)
- Python 3.8+
- PETSc 3.24.0+
- MPI (MPICH, OpenMPI, etc.)
- Pythia/Pyre 1.1.0+
- spatialdata 3.1.0+
- Proj 6+
- SWIG 4.0+ (if building from source with `--enable-swig`)

**Optional**:
- HDF5 4.0+ (for HDF5/Xdmf output)
- NetCDF 4+ (for Cubit EXODUS format)
- Catch2 (for C++ unit testing)

### Configuration Hierarchy

```
configure.ac (Autoconf input)
     │
     ├─> configure (generated script)
     │    └─> Makefile.in → Makefile (per directory)
     │
     └─> portinfo (config header)
```

## Autotools Configuration

### configure.ac Structure

**Location**: `/workspace/configure.ac`

**Key Sections**:

```bash
# 1. Initialization
AC_INIT([PyLith], [5.0.0dev], [https://geodynamics.org/resources/pylith])
AM_INIT_AUTOMAKE([foreign subdir-objects tar-pax])

# 2. Compiler checks
AC_PROG_CC
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(14)  # Require C++14

# 3. Python configuration
AM_PATH_PYTHON([3.8])
CIT_PYTHON_SYSCONFIG
CIT_CHECK_PYTHON_HEADER
CIT_PYTHON_MODULE([pythia],[1.1.0])
CIT_PYTHON_MODULE([numpy],[1.20.0])

# 4. MPI
AC_SEARCH_LIBS([MPI_Init], [mpi mpich])
AC_CHECK_HEADER([mpi.h])

# 5. PETSc
CIT_PATH_PETSC([3.24.0])
CIT_HEADER_PETSC
CIT_CHECK_LIB_PETSC

# 6. Optional features
AC_ARG_ENABLE([cubit], ...)  # Cubit support
AC_ARG_ENABLE([hdf5], ...)   # HDF5 support
AC_ARG_ENABLE([swig], ...)   # SWIG module generation
AC_ARG_ENABLE([testing], ...)  # Unit testing

# 7. Generate Makefiles
AC_CONFIG_FILES([
    Makefile
    pylith/Makefile
    libsrc/Makefile
    libsrc/pylith/Makefile
    # ... (many more)
])
```

### Configuration Options

Common configuration commands:

```bash
# Basic configuration (assumes dependencies in standard locations)
./configure

# Specify Python
./configure PYTHON=/usr/bin/python3

# Enable SWIG (regenerate bindings)
./configure --enable-swig

# Enable testing
./configure --enable-testing

# Enable HDF5 and Cubit support
./configure --enable-hdf5 --enable-cubit

# Specify PETSc location
./configure --with-petsc-dir=/path/to/petsc --with-petsc-arch=arch-name

# Enable sanitizers (debugging)
./configure --enable-sanitizer=address

# Disable optimization (debugging)
./configure CXXFLAGS="-g -O0"
```

### Custom Autoconf Macros

**Location**: `m4/` directory

PyLith defines custom macros for checking dependencies:

```m4
# Example: CIT_PATH_PETSC
# Checks for PETSc installation and sets variables
AC_DEFUN([CIT_PATH_PETSC], [
  AC_ARG_WITH([petsc-dir], ...)
  # Check for petsc.h
  # Check for libpetsc
  # Set PETSC_CPPFLAGS, PETSC_LDFLAGS, PETSC_LIBS
])
```

## Directory Structure

### Source Tree

```
pylith/
├── configure.ac              # Autoconf input
├── Makefile.am               # Top-level Automake input
│
├── libsrc/                   # C++ core implementation
│   ├── Makefile.am
│   └── pylith/
│       ├── Makefile.am
│       ├── bc/
│       │   ├── Makefile.am
│       │   ├── *.cc         # Implementation
│       │   └── *.hh         # Headers
│       ├── materials/
│       ├── problems/
│       └── ...
│
├── modulesrc/                # SWIG interfaces
│   ├── Makefile.am
│   └── bc/
│       ├── Makefile.am
│       └── *.i              # SWIG interface files
│
├── pylith/                   # Python package
│   ├── Makefile.am
│   ├── __init__.py
│   ├── apps/
│   ├── bc/
│   └── ...
│
├── tests/                    # Test suite
│   ├── libtests/            # C++ unit tests
│   ├── pytests/             # Python unit tests
│   ├── mmstests/            # Method of manufactured solutions
│   └── fullscale/           # Integration tests
│
└── examples/                 # Example simulations
```

### Makefile.am Example

**libsrc/pylith/materials/Makefile.am**:

```makefile
# Subdirectories (none for this directory)
subdir = libsrc/pylith/materials

# Headers to install
include_HEADERS = \
    Material.hh \
    Elasticity.hh \
    IsotropicLinearElasticity.hh \
    # ... more headers

# Library to build
noinst_LTLIBRARIES = libmaterials.la

# Library sources
libmaterials_la_SOURCES = \
    Material.cc \
    Elasticity.cc \
    IsotropicLinearElasticity.cc \
    # ... more sources

# Compiler flags
AM_CPPFLAGS = \
    -I$(top_srcdir)/libsrc \
    $(PETSC_CPPFLAGS) \
    $(PYTHON_CPPFLAGS)

# Linker flags
libmaterials_la_LIBADD = \
    $(top_builddir)/libsrc/pylith/feassemble/libfeassemble.la
```

## Compilation Process

### Full Build Sequence

```bash
# 1. Generate configure script (if building from git)
autoreconf --install --verbose --force

# 2. Configure
./configure [OPTIONS]

# 3. Compile (parallel build with 4 cores)
make -j4

# 4. Install (optional)
make install

# 5. Run tests (optional)
make check
```

### What Happens During Build

```
make
 │
 ├─> Build libsrc/pylith/*.la
 │   ├─ Compile *.cc → *.lo (libtool objects)
 │   ├─ Link into *.la (libtool libraries)
 │   │  Order: utils, topology, feassemble, materials, bc, faults, problems
 │   └─ Create libpylith.la (master library)
 │
 ├─> Build modulesrc/* (SWIG)
 │   ├─ Run SWIG: *.i → *_wrap.cxx + *.py
 │   ├─ Compile *_wrap.cxx → *_wrap.lo
 │   ├─ Link into *.la (Python extension modules)
 │   └─ Install *.la and *.py to site-packages
 │
 └─> Install Python files (pylith/*.py)
     └─ Copy to site-packages/pylith/
```

### Libtool Libraries

PyLith builds **libtool libraries** (`.la` files):

- **Shared libraries** (`.so` on Linux, `.dylib` on macOS)
- **Static libraries** (`.a`) - disabled by default
- **Python extension modules** (`.so` on Linux)

**Linking Order** (dependencies):

```
libpylith.la
├── libproblems.la
│   ├── libmaterials.la
│   │   └── libfeassemble.la
│   │       └── libtopology.la
│   │           └── libutils.la
│   ├── libbc.la
│   └── libfaults.la
├── libmeshio.la
└── external libs (PETSc, MPI, HDF5, etc.)
```

### Parallel Build

Automake supports parallel compilation:

```bash
make -j$(nproc)  # Use all CPU cores
make -j4         # Use 4 cores
```

**Build times** (approximate, depends on system):
- Fresh build: 5-15 minutes
- Incremental build: seconds to minutes

## SWIG Integration

### SWIG Purpose

**SWIG** (Simplified Wrapper and Interface Generator) generates Python bindings for C++ code.

**Input**: `.i` files (SWIG interface definitions)  
**Output**: `*_wrap.cxx` (C++ wrapper) + `*.py` (Python module)

### SWIG Interface Example

**modulesrc/materials/Material.i**:

```swig
// SWIG directives
%module material  // Python module name

// Include necessary headers
%{
#include "pylith/materials/Material.hh"
%}

// Import base classes
%import "pylith/problems/Physics.i"

// Namespace
%include "pylith/materials/materialsfwd.hh"

// Parse header to generate bindings
%include "pylith/materials/Material.hh"
```

### SWIG Build Process

```
Material.i
    │
    ├─> SWIG → material_wrap.cxx + material.py
    │            │
    │            ├─> g++ → material_wrap.o
    │            │          │
    │            │          └─> Link with libpylith.la
    │            │              → _material.so (Python extension)
    │            │
    │            └─> Install material.py
    │
    └─> Python import:
        from pylith.materials import material
        from pylith.materials.material import Material
```

### Makefile.am for SWIG

**modulesrc/materials/Makefile.am**:

```makefile
if ENABLE_SWIG
  # Generate wrapper with SWIG
  %_wrap.cxx %.py: %.i
      $(SWIG) -c++ -python -I$(top_srcdir)/modulesrc/include $<
endif

# Python extension module
pyexec_LTLIBRARIES = _materials.la

_materials_la_SOURCES = materials_wrap.cxx
_materials_la_LDFLAGS = -module -avoid-version
_materials_la_LIBADD = $(top_builddir)/libsrc/pylith/libpylith.la

# Install Python files
nobase_python_PYTHON = \
    __init__.py \
    materials.py
```

### Why SWIG?

**Advantages**:
- Automatic binding generation
- Handles complex C++ features (templates, inheritance, STL)
- Reduces boilerplate code
- Maintains type safety

**Alternative**: Could use pybind11, Cython, or manual CPython API (more work)

## Testing System

### Test Structure

```
tests/
├── libtests/              # C++ unit tests (Catch2)
│   ├── bc/
│   ├── materials/
│   ├── topology/
│   └── ...
│
├── pytests/               # Python unit tests (unittest)
│   ├── bc/
│   ├── materials/
│   └── ...
│
├── mmstests/              # Method of Manufactured Solutions
│   ├── linearelasticity/
│   ├── incompressibleelasticity/
│   └── poroelasticity/
│
└── fullscale/             # Full-scale integration tests
    ├── linearelasticity/
    ├── viscoelasticity/
    └── poroelasticity/
```

### Running Tests

```bash
# All tests
make check

# Specific test suite
make check -C tests/pytests
make check -C tests/fullscale/linearelasticity

# Individual test (pytest)
cd tests/pytests
pytest bc/TestDirichletTimeDependent.py

# C++ test (Catch2)
cd tests/libtests
./test_materials --list-tests
./test_materials "[Material]"
```

### Test Makefile Example

**tests/pytests/Makefile.am**:

```makefile
TESTS = test_pylith.py

test_pylith.py:
    @echo "#!/bin/bash" > $@
    @echo "pytest --verbose --junit-xml=junit.xml" >> $@
    @chmod +x $@

CLEANFILES = test_pylith.py

# Define test environment
TESTS_ENVIRONMENT = \
    PYTHONPATH=$(top_builddir):$(PYTHONPATH) \
    PATH=$(top_builddir)/applications:$(PATH)
```

### Test Coverage

Enable test coverage:

```bash
./configure --enable-test-coverage

make check

# Generate coverage report
make coverage-libtests  # C++ coverage
make coverage-pytests   # Python coverage
```

**Tools**:
- **C++**: lcov + genhtml
- **Python**: coverage.py

Coverage reports generated in:
- `tests/libtests/coverage/`
- `tests/pytests/htmlcov/`

### Continuous Integration

**GitHub Actions** (`.github/workflows/ci-main.yml`):

```yaml
name: CI

on: [push, pull_request]

jobs:
  build-and-test:
    runs-on: ubuntu-latest
    
    steps:
      - uses: actions/checkout@v2
      
      - name: Install dependencies
        run: |
          # Install PETSc, MPI, Python, etc.
      
      - name: Configure
        run: |
          ./configure --enable-testing
      
      - name: Build
        run: |
          make -j$(nproc)
      
      - name: Test
        run: |
          make check
```

## Build Optimization

### Speeding Up Builds

**1. Use ccache** (compiler cache):
```bash
./configure CC="ccache gcc" CXX="ccache g++"
```

**2. Parallel build**:
```bash
make -j$(nproc)
```

**3. Disable SWIG** (if not modifying C++):
```bash
./configure --disable-swig
# Uses pre-generated wrappers
```

**4. Incremental builds**:
```bash
# Only rebuild changed files
make
```

**5. Out-of-tree builds**:
```bash
mkdir build
cd build
../configure
make
# Keeps source tree clean
```

### Debug vs. Release Builds

**Debug** (default):
```bash
./configure CXXFLAGS="-g -O0"
# Symbols for gdb, no optimization
```

**Release** (optimized):
```bash
./configure CXXFLAGS="-O3 -DNDEBUG"
# Full optimization, no assertions
```

**With sanitizers** (detect memory errors):
```bash
./configure --enable-sanitizer=address
# AddressSanitizer for memory bugs
```

## Installation

### Standard Installation

```bash
./configure --prefix=/usr/local
make
make install

# Installs:
# /usr/local/bin/pylith
# /usr/local/lib/libpylith.so
# /usr/local/lib/python3.x/site-packages/pylith/
```

### Development Installation

For development, avoid `make install`:

```bash
# Build in-place
./configure
make

# Run from build directory
export PYTHONPATH=$PWD:$PYTHONPATH
export PATH=$PWD/applications:$PATH

pylith --version
```

### Uninstall

```bash
make uninstall
```

## Common Build Issues

### Issue: PETSc not found

```bash
# Solution: Specify PETSc location
./configure --with-petsc-dir=/path/to/petsc --with-petsc-arch=arch-name

# Or set environment variables
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-name
./configure
```

### Issue: Python module not found

```bash
# Solution: Check Python path
./configure PYTHON=/usr/bin/python3
```

### Issue: SWIG version too old

```bash
# Solution: Disable SWIG (use pre-generated wrappers)
./configure --disable-swig
```

### Issue: Parallel build fails

```bash
# Solution: Serial build
make -j1
```

## Summary

PyLith's build system:

✅ **Autotools**: Standard, portable build system  
✅ **SWIG**: Automatic Python binding generation  
✅ **Libtool**: Cross-platform shared library support  
✅ **Modular**: Clean separation of components  
✅ **Testable**: Comprehensive test suite  

Key commands:
```bash
./configure [OPTIONS]  # Configure
make -j4               # Build
make check             # Test
make install           # Install (optional)
```

Next: [05-Module-Structure.md](05-Module-Structure.md) - Python-C++ integration details
