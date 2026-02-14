# Gmsh Model Building Scripts

This directory contains scripts and examples for building Gmsh models compatible with PyLith, ranging from basic to very complex geometries.

## Structure

```
gmsh-models/
├── examples/
│   ├── basic/              # Most basic examples
│   │   └── box_2d_no_faults.py
│   ├── intermediate/       # Single fault examples
│   │   └── single_fault_2d.py
│   ├── complex/            # Multiple faults examples
│   │   └── multiple_faults_2d.py
│   └── very-complex/       # 3D with multiple faults
│       └── multiple_faults_3d.py
└── tests/                  # Test suite
    ├── test_mesh_generation.py
    └── test_pylith_integration.py
```

## Examples

### Basic: Box 2D (No Faults)

**File:** `examples/basic/box_2d_no_faults.py`

The most basic example - a simple rectangular domain with boundary conditions but no faults. This demonstrates:
- Domain creation
- Material marking
- Boundary condition marking
- Uniform mesh generation

**Usage:**
```bash
cd examples/basic
python box_2d_no_faults.py --write
```

### Intermediate: Single Fault 2D

**File:** `examples/intermediate/single_fault_2d.py`

A 2D domain with a single vertical strike-slip fault that splits the domain. This demonstrates:
- Domain split by a fault
- Multiple materials (one on each side)
- Fault marking
- Mesh refinement near fault

**Usage:**
```bash
cd examples/intermediate
python single_fault_2d.py --write
```

### Complex: Multiple Faults 2D

**File:** `examples/complex/multiple_faults_2d.py`

A 2D domain with multiple intersecting faults (main fault + two splay faults). This demonstrates:
- Multiple faults creating multiple material regions
- Fault intersection points
- Complex geometry handling
- Refinement near multiple faults

**Usage:**
```bash
cd examples/complex
python multiple_faults_2d.py --write
```

### Very Complex: Multiple Faults 3D

**File:** `examples/very-complex/multiple_faults_3d.py`

A 3D domain with multiple intersecting faults. This demonstrates:
- 3D domain with multiple faults
- Multiple materials in 3D
- 3D fault marking
- Complex 3D mesh refinement

**Usage:**
```bash
cd examples/very-complex
python multiple_faults_3d.py --write
```

## Running Examples

All examples use the `GenerateMesh` base class from `pylith.meshio.gmsh_utils`, which provides command-line options:

```bash
# Generate mesh (geometry + mark + generate + write)
python script.py --write

# Generate mesh and open in Gmsh GUI
python script.py --write --gui

# Create only geometry and open in Gmsh GUI
python script.py --geometry --gui

# See all options
python script.py --help
```

## Testing

Run the test suite to validate mesh generation and PyLith compatibility:

```bash
# Run all tests
cd tests
python -m pytest test_mesh_generation.py test_pylith_integration.py -v

# Run specific test class
python test_mesh_generation.py TestBasicMeshGeneration

# Run with unittest
python -m unittest discover -v
```

### Test Coverage

The test suite validates:

1. **Mesh Generation Tests** (`test_mesh_generation.py`):
   - Basic mesh generation (no faults)
   - Single fault mesh generation
   - Multiple faults mesh generation
   - PyLith compatibility requirements
   - Mesh quality metrics

2. **PyLith Integration Tests** (`test_pylith_integration.py`):
   - Mesh file format compatibility
   - Physical group naming conventions
   - PyLith mesh import capabilities (if PyLith is available)

## PyLith Compatibility

All examples follow PyLith requirements:

1. **Material Tags**: Must be positive integers, named as `material-id:N`
2. **Boundary Groups**: Named descriptively (e.g., `boundary_xneg`, `boundary_ypos`)
3. **Fault Groups**: Named descriptively (e.g., `fault`, `fault_main`, `fault_secondary`)
4. **File Format**: Mesh files must have `.msh` extension
5. **Physical Groups**: Correct dimensions (1 for curves in 2D, 2 for surfaces in 3D)

## Requirements

- Python 3.x
- Gmsh Python API
- PyLith (for integration tests)

## Notes

- All examples use the built-in Gmsh geometry engine except the 3D example which uses OpenCASCADE
- Mesh refinement uses geometric progression from faults
- Examples are designed to be educational and can be modified for specific use cases
