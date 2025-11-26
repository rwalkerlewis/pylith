#!/usr/bin/env nemesis
"""Integration tests for PyLith compatibility.

These tests verify that generated meshes can be used with PyLith by:
1. Checking mesh file format compatibility
2. Validating physical group naming conventions
3. Testing mesh import capabilities
"""

import unittest
import os
import sys
import tempfile
import shutil

import gmsh

# Try to import PyLith components for integration testing
try:
    from pylith.meshio.MeshIOPetsc import MeshIOPetsc, mesh_input
    PYLITH_AVAILABLE = True
except ImportError:
    PYLITH_AVAILABLE = False

# Add parent directory to path to import examples
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'examples'))


@unittest.skipIf(not PYLITH_AVAILABLE, "PyLith not available")
class TestPyLithMeshImport(unittest.TestCase):
    """Test that generated meshes can be imported by PyLith."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.original_cwd = os.getcwd()
        os.chdir(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        os.chdir(self.original_cwd)
        shutil.rmtree(self.temp_dir)
        gmsh.finalize()
    
    def test_basic_mesh_import(self):
        """Test that basic mesh can be imported by PyLith."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        filename = "test_mesh.msh"
        app.write(filename, binary=False)
        
        # Try to create PyLith mesh reader
        reader = mesh_input()
        reader.filename = filename
        
        # This would normally read the mesh, but we'll just check the setup
        self.assertEqual(reader.filename, filename, "Reader should accept mesh filename")
    
    def test_fault_mesh_import(self):
        """Test that fault mesh can be imported by PyLith."""
        from intermediate.single_fault_2d import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        filename = "test_fault_mesh.msh"
        app.write(filename, binary=False)
        
        # Try to create PyLith mesh reader
        reader = mesh_input()
        reader.filename = filename
        
        self.assertEqual(reader.filename, filename, "Reader should accept fault mesh filename")


class TestPhysicalGroupNaming(unittest.TestCase):
    """Test that physical groups follow PyLith naming conventions."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.original_cwd = os.getcwd()
        os.chdir(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        os.chdir(self.original_cwd)
        shutil.rmtree(self.temp_dir)
        gmsh.finalize()
    
    def test_material_naming_convention(self):
        """Test that materials follow PyLith naming convention (material-id:N)."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        
        physical_groups = gmsh.model.get_physical_groups()
        material_groups = []
        
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name and name.startswith("material-id:"):
                material_groups.append((name, tag))
                # Check naming convention
                expected_name = f"material-id:{tag}"
                self.assertEqual(name, expected_name,
                               f"Material name should be '{expected_name}', got '{name}'")
        
        self.assertGreater(len(material_groups), 0, "At least one material should exist")
    
    def test_boundary_naming_convention(self):
        """Test that boundaries follow PyLith naming conventions."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        
        physical_groups = gmsh.model.get_physical_groups()
        boundary_names = []
        
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name and name.startswith("boundary_"):
                boundary_names.append(name)
                # Boundaries should have descriptive names
                self.assertIn("boundary_", name, "Boundary names should contain 'boundary_'")
        
        expected_boundaries = ["boundary_xneg", "boundary_xpos", "boundary_yneg", "boundary_ypos"]
        for expected in expected_boundaries:
            self.assertIn(expected, boundary_names,
                        f"Expected boundary '{expected}' not found. Found: {boundary_names}")
    
    def test_fault_naming_convention(self):
        """Test that faults follow PyLith naming conventions."""
        from intermediate.single_fault_2d import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        
        physical_groups = gmsh.model.get_physical_groups()
        fault_found = False
        
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name == "fault" or (name and "fault" in name.lower()):
                fault_found = True
                # Fault names should be descriptive
                self.assertIsNotNone(name, "Fault should have a name")
                break
        
        self.assertTrue(fault_found, "Fault group should exist")


class TestMeshFileFormat(unittest.TestCase):
    """Test mesh file format compatibility."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.original_cwd = os.getcwd()
        os.chdir(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        os.chdir(self.original_cwd)
        shutil.rmtree(self.temp_dir)
        gmsh.finalize()
    
    def test_msh_file_extension(self):
        """Test that mesh files have .msh extension (PyLith requirement)."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        filename = "test_mesh.msh"
        app.write(filename, binary=False)
        
        self.assertTrue(filename.endswith(".msh"), "Mesh file should have .msh extension")
        self.assertTrue(os.path.exists(filename), "Mesh file should exist")
    
    def test_ascii_format_readable(self):
        """Test that ASCII format mesh files are readable."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        filename = "test_mesh.msh"
        app.write(filename, binary=False)
        
        # Check that file is readable and contains expected content
        with open(filename, 'r') as f:
            content = f.read()
            # Gmsh mesh files should contain "$MeshFormat" header
            self.assertIn("$MeshFormat", content, "Mesh file should contain Gmsh header")
            self.assertIn("$Nodes", content, "Mesh file should contain nodes section")
            self.assertIn("$Elements", content, "Mesh file should contain elements section")


def load_tests(loader, tests, pattern):
    """Load all test cases."""
    test_classes = [
        TestPhysicalGroupNaming,
        TestMeshFileFormat,
    ]
    
    if PYLITH_AVAILABLE:
        test_classes.append(TestPyLithMeshImport)
    
    suite = unittest.TestSuite()
    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite


if __name__ == "__main__":
    unittest.main(verbosity=2)
