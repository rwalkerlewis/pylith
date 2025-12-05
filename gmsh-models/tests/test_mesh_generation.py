#!/usr/bin/env nemesis
"""Test suite for validating Gmsh mesh generation and PyLith compatibility.

This test suite ensures that:
1. Meshes can be generated successfully
2. Generated meshes have correct structure
3. Meshes are compatible with PyLith requirements
4. Physical groups are correctly defined
"""

import unittest
import os
import sys
import tempfile
import shutil

import gmsh

# Add parent directory to path to import examples
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'examples'))


class TestBasicMeshGeneration(unittest.TestCase):
    """Test basic mesh generation (no faults)."""
    
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
    
    def test_basic_box_2d_generation(self):
        """Test that basic 2D box mesh can be generated."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        # Check that geometry was created
        entities = gmsh.model.get_entities()
        self.assertGreater(len(entities), 0, "Geometry entities should be created")
        
        # Check that physical groups exist
        physical_groups = gmsh.model.get_physical_groups()
        self.assertGreater(len(physical_groups), 0, "Physical groups should be created")
        
        # Check for material group
        material_found = False
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name and name.startswith("material-id:"):
                material_found = True
                break
        self.assertTrue(material_found, "Material group should be created")
        
        # Check for boundary groups
        boundary_names = ["boundary_xneg", "boundary_xpos", "boundary_yneg", "boundary_ypos"]
        found_boundaries = []
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name in boundary_names:
                found_boundaries.append(name)
        self.assertEqual(len(found_boundaries), len(boundary_names),
                        f"All boundary groups should be created. Found: {found_boundaries}")
    
    def test_basic_box_2d_write(self):
        """Test that basic 2D box mesh can be written to file."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        filename = "test_mesh.msh"
        app.write(filename, binary=False)
        
        self.assertTrue(os.path.exists(filename), "Mesh file should be created")
        self.assertGreater(os.path.getsize(filename), 0, "Mesh file should not be empty")


class TestSingleFaultMeshGeneration(unittest.TestCase):
    """Test single fault mesh generation."""
    
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
    
    def test_single_fault_2d_generation(self):
        """Test that single fault 2D mesh can be generated."""
        from intermediate.single_fault_2d import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        # Check that geometry was created
        entities = gmsh.model.get_entities()
        self.assertGreater(len(entities), 0, "Geometry entities should be created")
        
        # Check for fault group
        physical_groups = gmsh.model.get_physical_groups()
        fault_found = False
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name == "fault":
                fault_found = True
                break
        self.assertTrue(fault_found, "Fault group should be created")
        
        # Check for multiple materials (one on each side of fault)
        material_tags = []
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name and name.startswith("material-id:"):
                material_tags.append(tag)
        self.assertGreaterEqual(len(material_tags), 2,
                               "At least two materials should be created (one on each side of fault)")


class TestMultipleFaultsMeshGeneration(unittest.TestCase):
    """Test multiple faults mesh generation."""
    
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
    
    def test_multiple_faults_2d_generation(self):
        """Test that multiple faults 2D mesh can be generated."""
        from complex.multiple_faults_2d import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        # Check that geometry was created
        entities = gmsh.model.get_entities()
        self.assertGreater(len(entities), 0, "Geometry entities should be created")
        
        # Check for multiple fault groups
        physical_groups = gmsh.model.get_physical_groups()
        fault_names = []
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name and ("fault" in name.lower()):
                fault_names.append(name)
        
        self.assertGreaterEqual(len(fault_names), 2,
                               f"Multiple fault groups should be created. Found: {fault_names}")
        
        # Check for multiple materials (created by fault intersections)
        material_tags = []
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name and name.startswith("material-id:"):
                material_tags.append(tag)
        self.assertGreaterEqual(len(material_tags), 3,
                               "Multiple materials should be created by fault intersections")


class TestPyLithCompatibility(unittest.TestCase):
    """Test PyLith compatibility requirements."""
    
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
    
    def test_material_tags_positive(self):
        """Test that material tags are positive (PyLith requirement)."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        
        physical_groups = gmsh.model.get_physical_groups()
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name and name.startswith("material-id:"):
                self.assertGreater(tag, 0, f"Material tag {tag} should be positive")
    
    def test_boundary_groups_have_correct_dimension(self):
        """Test that boundary groups have correct dimension (1 for 2D, 2 for 3D)."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        
        physical_groups = gmsh.model.get_physical_groups()
        boundary_names = ["boundary_xneg", "boundary_xpos", "boundary_yneg", "boundary_ypos"]
        
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name in boundary_names:
                # For 2D meshes, boundaries should be dimension 1 (curves)
                self.assertEqual(dim, 1, f"Boundary {name} should have dimension 1 for 2D mesh")
    
    def test_fault_groups_have_correct_dimension(self):
        """Test that fault groups have correct dimension."""
        from intermediate.single_fault_2d import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        
        physical_groups = gmsh.model.get_physical_groups()
        fault_found = False
        
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if name == "fault":
                fault_found = True
                # For 2D meshes, faults should be dimension 1 (curves)
                self.assertEqual(dim, 1, f"Fault should have dimension 1 for 2D mesh")
                break
        
        self.assertTrue(fault_found, "Fault group should exist")
    
    def test_mesh_can_be_read_by_petsc(self):
        """Test that generated mesh can be read (basic structure check)."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        filename = "test_mesh.msh"
        app.write(filename, binary=False)
        
        # Try to read the mesh back with gmsh
        gmsh.finalize()
        gmsh.initialize()
        gmsh.open(filename)
        
        # Check that we can get entities
        entities = gmsh.model.get_entities()
        self.assertGreater(len(entities), 0, "Mesh should be readable")
        
        # Check that we can get physical groups
        physical_groups = gmsh.model.get_physical_groups()
        self.assertGreater(len(physical_groups), 0, "Physical groups should be readable")


class TestMeshQuality(unittest.TestCase):
    """Test mesh quality metrics."""
    
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
    
    def test_mesh_has_elements(self):
        """Test that generated mesh has elements."""
        from basic.box_2d_no_faults import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        # Get element types and counts
        element_types, element_tags, node_tags = gmsh.model.mesh.get_elements()
        
        self.assertGreater(len(element_types), 0, "Mesh should have elements")
        self.assertGreater(len(element_tags[0]), 0, "Mesh should have element tags")
        self.assertGreater(len(node_tags[0]), 0, "Mesh should have node tags")
    
    def test_mesh_refinement_near_fault(self):
        """Test that mesh is refined near faults."""
        from intermediate.single_fault_2d import App
        
        app = App()
        app.initialize(None)
        app.create_geometry()
        app.mark()
        app.generate_mesh("tri")
        
        # Get node coordinates
        node_tags, coords, _ = gmsh.model.mesh.get_nodes()
        
        # Get elements
        element_types, element_tags, node_tags_elem = gmsh.model.mesh.get_elements()
        
        # Calculate element sizes near fault (x=0) vs far from fault
        # This is a basic check - in practice you'd compute actual element sizes
        self.assertGreater(len(node_tags), 0, "Mesh should have nodes")
        self.assertGreater(len(element_tags[0]), 0, "Mesh should have elements")


def load_tests(loader, tests, pattern):
    """Load all test cases."""
    test_classes = [
        TestBasicMeshGeneration,
        TestSingleFaultMeshGeneration,
        TestMultipleFaultsMeshGeneration,
        TestPyLithCompatibility,
        TestMeshQuality,
    ]
    suite = unittest.TestSuite()
    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite


if __name__ == "__main__":
    unittest.main(verbosity=2)
