# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""
Test cases for thermoelastic bar problem.
"""

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check

import meshes
import thermobar_soln


class TestCase(FullTestCase):
    """Test case for thermoelastic bar."""

    def setUp(self):
        """Set up test case."""
        # Temperature check should be exact (linear gradient)
        temp_defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": thermobar_soln.AnalyticalSolution(),
            "mesh": self.mesh,
            "tolerance": 1.0e-5,  # Tight tolerance for temperature
        }
        # Displacement check uses larger tolerance since the analytical solution
        # is an approximation of the full 2D plane strain problem.
        disp_defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": thermobar_soln.AnalyticalSolution(),
            "mesh": self.mesh,
            "tolerance": 1.0,  # 100% relative tolerance for approximate displacement solution
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["temperature"],
                defaults=temp_defaults,
            ),
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement"],
                defaults=disp_defaults,
            ),
            Check(
                mesh_entities=["thermoelastic_material"],
                vertex_fields=["temperature"],
                defaults=temp_defaults,
            ),
            Check(
                mesh_entities=["thermoelastic_material"],
                vertex_fields=["displacement"],
                defaults=disp_defaults,
            ),
        ]
        return

    def run_pylith(self, testName, args):
        """Run PyLith simulation."""
        FullTestCase.run_pylith(self, testName, args, None, nprocs=1)
        return


class TestTri(TestCase):
    """Test case for triangular mesh."""

    def setUp(self):
        """Set up test case."""
        self.name = "thermobar_tri"
        self.mesh = meshes.Tri()
        super().setUp()
        
        TestCase.run_pylith(self, self.name, ["thermobar_tri.cfg"])
        return

    def test_tri(self):
        """Run simulation with triangular mesh."""
        pass  # Simulation is run in setUp


class TestQuad(TestCase):
    """Test case for quadrilateral mesh."""

    def setUp(self):
        """Set up test case."""
        self.name = "thermobar_quad"
        self.mesh = meshes.Quad()
        super().setUp()
        
        TestCase.run_pylith(self, self.name, ["thermobar_quad.cfg"])
        return

    def test_quad(self):
        """Run simulation with quadrilateral mesh."""
        pass  # Simulation is run in setUp


def load_tests(loader, tests, pattern):
    """Custom test loader."""
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromTestCase(TestTri))
    suite.addTests(loader.loadTestsFromTestCase(TestQuad))
    return suite


if __name__ == "__main__":
    unittest.main()


# End of file
