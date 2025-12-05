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
        FullTestCase.setUp(self)

        self.exactsoln = thermobar_soln.AnalyticalSolution()

        # Set field checks
        self.checks = [
            Check(
                mesh_entities=["domain"],
                filename="output/thermobar-domain.h5",
                vertex_fields=["displacement", "temperature"],
            ),
            Check(
                mesh_entities=["thermoelastic_material"],
                filename="output/thermobar-thermoelastic_material.h5",
                vertex_fields=["displacement", "temperature"],
            ),
        ]
        return

    def run_pylith(self, testName, args):
        """Run PyLith simulation."""
        FullTestCase.run_pylith(self, testName, args, generatedb=False, nprocs=1)
        return


class TestTri(TestCase):
    """Test case for triangular mesh."""

    def setUp(self):
        """Set up test case."""
        self.mesh = meshes.Tri()
        TestCase.setUp(self)
        return

    def test_tri(self):
        """Run simulation with triangular mesh."""
        self.run_pylith("tri", ["thermobar_tri.cfg"])
        return


class TestQuad(TestCase):
    """Test case for quadrilateral mesh."""

    def setUp(self):
        """Set up test case."""
        self.mesh = meshes.Quad()
        TestCase.setUp(self)
        return

    def test_quad(self):
        """Run simulation with quadrilateral mesh."""
        self.run_pylith("quad", ["thermobar_quad.cfg"])
        return


def load_tests(loader, tests, pattern):
    """Custom test loader."""
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromTestCase(TestTri))
    suite.addTests(loader.loadTestsFromTestCase(TestQuad))
    return suite


if __name__ == "__main__":
    unittest.main()


# End of file
