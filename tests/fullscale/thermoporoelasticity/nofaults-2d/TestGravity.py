#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Thermoporoelasticity fullscale test with gravity and temperature gradient.

NOTE: This test is currently disabled because the analytical solution expects a
steady-state equilibrium with both hydrostatic pressure and a geothermal temperature
gradient. However, achieving this steady-state requires either:
1. Fixing the temperature field (not simulating heat diffusion), or
2. Adding a heat source/flux boundary condition at depth.

Without these, the temperature diffuses to uniform T=300K at steady-state, and the
coupled thermoporoelastic solution differs significantly from the expected analytical
result. This test needs further development to properly model the intended physics.
"""

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check
from pylith.testing import TestCases

import meshes
import gravity_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Thermoporoelasticity fullscale test with gravity and temperature gradient."""

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": gravity_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["thermoporoelastic"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
@unittest.skip("Analytical solution requires fixed temperature gradient; needs further development")
class TestTriGmsh(TestCase):
    """Test case for triangular mesh."""

    def setUp(self):
        self.name = "gravity_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_tri.cfg"])


# -------------------------------------------------------------------------------------------------
@unittest.skip("Analytical solution requires fixed temperature gradient; needs further development")
class TestQuadGmsh(TestCase):
    """Test case for quadrilateral mesh."""

    def setUp(self):
        self.name = "gravity_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_quad.cfg"])


# -------------------------------------------------------------------------------------------------
def load_tests(loader, tests, pattern):
    TEST_CLASSES = (
        TestTriGmsh,
        TestQuadGmsh,
    )
    return TestCases.make_suite(test_classes=TEST_CLASSES, loader=loader)


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()
    unittest.main(verbosity=2, argv=["test"])


# End of file
