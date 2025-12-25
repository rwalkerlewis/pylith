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
"""Thermoporoelasticity fullscale test with gravity and temperature gradient (3D).

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
    """Thermoporoelasticity fullscale test with gravity and temperature gradient (3D)."""

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
                mesh_entities=["upper_crust"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["lower_crust"],
                vertex_fields=["displacement", "pressure", "temperature"],
                final_time_only=True,
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
@unittest.skip("Analytical solution requires fixed temperature gradient; needs further development")
class TestTet(TestCase):
    """Test case for tetrahedral mesh."""

    def setUp(self):
        self.name = "gravity_tet"
        self.mesh = meshes.Tet()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_tet.cfg"])


# -------------------------------------------------------------------------------------------------
@unittest.skip("Analytical solution requires fixed temperature gradient; needs further development")
class TestHex(TestCase):
    """Test case for hexahedral mesh."""

    def setUp(self):
        self.name = "gravity_hex"
        self.mesh = meshes.Hex()
        super().setUp()
        TestCase.run_pylith(self, self.name, ["gravity.cfg", "gravity_hex.cfg"])


# -------------------------------------------------------------------------------------------------
def load_tests(loader, tests, pattern):
    TEST_CLASSES = (
        TestTet,
        TestHex,
    )
    return TestCases.make_suite(test_classes=TEST_CLASSES, loader=loader)


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()
    unittest.main(verbosity=2, argv=["test"])


# End of file
