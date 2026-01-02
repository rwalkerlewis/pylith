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
"""Test cases for poroelasticity with prescribed nonzero fault slip.

This test verifies that PyLith correctly handles:
1. Coupled poroelastic physics with fault slip
2. Prescribed nonzero fault slip (0.5 m left-lateral)
3. Pore pressure evolution with drained boundaries
"""

import unittest

from pylith.testing.FullTestApp import FullTestCase, Check

import meshes
import faultslip_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):
    """Base class for fault slip tests."""

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": faultslip_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain"],
                vertex_fields=["displacement", "trace_strain", "pressure"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["poroelastic"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields=[
                    "solid_density",
                    "fluid_density",
                    "fluid_viscosity",
                    "porosity",
                    "shear_modulus",
                    "drained_bulk_modulus",
                    "biot_coefficient",
                    "biot_modulus",
                    "isotropic_permeability",
                ],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["poroelastic"],
                vertex_fields=[
                    "displacement",
                    "pressure",
                    "trace_strain",
                    "cauchy_strain",
                    "cauchy_stress",
                ],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=[
                    "bc_disp_xneg",
                    "bc_disp_xpos",
                    "bc_disp_yneg",
                    "bc_disp_ypos",
                ],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude", "normal_dir", "tangential_dir"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=[
                    "bc_disp_xneg",
                    "bc_disp_xpos",
                    "bc_disp_yneg",
                    "bc_disp_ypos",
                    "bc_press_yneg",
                    "bc_press_ypos",
                ],
                vertex_fields=["displacement", "trace_strain", "pressure"],
                final_time_only=True,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["normal_dir", "strike_dir"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                vertex_fields=["slip", "traction_change"],
                final_time_only=True,
                defaults=defaults,
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestQuadGmsh(TestCase):
    """Test with quadrilateral mesh from Gmsh."""

    def setUp(self):
        self.name = "faultslip_quad"
        self.mesh = meshes.QuadGmsh()
        super().setUp()

        TestCase.run_pylith(
            self, self.name, ["faultslip.cfg", "faultslip_quad.cfg"]
        )
        return


# -------------------------------------------------------------------------------------------------
class TestTriGmsh(TestCase):
    """Test with triangular mesh from Gmsh."""

    def setUp(self):
        self.name = "faultslip_tri"
        self.mesh = meshes.TriGmsh()
        super().setUp()

        TestCase.run_pylith(
            self, self.name, ["faultslip.cfg", "faultslip_tri.cfg"]
        )
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestQuadGmsh,
        TestTriGmsh,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
