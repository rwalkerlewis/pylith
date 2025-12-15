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

import unittest

from pylith.testing.FullTestApp import TestDriver, FullTestCase

import TestThermoporo


TEST_MODULES = (TestThermoporo,)


class TestApp(TestDriver):
    """Driver application for thermoporoelasticity full-scale tests (3D, no faults)."""

    def _suite(self):
        loader = unittest.defaultTestLoader
        suite = unittest.TestSuite()
        for mod in TEST_MODULES:
            suite.addTests(loader.loadTestsFromModule(mod))
        return suite


if __name__ == "__main__":
    FullTestCase.parse_args()
    TestApp().main()


# End of file
