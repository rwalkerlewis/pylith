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

from pylith.testing.FullTestApp import TestDriver, FullTestCase

import unittest


class TestApp(TestDriver):
    """Driver application for full-scale source tests."""

    def _suite(self):
        suite = unittest.TestSuite()

        import TestPointSourceRicker
        for test in TestPointSourceRicker.test_cases():
            suite.addTest(unittest.makeSuite(test))

        return suite


if __name__ == "__main__":
    FullTestCase.parse_args()
    TestApp().main()


# End of file

