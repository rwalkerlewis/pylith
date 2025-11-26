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

from pylith.testing.FullTestApp import TestDriver
import TestPointSource


class PointSourceTestDriver(TestDriver):
    """Driver for point source tests."""

    def __init__(self):
        TestDriver.__init__(self)
        return

    def _createSuite(self):
        """Create test suite."""
        suite = unittest.TestSuite()
        for test in TestPointSource.test_cases():
            suite.addTest(unittest.makeSuite(test))
        return suite


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    TestPointSource.FullTestCase.parse_args()
    driver = PointSourceTestDriver()
    driver.main()


# End of file
