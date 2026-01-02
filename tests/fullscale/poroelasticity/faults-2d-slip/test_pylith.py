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
"""Driver for running all tests in this directory."""

import unittest

from pylith.testing.FullTestApp import FullTestCase

import TestFaultSlip


# -------------------------------------------------------------------------------------------------
def test_cases():
    """Collect test cases from all test modules."""
    cases = []
    cases += TestFaultSlip.test_cases()
    return cases


# -------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
