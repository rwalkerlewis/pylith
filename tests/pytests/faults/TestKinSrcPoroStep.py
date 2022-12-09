#!/usr/bin/env nemesis
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ======================================================================
#
# @file tests/pytests/faults/TestKinSrcPoroStep.py
#
# @brief Unit testing of Python KinSrcPoroStep object.

import unittest

from pylith.faults.KinSrcPoroStep import KinSrcPoroStep
from pylith.tests.UnitTestApp import configureComponent


class TestKinSrcPoroStep(unittest.TestCase):
    """Unit testing of KinSrcPoroStep object.
    """

    def test_constructor(self):
        src = KinSrcPoroStep()

    def test_configure(self):
        src = KinSrcPoroStep()
        configureComponent(src)

    def test_factory(self):
        from pylith.faults.KinSrcPoroStep import eq_kinematic_src
        src = eq_kinematic_src()
        self.assertTrue(isinstance(src, KinSrcPoroStep))


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestKinSrcPoroStep))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file