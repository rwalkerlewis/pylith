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
# @file tests/pytests/faults/TestKinSrc.py
#
# @brief Unit testing of Python KinSrc object.

import unittest

from pylith.testing.UnitTestApp import (TestAbstractComponent, TestComponent)
import pylith.faults.KinSrcPoro
import pylith.faults.KinSrcPoroStep


class TestKinSrcPoro(TestAbstractComponent):
    """Unit testing of KinSrcPoro object.
    """
    _class = pylith.faults.KinSrcPoro.KinSrcPoro


class TestKinSrcPoroStep(TestComponent):
    """Unit testing of KinSrcPoroStep object.
    """
    _class = pylith.faults.KinSrcPoroStep.KinSrcPoroStep
    _factory = pylith.faults.KinSrcPoroStep.eq_kinematic_src


if __name__ == "__main__":
    suite = unittest.TestSuite()
    for cls in [
        TestKinSrcPoro,
        TestKinSrcPoroStep,
    ]:
        suite.addTest(unittest.makeSuite(cls))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
