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
# @file tests/pytests/sources/TestWellboreSource.py
#
# @brief Unit testing of Python TimeHistoryWavelet source time function.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.sources.TimeHistoryWavelet import (TimeHistoryWavelet, momenttensorforce_sourcetimefunction)


class TestTimeHistoryWavelet(TestComponent):
    """Unit testing of TimeHistoryWavelet source time function."""
    _class = TimeHistoryWavelet
    _factory = momenttensorforce_sourcetimefunction


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTimeHistoryWavelet))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
