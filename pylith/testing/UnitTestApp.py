# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
#
# @file pylith/testing/UnitTestApp.py
#
# @brief Backward-compatible imports for unit test helpers.

from .TestCases import (  # noqa: F401
    TestAbstractComponent,
    TestComponent,
    configureComponent,
    make_suite,
)

# End of file

