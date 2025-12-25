# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

# NOTE: These tests are currently disabled because the test solution (rigid block motion
# with uniform pressure, temperature, and trace strain) is not a valid analytical solution
# for the coupled thermoporoelasticity equations. The test case needs to be redesigned
# with a proper manufactured solution that satisfies all 4 governing equations:
# 1. Momentum balance
# 2. Mass balance (fluid flow)
# 3. Trace strain constraint (div(u) = epsilon_v)
# 4. Energy balance (heat equation)

include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_thermoporoelasticity_faults2d

# Tests are disabled until the test case is redesigned
TESTS = 

check_SCRIPTS =

check_PROGRAMS =  $(MMS_DRIVER)


mmstest_thermoporoelasticity_faults2d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestFaultKinThermoporoelasticity.cc \
	TestCases.cc \
	TwoBlocksStaticTHM.cc

run_%.sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
	chmod +x $@

dist_noinst_HEADERS = \
	TestFaultKinThermoporoelasticity.hh \
	TwoBlocksStaticTHM.hh


dist_noinst_DATA = \
	data/tri.mesh \
	data/quad.mesh


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
