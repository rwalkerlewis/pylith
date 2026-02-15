# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================


include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_thermoporoelasticity_faults2d

TESTS = \
	run_TwoBlocksStaticTHM.sh

check_SCRIPTS = $(TESTS)

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
