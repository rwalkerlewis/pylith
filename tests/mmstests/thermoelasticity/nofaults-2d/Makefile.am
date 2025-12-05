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

MMS_DRIVER = mmstest_thermoelasticity_nofaults2d

TESTS = $(MMS_DRIVER)

check_PROGRAMS = $(MMS_DRIVER)

mmstest_thermoelasticity_nofaults2d_SOURCES = \
	TestThermoelasticity.cc \
	UniformThermoelasticity2D.cc \
	TestCases.cc \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc

dist_noinst_HEADERS = \
	TestThermoelasticity.hh \
	UniformThermoelasticity2D.hh

dist_noinst_DATA = \
	data/tri.mesh \
	data/quad.mesh

noinst_TMP =


export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
