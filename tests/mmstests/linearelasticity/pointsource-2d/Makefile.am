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

MMS_DRIVER = mmstest_linearelasticity_pointsource2d

TESTS = run_$(MMS_DRIVER).sh

check_PROGRAMS = $(MMS_DRIVER)

mmstest_linearelasticity_pointsource2d_SOURCES = \
	TestPointForce.cc \
	TestPointForce.hh \
	PointForceExplosion2D.cc \
	PointForceExplosion2D.hh \
	TestCases.cc

run_$(MMS_DRIVER).sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) -v yes [TriP1]" >> $@
	chmod +x $@

dist_noinst_DATA = \
	data/tri.mesh \
	data/quad.mesh


CLEANFILES = run_$(MMS_DRIVER).sh


# End of file
