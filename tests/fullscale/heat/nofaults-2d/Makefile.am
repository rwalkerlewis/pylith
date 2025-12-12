# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

include $(top_srcdir)/tests/check_pylith.am

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	TestHeatBar.py \
	heatbar_soln.py

dist_noinst_DATA = \
	pylithapp.cfg \
	heatbar_tri.cfg \
	heatbar_quad.cfg \
	mesh_tri.mesh \
	mesh_quad.mesh

noinst_TMP =


export_datadir = $(abs_top_builddir)/tests/fullscale/heat/nofaults-2d
include $(top_srcdir)/tests/data.am


clean-local: clean-local-tmp clean-data
CLEANFILES = $(noinst_TMP)


# End of file
