# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	thermoporo_soln.py \
	TestThermoporo.py \
	gravity_soln.py \
	TestGravity.py

dist_noinst_DATA = \
	mesh_tri.msh \
	mesh_quad.msh \
	pylithapp.cfg \
	thermoporo.cfg \
	thermoporo_tri.cfg \
	thermoporo_quad.cfg \
	gravity.cfg \
	gravity_tri.cfg \
	gravity_quad.cfg \
	gravity_temperature_ic.spatialdb \
	gravity_pressure_ic.spatialdb

noinst_TMP =

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	$(RM) $(RM_FLAGS) -r output __pycache__


# End of file
