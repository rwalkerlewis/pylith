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
	blackoil_soln.py \
	TestBlackOil.py \
	generate_gmsh.py

dist_noinst_DATA = \
	pylithapp.cfg \
	blackoil.cfg \
	blackoil_quad.cfg \
	blackoil_tri.cfg \
	blackoil_properties.spatialdb \
	mesh_quad.msh \
	mesh_tri.msh

noinst_TMP =

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# 'export' the input files by performing a mock install.
export-data: $(dist_noinst_DATA)
	(mkdir -p $(export_datadir) && cp -f $(dist_noinst_DATA) $(export_datadir) && cp -f $(dist_noinst_PYTHON) $(export_datadir))


# End of file
