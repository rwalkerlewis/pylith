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
include $(top_srcdir)/tests/check_cppunit.am

TESTS = test_pylith.py
dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	generate_gmsh.py \
	TestFaultSlip.py \
	faultslip_soln.py

dist_noinst_DATA = \
	pylithapp.cfg \
	faultslip.cfg \
	faultslip_tri.cfg \
	faultslip_quad.cfg \
	mesh_tri.msh \
	mesh_quad.msh

noinst_TMP = \
	output

export_datadir = $(abs_top_builddir)/tests/fullscale/poroelasticity/faults-2d-slip
include $(top_srcdir)/tests/data.am


# End of file
