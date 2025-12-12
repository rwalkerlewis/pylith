# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

include $(top_srcdir)/tests/check_cppunit.am

TESTS = test_pylith.py

dist_check_SCRIPTS = test_pylith.py

dist_noinst_PYTHON = \
	meshes.py \
	thermobar_soln.py \
	TestThermoBar.py

dist_noinst_DATA = \
	mesh_tri.mesh \
	mesh_quad.mesh \
	pylithapp.cfg \
	thermobar_tri.cfg \
	thermobar_quad.cfg


# End of file
