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
	TestPointSourceRicker.py

dist_noinst_DATA = \
	pylithapp.cfg \
	pointsource_ricker.cfg

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am

clean-local: clean-data
	$(RM) $(RM_FLAGS) -r output __pycache__

# End of file

