# -*- Makefile -*-
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

include $(top_srcdir)/check.am

TESTS = testpylith.py

dist_check_SCRIPTS = testpylith.py

dist_noinst_PYTHON = \
	TestQuad4.py \
	TestRigidRotate.py \
	rigidrotate_soln.py \
	rigidrotate_gendb.py \
	TestCompressRotate.py \
	compressrotate_soln.py \
	compressrotate_gendb.py \
	TestShearRotate.py \
	shearrotate_soln.py \
	shearrotate_gendb.py

dist_noinst_DATA = \
	mesh.jou \
	mesh.exo \
	matprops.spatialdb \
	rigidrotate.cfg \
	compressrotate.cfg \
	shearrotate.cfg

noinst_TMP = \
	rigidrotate_disp.spatialdb \
	compressrotate_disp.spatialdb \
	shearrotate_disp.spatialdb


# 'export' the input files by performing a mock install
export_datadir = $(top_builddir)/tests_auto/2d/finitestrain
export-data: $(dist_noinst_PYTHON) $(dist_noinst_DATA)
	for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA); do $(install_sh_DATA) $(srcdir)/$$f $(export_datadir); done

clean-data:
	if [ "X$(top_srcdir)" != "X$(top_builddir)" ]; then for f in $(dist_noinst_PYTHON) $(dist_noinst_DATA) $(noinst_TMP); do $(RM) $(RM_FLAGS) $(export_datadir)/$$f; done; fi


BUILT_SOURCES = export-data

clean-local: clean-local-tmp clean-data
.PHONY: clean-local-tmp
clean-local-tmp:
	-rm *.h5 *.xmf *.pyc


# End of file 
