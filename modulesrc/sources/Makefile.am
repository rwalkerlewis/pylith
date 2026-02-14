# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

subpackage = sources
include $(top_srcdir)/subpackage.am
include $(top_srcdir)/modulesrc/module.am

subpkgpython_LTLIBRARIES = _sources.la
subpkgpython_PYTHON = sources.py

swig_sources = \
	sources.i \
	PointForce.i

swig_generated = \
	sources_wrap.cxx \
	sources.py

_sources_la_LDFLAGS = -module -avoid-version \
	$(AM_LDFLAGS) $(PYTHON_LA_LDFLAGS)

dist__sources_la_SOURCES = $(swig_sources) $(swig_generated)

_sources_la_LIBADD = \
	$(top_builddir)/libsrc/pylith/libpylith.la \
	$(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) $(PYTHON_SYSLIBS)

if ENABLE_SWIG
$(srcdir)/sources_wrap.cxx $(srcdir)/sources.py: $(swig_sources)
	$(SWIG) -Wall -c++ -python $<
else
$(srcdir)/sources_wrap.cxx $(srcdir)/sources.py:
	$(error Missing SWIG generated files. Make sure SWIG is installed and reconfigure with --enable-swig)
endif


MAINTAINERCLEANFILES = \
	$(srcdir)/sources_wrap.cxx \
	$(srcdir)/sources.py


# End of file
