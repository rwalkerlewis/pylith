# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pylith.testing.FullTestApp import MeshEntity


class TriGmsh(object):
    """Mesh information for tri mesh using Gmsh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=42, ncorners=3, nvertices=30),
        # Materials
        "thermoporoelastic": MeshEntity(ncells=42, ncorners=3, nvertices=30),
    }


class QuadGmsh(object):
    """Mesh information for quad mesh using Gmsh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=16, ncorners=4, nvertices=25),
        # Materials
        "thermoporoelastic": MeshEntity(ncells=16, ncorners=3, nvertices=25),
    }


# End of file
