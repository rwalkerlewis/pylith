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


class TetGmsh(object):
    """Mesh information for tet mesh using Gmsh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=571, ncorners=4, nvertices=182 + 38),
        "thermoporoelastic": MeshEntity(ncells=571, ncorners=4, nvertices=182 + 38),
        "fault": MeshEntity(ncells=56, ncorners=3, nvertices=38),
    }


class HexGmsh(object):
    """Mesh information for hex mesh using Gmsh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=150, ncorners=8, nvertices=252 + 36),
        "thermoporoelastic": MeshEntity(ncells=150, ncorners=8, nvertices=252 + 36),
        "fault": MeshEntity(ncells=25, ncorners=4, nvertices=36),
    }


# End of file
