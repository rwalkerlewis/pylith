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


class Tet(object):
    """Mesh information for tet mesh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=374, ncorners=4, nvertices=112),
        "upper_crust": MeshEntity(ncells=140, ncorners=4, nvertices=57),
        "lower_crust": MeshEntity(ncells=234, ncorners=4, nvertices=80),
    }


class Hex(object):
    """Mesh information for hex mesh."""

    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=8, nvertices=100),
        "upper_crust": MeshEntity(ncells=16, ncorners=8, nvertices=50),
        "lower_crust": MeshEntity(ncells=32, ncorners=8, nvertices=75),
    }


# End of file
