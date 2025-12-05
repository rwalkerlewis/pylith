# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
"""
Mesh information for thermoelasticity fullscale tests.
"""

from pylith.testing.FullTestApp import MeshEntity


class Tri:
    """Mesh information for triangular mesh."""
    ENTITIES = {
        "domain": MeshEntity(ncells=64, ncorners=3, nvertices=45),
        "thermoelastic_material": MeshEntity(ncells=64, ncorners=3, nvertices=45),
        "bc_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }


class Quad:
    """Mesh information for quadrilateral mesh."""
    ENTITIES = {
        "domain": MeshEntity(ncells=32, ncorners=4, nvertices=45),
        "thermoelastic_material": MeshEntity(ncells=32, ncorners=4, nvertices=45),
        "bc_xneg": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_xpos": MeshEntity(ncells=4, ncorners=2, nvertices=5),
        "bc_yneg": MeshEntity(ncells=8, ncorners=2, nvertices=9),
        "bc_ypos": MeshEntity(ncells=8, ncorners=2, nvertices=9),
    }


# End of file
