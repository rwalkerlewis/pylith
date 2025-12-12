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


class Tri(object):
    """Mesh information for tri mesh.

    The mesh has 12 vertices and 14 triangular cells.
    - boundary_xneg has 2 edges (cells 0,2 which touch vertices 0,1,2)
    - boundary_xpos has 2 edges (cells 11,13 which touch vertices 9,10,11)
    """

    ENTITIES = {
        "domain": MeshEntity(ncells=14, ncorners=3, nvertices=12),
        # Materials
        "heat_material": MeshEntity(ncells=14, ncorners=3, nvertices=12),
        # Boundaries
        "bc_xneg": MeshEntity(ncells=2, ncorners=2, nvertices=3),
        "bc_xpos": MeshEntity(ncells=2, ncorners=2, nvertices=3),
    }


class Quad(object):
    """Mesh information for quad mesh.

    The mesh has 16 vertices and 9 quadrilateral cells.
    - boundary_xneg has 3 edges (vertices 0,1,2,3)
    - boundary_xpos has 3 edges (vertices 12,13,14,15)
    """

    ENTITIES = {
        "domain": MeshEntity(ncells=9, ncorners=4, nvertices=16),
        # Materials
        "heat_material": MeshEntity(ncells=9, ncorners=4, nvertices=16),
        # Boundaries
        "bc_xneg": MeshEntity(ncells=3, ncorners=2, nvertices=4),
        "bc_xpos": MeshEntity(ncells=3, ncorners=2, nvertices=4),
    }


# End of file
