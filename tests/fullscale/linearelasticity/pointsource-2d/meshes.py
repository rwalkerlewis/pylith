#!/usr/bin/env nemesis
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


class TriGmsh(MeshEntity):
    """Mesh information for tri mesh using Gmsh."""
    ENTITIES = {
        "domain": MeshEntity(ncells=512, ncorners=3, nvertices=289),
        "bc_xneg": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "bc_xpos": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "bc_yneg": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "bc_ypos": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "elastic": MeshEntity(ncells=512, ncorners=3, nvertices=289),
    }

    def __init__(self):
        MeshEntity.__init__(self)
        for entity, info in self.ENTITIES.items():
            setattr(self, entity, info)


class QuadGmsh(MeshEntity):
    """Mesh information for quad mesh using Gmsh."""
    ENTITIES = {
        "domain": MeshEntity(ncells=256, ncorners=4, nvertices=289),
        "bc_xneg": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "bc_xpos": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "bc_yneg": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "bc_ypos": MeshEntity(ncells=16, ncorners=2, nvertices=17),
        "elastic": MeshEntity(ncells=256, ncorners=4, nvertices=289),
    }

    def __init__(self):
        MeshEntity.__init__(self)
        for entity, info in self.ENTITIES.items():
            setattr(self, entity, info)


# End of file
