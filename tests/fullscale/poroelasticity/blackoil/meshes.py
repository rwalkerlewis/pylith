# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""Mesh information for tests."""


class QuadGmsh(object):
    """Mesh information for quad mesh using Gmsh."""

    NCELLS = 64
    NVERTICES = 81
    CELL_DIM = 2

    ENTITIES = {
        "domain": {
            "ncells": NCELLS,
            "nvertices": NVERTICES,
        },
        "poroelastic": {
            "ncells": NCELLS,
            "nvertices": NVERTICES,
        },
        "bc_xneg": {
            "ncells": 8,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 8,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 8,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 8,
            "nvertices": 9,
        },
    }


class TriGmsh(object):
    """Mesh information for tri mesh using Gmsh."""

    NCELLS = 128
    NVERTICES = 81
    CELL_DIM = 2

    ENTITIES = {
        "domain": {
            "ncells": NCELLS,
            "nvertices": NVERTICES,
        },
        "poroelastic": {
            "ncells": NCELLS,
            "nvertices": NVERTICES,
        },
        "bc_xneg": {
            "ncells": 8,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 8,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 8,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 8,
            "nvertices": 9,
        },
    }


# End of file
