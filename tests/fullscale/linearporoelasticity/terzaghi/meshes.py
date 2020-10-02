#!/usr/bin/env python
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
# @file tests/fullscale/poroelasticty/terzaghi/meshes.py
#
# @brief Mesh information for test cases.


class Tri(object):
    """
    Mesh information for tri mesh.
    """
    DOMAIN = {
        "ncells": 234,
        "ncorners": 3,
        "nvertices": 138,
    }
    MATERIALS = {
        "poroelastic": {
            "ncells": 234,
            "ncorners": 3,
            "nvertices": 138,
        }
    }
    BOUNDARIES = {
        "edge_xneg": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
        "edge_xpos": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
        "edge_yneg": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
        "edge_ypos": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
    }


class Quad(object):
    """
    Mesh information for quad mesh.
    """
    DOMAIN = {
        "ncells": 100,
        "ncorners": 4,
        "nvertices": 121,
    }
    MATERIALS = {
        "poroelastic_xpos": {
            "ncells": 100,
            "ncorners": 4,
            "nvertices": 121,
        }
    }
    BOUNDARIES = {
        "edge_xneg": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
        "edge_xpos": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
        "edge_yneg": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
        "edge_ypos": {
            "ncells": 10,
            "ncorners": 2,
            "nvertices": 11,
        },
    }


# End of file
