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
# @file tests/fullscale/linearelasticity/nofaults/2d/meshes.py
#
# @brief Mesh information for test cases.


class Tri(object):
    """Mesh information for tri mesh.
    """
    DOMAIN = {
        "ncells": 124,
        "ncorners": 3,
        "nvertices": 88,
    }
    MATERIALS = {
        "mat_xneg": {
            "ncells": 30,
            "ncorners": 3,
            "nvertices": 26,
        },
        "mat_xmid": {
            "ncells": 30,
            "ncorners": 3,
            "nvertices": 26,
        },
        "mat_xposypos": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
        "mat_xposyneg": {
            "ncells": 32,
            "ncorners": 3,
            "nvertices": 25,
        },
    }
    FAULTS = {
        "fault": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        }
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 10,
        },
    }


class Quad(object):
    """Mesh information for tri mesh.
    """
    DOMAIN = {
        "ncells": 64,
        "ncorners": 4,
        "nvertices": 90,
    }
    MATERIALS = {
        "mat_xneg": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 27,
        },
        "mat_xmid": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 27,
        },
        "mat_xposypos": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
        },
        "mat_xposyneg": {
            "ncells": 16,
            "ncorners": 4,
            "nvertices": 25,
        },
    }
    FAULTS = {
        "fault": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        }
    }
    BOUNDARIES = {
        "bc_xneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_xpos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_yneg": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 9,
        },
        "bc_ypos": {
            "ncells": 8,
            "ncorners": 2,
            "nvertices": 10,
        },
    }


# End of file
