# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/faults/SingleRupturePoro.py
#
# @brief Python kinematic rupture container with one poroelastic rupture.

from pylith.utils.PetscComponent import PetscComponent


class SingleRupturePoro(PetscComponent):
    """Python kinematic rupture container with one rupture.
    """

    import pythia.pyre.inventory

    from .KinSrcPoroStep import KinSrcPoroStep
    rupture = pythia.pyre.inventory.facility("rupture", family="eq_kinematic_src", factory=KinSrcPoroStep)
    rupture.meta['tip'] = "Kinematic, poroelastic earthquake rupture in problem."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="singleruptureporo"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="rupture")
        return


# End of file
