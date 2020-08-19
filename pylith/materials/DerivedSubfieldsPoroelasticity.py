# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/materials/AuxSubieldsPoroelasticity.py
#
# @brief Python container for elasticity equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class DerivedSubfieldsPoroelasticity(PetscComponent):
    """
    Python container for derived subfields for poroelasticity.

    INVENTORY

    Properties
      - None

    Facilities
      - *cauchy_stress* Cauchy stress subfield.
      - *cauchy_strain* Cauchy strain subfield.
    """

    import pyre.inventory

    from pylith.topology.Subfield import Subfield

    cauchyStress = pyre.inventory.facility("cauchy_stress", family="subfield", factory=Subfield)
    cauchyStress.meta['tip'] = "Cauchy stress subfield."

    cauchyStrain = pyre.inventory.facility("cauchy_strain", family="subfield", factory=Subfield)
    cauchyStrain.meta['tip'] = "Cauchy strain subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="derivedsubfieldsporoelasticity"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="derived_subfields")
        return


# FACTORIES ////////////////////////////////////////////////////////////

def derived_subfields():
    """
    Factory associated with DerivedSubfieldsPoroelasticity.
    """
    return DerivedSubfieldsPoroelasticity()


# End of file
