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
# @brief Python container for poroelasticity equation subfields.

from pylith.utils.PetscComponent import PetscComponent


class AuxSubfieldsPoroelasticity(PetscComponent):
    """
    Python container for poroelasticity equation subfields.

    INVENTORY

    Properties
      - None

    Facilities
      - *porosity* Porosity subfield.
      - *density* Rock Density subfield.
      - *fluid_density* Fluid Density subfield.
      - *fluid_viscosity* Fluid viscosity subfield.
      - *body_force* Body force.
      - *gravitational_acceleration* Gravitational acceleration subfield.
    """

    import pyre.inventory

    from pylith.topology.Subfield import Subfield

    porosity = pyre.inventory.facility("porosity", family="auxiliary_subfield", factory=Subfield)
    porosity.meta['tip'] = "Porosity subfield."

    density = pyre.inventory.facility("density", family="auxiliary_subfield", factory=Subfield)
    density.meta['tip'] = "Density subfield."

    fluidDensity = pyre.inventory.facility("fluid_density", family="auxiliary_subfield", factory=Subfield)
    fluidDensity.meta['tip'] = "Fluid density subfield."

    fluidViscosity = pyre.inventory.facility("fluid_viscosity", family="auxiliary_subfield", factory=Subfield)
    fluidViscosity.meta['tip'] = "Fluid viscosity subfield."

    bodyForce = pyre.inventory.facility("body_force", family="auxiliary_subfield", factory=Subfield)
    bodyForce.meta['tip'] = "Body force subfield."

    gravitationalAcceleration = pyre.inventory.facility(
        "gravitational_acceleration", family="auxiliary_subfield", factory=Subfield)
    gravitationalAcceleration.meta['tip'] = "Gravitational acceleration subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxsubfieldsporoelasticity"):
        """
        Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """
    Factory associated with AuxSubfieldsPoroelasticity.
    """
    return AuxSubfieldsPoroelasticity()


# End of file
