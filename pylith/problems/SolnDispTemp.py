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
# @file pylith/problems/SolnDispTemp.py
#
# @brief Python subfields container with displacement and temperature subfields.

from pylith.utils.PetscComponent import PetscComponent
from .Solution import Solution as SolutionBase


class SolnDispTemp(PetscComponent):
    """Python subfields container with displacement and temperature subfields.

    IMPORTANT: Use the Solution class (below) to set this object as the default facilities array for the solution
    subfields.
    """

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pythia.pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldTemperature import SubfieldTemperature
    temperature = pythia.pyre.inventory.facility("temperature", family="soln_subfield", factory=SubfieldTemperature)
    temperature.meta['tip'] = "Temperature subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="solndisptemp"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")
        return

    def _configure(self):
        PetscComponent._configure(self)
        return

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement, temperature].

        """
        return [self.displacement, self.temperature]


class Solution(SolutionBase):
    """Python solution field with displacement, temperature, and Lagrange multiplier subfields.
    """

    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray("subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispTemp)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
