# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.utils.PetscComponent import PetscComponent
from .Solution import Solution as SolutionBase


class SolnDispTempLagrange(PetscComponent):
    """
    Container for solution subfields with displacement, temperature, and fault Lagrange multiplier subfields.

    This solution type is used for thermoelasticity problems with faults, enabling:
    - Coupled thermal-mechanical deformation
    - Prescribed slip on faults (earthquake rupture)
    - Heat generation from fault slip (shear heating)
    - Heat conduction throughout the domain
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem]
            solution = pylith.problems.SolnDispTempLagrange
        """
    }

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pythia.pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldTemperature import SubfieldTemperature
    temperature = pythia.pyre.inventory.facility("temperature", family="soln_subfield", factory=SubfieldTemperature)
    temperature.meta['tip'] = "Temperature subfield."

    from .SubfieldLagrangeFault import SubfieldLagrangeFault
    lagrangeFault = pythia.pyre.inventory.facility("lagrange_multiplier_fault", family="soln_subfield", factory=SubfieldLagrangeFault)
    lagrangeFault.meta['tip'] = "Fault Lagrange multiplier subfield."

    def __init__(self, name="solndisptemplagrange"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")

    def _configure(self):
        PetscComponent._configure(self)

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement, temperature, lagrange_multiplier_fault].

        """
        return [self.displacement, self.temperature, self.lagrangeFault]


class Solution(SolutionBase):
    """Python solution field with displacement, temperature, and Lagrange multiplier subfields.
    """
    import pythia.pyre.inventory

    from .SolutionSubfield import subfieldFactory
    subfields = pythia.pyre.inventory.facilityArray("subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispTempLagrange)
    subfields.meta['tip'] = "Subfields in solution."


# FACTORIES ////////////////////////////////////////////////////////////
def solution():
    """Factory associated with Solution.
    """
    return Solution()


# End of file
