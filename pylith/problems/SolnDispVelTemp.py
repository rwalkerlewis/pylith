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
from .SolutionSubfield import subfieldFactory


class SolnDispVelTemp(PetscComponent):
    """
    Solution field with displacement, velocity, and temperature subfields.

    Used for dynamic thermoelasticity with two-way coupling.

    :::{important}
    Use this solution for dynamic thermoelasticity simulations where
    thermoelastic heating (adiabatic heating/cooling from volumetric strain rate)
    is important. For quasistatic problems, use SolnDispTemp instead.
    :::
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem]
            solution = pylith.problems.SolnDispVelTemp
        """
    }

    import pythia.pyre.inventory

    from .SubfieldDisplacement import SubfieldDisplacement
    displacement = pythia.pyre.inventory.facility("displacement", family="soln_subfield", factory=SubfieldDisplacement)
    displacement.meta['tip'] = "Displacement subfield."

    from .SubfieldVelocity import SubfieldVelocity
    velocity = pythia.pyre.inventory.facility("velocity", family="soln_subfield", factory=SubfieldVelocity)
    velocity.meta['tip'] = "Velocity subfield."

    from .SubfieldTemperature import SubfieldTemperature
    temperature = pythia.pyre.inventory.facility("temperature", family="soln_subfield", factory=SubfieldTemperature)
    temperature.meta['tip'] = "Temperature subfield."

    subfields = pythia.pyre.inventory.facilityArray(
        "subfields", family="soln_subfields", itemFactory=subfieldFactory, factory=SolnDispVelTemp)
    subfields.meta['tip'] = "Solution subfields."

    def __init__(self, name="solndispveltemp"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="soln_subfields")

    def _defaults(self):
        self.subfields = ["displacement", "velocity", "temperature"]

    def components(self):
        """Get list of solution field components.
        """
        return [self.displacement, self.velocity, self.temperature]


# Factories

def soln_subfields():
    """Factory associated with SolnDispVelTemp.
    """
    return SolnDispVelTemp()


# End of file
