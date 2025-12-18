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


class AuxSubfieldsHeat(PetscComponent):
    """
    Auxiliary subfields associated with the heat equation.

    Setting the parameters for a subfield does not turn on its use.
    The [`Heat` Component](Heat.md) has flags for including or excluding terms in the heat equation.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_heat.auxiliary_fields]
            density.basis_order = 0
            specific_heat.basis_order = 0
            thermal_conductivity.basis_order = 0
            heat_source.basis_order = 0
        """
    }
    """Python container for heat equation subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    density = pythia.pyre.inventory.facility("density", family="auxiliary_subfield", factory=Subfield)
    density.meta['tip'] = "Density subfield."

    specificHeat = pythia.pyre.inventory.facility("specific_heat", family="auxiliary_subfield", factory=Subfield)
    specificHeat.meta['tip'] = "Specific heat subfield."

    thermalConductivity = pythia.pyre.inventory.facility("thermal_conductivity", family="auxiliary_subfield", factory=Subfield)
    thermalConductivity.meta['tip'] = "Thermal conductivity subfield."

    heatSource = pythia.pyre.inventory.facility("heat_source", family="auxiliary_subfield", factory=Subfield)
    heatSource.meta['tip'] = "Heat source subfield."

    def __init__(self, name="auxsubfieldsheat"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_subfields")

    def _configure(self):
        PetscComponent._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////

def auxiliary_subfields():
    """Factory associated with AuxSubfieldsHeat.
    """
    return AuxSubfieldsHeat()


# End of file
