# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from .Material import Material
from .materials import Heat as ModuleHeat

from .IsotropicHeat import IsotropicHeat


class Heat(Material, ModuleHeat):
    """
    Material behavior governed by the heat equation.

    Implements `Material`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_heat]
            label_value = 3
            use_heat_source = False
            bulk_rheology = pylith.materials.IsotropicHeat

            auxiliary_subfields.density.basis_order = 0
            auxiliary_subfields.specific_heat.basis_order = 0
            auxiliary_subfields.thermal_conductivity.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useHeatSource = pythia.pyre.inventory.bool("use_heat_source", default=False)
    useHeatSource.meta['tip'] = "Include heat source term in heat equation."

    rheology = pythia.pyre.inventory.facility("bulk_rheology", family="heat_rheology", factory=IsotropicHeat)
    rheology.meta['tip'] = "Bulk rheology for heat conduction."

    def __init__(self, name="heat"):
        """Constructor.
        """
        Material.__init__(self, name)

    def _defaults(self):
        from .AuxSubfieldsHeat import AuxSubfieldsHeat
        self.auxiliarySubfields = AuxSubfieldsHeat("auxiliary_subfields")

        from .DerivedSubfieldsHeat import DerivedSubfieldsHeat
        self.derivedSubfields = DerivedSubfieldsHeat("derived_subfields")

    def preinitialize(self, problem):
        """Setup material.
        """
        self.rheology.preinitialize(problem)
        Material.preinitialize(self, problem)

        self.rheology.addAuxiliarySubfields(self, problem)

        ModuleHeat.useHeatSource(self, self.useHeatSource)

    def _createModuleObj(self):
        """Create handle to C++ Heat.
        """
        ModuleHeat.__init__(self)
        ModuleHeat.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.


# Factories

def material():
    """Factory associated with Heat.
    """
    return Heat()


# End of file
