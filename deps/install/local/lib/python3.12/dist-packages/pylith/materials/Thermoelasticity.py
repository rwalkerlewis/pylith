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
from .materials import Thermoelasticity as ModuleThermoelasticity

from .IsotropicLinearThermoelasticity import IsotropicLinearThermoelasticity


class Thermoelasticity(Material, ModuleThermoelasticity):
    """
    Material behavior governed by coupled thermoelasticity with two-way coupling.

    This material provides full two-way thermoelastic coupling:

    **Quasistatic formulation** (solution = SolnDispTemp):
    - Temperature → Stress: thermal strain σ = C:(ε - α(T-Tref)I)
    - One-way coupling (thermoelastic heating negligible for slow processes)

    **Dynamic formulation** (solution = SolnDispVelTemp):
    - Temperature → Stress: thermal strain
    - Velocity → Temperature: thermoelastic heating ρc Ṫ = k∇²T + Tα3K∇·v
    - Full two-way coupling for adiabatic heating/cooling

    The thermoelastic heating term represents:
    - Compression (∇·v < 0): material heats up
    - Expansion (∇·v > 0): material cools down

    Implements `Material`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_thermoelastic]
            label_value = 3
            use_body_force = False
            use_heat_source = False
            bulk_rheology = pylith.materials.IsotropicLinearThermoelasticity

            auxiliary_subfields.density.basis_order = 0
            auxiliary_subfields.specific_heat.basis_order = 0
            auxiliary_subfields.thermal_conductivity.basis_order = 0
            auxiliary_subfields.reference_temperature.basis_order = 0
            auxiliary_subfields.thermal_expansion_coefficient.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force in momentum equation."

    useHeatSource = pythia.pyre.inventory.bool("use_heat_source", default=False)
    useHeatSource.meta['tip'] = "Include heat source term in heat equation."

    rheology = pythia.pyre.inventory.facility("bulk_rheology", family="thermoelasticity_rheology", factory=IsotropicLinearThermoelasticity)
    rheology.meta['tip'] = "Bulk rheology for thermoelasticity."

    def __init__(self, name="thermoelasticity"):
        """Constructor.
        """
        Material.__init__(self, name)

    def _defaults(self):
        from .AuxSubfieldsThermoelasticity import AuxSubfieldsThermoelasticity
        self.auxiliarySubfields = AuxSubfieldsThermoelasticity("auxiliary_subfields")

        from .DerivedSubfieldsElasticity import DerivedSubfieldsElasticity
        self.derivedSubfields = DerivedSubfieldsElasticity("derived_subfields")

    def preinitialize(self, problem):
        """Setup material.
        """
        self.rheology.preinitialize(problem)
        Material.preinitialize(self, problem)

        self.rheology.addAuxiliarySubfields(self, problem)

        ModuleThermoelasticity.useBodyForce(self, self.useBodyForce)
        ModuleThermoelasticity.useHeatSource(self, self.useHeatSource)

    def _createModuleObj(self):
        """Create handle to C++ Thermoelasticity.
        """
        ModuleThermoelasticity.__init__(self)
        ModuleThermoelasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.


# Factories

def material():
    """Factory associated with Thermoelasticity.
    """
    return Thermoelasticity()


# End of file
