# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

from pylith.materials.Material import Material
from .materials import Thermoporoelasticity as ModuleThermoporoelasticity


class Thermoporoelasticity(Material, ModuleThermoporoelasticity):
    """
    Thermoporoelasticity material for fully coupled thermal-hydraulic-mechanical (THM) simulations.
    
    This material couples:
    1. Momentum balance (elasticity with thermal and pore pressure effects)
    2. Fluid mass balance (Darcy flow with thermal effects)
    3. Heat equation (conduction with optional coupling)

    Implements `Material`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_thermoporoelastic]
            description = Thermoporoelastic material for THM coupling
            label_value = 1
            
            db_auxiliary_field = spatialdata.spatialdb.SimpleDB
            db_auxiliary_field.description = Thermoporoelastic properties
            db_auxiliary_field.iohandler.filename = mat_thermoporoelastic.spatialdb
            
            bulk_rheology = pylith.materials.IsotropicLinearThermoporoelasticity
        """
    }

    import pythia.pyre.inventory
    
    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in momentum equation."

    useSourceDensity = pythia.pyre.inventory.bool("use_source_density", default=False)
    useSourceDensity.meta['tip'] = "Include source density term (fluid source) in pressure equation."

    useHeatSource = pythia.pyre.inventory.bool("use_heat_source", default=False)
    useHeatSource.meta['tip'] = "Include heat source term in temperature equation."

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress and strain in computation of stress and strain."

    from .IsotropicLinearThermoporoelasticity import IsotropicLinearThermoporoelasticity
    bulkRheology = pythia.pyre.inventory.facility("bulk_rheology", family="thermoporoelasticity_rheology", 
                                                   factory=IsotropicLinearThermoporoelasticity)
    bulkRheology.meta['tip'] = "Bulk rheology for thermoporoelastic material."

    def __init__(self, name="thermoporoelasticity"):
        """Constructor.
        """
        Material.__init__(self, name)

    def _defaults(self):
        self.auxiliarySubfields.subfields.solid_density.basis_order = 0
        self.auxiliarySubfields.subfields.fluid_density.basis_order = 0
        self.auxiliarySubfields.subfields.fluid_viscosity.basis_order = 0
        self.auxiliarySubfields.subfields.porosity.basis_order = 0

    def preinitialize(self, problem):
        """Setup material.
        """
        self.bulkRheology.preinitialize(problem)

        Material.preinitialize(self, problem)

        ModuleThermoporoelasticity.useBodyForce(self, self.useBodyForce)
        ModuleThermoporoelasticity.useSourceDensity(self, self.useSourceDensity)
        ModuleThermoporoelasticity.useHeatSource(self, self.useHeatSource)
        ModuleThermoporoelasticity.useReferenceState(self, self.useReferenceState)
        ModuleThermoporoelasticity.setBulkRheology(self, self.bulkRheology)

    def _createModuleObj(self):
        """Create handle to C++ object.
        """
        ModuleThermoporoelasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////
def material():
    """Factory associated with Thermoporoelasticity.
    """
    return Thermoporoelasticity()


# End of file
