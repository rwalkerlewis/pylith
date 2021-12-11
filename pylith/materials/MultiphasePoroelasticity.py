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
# @file pylith/materials/MultiphasePoroelasticity.py
#
# @brief Python object for solving the MultiphasePoroelasticity equation.
#
# Factory: material

from .Material import Material
from .materials import MultiphasePoroelasticity as ModuleMultiphasePoroelasticity

from .IsotropicLinearMultiphasePoroelasticity import IsotropicLinearMultiphasePoroelasticity


class MultiphasePoroelasticity(Material, ModuleMultiphasePoroelasticity):
    """Python material property manager.

    FACTORY: material
    """

    import pythia.pyre.inventory

    useBodyForce = pythia.pyre.inventory.bool("use_body_force", default=False)
    useBodyForce.meta['tip'] = "Include body force term in MultiphasePoroelasticity equation."

    useSourceDensity = pythia.pyre.inventory.bool("use_source_density", default=False)
    useSourceDensity.meta['tip'] = "Include source_density term in MultiphasePoroelasticity equation."

    useStateVars = pythia.pyre.inventory.bool(
        "use_state_variables", default=False)
    useStateVars.meta['tip'] = "Update auxiliary field terms with run."

    rheology = pythia.pyre.inventory.facility(
        "bulk_rheology", family="multiphaseporoelasticity_rheology", factory=IsotropicLinearMultiphasePoroelasticity)
    rheology.meta['tip'] = "Bulk rheology for poroelastic material."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="multiphaseporoelasticity"):
        """Constructor.
        """
        Material.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsMultiphasePoroelasticity import AuxSubfieldsMultiphasePoroelasticity
        self.auxiliarySubfields = AuxSubfieldsMultiphasePoroelasticity("auxiliary_subfields")

        from .DerivedSubfieldsElasticity import DerivedSubfieldsElasticity
        self.derivedSubfields = DerivedSubfieldsElasticity("derived_subfields")

    def preinitialize(self, problem):
        """Setup material.
        """
        self.rheology.preinitialize(problem)
        Material.preinitialize(self, problem)

        self.rheology.addAuxiliarySubfields(self, problem)

        ModuleMultiphasePoroelasticity.useBodyForce(self, self.useBodyForce)
        ModuleMultiphasePoroelasticity.useSourceDensity(self, self.useSourceDensity)
        ModuleMultiphasePoroelasticity.useStateVars(self, self.useStateVars)        
        return

    def _createModuleObj(self):
        """Create handle to C++ MultiphasePoroelasticity.
        """
        ModuleMultiphasePoroelasticity.__init__(self)
        ModuleMultiphasePoroelasticity.setBulkRheology(self, self.rheology)  # Material sets auxiliary db in rheology.
        return


# Factories

def material():
    """Factory associated with MultiphasePoroelasticity.
    """
    return MultiphasePoroelasticity()


# End of file
