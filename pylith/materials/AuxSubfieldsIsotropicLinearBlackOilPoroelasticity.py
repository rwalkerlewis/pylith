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


class AuxSubfieldsIsotropicLinearBlackOilPoroelasticity(PetscComponent):
    """
    Auxiliary subfields associated with the isotropic linear black oil poroelastic bulk rheology.

    The black oil model extends standard poroelasticity with pressure-dependent
    fluid properties. Additional fields include:
    - reference_pressure: Reference pressure for fluid property correlations
    - fluid_compressibility: Fluid compressibility at reference pressure
    - fluid_compressibility_coefficient: Rate of change of compressibility with pressure
    - viscosity_coefficient: Rate of change of viscosity with pressure
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_blackoil.rheology.auxiliary_fields]
            shear_modulus.basis_order = 1
            biot_coefficient.basis_order = 0
            isotropic_permeability.basis_order = 0
            drained_bulk_modulus.basis_order = 1
            biot_modulus.basis_order = 1
            reference_pressure.basis_order = 0
            fluid_compressibility.basis_order = 0
            fluid_compressibility_coefficient.basis_order = 0
            viscosity_coefficient.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    shearModulus = pythia.pyre.inventory.facility("shear_modulus", family="auxiliary_subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    biotCoefficient = pythia.pyre.inventory.facility("biot_coefficient", family="auxiliary_subfield", factory=Subfield)
    biotCoefficient.meta['tip'] = "Biot coefficient subfield."

    isotropicPermeability = pythia.pyre.inventory.facility("isotropic_permeability", family="auxiliary_subfield", factory=Subfield)
    isotropicPermeability.meta['tip'] = "Isotropic permeability subfield."

    tensorPermeability = pythia.pyre.inventory.facility("tensor_permeability", family="auxiliary_subfield", factory=Subfield)
    tensorPermeability.meta['tip'] = "Tensor permeability subfield."

    drainedBulkModulus = pythia.pyre.inventory.facility("drained_bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    drainedBulkModulus.meta['tip'] = "Drained bulk modulus subfield."

    biotModulus = pythia.pyre.inventory.facility("biot_modulus", family="auxiliary_subfield", factory=Subfield)
    biotModulus.meta['tip'] = "Biot modulus subfield."

    referenceStress = pythia.pyre.inventory.facility("reference_stress", family="auxiliary_subfield", factory=Subfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pythia.pyre.inventory.facility("reference_strain", family="auxiliary_subfield", factory=Subfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    # Black oil specific fields
    referencePressure = pythia.pyre.inventory.facility("reference_pressure", family="auxiliary_subfield", factory=Subfield)
    referencePressure.meta['tip'] = "Reference pressure for fluid property correlations."

    fluidCompressibility = pythia.pyre.inventory.facility("fluid_compressibility", family="auxiliary_subfield", factory=Subfield)
    fluidCompressibility.meta['tip'] = "Fluid compressibility at reference pressure."

    fluidCompressibilityCoefficient = pythia.pyre.inventory.facility("fluid_compressibility_coefficient", family="auxiliary_subfield", factory=Subfield)
    fluidCompressibilityCoefficient.meta['tip'] = "Rate of change of fluid compressibility with pressure."

    viscosityCoefficient = pythia.pyre.inventory.facility("viscosity_coefficient", family="auxiliary_subfield", factory=Subfield)
    viscosityCoefficient.meta['tip'] = "Rate of change of fluid viscosity with pressure."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxfieldsisotropiclinearblackoilporoelasticity"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_fields")

    def _configure(self):
        PetscComponent._configure(self)

# FACTORIES ////////////////////////////////////////////////////////////


def auxiliary_subfields():
    """Factory associated with AuxSubfieldsIsotropicLinearBlackOilPoroelasticity.
    """
    return AuxSubfieldsIsotropicLinearBlackOilPoroelasticity()


# End of file
