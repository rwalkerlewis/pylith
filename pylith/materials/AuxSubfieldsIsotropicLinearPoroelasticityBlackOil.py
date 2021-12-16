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

# @file pylith/materials/AuxFieldsIsotropicLinearPoroelasticityBlackOil.py
##
# @brief Python subfields container for isotropic, linear poroelasticityblackoil
# subfields.

from pylith.utils.PetscComponent import PetscComponent

# AuxFieldsIsotropicLinearPoroelasticityBlackOil class


class AuxSubfieldsIsotropicLinearPoroelasticityBlackOil(PetscComponent):
    """Python subfields container for isotropic, linear poroelasticityblackoil subfields.

    FACTORY: auxiliary_subfields
    """

    import pythia.pyre.inventory

    from pylith.topology.Subfield import Subfield

    referenceStress = pythia.pyre.inventory.facility(
        "reference_stress", family="auxiliary_subfield", factory=Subfield)
    referenceStress.meta['tip'] = "Reference stress subfield."

    referenceStrain = pythia.pyre.inventory.facility(
        "reference_strain", family="auxiliary_subfield", factory=Subfield)
    referenceStrain.meta['tip'] = "Reference strain subfield."

    shearModulus = pythia.pyre.inventory.facility(
        "shear_modulus", family="auxiliary_subfield", factory=Subfield)
    shearModulus.meta['tip'] = "Shear modulus subfield."

    biotCoefficient = pythia.pyre.inventory.facility(
        "biot_coefficient", family="auxiliary_subfield", factory=Subfield)
    biotCoefficient.meta['tip'] = "Biot coefficient subfield."

    isotropicPermeability = pythia.pyre.inventory.facility(
        "isotropic_permeability", family="auxiliary_subfield", factory=Subfield)
    isotropicPermeability.meta['tip'] = "Isotropic permeability subfield."

    tensorPermeability = pythia.pyre.inventory.facility(
        "tensor_permeability", family="auxiliary_subfield", factory=Subfield)
    tensorPermeability.meta['tip'] = "Tensor permeability subfield."

    solidBulkModulus = pythia.pyre.inventory.facility(
        "solid_bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    solidBulkModulus.meta['tip'] = "Solid bulk modulus subfield."

    fluidBulkModulus = pythia.pyre.inventory.facility(
        "fluid_bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    fluidBulkModulus.meta['tip'] = "Fluid bulk modulus subfield."

    drainedBulkModulus = pythia.pyre.inventory.facility(
        "drained_bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    drainedBulkModulus.meta['tip'] = "Drained bulk modulus subfield."

    undrainedBulkModulus = pythia.pyre.inventory.facility(
        "undrained_bulk_modulus", family="auxiliary_subfield", factory=Subfield)
    undrainedBulkModulus.meta['tip'] = "Undrained bulk modulus subfield."

    fluidViscosity = pythia.pyre.inventory.facility("fluid_viscosity", family="auxiliary_subfield", factory=Subfield)
    fluidViscosity.meta['tip'] = "Fluid viscosity subfield."

    fluidSaturation = pythia.pyre.inventory.facility(
        "fluid_saturation", family="auxiliary_subfield", factory=Subfield)
    fluidSaturation.meta['tip'] = "Fluid saturation subfield."

    relativePermeability = pythia.pyre.inventory.facility(
        "relative_permeability", family="auxiliary_subfield", factory=Subfield)
    relativePermeability.meta['tip'] = "Relative permeability subfield."

    formationVolumeFactor = pythia.pyre.inventory.facility(
        "formation_volume_factor", family="auxiliary_subfield", factory=Subfield)
    formationVolumeFactor.meta['tip'] = "Formation volume factor subfield."

    solutionGasOilRatio = pythia.pyre.inventory.facility(
        "solution_gas_oil_ratio", family="auxiliary_subfield", factory=Subfield)
    solutionGasOilRatio.meta['tip'] = "Solution gas oil ratio."

    solutionOilGasRatio = pythia.pyre.inventory.facility(
        "solution_oil_gas_ratio", family="auxiliary_subfield", factory=Subfield)
    solutionOilGasRatio.meta['tip'] = "Solution oil gas ratio."

    fluidDensity = pythia.pyre.inventory.facility("fluid_density", family="auxiliary_subfield", factory=Subfield)
    fluidDensity.meta['tip'] = "Fluid density subfield."

    youngsModulus = pythia.pyre.inventory.facility(
        "youngs_modulus", family="auxiliary_subfield", factory=Subfield)
    youngsModulus.meta['tip'] = "Young's modulus subfield."

    poissonsRatio = pythia.pyre.inventory.facility(
        "poissons_ratio", family="auxiliary_subfield", factory=Subfield)
    poissonsRatio.meta['tip'] = "Poisson's ratio subfield."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="auxfieldsisotropiclinearporoelasticityblackoil"):
        """Constructor.
        """
        PetscComponent.__init__(self, name, facility="auxiliary_fields")
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _configure(self):
        PetscComponent._configure(self)
        return

# FACTORIES ////////////////////////////////////////////////////////////


def auxiliary_subfields():
    """Factory associated with AuxSubfieldsIsotropicLinearPoroelasticityBlackOil.
    """
    return AuxSubfieldsIsotropicLinearPoroelasticityBlackOil()


# End of file
