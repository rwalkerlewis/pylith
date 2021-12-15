// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AuxiliaryFactoryPoroelasticBlackOil.hh" // implementation of object methods

#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::AuxiliaryFactoryPoroelasticBlackOil(void) {
    GenericComponent::setName("AuxiliaryFactoryPoroelasticBlackOil");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::~AuxiliaryFactoryPoroelasticBlackOil(void) {}


// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addIsotropicPermeability(void) { // isotropicPermeablity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addIsotropicPermeability(void)");

    const char* subfieldName = "isotropic_permeability";

    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal permeabilityScale = lengthScale*lengthScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = permeabilityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addIsotropicPermeability


// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addTensorPermeability(void) { // isotropicPermeablity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTensorPermeability(void)");

    const char* subfieldName = "tensor_permeability";
    const char* componentNames[6] = {
        "permeability_xx",
        "permeability_yy",
        "permeability_zz",
        "permeability_xy",
        "permeability_yz",
        "permeability_xz"
    };
    const int tensorSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal permeabilityScale = lengthScale*lengthScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = tensorSize;
    description.componentNames.resize(tensorSize);
    for (int i = 0; i < tensorSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = permeabilityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addTensorPermeability


// --------------------------------------------------------------------
// Add solid bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addSolidBulkModulus(void) { // solidBulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSolidBulkModulus(void)");

    const char* subfieldName = "solid_bulk_modulus";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSolidBulkModulus


// ---------------------------------------------------------------------
// Add biot coefficient subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addBiotCoefficient(void) { // biotCoefficient
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBiotCoefficient(void)");

    const char* subfieldName = "biot_coefficient";

    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBiotCoefficient


// --------------------------------------------------------------------
// Add drained bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addDrainedBulkModulus(void) { // DrainedBulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDrainedBulkModulus(void)");

    const char* subfieldName = "drained_bulk_modulus";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addDrainedBulkModulus


// --------------------------------------------------------------------
// Add undrained bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addUndrainedBulkModulus(void) { // UndrainedBulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addUndrainedBulkModulus(void)");

    const char* subfieldName = "undrained_bulk_modulus";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addUndrainedBulkModulus


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* subfieldName = "shear_modulus";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addShearModulus


// --------------------------------------------------------------------
// Add fluid bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addThreePhaseFluidModulus(void) { // threePhaseFluidModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThreePhaseFluidModulus(void)");

    const char* subfieldName = "fluid_bulk_modulus";
    const char* componentNames[3] = {
        "water_bulk_modulus",
        "oil_bulk_modulus",
        "gas_bulk_modulus"
    };
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 3;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThreePhaseFluidBulkModulus


// --------------------------------------------------------------------
// Add fluid saturation values subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addThreePhaseSaturation(void) { // threePhaseSaturation
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThreePhaseSaturation(void)");

    const char* subfieldName = "fluid_saturation";
    const PylithReal noScale = 1;
    const char* componentNames[3] = {
        "water_saturation",
        "oil_saturation",
        "gas_saturation"
    };

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 3;
    description.hasHistory = true;
    description.historySize = 1;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThreePhaseSaturation


// --------------------------------------------------------------------
// Add fluid viscosity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addThreePhaseFluidViscosity(void) { // threePhaseFluidViscosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThreePhaseFluidViscosity(void)");

    const char* subfieldName = "fluid_viscosity";
    const char* componentNames[3] = {
        "water_viscosity",
        "oil_viscosity",
        "gas_viscosity"
    };
    const PylithReal pressureScale = _normalizer->getPressureScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal viscosityScale = pressureScale*timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 3;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = viscosityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThreePhaseFluidViscosity


// --------------------------------------------------------------------
// Add relative permeability values subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addRelativePermeability(void) { // relativePermeability
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRelativePermeability(void)");

    const char* subfieldName = "relative_permeability";
    const PylithReal noScale = 1;
    const char* componentNames[3] = {
        "relative_water",
        "relative_oil",
        "relative_gas"
    };

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 3;
    description.hasHistory = true;
    description.historySize = 1;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addRelativePermeability


// --------------------------------------------------------------------
// Add formation volume factor values subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addFormationVolumeFactors(void) { // formationVolumeFactors
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFormationVolumeFactors(void)");

    const char* subfieldName = "formation_volume_factor";
    const PylithReal noScale = 1;
    const char* componentNames[3] = {
        "B_w",
        "B_o",
        "B_g"
    };

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 3;
    description.hasHistory = true;
    description.historySize = 1;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFormationVolumeFactors


// ---------------------------------------------------------------------
// Add solution gas oil ratio subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addSolutionGasOilRatio(void) { // solutionGasOilRatio
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("solutionGasOilRatio(void)");

    const char* subfieldName = "solution_gas_oil_ratio";
    const PylithReal noScale = 1;
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSolutionGasOilRatio


// ---------------------------------------------------------------------
// Add solution oil gas ratio subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addSolutionOilGasRatio(void) { // solutionGasOilRatio
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("solutionGasOilRatio(void)");

    const char* subfieldName = "solution_oil_gas_ratio";
    const PylithReal noScale = 1;
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSolutionOilGasRatio


// ----------------------------------------------------------------------
// Add fluid density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addFluidDensities(void) { // fluidDensities
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidDensities(void)");

    const char* subfieldName = "fluid_density";
    const char* componentNames[3] = {
        "water_density",
        "oil_density",
        "gas_density"
    };
    const PylithReal densityScale = _normalizer->getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 3;
    description.hasHistory = true;
    description.historySize = 1;
    description.componentNames.resize(3);
    for (int i = 0; i < 3; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidDensities


// ----------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addReferenceStress(void) { // referenceStress
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addReferenceStress(void)");

    const char* subfieldName = "reference_stress";
    const char* componentNames[6] = {
        "reference_stress_xx",
        "reference_stress_yy",
        "reference_stress_zz",
        "reference_stress_xy",
        "reference_stress_yz",
        "reference_stress_xz"
    };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferenceStress


// ----------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addReferenceStrain(void) { // addReferenceStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRefrenceStrain(void)");

    const char* subfieldName = "reference_strain";
    const char* componentNames[6] = {
        "reference_strain_xx",
        "reference_strain_yy",
        "reference_strain_zz",
        "reference_strain_xy",
        "reference_strain_yz",
        "reference_strain_xz"
    };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferenceStrain


// ---------------------------------------------------------------------------------------------------------------------
// Add young's modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addYoungsModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addYoungsModulus(void)");

    const char* subfieldName = "youngs_modulus";
    const PylithReal pressureScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addYoungsModulus


// ----------------------------------------------------------------------
// Add poisson's ratio subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticBlackOil::addPoissonsRatio(void) { // porosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPoissonsRatio(void)");

    const char* subfieldName = "poissons_ratio";

    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = noScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPoissonsRatio


// End of file
