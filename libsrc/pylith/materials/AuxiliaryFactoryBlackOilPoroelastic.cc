// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/materials/AuxiliaryFactoryBlackOilPoroelastic.hh" // implementation of object methods

#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "pylith/scales/Scales.hh" // USES Scales
#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::AuxiliaryFactoryBlackOilPoroelastic(void) {
    GenericComponent::setName("AuxiliaryFactoryBlackOilPoroelastic");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::~AuxiliaryFactoryBlackOilPoroelastic(void) {}


// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addIsotropicPermeability(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addIsotropicPermeability(void)");

    const char* subfieldName = "isotropic_permeability";
    const PylithReal permeabilityScale = pylith::scales::ElasticityScales::getPermeabilityScale(*_scales);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = permeabilityScale;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addIsotropicPermeability


// ----------------------------------------------------------------------------
// Add tensor permeability subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addTensorPermeability(void) {
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
    const PylithReal permeabilityScale = pylith::scales::ElasticityScales::getPermeabilityScale(*_scales);

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

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addTensorPermeability


// --------------------------------------------------------------------
// Add drained bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addDrainedBulkModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDrainedBulkModulus(void)");

    const char* subfieldName = "drained_bulk_modulus";
    const PylithReal rigidityScale = _scales->getRigidityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = rigidityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addDrainedBulkModulus


// ---------------------------------------------------------------------
// Add biot coefficient subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addBiotCoefficient(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBiotCoefficient(void)");

    const char* subfieldName = "biot_coefficient";
    const PylithReal scale = 1.0;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = scale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBiotCoefficient


// ---------------------------------------------------------------------
// Add biot modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addBiotModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBiotModulus(void)");

    const char* subfieldName = "biot_modulus";
    const PylithReal rigidityScale = _scales->getRigidityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = rigidityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::biotModulusFromInput(subfieldName, this);

    PYLITH_METHOD_END;
} // addBiotModulus


// ----------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addReferenceStress(void) {
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
    const PylithReal stressScale = pylith::scales::ElasticityScales::getStressScale(*_scales);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = stressScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferenceStress


// ----------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addReferenceStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addReferenceStrain(void)");

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
    const PylithReal strainScale = pylith::scales::ElasticityScales::getStrainScale(*_scales);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = strainScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferenceStrain


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* subfieldName = "shear_modulus";
    const PylithReal rigidityScale = _scales->getRigidityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = rigidityScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addShearModulus


// ---------------------------------------------------------------------------------------------------------------------
// Add reference pressure subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addReferencePressure(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addReferencePressure(void)");

    const char* subfieldName = "reference_pressure";
    const PylithReal pressureScale = _scales->getRigidityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = NULL; // Can be negative for gauge pressure

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferencePressure


// ---------------------------------------------------------------------------------------------------------------------
// Add fluid compressibility subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addFluidCompressibility(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidCompressibility(void)");

    const char* subfieldName = "fluid_compressibility";
    const PylithReal rigidityScale = _scales->getRigidityScale();
    const PylithReal compressibilityScale = 1.0 / rigidityScale; // 1/Pa

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = compressibilityScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidCompressibility


// ---------------------------------------------------------------------------------------------------------------------
// Add fluid compressibility coefficient (pressure dependence) subfield.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addFluidCompressibilityCoefficient(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidCompressibilityCoefficient(void)");

    const char* subfieldName = "fluid_compressibility_coefficient";
    const PylithReal rigidityScale = _scales->getRigidityScale();
    const PylithReal coefficientScale = 1.0 / (rigidityScale * rigidityScale); // 1/Pa^2

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = coefficientScale;
    description.validator = NULL; // Can be any value

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidCompressibilityCoefficient


// ---------------------------------------------------------------------------------------------------------------------
// Add viscosity coefficient (pressure dependence) subfield.
void
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::addViscosityCoefficient(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addViscosityCoefficient(void)");

    const char* subfieldName = "viscosity_coefficient";
    const PylithReal rigidityScale = _scales->getRigidityScale();
    const PylithReal coefficientScale = 1.0 / rigidityScale; // 1/Pa

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = coefficientScale;
    description.validator = NULL; // Can be any value

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addViscosityCoefficient


// End of file
