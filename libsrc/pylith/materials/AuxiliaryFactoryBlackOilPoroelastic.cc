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
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::AuxiliaryFactoryBlackOilPoroelastic(void) :
    AuxiliaryFactoryPoroelastic() {
    GenericComponent::setName("AuxiliaryFactoryBlackOilPoroelastic");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryBlackOilPoroelastic::~AuxiliaryFactoryBlackOilPoroelastic(void) {}


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
