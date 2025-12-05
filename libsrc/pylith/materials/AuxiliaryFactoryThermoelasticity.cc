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

#include "pylith/materials/AuxiliaryFactoryThermoelasticity.hh" // implementation of object methods

#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryThermoelasticity::AuxiliaryFactoryThermoelasticity(void) {
    GenericComponent::setName("auxiliaryfactorythermoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryThermoelasticity::~AuxiliaryFactoryThermoelasticity(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDensity(void)");

    const char* subfieldName = "density";
    const PylithReal densityScale = _scales->getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add body force subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addBodyForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBodyForce(void)");

    const char* subfieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    const PylithReal forceScale = _scales->getDensityScale() * _scales->getLengthScale() /
                                  (_scales->getTimeScale() * _scales->getTimeScale());

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = forceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBodyForce


// ---------------------------------------------------------------------------------------------------------------------
// Add gravity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addGravityField(spatialdata::spatialdb::GravityField* gf) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addGravityField(gf="<<gf<<")");

    const char* subfieldName = "gravitational_acceleration";
    const char* componentNames[3] = { "gravitational_acceleration_x", "gravitational_acceleration_y", "gravitational_acceleration_z" };

    const PylithReal accelerationScale = _scales->getLengthScale() / (_scales->getTimeScale() * _scales->getTimeScale());

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = accelerationScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    pylith::materials::Query::gravityFieldFromDB(subfieldName, this, gf);

    PYLITH_METHOD_END;
} // addGravityField


// ---------------------------------------------------------------------------------------------------------------------
// Add specific heat subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addSpecificHeat(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSpecificHeat(void)");

    const char* subfieldName = "specific_heat";
    const PylithReal pressureScale = _scales->getPressureScale();
    const PylithReal temperatureScale = _scales->getTemperatureScale();
    const PylithReal specificHeatScale = pressureScale / (_scales->getDensityScale() * temperatureScale);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = specificHeatScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSpecificHeat


// ---------------------------------------------------------------------------------------------------------------------
// Add thermal conductivity subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addThermalConductivity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThermalConductivity(void)");

    const char* subfieldName = "thermal_conductivity";
    const PylithReal pressureScale = _scales->getPressureScale();
    const PylithReal lengthScale = _scales->getLengthScale();
    const PylithReal timeScale = _scales->getTimeScale();
    const PylithReal temperatureScale = _scales->getTemperatureScale();
    const PylithReal thermalConductivityScale = pressureScale * lengthScale / (timeScale * temperatureScale);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = thermalConductivityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThermalConductivity


// ---------------------------------------------------------------------------------------------------------------------
// Add heat source subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addHeatSource(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addHeatSource(void)");

    const char* subfieldName = "heat_source";
    const PylithReal pressureScale = _scales->getPressureScale();
    const PylithReal timeScale = _scales->getTimeScale();
    const PylithReal heatSourceScale = pressureScale / timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = heatSourceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addHeatSource


// ---------------------------------------------------------------------------------------------------------------------
// Add reference temperature subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addReferenceTemperature(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addReferenceTemperature(void)");

    const char* subfieldName = "reference_temperature";
    const PylithReal temperatureScale = _scales->getTemperatureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = temperatureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addReferenceTemperature


// ---------------------------------------------------------------------------------------------------------------------
// Add thermal expansion coefficient subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addThermalExpansionCoeff(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThermalExpansionCoeff(void)");

    const char* subfieldName = "thermal_expansion_coefficient";
    // Thermal expansion coefficient has units of 1/K (per Kelvin)
    const PylithReal temperatureScale = _scales->getTemperatureScale();
    const PylithReal thermalExpansionScale = 1.0 / temperatureScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = thermalExpansionScale;
    description.validator = NULL; // Can be positive or negative (though usually positive)

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThermalExpansionCoeff


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* subfieldName = "shear_modulus";
    const PylithReal pressureScale = _scales->getPressureScale();

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
    pylith::materials::Query::dbQueryVs(subfieldName, this);

    PYLITH_METHOD_END;
} // addShearModulus


// ---------------------------------------------------------------------------------------------------------------------
// Add bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryThermoelasticity::addBulkModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulus(void)");

    const char* subfieldName = "bulk_modulus";
    const PylithReal pressureScale = _scales->getPressureScale();

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
    pylith::materials::Query::dbQueryBulkModulus(subfieldName, this);

    PYLITH_METHOD_END;
} // addBulkModulus


// End of file
