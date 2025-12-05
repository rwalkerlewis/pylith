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

#include "pylith/materials/AuxiliaryFactoryHeat.hh" // implementation of object methods

#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryHeat::AuxiliaryFactoryHeat(void) {
    GenericComponent::setName("auxiliaryfactoryheat");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryHeat::~AuxiliaryFactoryHeat(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add density subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryHeat::addDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDensity(void)");

    const char* subfieldName = "density";
    // Density scale: M / L^3
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
// Add specific heat subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryHeat::addSpecificHeat(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSpecificHeat(void)");

    const char* subfieldName = "specific_heat";
    // Specific heat scale: L^2 / (T^2 * Theta), where Theta is temperature
    // Using pressure scale / (density * temperature) = L^2 / (t^2 * Theta)
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
pylith::materials::AuxiliaryFactoryHeat::addThermalConductivity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThermalConductivity(void)");

    const char* subfieldName = "thermal_conductivity";
    // Thermal conductivity scale: M * L / (T^3 * Theta)
    // Using pressure_scale * length_scale / (time_scale * temperature_scale)
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
pylith::materials::AuxiliaryFactoryHeat::addHeatSource(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addHeatSource(void)");

    const char* subfieldName = "heat_source";
    // Heat source scale: Power / Volume = M / (L * T^3)
    // Using pressure_scale / time_scale
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
    description.validator = NULL; // Heat source can be positive or negative

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addHeatSource


// End of file
