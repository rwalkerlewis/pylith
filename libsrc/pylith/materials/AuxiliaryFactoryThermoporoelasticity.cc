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

#include "pylith/materials/AuxiliaryFactoryThermoporoelasticity.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryThermoporoelasticity::AuxiliaryFactoryThermoporoelasticity(void) {
    GenericComponent::setName("auxiliaryfactorythermoporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryThermoporoelasticity::~AuxiliaryFactoryThermoporoelasticity(void) {}


// ============================= Base poroelastic fields =============================

// ---------------------------------------------------------------------------------------------------------------------
// Add solid density subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addSolidDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addSolidDensity(void)");

    const char* subfieldName = "solid_density";
    const PylithReal densityScale = _normalizer->getDensityScale();

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
} // addSolidDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add fluid density subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addFluidDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addFluidDensity(void)");

    const char* subfieldName = "fluid_density";
    const PylithReal densityScale = _normalizer->getDensityScale();

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
} // addFluidDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add fluid viscosity subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addFluidViscosity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addFluidViscosity(void)");

    const char* subfieldName = "fluid_viscosity";
    const PylithReal pressureScale = _normalizer->getPressureScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal viscosityScale = pressureScale * timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = viscosityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidViscosity


// ---------------------------------------------------------------------------------------------------------------------
// Add porosity subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addPorosity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addPorosity(void)");

    const char* subfieldName = "porosity";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPorosity


// ============================= Optional fields =============================

// ---------------------------------------------------------------------------------------------------------------------
// Add body force subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addBodyForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addBodyForce(void)");

    const char* subfieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    const PylithReal forceScale = _normalizer->getPressureScale() / _normalizer->getLengthScale();

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
// Add gravity subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addGravityField(spatialdata::spatialdb::GravityField* gf) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addGravityField(gf="<<gf<<")");

    const char* subfieldName = "gravitational_acceleration";
    const char* componentNames[3] = { "gravitational_acceleration_x", "gravitational_acceleration_y", "gravitational_acceleration_z" };

    const PylithReal accelerationScale = _normalizer->getLengthScale() / (_normalizer->getTimeScale() * _normalizer->getTimeScale());

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
    this->setSubfieldQuery(subfieldName, gf);

    PYLITH_METHOD_END;
} // addGravityField


// ---------------------------------------------------------------------------------------------------------------------
// Add source density subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addSourceDensity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addSourceDensity(void)");

    const char* subfieldName = "source_density";
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0 / timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSourceDensity


// ---------------------------------------------------------------------------------------------------------------------
// Add heat source subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addHeatSource(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addHeatSource(void)");

    const char* subfieldName = "heat_source";
    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal pressureScale = _normalizer->getPressureScale();
    // Heat source has units of power per volume: W/m^3 = kg/(m*s^3) = pressure / time
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


// ============================= Rheology fields =============================

// ---------------------------------------------------------------------------------------------------------------------
// Add Biot coefficient subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addBiotCoefficient(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addBiotCoefficient(void)");

    const char* subfieldName = "biot_coefficient";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBiotCoefficient


// ---------------------------------------------------------------------------------------------------------------------
// Add Biot modulus subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addBiotModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addBiotModulus(void)");

    const char* subfieldName = "biot_modulus";
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
} // addBiotModulus


// ---------------------------------------------------------------------------------------------------------------------
// Add drained bulk modulus subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addDrainedBulkModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addDrainedBulkModulus(void)");

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


// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addShearModulus(void)");

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


// ---------------------------------------------------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addIsotropicPermeability(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addIsotropicPermeability(void)");

    const char* subfieldName = "isotropic_permeability";
    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal permeabilityScale = lengthScale * lengthScale; // m^2

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = permeabilityScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addIsotropicPermeability


// ---------------------------------------------------------------------------------------------------------------------
// Add reference temperature subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addReferenceTemperature(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addReferenceTemperature(void)");

    const char* subfieldName = "reference_temperature";
    const PylithReal temperatureScale = _normalizer->getTemperatureScale();

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
// Add thermal expansion coefficient subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addThermalExpansionCoeff(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addThermalExpansionCoeff(void)");

    const char* subfieldName = "thermal_expansion_coefficient";
    const PylithReal temperatureScale = _normalizer->getTemperatureScale();
    const PylithReal thermalExpansionScale = 1.0 / temperatureScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = thermalExpansionScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThermalExpansionCoeff


// ---------------------------------------------------------------------------------------------------------------------
// Add fluid thermal expansion subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addFluidThermalExpansion(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addFluidThermalExpansion(void)");

    const char* subfieldName = "fluid_thermal_expansion";
    const PylithReal temperatureScale = _normalizer->getTemperatureScale();
    const PylithReal fluidThermalExpScale = 1.0 / temperatureScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = fluidThermalExpScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidThermalExpansion


// ---------------------------------------------------------------------------------------------------------------------
// Add thermal conductivity subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addThermalConductivity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addThermalConductivity(void)");

    const char* subfieldName = "thermal_conductivity";
    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal pressureScale = _normalizer->getPressureScale();
    const PylithReal temperatureScale = _normalizer->getTemperatureScale();
    // Thermal conductivity: W/(m*K) = kg*m/(s^3*K) = pressure * length / (time * temperature)
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
// Add specific heat subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addSpecificHeat(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addSpecificHeat(void)");

    const char* subfieldName = "specific_heat";
    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal timeScale = _normalizer->getTimeScale();
    const PylithReal temperatureScale = _normalizer->getTemperatureScale();
    // Specific heat: J/(kg*K) = m^2/(s^2*K) = length^2 / (time^2 * temperature)
    const PylithReal specificHeatScale = (lengthScale * lengthScale) / (timeScale * timeScale * temperatureScale);

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
// Add reference stress subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addReferenceStress(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addReferenceStress(void)");

    const char* subfieldName = "reference_stress";
    const char* componentNames[6] = {
        "reference_stress_xx",
        "reference_stress_yy",
        "reference_stress_zz",
        "reference_stress_xy",
        "reference_stress_yz",
        "reference_stress_xz"
    };
    const int stressSize = (_spaceDim == 2) ? 4 : 6;
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


// ---------------------------------------------------------------------------------------------------------------------
// Add reference strain subfield to auxiliary subfields.
void
pylith::materials::AuxiliaryFactoryThermoporoelasticity::addReferenceStrain(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addReferenceStrain(void)");

    const char* subfieldName = "reference_strain";
    const char* componentNames[6] = {
        "reference_strain_xx",
        "reference_strain_yy",
        "reference_strain_zz",
        "reference_strain_xy",
        "reference_strain_yz",
        "reference_strain_xz"
    };
    const int strainSize = (_spaceDim == 2) ? 4 : 6;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
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


// End of file
