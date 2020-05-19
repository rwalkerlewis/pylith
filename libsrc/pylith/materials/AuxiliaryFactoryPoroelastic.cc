// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "AuxiliaryFactoryPoroelastic.hh" // implementation of object methods

#include "Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::AuxiliaryFactoryPoroelastic::AuxiliaryFactoryPoroelastic(void) {
    GenericComponent::setName("AuxiliaryFactoryPoroelastic");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::AuxiliaryFactoryPoroelastic::~AuxiliaryFactoryPoroelastic(void) {}

//JS
// ----------------------------------------------------------------------------
// Add isotropic permeability subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addIsotropicPermeability(void)
{ // isotropicPermeablity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("isotropicPermeability(void)");

    const char* fieldName = "isotropic_permeability";

    const PylithReal lengthScale = _normalizer->lengthScale();
    const PylithReal permeabilityScale = lengthScale*lengthScale;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = permeabilityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addIsotropicPermeability

// --------------------------------------------------------------------
// Add fluid bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addFluidBulkModulus(void)
{ // fluidBulkModulus
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("fluidBulkModulus(void)");

    const char* fieldName = "fluid_bulk_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addFluidBulkModulus

// ---------------------------------------------------------------------
// Add biot coefficient subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addBiotCoefficient(void)
{ // biotCoefficient
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("biotCoefficient(void)");

    const char* fieldName = "biot_coefficient";

    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = noScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addBiotCoefficient

// ----------------------------------------------------------------------
// Add reference stress subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addReferenceStress(void)
{ // referenceStress
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("referenceStress(void)");

    const char* fieldName = "reference_stress";
    const char* componentNames[6] = { "reference_stress_xx", "reference_stress_yy", "reference_stress_zz", "reference_stress_xy", "reference_stress_yz", "reference_stress_xz" };
    const int stressSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = stressSize;
    description.componentNames.resize(stressSize);
    for (int i = 0; i < stressSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addReferenceStress


// ----------------------------------------------------------------------
// Add reference strain subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addReferenceStrain(void)
{ // addReferenceStrain
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("refrenceStrain(void)");

    const char* fieldName = "reference_strain";
    const char* componentNames[6] = { "reference_strain_xx", "reference_strain_yy", "reference_strain_zz", "reference_strain_xy", "reference_strain_yz", "reference_strain_xz" };
    const int strainSize = (3 == _spaceDim) ? 6 : (2 == _spaceDim) ? 4 : 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = (3 == _spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = strainSize;
    description.componentNames.resize(strainSize);
    for (int i = 0; i < strainSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addReferenceStrain

// ---------------------------------------------------------------------------------------------------------------------
// Add shear modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addShearModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulus(void)");

    const char* fieldName = "shear_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);
    //_setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryShearModulus);

    PYLITH_METHOD_END;
} // addShearModulus

// ---------------------------------------------------------------------------------------------------------------------
// Add solid bulk modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addSolidBulkModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSolidBulkModulus(void)");

    const char* fieldName = "solid_bulk_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);
    //_setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQuerySolidBulkModulus);

    PYLITH_METHOD_END;
} // addSolidBulkModulus

// ---------------------------------------------------------------------------------------------------------------------
// Add young's modulus subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelastic::addYoungsModulus(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addYoungsModulus(void)");

    const char* fieldName = "youngs_modulus";
    const PylithReal pressureScale = _normalizer->pressureScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);
    //_setSubfieldQueryFn(fieldName, pylith::materials::Query::dbQueryYoungsModulus);

    PYLITH_METHOD_END;
} // addYoungsModulus

// ----------------------------------------------------------------------
// Add poisson's ratio subfield to auxiliary fields.
void
pylith::materials::AuxiliaryFactoryPoroelasticity::addPoissonsRatio(void)
{ // porosity
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPoissonsRatio(void)");

    const char* fieldName = "poissons_ratio";

    const PylithReal noScale = 1;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = fieldName;
    description.scale = noScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addPoissonsRatio


// End of file
