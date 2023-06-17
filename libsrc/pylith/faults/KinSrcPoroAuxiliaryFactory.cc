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

#include "KinSrcPoroAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcPoroAuxiliaryFactory::KinSrcPoroAuxiliaryFactory(void) {
    GenericComponent::setName("kinsrcporoauxiliaryfactory");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcPoroAuxiliaryFactory::~KinSrcPoroAuxiliaryFactory(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add slip initiation time (relative to origin time) subfield to auxiliary field.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addInitiationTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addInitiationTime(void)");

    const char* subfieldName = "initiation_time";
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addInitiationTime


// ---------------------------------------------------------------------------------------------------------------------
// Add riseTime subfield to auxiliary field.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addRiseTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRiseTime(void)");

    const char* subfieldName = "rise_time";
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addRiseTime


// ---------------------------------------------------------------------------------------------------------------------
// Add final slip subfield to auxiliary field.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addFinalSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFinalSlip(void)");

    const char* subfieldName = "slip";
    const char* componentNames[3] = { "final_slip_opening", "final_slip_left_lateral", "final_slip_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFinalSlip


// ---------------------------------------------------------------------------------------------------------------------
// Add  thickness subfield to auxiliary field.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addThickness(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addThickness(void)");

    const char* subfieldName = "thickness";

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = lengthScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    // No query; populated by kinematic fault at beginning of time step.

    // this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addThickness


// ---------------------------------------------------------------------
// Add layer porosity to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addPorosity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPorosity(void)");

    const char* subfieldName = "porosity";

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
    // No query; populated by kinematic fault at beginning of time step.
    // ** TO DO **
    // Is this step correct?
    // this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPorosity


// ---------------------------------------------------------------------
// Add beta p to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addBetaP(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBetaP(void)");

    const char* subfieldName = "beta_p";

    // ** TO DO **
    // This works? The scale should be 1/Pa
    const PylithReal betaScale = 1. / _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = betaScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    // No query; populated by kinematic fault at beginning of time step.
    // ** TO DO **
    // Is this step correct?
    // this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBetaP


// Add _beta_sigma to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addBetaSigma(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBetaSigma(void)");

    const char* subfieldName = "beta_sigma";

    const PylithReal betaScale = 1. / _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = betaScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    // No query; populated by kinematic fault at beginning of time step.
    // ** TO DO **
    // Is this step correct?
    // this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBetaSigma


// Add tangential permeability to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addPermeabilityTangential(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPermeabilityTangential(void)");

    const char* subfieldName = "permeability_tangential";

    // ** TO DO **
    // Please verify this following line
    // Permeability scale should be m^2
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
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    // No query; populated by kinematic fault at beginning of time step.
    // ** TO DO **
    // Is this step correct?
    // this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPermeabilityTangential


// Add  normal permeability to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addPermeabilityNormal(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addPermeabilityNormal(void)");

    const char* subfieldName = "permeability_normal";

    // ** TO DO **
    // Please verify this following line
    // Permeability scale should be m^2
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
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    // No query; populated by kinematic fault at beginning of time step.
    // ** TO DO **
    // Is this step correct?
    // this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addPermeabilityNormal


// ---------------------------------------------------------------------
// Add  fluid viscosity to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addFluidViscosity(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addFluidViscosity(void)");

    const char* subfieldName = "fluid_viscosity";

    // ** TO DO **
    // Please verify this following line
    // viscosity scale should be pa s
    const PylithReal fluidViscosityScale = _normalizer->getPressureScale()
                                           * _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = fluidViscosityScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    // No query; populated by kinematic fault at beginning of time step.
    // ** TO DO **
    // Is this step correct?
    // this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addFluidViscosity


// ---------------------------------------------------------------------
// Add negative side  bulk modulus to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addBulkModulusNegative(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulusNegative(void)");

    const char* subfieldName = "bulk_modulus_negative";

    const PylithReal bulkModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = bulkModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBulkModulusNegative


// ---------------------------------------------------------------------
// Add positive side bulk modulus to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addBulkModulusPositive(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBulkModulusPositive(void)");

    const char* subfieldName = "bulk_modulus_positive";

    const PylithReal bulkModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = bulkModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addBulkModulusPositive


// ---------------------------------------------------------------------
// Add negative side shear modulus to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addShearModulusNegative(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulusNegative(void)");

    const char* subfieldName = "shear_modulus_negative";

    const PylithReal shearModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = shearModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addShearModulusNegative


// ---------------------------------------------------------------------
// Add positive side shear modulus to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addShearModulusPositive(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addShearModulusPositive(void)");

    const char* subfieldName = "shear_modulus_positive";

    const PylithReal shearModulusScale = _normalizer->getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = shearModulusScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addShearModulusPositive


// ---------------------------------------------------------------------
// Add body force subfield to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addBodyForce(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addBodyForce(void)");

    const char* fieldName = "body_force";
    const char* componentNames[3] = { "body_force_x", "body_force_y", "body_force_z" };

    // ** TO DO **
    // Verify this line
    // The scale should be pa / m
    const PylithReal bodyForceScale = _normalizer->getPressureScale() / _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = bodyForceScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    // ** TO DO **
    // Is this step correct?
    // this->setSubfieldQuery(fieldName);

    PYLITH_METHOD_END;
} // addBodyForce


// ---------------------------------------------------------------------
// Add source to auxiliary fields.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addSource(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSource(void)");

    const char* subfieldName = "source";
    // ** TO DO **
    // Verify this line
    // The scale for source is 1 / s
    const PylithReal sourceScale = 1. / _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = sourceScale;
    description.validator = pylith::topology::FieldQuery::validatorNonnegative;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // ** TO DO **
    // Is this step correct?
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSource


// ---------------------------------------------------------------------------------------------------------------------
// Add slip rate subfield to auxiliary field.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addSlipRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlipRate(void)");

    const char* subfieldName = "slip_rate";
    const char* componentNames[3] = { "slip_rate_opening", "slip_rate_left_lateral", "slip_rate_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();
    const PylithReal timeScale = _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale / timeScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addSlipRate


// ---------------------------------------------------------------------------------------------------------------------
// Add time history value subfield to auxiliary field.
void
pylith::faults::KinSrcPoroAuxiliaryFactory::addTimeHistoryValue(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryValue(void)");

    const char* subfieldName = "time_history_value";

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization("_slip"));
    // No subfield query; populated at begining of time step.

    PYLITH_METHOD_END;
} // addTimeHistoryValue


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::faults::KinSrcPoroAuxiliaryFactory::updateTimeHistoryValue(pylith::topology::Field* auxiliaryField,
                                                                   const PylithReal t,
                                                                   const PylithReal timeScale,
                                                                   spatialdata::spatialdb::TimeHistory* const dbTimeHistory) {
    PYLITH_METHOD_BEGIN;
    pythia::journal::debug_t debug("kinsrcporoauxiliaryfactory");
    debug << pythia::journal::at(__HERE__)
          << "KinSrcPoroAuxiliaryFactory::updateTimeHistoryValue(auxiliaryField="<<auxiliaryField<<", t="<<t
          <<", timeScale="<<timeScale<<", dbTimeHistory="<<dbTimeHistory<<")"
          << pythia::journal::endl;

    assert(auxiliaryField);
    assert(dbTimeHistory);

    PetscErrorCode err = 0;

    PetscSection auxiliaryFieldSection = auxiliaryField->getLocalSection();assert(auxiliaryFieldSection);
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryFieldSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    pylith::topology::VecVisitorMesh auxiliaryFieldVisitor(*auxiliaryField);
    PetscScalar* auxiliaryFieldArray = auxiliaryFieldVisitor.localArray();assert(auxiliaryFieldArray);

    // Compute offset of time history subfields in auxiliary field.
    const PetscInt i_startTime = auxiliaryField->getSubfieldInfo("initiation_time").index;
    const PetscInt i_value = auxiliaryField->getSubfieldInfo("time_history_value").index;

    // Loop over all points in section.
    for (PetscInt p = pStart; p < pEnd; ++p) {
        // Skip points without values in section.
        if (!auxiliaryFieldVisitor.sectionDof(p)) {continue;}

        // Get starting time and compute relative time for point.
        const PetscInt offStartTime = auxiliaryFieldVisitor.sectionSubfieldOffset(i_startTime, p);
        const PylithScalar tStart = auxiliaryFieldArray[offStartTime];
        const PylithScalar tRel = t - tStart;

        // Query time history for value (normalized amplitude).
        PylithScalar value = 0.0;
        if (tRel >= 0.0) {
            PylithScalar tDim = tRel * timeScale;
            const int err = dbTimeHistory->query(&value, tDim);
            if (err) {
                std::ostringstream msg;
                msg << "Error querying for time '" << tDim << "' in time history database '" << dbTimeHistory->getDescription() << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } // if

        // Update value (normalized amplitude) in auxiliary field.
        const PetscInt offValue = auxiliaryFieldVisitor.sectionSubfieldOffset(i_value, p);
        auxiliaryFieldArray[offValue] = value;
    } // for

    PYLITH_METHOD_END;
} // updateAuxiliaryField


// End of file
