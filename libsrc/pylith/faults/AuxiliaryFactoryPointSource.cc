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

#include "AuxiliaryFactoryPointSource.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::AuxiliaryFactoryPointSource::AuxiliaryFactoryPointSource(void) {
    GenericComponent::setName("auxiliaryfactorypointsource");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::AuxiliaryFactoryPointSource::~AuxiliaryFactoryPointSource(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Add fault strike direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryPointSource::addMomentTensor(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addMomentTensor(void)");

    const char* subfieldName = "moment_tensor";
    const char* componentNames[6] = {
        "m_rr",
        "m_tt",
        "m_pp",
        "m_rt",
        "m_rp",
        "m_tp"
        };
    const int tensorSize = 6;

    const PylithReal pressureScale = _normalizer->getPressureScale();
    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = tensorSize;
    description.componentNames.resize(tensorSize);
    for (int i = 0; i < tensorSize; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = pressureScale * lengthScale * lengthScale * lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
} // addMomentTensor

// ---------------------------------------------------------------------------------------------------------------------
// Add origin time for point sources.
void
pylith::faults::AuxiliaryFactoryKinematic::addOriginTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addOriginTime(void)");

    const char* subfieldName = "origin_time";
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

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    this->setSubfieldQuery(subfieldName);

    PYLITH_METHOD_END;
    } // addOriginTime





// ---------------------------------------------------------------------------------------------------------------------
// Add fault normal direction subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addNormalDir(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addNormalDir(void)");

    const char* fieldName = "normal_dir";
    const char* componentNames[3] = { "normal_dir_x", "normal_dir_y", "normal_dir_z" };

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; computed during initialization.

    PYLITH_METHOD_END;
} // addNormalDir


// ---------------------------------------------------------------------------------------------------------------------
// Add fault slip subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addSlip(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlip(void)");

    const char* fieldName = "slip";
    const char* componentNames[3] = { "slip_opening", "slip_left_lateral", "slip_reverse" };

    const PylithReal lengthScale = _normalizer->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlip


// ---------------------------------------------------------------------------------------------------------------------
// Add fault slip rate subfield to auxiliary fields.
void
pylith::faults::AuxiliaryFactoryKinematic::addSlipRate(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSlipRate(void)");

    const char* fieldName = "slip_rate";
    const char* componentNames[3] = { "slip_rate_opening", "slip_rate_left_lateral", "slip_rate_reverse" };

    const PylithReal velocityScale = _normalizer->getLengthScale() / _normalizer->getTimeScale();

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = velocityScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));
    // No query; populated by kinematic source at beginning of time step.

    PYLITH_METHOD_END;
} // addSlipRate


// End of file
