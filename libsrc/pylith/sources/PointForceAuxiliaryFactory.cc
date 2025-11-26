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

#include "pylith/sources/PointForceAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "pylith/scales/ElasticityScales.hh" // USES ElasticityScales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::sources::PointForceAuxiliaryFactory::PointForceAuxiliaryFactory(void) {
    GenericComponent::setName("pointforceauxiliaryfactory");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::sources::PointForceAuxiliaryFactory::~PointForceAuxiliaryFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Add source location subfield to auxiliary subfields.
void
pylith::sources::PointForceAuxiliaryFactory::addSourceLocation(const PylithReal* location,
                                                               const int size) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSourceLocation(location="<<location<<", size="<<size<<")");

    assert(location);
    assert(size > 0 && size <= 3);

    const char* subfieldName = "source_location";
    const char* componentNames[3] = {
        "source_location_x",
        "source_location_y",
        "source_location_z",
    };

    const PylithReal lengthScale = _scales->getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = size;
    description.componentNames.resize(size);
    for (int i = 0; i < size; ++i) {
        description.componentNames[i] = componentNames[i];
    }
    description.scale = lengthScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));

    // Set values directly (uniform field)
    // Note: In practice, this would be done via spatial database
    // For simplicity, we use constant values
    this->setSubfieldQueryFn(subfieldName,
                             pylith::topology::FieldQuery::dbQueryGeneric,
                             location);

    PYLITH_METHOD_END;
} // addSourceLocation


// ------------------------------------------------------------------------------------------------
// Add moment tensor subfield to auxiliary subfields.
void
pylith::sources::PointForceAuxiliaryFactory::addMomentTensor(const PylithReal* momentTensor,
                                                             const int size) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addMomentTensor(momentTensor="<<momentTensor<<", size="<<size<<")");

    assert(momentTensor);
    assert(size == 3 || size == 4 || size == 6);

    const char* subfieldName = "moment_tensor";
    const char* componentNames2D[3] = {
        "moment_tensor_xx",
        "moment_tensor_yy",
        "moment_tensor_xy",
    };
    const char* componentNames3D[6] = {
        "moment_tensor_xx",
        "moment_tensor_yy",
        "moment_tensor_zz",
        "moment_tensor_xy",
        "moment_tensor_xz",
        "moment_tensor_yz",
    };

    // Moment tensor is dimensionless (normalized direction)
    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (size == 6) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = size;
    description.componentNames.resize(size);
    const char** names = (size <= 3) ? componentNames2D : componentNames3D;
    for (int i = 0; i < size; ++i) {
        description.componentNames[i] = names[i];
    }
    description.scale = 1.0;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQueryFn(subfieldName,
                             pylith::topology::FieldQuery::dbQueryGeneric,
                             momentTensor);

    PYLITH_METHOD_END;
} // addMomentTensor


// ------------------------------------------------------------------------------------------------
// Add magnitude subfield to auxiliary subfields.
void
pylith::sources::PointForceAuxiliaryFactory::addMagnitude(const PylithReal value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addMagnitude(value="<<value<<")");

    const char* subfieldName = "magnitude";

    // Magnitude has units of [Force * Length] = [Rigidity * Length^3]
    const PylithReal rigidityScale = _scales->getRigidityScale();
    const PylithReal lengthScale = _scales->getLengthScale();
    const PylithReal momentScale = rigidityScale * lengthScale * lengthScale * lengthScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = momentScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQueryFn(subfieldName,
                             pylith::topology::FieldQuery::dbQueryGeneric,
                             &value);

    PYLITH_METHOD_END;
} // addMagnitude


// ------------------------------------------------------------------------------------------------
// Add origin time subfield to auxiliary subfields.
void
pylith::sources::PointForceAuxiliaryFactory::addOriginTime(const PylithReal value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addOriginTime(value="<<value<<")");

    const char* subfieldName = "origin_time";
    const PylithReal timeScale = _scales->getTimeScale();

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
    this->setSubfieldQueryFn(subfieldName,
                             pylith::topology::FieldQuery::dbQueryGeneric,
                             &value);

    PYLITH_METHOD_END;
} // addOriginTime


// ------------------------------------------------------------------------------------------------
// Add dominant frequency subfield to auxiliary subfields.
void
pylith::sources::PointForceAuxiliaryFactory::addDominantFrequency(const PylithReal value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addDominantFrequency(value="<<value<<")");

    const char* subfieldName = "dominant_frequency";
    const PylithReal timeScale = _scales->getTimeScale();
    const PylithReal frequencyScale = 1.0 / timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = frequencyScale;
    description.validator = pylith::topology::FieldQuery::validatorPositive;

    _field->subfieldAdd(description, getSubfieldDiscretization(subfieldName));
    this->setSubfieldQueryFn(subfieldName,
                             pylith::topology::FieldQuery::dbQueryGeneric,
                             &value);

    PYLITH_METHOD_END;
} // addDominantFrequency


// ------------------------------------------------------------------------------------------------
// Add time delay subfield to auxiliary subfields.
void
pylith::sources::PointForceAuxiliaryFactory::addTimeDelay(const PylithReal value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeDelay(value="<<value<<")");

    const char* subfieldName = "time_delay";
    const PylithReal timeScale = _scales->getTimeScale();

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
    this->setSubfieldQueryFn(subfieldName,
                             pylith::topology::FieldQuery::dbQueryGeneric,
                             &value);

    PYLITH_METHOD_END;
} // addTimeDelay


// End of file
