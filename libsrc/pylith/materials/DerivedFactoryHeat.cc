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

#include "pylith/materials/DerivedFactoryHeat.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "pylith/scales/Scales.hh" // USES Scales

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::DerivedFactoryHeat::DerivedFactoryHeat(void) {
    GenericComponent::setName("derivedfactoryheat");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::DerivedFactoryHeat::~DerivedFactoryHeat(void) {}


// ------------------------------------------------------------------------------------------------
// Add heat flux subfield to derived fields.
void
pylith::materials::DerivedFactoryHeat::addHeatFlux(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addHeatFlux(void)");

    const char* fieldName = "heat_flux";
    const char* componentNames[3] = {
        "heat_flux_x",
        "heat_flux_y",
        "heat_flux_z"
    };

    // Heat flux scale: Power / Area = M / T^3
    // Using pressure_scale * length_scale / time_scale
    const PylithReal pressureScale = _scales->getPressureScale();
    const PylithReal lengthScale = _scales->getLengthScale();
    const PylithReal timeScale = _scales->getTimeScale();
    const PylithReal heatFluxScale = pressureScale * lengthScale / timeScale;

    pylith::topology::Field::Description description;
    description.label = fieldName;
    description.alias = fieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = _spaceDim;
    description.componentNames.resize(_spaceDim);
    for (int i = 0; i < _spaceDim; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = heatFluxScale;
    description.validator = NULL;

    _field->subfieldAdd(description, getSubfieldDiscretization(fieldName));

    PYLITH_METHOD_END;
} // addHeatFlux


// ------------------------------------------------------------------------------------------------
// Add subfields using discretizations provided.
void
pylith::materials::DerivedFactoryHeat::addSubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfields(void)");

    if (_subfieldDiscretizations.find("heat_flux") != _subfieldDiscretizations.end()) {
        addHeatFlux();
    } // if

    PYLITH_METHOD_END;
} // addSubfields


// End of file
