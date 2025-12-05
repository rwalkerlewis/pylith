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

#include "pylith/materials/RheologyThermoporoelasticity.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT*

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::RheologyThermoporoelasticity::RheologyThermoporoelasticity(void) {
    PyreComponent::setName("rheologythermoporoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::RheologyThermoporoelasticity::~RheologyThermoporoelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::RheologyThermoporoelasticity::deallocate(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::RheologyThermoporoelasticity::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                       const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    // Default implementation does nothing.
    PYLITH_METHOD_END;
} // updateKernelConstants


// End of file
