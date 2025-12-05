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

#include "pylith/materials/RheologyThermoelasticity.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::RheologyThermoelasticity::RheologyThermoelasticity(void) {
    PyreComponent::setName("rheologythermoelasticity");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::RheologyThermoelasticity::~RheologyThermoelasticity(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::RheologyThermoelasticity::deallocate(void) {
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::RheologyThermoelasticity::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                   const PylithReal dt) const {
    // Default is to do nothing.
} // updateKernelConstants


// End of file
