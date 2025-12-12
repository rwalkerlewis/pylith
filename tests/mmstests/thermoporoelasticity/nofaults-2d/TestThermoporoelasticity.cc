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

#include "TestThermoporoelasticity.hh" // Implementation of class methods

#include "pylith/materials/Thermoporoelasticity.hh" // USES Thermoporoelasticity
#include "pylith/materials/IsotropicLinearThermoporoelasticity.hh" // USES IsotropicLinearThermoporoelasticity

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/problems/TimeDependent.hh" // USES TimeDependent

#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::TestThermoporoelasticity_Data::TestThermoporoelasticity_Data(void) :
    meshFilename(NULL),
    lengthScale(1.0e+3),
    pressureScale(2.25e+10),
    timeScale(2.0),
    densityScale(3.0e+3),
    temperatureScale(1.0),
    solidDensity(2500.0),
    fluidDensity(1000.0),
    fluidViscosity(1.0e-3),
    porosity(0.1),
    biotCoefficient(0.8),
    biotModulus(1.0e+10),
    drainedBulkModulus(5.0e+9),
    shearModulus(3.0e+9),
    isotropicPermeability(1.0e-15),
    referenceTemperature(293.0),
    thermalExpansionCoeff(1.0e-5),
    fluidThermalExpansion(2.1e-4),
    thermalConductivity(3.0),
    specificHeat(800.0),
    bcLabel("boundary"),
    bcLabelId(1),
    solnExactDisp(NULL),
    solnExactPres(NULL),
    solnExactTemp(NULL),
    bodyForceFn(NULL),
    sourceDensityFn(NULL),
    heatSourceFn(NULL) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::TestThermoporoelasticity_Data::~TestThermoporoelasticity_Data(void) {}


// ------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::TestThermoporoelasticity::setUp(void) {
    MMSTest::setUp();

    _data = NULL;
} // setUp


// ------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::TestThermoporoelasticity::tearDown(void) {
    delete _data;_data = NULL;

    MMSTest::tearDown();
} // tearDown


// ------------------------------------------------------------------------------------------------
// Set exact solution in domain.
void
pylith::TestThermoporoelasticity::setExactSolution(void) {
    // Placeholder for exact solution setup
    // Will be implemented in specific test cases
} // setExactSolution


// End of file
