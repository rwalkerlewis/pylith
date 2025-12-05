// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/testing/MMSTest.hh" // ISA MMSTest

namespace pylith {
    class TestThermoporoelasticity;
    class TestThermoporoelasticity_Data;
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::TestThermoporoelasticity : public pylith::testing::MMSTest {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Setup testing data.
     */
    virtual void setUp(void);

    /// Deallocate testing data.
    virtual void tearDown(void);

    /// Set exact solution in domain.
    void setExactSolution(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestThermoporoelasticity_Data* _data; ///< Data for testing.

}; // class TestThermoporoelasticity

// ------------------------------------------------------------------------------------------------
class pylith::TestThermoporoelasticity_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestThermoporoelasticity_Data(void);

    /// Destructor
    ~TestThermoporoelasticity_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    const char* meshFilename; ///< Name of file with ASCII mesh.

    // Scales
    PylithReal lengthScale; ///< Length scale.
    PylithReal pressureScale; ///< Pressure scale.
    PylithReal timeScale; ///< Time scale.
    PylithReal densityScale; ///< Density scale.
    PylithReal temperatureScale; ///< Temperature scale.

    // Material properties
    PylithReal solidDensity; ///< Solid density.
    PylithReal fluidDensity; ///< Fluid density.
    PylithReal fluidViscosity; ///< Fluid viscosity.
    PylithReal porosity; ///< Porosity.
    PylithReal biotCoefficient; ///< Biot coefficient.
    PylithReal biotModulus; ///< Biot modulus.
    PylithReal drainedBulkModulus; ///< Drained bulk modulus.
    PylithReal shearModulus; ///< Shear modulus.
    PylithReal isotropicPermeability; ///< Isotropic permeability.
    PylithReal referenceTemperature; ///< Reference temperature.
    PylithReal thermalExpansionCoeff; ///< Solid thermal expansion coefficient.
    PylithReal fluidThermalExpansion; ///< Fluid thermal expansion.
    PylithReal thermalConductivity; ///< Thermal conductivity.
    PylithReal specificHeat; ///< Specific heat.

    // Boundary conditions
    const char* bcLabel; ///< Label for boundary conditions.
    PylithInt bcLabelId; ///< Label ID for boundary conditions.

    // Exact solution functions
    PetscPointFunc solnExactDisp; ///< Displacement exact solution.
    PetscPointFunc solnExactPres; ///< Pressure exact solution.
    PetscPointFunc solnExactTemp; ///< Temperature exact solution.

    // Auxiliary fields
    PetscPointFunc bodyForceFn; ///< Body force function (if any).
    PetscPointFunc sourceDensityFn; ///< Source density function (if any).
    PetscPointFunc heatSourceFn; ///< Heat source function (if any).

}; // class TestThermoporoelasticity_Data

// End of file
