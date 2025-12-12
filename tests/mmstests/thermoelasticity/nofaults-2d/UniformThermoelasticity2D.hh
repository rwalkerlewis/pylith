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

#include "TestThermoelasticity.hh"

namespace pylith {
    class UniformThermoelasticity2D;
}

// ------------------------------------------------------------------------------------------------
class pylith::UniformThermoelasticity2D {
    // Uniform thermoelasticity test case with:
    // - Linear displacement field
    // - Constant temperature field (no thermal strain)
    // - Zero body force (residual = 0 for linear solution)

    // Coefficient for nondimensionalization
    static const double LENGTHSCALE;
    static const double TIMESCALE;
    static const double PRESSURESCALE;
    static const double TEMPERATURESCALE;

    // Auxiliary field values
    static const double DENSITY;
    static const double VS;
    static const double VP;
    static const double SPECIFIC_HEAT;
    static const double THERMAL_CONDUCTIVITY;
    static const double REFERENCE_TEMPERATURE;
    static const double THERMAL_EXPANSION_COEFF;

    // Solution field parameters
    static const double DISP_GRADIENT; // Linear displacement gradient
    static const double TEMPERATURE;   // Constant temperature

protected:

    // Density
    static double density(const double x, const double y) {
        return DENSITY;
    }

    static const char* density_units(void) {
        return "kg/m**3";
    }

    // Vs
    static double vs(const double x, const double y) {
        return VS;
    }

    static const char* vs_units(void) {
        return "m/s";
    }

    // Vp
    static double vp(const double x, const double y) {
        return VP;
    }

    static const char* vp_units(void) {
        return "m/s";
    }

    // Specific heat
    static double specific_heat(const double x, const double y) {
        return SPECIFIC_HEAT;
    }

    static const char* specific_heat_units(void) {
        return "J/(kg*K)";
    }

    // Thermal conductivity
    static double thermal_conductivity(const double x, const double y) {
        return THERMAL_CONDUCTIVITY;
    }

    static const char* thermal_conductivity_units(void) {
        return "W/(m*K)";
    }

    // Reference temperature
    static double reference_temperature(const double x, const double y) {
        return REFERENCE_TEMPERATURE;
    }

    static const char* reference_temperature_units(void) {
        return "K";
    }

    // Thermal expansion coefficient
    static double thermal_expansion_coefficient(const double x, const double y) {
        return THERMAL_EXPANSION_COEFF;
    }

    static const char* thermal_expansion_coefficient_units(void) {
        return "1/K";
    }

    // Solution: Displacement (linear field, no body force needed)
    static PetscErrorCode solnkernel_disp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* f,
                                          void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(x);
        assert(f);

        // Linear displacement field: u_i = a * x_i
        f[0] = DISP_GRADIENT * x[0];
        f[1] = DISP_GRADIENT * x[1];

        return PETSC_SUCCESS;
    }

    // Solution: Temperature (constant)
    static PetscErrorCode solnkernel_temp(PetscInt spaceDim,
                                          PetscReal t,
                                          const PetscReal x[],
                                          PetscInt numComponents,
                                          PetscScalar* f,
                                          void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(x);
        assert(f);

        // Constant temperature field (at reference, so no thermal strain)
        f[0] = TEMPERATURE;

        return PETSC_SUCCESS;
    }

    // Solution: Displacement time derivative (zero for steady-state)
    static PetscErrorCode solnkernel_disp_dot(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar* f,
                                              void* context) {
        assert(2 == spaceDim);
        assert(2 == numComponents);
        assert(f);

        f[0] = 0.0;
        f[1] = 0.0;

        return PETSC_SUCCESS;
    }

    // Solution: Temperature time derivative (zero for steady-state)
    static PetscErrorCode solnkernel_temp_dot(PetscInt spaceDim,
                                              PetscReal t,
                                              const PetscReal x[],
                                              PetscInt numComponents,
                                              PetscScalar* f,
                                              void* context) {
        assert(2 == spaceDim);
        assert(1 == numComponents);
        assert(f);

        f[0] = 0.0;

        return PETSC_SUCCESS;
    }

public:

    static TestThermoelasticity_Data* createData(void);

}; // class UniformThermoelasticity2D

// End of file
