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

#include "UniformThermoelasticity2D.hh" // Implementation of class methods

#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory

namespace pylith {
    class _UniformThermoelasticity2D;
}

// ------------------------------------------------------------------------------------------------
class pylith::_UniformThermoelasticity2D {
public:

    // Spatial database user functions for auxiliary fields
    static double density(const double x, const double y) {
        return UniformThermoelasticity2D::density(x, y);
    }

    static double vs(const double x, const double y) {
        return UniformThermoelasticity2D::vs(x, y);
    }

    static double vp(const double x, const double y) {
        return UniformThermoelasticity2D::vp(x, y);
    }

    static double specific_heat(const double x, const double y) {
        return UniformThermoelasticity2D::specific_heat(x, y);
    }

    static double thermal_conductivity(const double x, const double y) {
        return UniformThermoelasticity2D::thermal_conductivity(x, y);
    }

    static double reference_temperature(const double x, const double y) {
        return UniformThermoelasticity2D::reference_temperature(x, y);
    }

    static double thermal_expansion_coefficient(const double x, const double y) {
        return UniformThermoelasticity2D::thermal_expansion_coefficient(x, y);
    }

    static const char* density_units(void) {
        return UniformThermoelasticity2D::density_units();
    }

    static const char* vs_units(void) {
        return UniformThermoelasticity2D::vs_units();
    }

    static const char* vp_units(void) {
        return UniformThermoelasticity2D::vp_units();
    }

    static const char* specific_heat_units(void) {
        return UniformThermoelasticity2D::specific_heat_units();
    }

    static const char* thermal_conductivity_units(void) {
        return UniformThermoelasticity2D::thermal_conductivity_units();
    }

    static const char* reference_temperature_units(void) {
        return UniformThermoelasticity2D::reference_temperature_units();
    }

    static const char* thermal_expansion_coefficient_units(void) {
        return UniformThermoelasticity2D::thermal_expansion_coefficient_units();
    }

}; // _UniformThermoelasticity2D

// ------------------------------------------------------------------------------------------------
// Static coefficient values
const double pylith::UniformThermoelasticity2D::LENGTHSCALE = 1.0e+3;
const double pylith::UniformThermoelasticity2D::TIMESCALE = 1.0e+3;
const double pylith::UniformThermoelasticity2D::PRESSURESCALE = 1.0e+9;
const double pylith::UniformThermoelasticity2D::TEMPERATURESCALE = 1.0e+3;

const double pylith::UniformThermoelasticity2D::DENSITY = 2500.0;
const double pylith::UniformThermoelasticity2D::VS = 3000.0;
const double pylith::UniformThermoelasticity2D::VP = 5196.0;
const double pylith::UniformThermoelasticity2D::SPECIFIC_HEAT = 1000.0;
const double pylith::UniformThermoelasticity2D::THERMAL_CONDUCTIVITY = 3.0;
const double pylith::UniformThermoelasticity2D::REFERENCE_TEMPERATURE = 300.0;
const double pylith::UniformThermoelasticity2D::THERMAL_EXPANSION_COEFF = 1.0e-5;

const double pylith::UniformThermoelasticity2D::DISP_GRADIENT = 1.0e-4;
const double pylith::UniformThermoelasticity2D::TEMPERATURE = 300.0; // At reference temp, no thermal strain

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
pylith::UniformThermoelasticity2D::createData(void) {
    TestThermoelasticity_Data* data = new TestThermoelasticity_Data();

    data->journalName = "UniformThermoelasticity2D";
    data->spaceDim = 2;
    data->boundaryLabel = "boundary";
    data->useAsciiMesh = true;

    // Scales
    data->scales.setLengthScale(LENGTHSCALE);
    data->scales.setTimeScale(TIMESCALE);
    data->scales.setPressureScale(PRESSURESCALE);
    data->scales.setTemperatureScale(TEMPERATURESCALE);
    data->scales.computeDensityScale();

    // Test parameters
    data->t = 0.0;
    data->dt = 0.05;
    data->tolerance = 1.0e-9;
    data->isJacobianLinear = true;
    data->allowZeroResidual = true;
    data->jacobianConvergenceRate = 1.0;
    data->formulation = pylith::problems::Physics::QUASISTATIC;

    // Solution discretizations: displacement (P1), temperature (P1)
    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // displacement
        pylith::topology::Field::Discretization(1, 1), // temperature
    };
    data->numSolnSubfields = 2;
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    // Solution functions
    static pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
        pylith::UniformThermoelasticity2D::solnkernel_disp,
        pylith::UniformThermoelasticity2D::solnkernel_temp,
    };
    data->exactSolnFns = _exactSolnFns;

    static pylith::testing::MMSTest::solution_fn _exactSolnDotFns[2] = {
        pylith::UniformThermoelasticity2D::solnkernel_disp_dot,
        pylith::UniformThermoelasticity2D::solnkernel_temp_dot,
    };
    data->exactSolnDotFns = _exactSolnDotFns;

    // Auxiliary fields
    static const char* _auxSubfields[7] = {
        "density",
        "specific_heat",
        "thermal_conductivity",
        "reference_temperature",
        "thermal_expansion_coefficient",
        "shear_modulus",
        "bulk_modulus",
    };
    data->numAuxSubfields = 7;
    data->auxSubfields = _auxSubfields;

    static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
        pylith::topology::Field::Discretization(0, 1), // density
        pylith::topology::Field::Discretization(0, 1), // specific_heat
        pylith::topology::Field::Discretization(0, 1), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 1), // reference_temperature
        pylith::topology::Field::Discretization(0, 1), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 1), // shear_modulus
        pylith::topology::Field::Discretization(0, 1), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    // Spatial database for auxiliary fields
    data->auxDB.addValue("density", _UniformThermoelasticity2D::density, _UniformThermoelasticity2D::density_units());
    data->auxDB.addValue("specific_heat", _UniformThermoelasticity2D::specific_heat, _UniformThermoelasticity2D::specific_heat_units());
    data->auxDB.addValue("thermal_conductivity", _UniformThermoelasticity2D::thermal_conductivity, _UniformThermoelasticity2D::thermal_conductivity_units());
    data->auxDB.addValue("reference_temperature", _UniformThermoelasticity2D::reference_temperature, _UniformThermoelasticity2D::reference_temperature_units());
    data->auxDB.addValue("thermal_expansion_coefficient", _UniformThermoelasticity2D::thermal_expansion_coefficient, _UniformThermoelasticity2D::thermal_expansion_coefficient_units());
    data->auxDB.addValue("vs", _UniformThermoelasticity2D::vs, _UniformThermoelasticity2D::vs_units());
    data->auxDB.addValue("vp", _UniformThermoelasticity2D::vp, _UniformThermoelasticity2D::vp_units());

    return data;
} // createData


// End of file
