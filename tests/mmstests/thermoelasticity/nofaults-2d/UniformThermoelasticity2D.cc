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
#include "pylith/bc/DirichletUserFn.hh" // USES DirichletUserFn

namespace pylith {
    class _UniformThermoelasticity2D;
}

// ------------------------------------------------------------------------------------------------
class pylith::_UniformThermoelasticity2D {
public:

    // Spatial database user functions for auxiliary fields
    static double density(const double x, const double y) {
        return UniformThermoelasticity2D::DENSITY;
    }

    static double vs(const double x, const double y) {
        return UniformThermoelasticity2D::VS;
    }

    static double vp(const double x, const double y) {
        return UniformThermoelasticity2D::VP;
    }

    static double specific_heat(const double x, const double y) {
        return UniformThermoelasticity2D::SPECIFIC_HEAT;
    }

    static double thermal_conductivity(const double x, const double y) {
        return UniformThermoelasticity2D::THERMAL_CONDUCTIVITY;
    }

    static double reference_temperature(const double x, const double y) {
        return UniformThermoelasticity2D::REFERENCE_TEMPERATURE;
    }

    static double thermal_expansion_coefficient(const double x, const double y) {
        return UniformThermoelasticity2D::THERMAL_EXPANSION_COEFF;
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

    // Test parameters
    data->t = 0.0;
    data->dt = 0.05;
    data->tolerance = 1.0e-4;
    data->isJacobianLinear = true;
    data->allowZeroResidual = true;
    data->jacobianConvergenceRate = 1.0;
    data->formulation = pylith::problems::Physics::QUASISTATIC;

    // Material settings
    data->material.setFormulation(pylith::problems::Physics::QUASISTATIC);
    data->material.useBodyForce(false);
    data->material.setIdentifier("thermoelasticity");
    data->material.setName("material-id=24");
    data->material.setLabelValue(24);

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
    data->auxDB.addValue("density", _UniformThermoelasticity2D::density, "kg/m**3");
    data->auxDB.addValue("specific_heat", _UniformThermoelasticity2D::specific_heat, "J/(kg*K)");
    data->auxDB.addValue("thermal_conductivity", _UniformThermoelasticity2D::thermal_conductivity, "watt/(meter*kelvin)");
    data->auxDB.addValue("reference_temperature", _UniformThermoelasticity2D::reference_temperature, "K");
    data->auxDB.addValue("thermal_expansion_coefficient", _UniformThermoelasticity2D::thermal_expansion_coefficient, "1/K");
    data->auxDB.addValue("vs", _UniformThermoelasticity2D::vs, "m/s");
    data->auxDB.addValue("vp", _UniformThermoelasticity2D::vp, "m/s");
    data->auxDB.setCoordSys(data->cs);

    // Boundary conditions (Dirichlet on all boundaries for displacement and temperature)
    static const PylithInt constrainedDispDOF[2] = {0, 1};
    static const PylithInt numConstrainedDisp = 2;
    static const PylithInt constrainedTempDOF[1] = {0};
    static const PylithInt numConstrainedTemp = 1;

    data->bcs.resize(2);

    pylith::bc::DirichletUserFn* bcDisp = new pylith::bc::DirichletUserFn();
    bcDisp->setSubfieldName("displacement");
    bcDisp->setLabelName("boundary");
    bcDisp->setLabelValue(1);
    bcDisp->setConstrainedDOF(constrainedDispDOF, numConstrainedDisp);
    bcDisp->setUserFn(solnkernel_disp);
    data->bcs[0] = bcDisp;

    pylith::bc::DirichletUserFn* bcTemp = new pylith::bc::DirichletUserFn();
    bcTemp->setSubfieldName("temperature");
    bcTemp->setLabelName("boundary");
    bcTemp->setLabelValue(1);
    bcTemp->setConstrainedDOF(constrainedTempDOF, numConstrainedTemp);
    bcTemp->setUserFn(solnkernel_temp);
    data->bcs[1] = bcTemp;

    // Solution functions
    static pylith::testing::MMSTest::solution_fn _exactSolnFns[2] = {
        pylith::UniformThermoelasticity2D::solnkernel_disp,
        pylith::UniformThermoelasticity2D::solnkernel_temp,
    };
    data->exactSolnFns = _exactSolnFns;
    data->exactSolnDotFns = NULL;  // Quasistatic, no time derivatives

    return data;
} // createData


// End of file
