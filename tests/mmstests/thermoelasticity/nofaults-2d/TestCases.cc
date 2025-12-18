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

#include "TestThermoelasticity.hh" // USES TestThermoelasticity_Data
#include "UniformThermoelasticity2D.hh" // USES UniformThermoelasticity2D

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
createTriP1(void) {
    pylith::TestThermoelasticity_Data* data = pylith::UniformThermoelasticity2D::createData();
    data->meshFilename = "data/tri.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // displacement
        pylith::topology::Field::Discretization(1, 1), // temperature
    };
    data->numSolnSubfields = 2;
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
createTriP2(void) {
    pylith::TestThermoelasticity_Data* data = pylith::UniformThermoelasticity2D::createData();
    data->meshFilename = "data/tri.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // specific_heat
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(2, 2), // temperature
    };
    data->numSolnSubfields = 2;
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
createQuadQ1(void) {
    pylith::TestThermoelasticity_Data* data = pylith::UniformThermoelasticity2D::createData();
    data->meshFilename = "data/quad.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // displacement
        pylith::topology::Field::Discretization(1, 1), // temperature
    };
    data->numSolnSubfields = 2;
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
createQuadQ2(void) {
    pylith::TestThermoelasticity_Data* data = pylith::UniformThermoelasticity2D::createData();
    data->meshFilename = "data/quad.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _auxDiscretizations[7] = {
        pylith::topology::Field::Discretization(0, 2), // density
        pylith::topology::Field::Discretization(0, 2), // specific_heat
        pylith::topology::Field::Discretization(0, 2), // thermal_conductivity
        pylith::topology::Field::Discretization(0, 2), // reference_temperature
        pylith::topology::Field::Discretization(0, 2), // thermal_expansion_coefficient
        pylith::topology::Field::Discretization(0, 2), // shear_modulus
        pylith::topology::Field::Discretization(0, 2), // bulk_modulus
    };
    data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_auxDiscretizations);

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(2, 2), // temperature
    };
    data->numSolnSubfields = 2;
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
// Test cases
// ------------------------------------------------------------------------------------------------

// TriP1
TEST_CASE("UniformThermoelasticity2D::TriP1::testDiscretization", "[UniformThermoelasticity2D][TriP1][testDiscretization]") {
    pylith::TestThermoelasticity(createTriP1()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::TriP1::testResidual", "[UniformThermoelasticity2D][TriP1][testResidual]") {
    pylith::TestThermoelasticity(createTriP1()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::TriP1::testJacobianTaylorSeries", "[UniformThermoelasticity2D][TriP1][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(createTriP1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::TriP1::testJacobianFiniteDiff", "[UniformThermoelasticity2D][TriP1][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(createTriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("UniformThermoelasticity2D::TriP2::testDiscretization", "[UniformThermoelasticity2D][TriP2][testDiscretization]") {
    pylith::TestThermoelasticity(createTriP2()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::TriP2::testResidual", "[UniformThermoelasticity2D][TriP2][testResidual]") {
    pylith::TestThermoelasticity(createTriP2()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::TriP2::testJacobianTaylorSeries", "[UniformThermoelasticity2D][TriP2][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(createTriP2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::TriP2::testJacobianFiniteDiff", "[UniformThermoelasticity2D][TriP2][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(createTriP2()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testDiscretization", "[UniformThermoelasticity2D][QuadQ1][testDiscretization]") {
    pylith::TestThermoelasticity(createQuadQ1()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testResidual", "[UniformThermoelasticity2D][QuadQ1][testResidual]") {
    pylith::TestThermoelasticity(createQuadQ1()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testJacobianTaylorSeries", "[UniformThermoelasticity2D][QuadQ1][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(createQuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testJacobianFiniteDiff", "[UniformThermoelasticity2D][QuadQ1][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(createQuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testDiscretization", "[UniformThermoelasticity2D][QuadQ2][testDiscretization]") {
    pylith::TestThermoelasticity(createQuadQ2()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testResidual", "[UniformThermoelasticity2D][QuadQ2][testResidual]") {
    pylith::TestThermoelasticity(createQuadQ2()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testJacobianTaylorSeries", "[UniformThermoelasticity2D][QuadQ2][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(createQuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testJacobianFiniteDiff", "[UniformThermoelasticity2D][QuadQ2][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(createQuadQ2()).testJacobianFiniteDiff();
}

// End of file
