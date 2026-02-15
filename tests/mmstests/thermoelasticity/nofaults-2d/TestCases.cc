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
// Helper class for setting up test data with mesh configuration.
class pylith::_TestThermoelasticity {
public:

    static TestThermoelasticity_Data* TriP1(void);
    static TestThermoelasticity_Data* TriP2(void);
    static TestThermoelasticity_Data* QuadQ1(void);
    static TestThermoelasticity_Data* QuadQ2(void);

};

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
pylith::_TestThermoelasticity::TriP1(void) {
    TestThermoelasticity_Data* data = UniformThermoelasticity2D::createData();
    data->meshFilename = "data/tri.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // displacement
        pylith::topology::Field::Discretization(1, 1), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
pylith::_TestThermoelasticity::TriP2(void) {
    TestThermoelasticity_Data* data = UniformThermoelasticity2D::createData();
    data->meshFilename = "data/tri.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(2, 2), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
pylith::_TestThermoelasticity::QuadQ1(void) {
    TestThermoelasticity_Data* data = UniformThermoelasticity2D::createData();
    data->meshFilename = "data/quad.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(1, 1), // displacement
        pylith::topology::Field::Discretization(1, 1), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
pylith::TestThermoelasticity_Data*
pylith::_TestThermoelasticity::QuadQ2(void) {
    TestThermoelasticity_Data* data = UniformThermoelasticity2D::createData();
    data->meshFilename = "data/quad.mesh";
    data->useAsciiMesh = true;

    static const pylith::topology::Field::Discretization _solnDiscretizations[2] = {
        pylith::topology::Field::Discretization(2, 2), // displacement
        pylith::topology::Field::Discretization(2, 2), // temperature
    };
    data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(_solnDiscretizations);

    return data;
}

// ------------------------------------------------------------------------------------------------
// Test cases
// ------------------------------------------------------------------------------------------------

// TriP1
TEST_CASE("UniformThermoelasticity2D::TriP1::testDiscretization", "[UniformThermoelasticity2D][TriP1][testDiscretization]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP1()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::TriP1::testResidual", "[UniformThermoelasticity2D][TriP1][testResidual]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP1()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::TriP1::testJacobianTaylorSeries", "[UniformThermoelasticity2D][TriP1][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::TriP1::testJacobianFiniteDiff", "[UniformThermoelasticity2D][TriP1][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("UniformThermoelasticity2D::TriP2::testDiscretization", "[UniformThermoelasticity2D][TriP2][testDiscretization]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP2()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::TriP2::testResidual", "[UniformThermoelasticity2D][TriP2][testResidual]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP2()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::TriP2::testJacobianTaylorSeries", "[UniformThermoelasticity2D][TriP2][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::TriP2::testJacobianFiniteDiff", "[UniformThermoelasticity2D][TriP2][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::TriP2()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testDiscretization", "[UniformThermoelasticity2D][QuadQ1][testDiscretization]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ1()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testResidual", "[UniformThermoelasticity2D][QuadQ1][testResidual]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ1()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testJacobianTaylorSeries", "[UniformThermoelasticity2D][QuadQ1][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ1::testJacobianFiniteDiff", "[UniformThermoelasticity2D][QuadQ1][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testDiscretization", "[UniformThermoelasticity2D][QuadQ2][testDiscretization]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ2()).testDiscretization();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testResidual", "[UniformThermoelasticity2D][QuadQ2][testResidual]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ2()).testResidual();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testJacobianTaylorSeries", "[UniformThermoelasticity2D][QuadQ2][testJacobianTaylorSeries]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformThermoelasticity2D::QuadQ2::testJacobianFiniteDiff", "[UniformThermoelasticity2D][QuadQ2][testJacobianFiniteDiff]") {
    pylith::TestThermoelasticity(pylith::_TestThermoelasticity::QuadQ2()).testJacobianFiniteDiff();
}

// End of file
