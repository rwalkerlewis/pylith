// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestFaultKinThermoelasticity
 */

#include "TestFaultKinThermoelasticity.hh" // USES TestFaultKinThermoelasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "TwoBlocksStaticThermal.hh"

// TriP1
TEST_CASE("TwoBlocksStaticThermal::TriP1::testDiscretization", "[TwoBlocksStaticThermal][TriP1][discretization]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP1()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticThermal::TriP1::testResidual", "[TwoBlocksStaticThermal][TriP1][residual]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP1()).testResidual();
}
TEST_CASE("TwoBlocksStaticThermal::TriP1::testJacobianTaylorSeries", "[TwoBlocksStaticThermal][TriP1][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticThermal::TriP1::testJacobianFiniteDiff", "[TwoBlocksStaticThermal][TriP1][Jacobian finite difference]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("TwoBlocksStaticThermal::TriP2::testDiscretization", "[TwoBlocksStaticThermal][TriP2][discretization]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP2()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticThermal::TriP2::testResidual", "[TwoBlocksStaticThermal][TriP2][residual]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP2()).testResidual();
}
TEST_CASE("TwoBlocksStaticThermal::TriP2::testJacobianTaylorSeries", "[TwoBlocksStaticThermal][TriP2][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticThermal::TriP2::testJacobianFiniteDiff", "[TwoBlocksStaticThermal][TriP2][Jacobian finite difference]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("TwoBlocksStaticThermal::TriP3::testDiscretization", "[TwoBlocksStaticThermal][TriP3][discretization]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP3()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticThermal::TriP3::testResidual", "[TwoBlocksStaticThermal][TriP3][residual]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP3()).testResidual();
}
TEST_CASE("TwoBlocksStaticThermal::TriP3::testJacobianTaylorSeries", "[TwoBlocksStaticThermal][TriP3][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticThermal::TriP3::testJacobianFiniteDiff", "[TwoBlocksStaticThermal][TriP3][Jacobian finite difference]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("TwoBlocksStaticThermal::QuadQ1::testDiscretization", "[TwoBlocksStaticThermal][QuadQ1][discretization]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ1()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ1::testResidual", "[TwoBlocksStaticThermal][QuadQ1][residual]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ1()).testResidual();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ1::testJacobianTaylorSeries", "[TwoBlocksStaticThermal][QuadQ1][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ1::testJacobianFiniteDiff", "[TwoBlocksStaticThermal][QuadQ1][Jacobian finite difference]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("TwoBlocksStaticThermal::QuadQ2::testDiscretization", "[TwoBlocksStaticThermal][QuadQ2][discretization]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ2()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ2::testResidual", "[TwoBlocksStaticThermal][QuadQ2][residual]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ2()).testResidual();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ2::testJacobianTaylorSeries", "[TwoBlocksStaticThermal][QuadQ2][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ2::testJacobianFiniteDiff", "[TwoBlocksStaticThermal][QuadQ2][Jacobian finite difference]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("TwoBlocksStaticThermal::QuadQ3::testDiscretization", "[TwoBlocksStaticThermal][QuadQ3][discretization]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ3()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ3::testResidual", "[TwoBlocksStaticThermal][QuadQ3][residual]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ3()).testResidual();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ3::testJacobianTaylorSeries", "[TwoBlocksStaticThermal][QuadQ3][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticThermal::QuadQ3::testJacobianFiniteDiff", "[TwoBlocksStaticThermal][QuadQ3][Jacobian finite difference]") {
    pylith::TestFaultKinThermoelasticity(pylith::TwoBlocksStaticThermal::QuadQ3()).testJacobianFiniteDiff();
}

// End of file
