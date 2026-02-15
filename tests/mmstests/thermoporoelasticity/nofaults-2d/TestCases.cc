// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestThermoporoelasticity
 */

#include "TestThermoporoelasticity.hh" // USES TestThermoporoelasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "PressureGradientTemp.hh"

// TriP2P1P1P1
TEST_CASE("PressureGradientTemp::TriP2P1P1P1::testDiscretization", "[PressureGradientTemp][TriP2P1P1P1][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::TriP2P1P1P1::testResidual", "[PressureGradientTemp][TriP2P1P1P1][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1()).testResidual();
}
TEST_CASE("PressureGradientTemp::TriP2P1P1P1::testJacobianTaylorSeries", "[PressureGradientTemp][TriP2P1P1P1][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::TriP2P1P1P1::testJacobianFiniteDiff", "[PressureGradientTemp][TriP2P1P1P1][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1()).testJacobianFiniteDiff();
}

// TriP3P2P2P2
TEST_CASE("PressureGradientTemp::TriP3P2P2P2::testDiscretization", "[PressureGradientTemp][TriP3P2P2P2][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::TriP3P2P2P2::testResidual", "[PressureGradientTemp][TriP3P2P2P2][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2()).testResidual();
}
TEST_CASE("PressureGradientTemp::TriP3P2P2P2::testJacobianTaylorSeries", "[PressureGradientTemp][TriP3P2P2P2][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::TriP3P2P2P2::testJacobianFiniteDiff", "[PressureGradientTemp][TriP3P2P2P2][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1Q1
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1::testDiscretization", "[PressureGradientTemp][QuadQ2Q1Q1Q1][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1::testResidual", "[PressureGradientTemp][QuadQ2Q1Q1Q1][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1()).testResidual();
}
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1::testJacobianTaylorSeries", "[PressureGradientTemp][QuadQ2Q1Q1Q1][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1::testJacobianFiniteDiff", "[PressureGradientTemp][QuadQ2Q1Q1Q1][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2Q2
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2::testDiscretization", "[PressureGradientTemp][QuadQ3Q2Q2Q2][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2::testResidual", "[PressureGradientTemp][QuadQ3Q2Q2Q2][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2()).testResidual();
}
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2::testJacobianTaylorSeries", "[PressureGradientTemp][QuadQ3Q2Q2Q2][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2::testJacobianFiniteDiff", "[PressureGradientTemp][QuadQ3Q2Q2Q2][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2()).testJacobianFiniteDiff();
}

// TriP2P1P1P1 w/state variables
TEST_CASE("PressureGradientTemp::TriP2P1P1P1_StateVars::testDiscretization", "[PressureGradientTemp][TriP2P1P1P1_StateVars][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::TriP2P1P1P1_StateVars::testResidual", "[PressureGradientTemp][TriP2P1P1P1_StateVars][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1_StateVars()).testResidual();
}
TEST_CASE("PressureGradientTemp::TriP2P1P1P1_StateVars::testJacobianTaylorSeries", "[PressureGradientTemp][TriP2P1P1P1_StateVars][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::TriP2P1P1P1_StateVars::testJacobianFiniteDiff", "[PressureGradientTemp][TriP2P1P1P1_StateVars][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP2P1P1P1_StateVars()).testJacobianFiniteDiff();
}

// TriP3P2P2P2 with state variables
TEST_CASE("PressureGradientTemp::TriP3P2P2P2_StateVars::testDiscretization", "[PressureGradientTemp][TriP3P2P2P2_StateVars][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::TriP3P2P2P2_StateVars::testResidual", "[PressureGradientTemp][TriP3P2P2P2_StateVars][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2_StateVars()).testResidual();
}
TEST_CASE("PressureGradientTemp::TriP3P2P2P2_StateVars::testJacobianTaylorSeries", "[PressureGradientTemp][TriP3P2P2P2_StateVars][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::TriP3P2P2P2_StateVars::testJacobianFiniteDiff", "[PressureGradientTemp][TriP3P2P2P2_StateVars][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::TriP3P2P2P2_StateVars()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1Q1 with state variables
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars::testDiscretization", "[PressureGradientTemp][QuadQ2Q1Q1Q1_StateVars][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars::testResidual", "[PressureGradientTemp][QuadQ2Q1Q1Q1_StateVars][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars()).testResidual();
}
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars::testJacobianTaylorSeries", "[PressureGradientTemp][QuadQ2Q1Q1Q1_StateVars][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars::testJacobianFiniteDiff", "[PressureGradientTemp][QuadQ2Q1Q1Q1_StateVars][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ2Q1Q1Q1_StateVars()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2Q2 with state variables
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars::testDiscretization", "[PressureGradientTemp][QuadQ3Q2Q2Q2_StateVars][discretization]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars::testResidual", "[PressureGradientTemp][QuadQ3Q2Q2Q2_StateVars][residual]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars()).testResidual();
}
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars::testJacobianTaylorSeries", "[PressureGradientTemp][QuadQ3Q2Q2Q2_StateVars][Jacobian Taylor series]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars::testJacobianFiniteDiff", "[PressureGradientTemp][QuadQ3Q2Q2Q2_StateVars][Jacobian finite difference]") {
    pylith::TestThermoporoelasticity(pylith::PressureGradientTemp::QuadQ3Q2Q2Q2_StateVars()).testJacobianFiniteDiff();
}

// End of file
