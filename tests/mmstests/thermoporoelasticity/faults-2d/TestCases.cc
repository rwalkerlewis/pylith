// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestFaultKinThermoporoelasticity
 */

#include "TestFaultKinThermoporoelasticity.hh" // USES TestFaultKinThermoporoelasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "TwoBlocksStaticTHM.hh"

// TriP2P1P1P1
TEST_CASE("TwoBlocksStaticTHM::TriP2P1P1P1::testDiscretization", "[TwoBlocksStaticTHM][TriP2P1P1P1][discretization]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP2P1P1P1()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticTHM::TriP2P1P1P1::testResidual", "[TwoBlocksStaticTHM][TriP2P1P1P1][residual]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP2P1P1P1()).testResidual();
}
TEST_CASE("TwoBlocksStaticTHM::TriP2P1P1P1::testJacobianTaylorSeries", "[TwoBlocksStaticTHM][TriP2P1P1P1][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP2P1P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticTHM::TriP2P1P1P1::testJacobianFiniteDiff", "[TwoBlocksStaticTHM][TriP2P1P1P1][Jacobian finite difference]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP2P1P1P1()).testJacobianFiniteDiff();
}

// TriP3P2P2P2
TEST_CASE("TwoBlocksStaticTHM::TriP3P2P2P2::testDiscretization", "[TwoBlocksStaticTHM][TriP3P2P2P2][discretization]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP3P2P2P2()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticTHM::TriP3P2P2P2::testResidual", "[TwoBlocksStaticTHM][TriP3P2P2P2][residual]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP3P2P2P2()).testResidual();
}
TEST_CASE("TwoBlocksStaticTHM::TriP3P2P2P2::testJacobianTaylorSeries", "[TwoBlocksStaticTHM][TriP3P2P2P2][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP3P2P2P2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticTHM::TriP3P2P2P2::testJacobianFiniteDiff", "[TwoBlocksStaticTHM][TriP3P2P2P2][Jacobian finite difference]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::TriP3P2P2P2()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1Q1
TEST_CASE("TwoBlocksStaticTHM::QuadQ2Q1Q1Q1::testDiscretization", "[TwoBlocksStaticTHM][QuadQ2Q1Q1Q1][discretization]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ2Q1Q1Q1()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticTHM::QuadQ2Q1Q1Q1::testResidual", "[TwoBlocksStaticTHM][QuadQ2Q1Q1Q1][residual]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ2Q1Q1Q1()).testResidual();
}
TEST_CASE("TwoBlocksStaticTHM::QuadQ2Q1Q1Q1::testJacobianTaylorSeries", "[TwoBlocksStaticTHM][QuadQ2Q1Q1Q1][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ2Q1Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticTHM::QuadQ2Q1Q1Q1::testJacobianFiniteDiff", "[TwoBlocksStaticTHM][QuadQ2Q1Q1Q1][Jacobian finite difference]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ2Q1Q1Q1()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2Q2
TEST_CASE("TwoBlocksStaticTHM::QuadQ3Q2Q2Q2::testDiscretization", "[TwoBlocksStaticTHM][QuadQ3Q2Q2Q2][discretization]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ3Q2Q2Q2()).testDiscretization();
}
TEST_CASE("TwoBlocksStaticTHM::QuadQ3Q2Q2Q2::testResidual", "[TwoBlocksStaticTHM][QuadQ3Q2Q2Q2][residual]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ3Q2Q2Q2()).testResidual();
}
TEST_CASE("TwoBlocksStaticTHM::QuadQ3Q2Q2Q2::testJacobianTaylorSeries", "[TwoBlocksStaticTHM][QuadQ3Q2Q2Q2][Jacobian Taylor series]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ3Q2Q2Q2()).testJacobianTaylorSeries();
}
TEST_CASE("TwoBlocksStaticTHM::QuadQ3Q2Q2Q2::testJacobianFiniteDiff", "[TwoBlocksStaticTHM][QuadQ3Q2Q2Q2][Jacobian finite difference]") {
    pylith::TestFaultKinThermoporoelasticity(pylith::TwoBlocksStaticTHM::QuadQ3Q2Q2Q2()).testJacobianFiniteDiff();
}

// End of file
