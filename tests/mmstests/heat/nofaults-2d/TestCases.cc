// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestHeat
 */

#include "TestHeat.hh" // USES TestHeat

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "UniformHeat2D.hh"
// TriP1
TEST_CASE("UniformHeat2D::TriP1::testDiscretization", "[UniformHeat2D][TriP1][discretization]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP1()).testDiscretization();
}
TEST_CASE("UniformHeat2D::TriP1::testResidual", "[UniformHeat2D][TriP1][residual]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP1()).testResidual();
}
TEST_CASE("UniformHeat2D::TriP1::testJacobianTaylorSeries", "[UniformHeat2D][TriP1][Jacobian Taylor series]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformHeat2D::TriP1::testJacobianFiniteDiff", "[UniformHeat2D][TriP1][Jacobian finite difference]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP1()).testJacobianFiniteDiff();
}

// TriP2
TEST_CASE("UniformHeat2D::TriP2::testDiscretization", "[UniformHeat2D][TriP2][discretization]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP2()).testDiscretization();
}
TEST_CASE("UniformHeat2D::TriP2::testResidual", "[UniformHeat2D][TriP2][residual]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP2()).testResidual();
}
TEST_CASE("UniformHeat2D::TriP2::testJacobianTaylorSeries", "[UniformHeat2D][TriP2][Jacobian Taylor series]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformHeat2D::TriP2::testJacobianFiniteDiff", "[UniformHeat2D][TriP2][Jacobian finite difference]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP2()).testJacobianFiniteDiff();
}

// TriP3
TEST_CASE("UniformHeat2D::TriP3::testDiscretization", "[UniformHeat2D][TriP3][discretization]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP3()).testDiscretization();
}
TEST_CASE("UniformHeat2D::TriP3::testResidual", "[UniformHeat2D][TriP3][residual]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP3()).testResidual();
}
TEST_CASE("UniformHeat2D::TriP3::testJacobianTaylorSeries", "[UniformHeat2D][TriP3][Jacobian Taylor series]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformHeat2D::TriP3::testJacobianFiniteDiff", "[UniformHeat2D][TriP3][Jacobian finite difference]") {
    pylith::TestHeat(pylith::UniformHeat2D::TriP3()).testJacobianFiniteDiff();
}

// QuadQ1
TEST_CASE("UniformHeat2D::QuadQ1::testDiscretization", "[UniformHeat2D][QuadQ1][discretization]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ1()).testDiscretization();
}
TEST_CASE("UniformHeat2D::QuadQ1::testResidual", "[UniformHeat2D][QuadQ1][residual]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ1()).testResidual();
}
TEST_CASE("UniformHeat2D::QuadQ1::testJacobianTaylorSeries", "[UniformHeat2D][QuadQ1][Jacobian Taylor series]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ1()).testJacobianTaylorSeries();
}
TEST_CASE("UniformHeat2D::QuadQ1::testJacobianFiniteDiff", "[UniformHeat2D][QuadQ1][Jacobian finite difference]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ1()).testJacobianFiniteDiff();
}

// QuadQ2
TEST_CASE("UniformHeat2D::QuadQ2::testDiscretization", "[UniformHeat2D][QuadQ2][discretization]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ2()).testDiscretization();
}
TEST_CASE("UniformHeat2D::QuadQ2::testResidual", "[UniformHeat2D][QuadQ2][residual]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ2()).testResidual();
}
TEST_CASE("UniformHeat2D::QuadQ2::testJacobianTaylorSeries", "[UniformHeat2D][QuadQ2][Jacobian Taylor series]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ2()).testJacobianTaylorSeries();
}
TEST_CASE("UniformHeat2D::QuadQ2::testJacobianFiniteDiff", "[UniformHeat2D][QuadQ2][Jacobian finite difference]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ2()).testJacobianFiniteDiff();
}

// QuadQ3
TEST_CASE("UniformHeat2D::QuadQ3::testDiscretization", "[UniformHeat2D][QuadQ3][discretization]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ3()).testDiscretization();
}
TEST_CASE("UniformHeat2D::QuadQ3::testResidual", "[UniformHeat2D][QuadQ3][residual]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ3()).testResidual();
}
TEST_CASE("UniformHeat2D::QuadQ3::testJacobianTaylorSeries", "[UniformHeat2D][QuadQ3][Jacobian Taylor series]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ3()).testJacobianTaylorSeries();
}
TEST_CASE("UniformHeat2D::QuadQ3::testJacobianFiniteDiff", "[UniformHeat2D][QuadQ3][Jacobian finite difference]") {
    pylith::TestHeat(pylith::UniformHeat2D::QuadQ3()).testJacobianFiniteDiff();
}

// End of file
