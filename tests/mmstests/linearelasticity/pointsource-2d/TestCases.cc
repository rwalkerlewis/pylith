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

#include "TestPointForce.hh" // USES TestPointForce
#include "PointForceExplosion2D.hh" // USES PointForceExplosion2D

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Explosion source tests
// ------------------------------------------------------------------------------------------------

// Triangle P1
TEST_CASE("PointForceExplosion2D::TriP1::testDiscretization", "[PointForceExplosion2D][TriP1][testDiscretization]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::TriP1());
    test.testDiscretization();
}
TEST_CASE("PointForceExplosion2D::TriP1::testResidual", "[PointForceExplosion2D][TriP1][testResidual]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::TriP1());
    test.testResidual();
}
TEST_CASE("PointForceExplosion2D::TriP1::testJacobianTaylorSeries", "[PointForceExplosion2D][TriP1][testJacobianTaylorSeries]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::TriP1());
    test.testJacobianTaylorSeries();
}

// Triangle P2
TEST_CASE("PointForceExplosion2D::TriP2::testDiscretization", "[PointForceExplosion2D][TriP2][testDiscretization]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::TriP2());
    test.testDiscretization();
}
TEST_CASE("PointForceExplosion2D::TriP2::testResidual", "[PointForceExplosion2D][TriP2][testResidual]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::TriP2());
    test.testResidual();
}
TEST_CASE("PointForceExplosion2D::TriP2::testJacobianTaylorSeries", "[PointForceExplosion2D][TriP2][testJacobianTaylorSeries]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::TriP2());
    test.testJacobianTaylorSeries();
}

// Quad Q1
TEST_CASE("PointForceExplosion2D::QuadQ1::testDiscretization", "[PointForceExplosion2D][QuadQ1][testDiscretization]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::QuadQ1());
    test.testDiscretization();
}
TEST_CASE("PointForceExplosion2D::QuadQ1::testResidual", "[PointForceExplosion2D][QuadQ1][testResidual]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::QuadQ1());
    test.testResidual();
}
TEST_CASE("PointForceExplosion2D::QuadQ1::testJacobianTaylorSeries", "[PointForceExplosion2D][QuadQ1][testJacobianTaylorSeries]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::QuadQ1());
    test.testJacobianTaylorSeries();
}

// Quad Q2
TEST_CASE("PointForceExplosion2D::QuadQ2::testDiscretization", "[PointForceExplosion2D][QuadQ2][testDiscretization]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::QuadQ2());
    test.testDiscretization();
}
TEST_CASE("PointForceExplosion2D::QuadQ2::testResidual", "[PointForceExplosion2D][QuadQ2][testResidual]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::QuadQ2());
    test.testResidual();
}
TEST_CASE("PointForceExplosion2D::QuadQ2::testJacobianTaylorSeries", "[PointForceExplosion2D][QuadQ2][testJacobianTaylorSeries]") {
    pylith::TestPointForce test;
    test.setData(pylith::PointForceExplosion2D::QuadQ2());
    test.testJacobianTaylorSeries();
}


// End of file
