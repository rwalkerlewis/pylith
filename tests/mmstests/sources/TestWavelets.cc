// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include "catch2/catch_test_macros.hpp"

#include "pylith/fekernels/RickerWavelet.hh"
#include "pylith/fekernels/SquareWavelet.hh"
#include "pylith/fekernels/GaussianWavelet.hh"

#include <cmath>

namespace {

static void _initOffsets(PylithInt aOff[], PylithInt aOff_x[]) {
    aOff[0] = 0; // moment_tensor (dim*dim)
    aOff[1] = 4; // time_delay
    aOff[2] = 5; // center_frequency (numA-1)

    aOff_x[0] = 0;
    aOff_x[1] = 0;
    aOff_x[2] = 0;
}

static void _initAux(PylithScalar a[],
                     const PylithScalar mt[],
                     const PylithScalar timeDelay,
                     const PylithScalar centerFrequency) {
    for (int i = 0; i < 4; ++i) { a[i] = mt[i]; }
    a[4] = timeDelay;
    a[5] = centerFrequency;
}

} // anonymous

TEST_CASE("Source time function kernels: 2D g1v", "[wavelets]") {
    const PylithInt dim = 2;
    const PylithInt numS = 1;
    const PylithInt numA = 3;
    const PylithInt numConstants = 0;

    const PylithInt sOff[numS] = {0};
    const PylithInt sOff_x[numS] = {0};
    const PylithScalar s[1] = {0.0};
    const PylithScalar s_t[1] = {0.0};
    const PylithScalar s_x[1] = {0.0};

    PylithInt aOff[numA];
    PylithInt aOff_x[numA];
    _initOffsets(aOff, aOff_x);

    PylithScalar a[6] = {0.0};
    const PylithScalar a_t[1] = {0.0};
    const PylithScalar a_x[1] = {0.0};
    const PylithScalar x[2] = {0.0, 0.0};
    const PylithScalar* constants = nullptr;

    const PylithScalar mt[4] = {1.0, 2.0, 3.0, 4.0};
    const PylithScalar timeDelay = 0.25;
    const PylithScalar f0 = 2.0;
    _initAux(a, mt, timeDelay, f0);

    SECTION("RickerWaveletPlaneStrain::g1v matches formula") {
        const PylithReal t = 0.5;
        PylithScalar g1[4] = {0.0, 0.0, 0.0, 0.0};
        pylith::fekernels::RickerWaveletPlaneStrain::g1v(dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                                        aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, g1);

        const double rt = static_cast<double>(t - timeDelay);
        const double pi = std::acos(-1.0);
        const double ricker = (1.0 - 2.0*pi*pi*f0*f0*rt*rt) * std::exp(-pi*pi*f0*f0*rt*rt);
        for (int i = 0; i < 4; ++i) {
            REQUIRE(g1[i] == Catch::Approx(-mt[i] * ricker));
        }
    }

    SECTION("SquareWaveletPlaneStrain::g1v matches step function") {
        // Before time delay => 0
        {
            const PylithReal t = 0.1;
            PylithScalar g1[4] = {0.0, 0.0, 0.0, 0.0};
            pylith::fekernels::SquareWaveletPlaneStrain::g1v(dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, g1);
            for (int i = 0; i < 4; ++i) {
                REQUIRE(g1[i] == Catch::Approx(0.0));
            }
        }

        // After time delay => 1
        {
            const PylithReal t = 0.5;
            PylithScalar g1[4] = {0.0, 0.0, 0.0, 0.0};
            pylith::fekernels::SquareWaveletPlaneStrain::g1v(dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                                            aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, g1);
            for (int i = 0; i < 4; ++i) {
                REQUIRE(g1[i] == Catch::Approx(-mt[i]));
            }
        }
    }

    SECTION("GaussianWaveletPlaneStrain::g1v matches formula") {
        const PylithReal t = 0.5;
        PylithScalar g1[4] = {0.0, 0.0, 0.0, 0.0};
        pylith::fekernels::GaussianWaveletPlaneStrain::g1v(dim, numS, numA, sOff, sOff_x, s, s_t, s_x,
                                                          aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, g1);

        const double rt = static_cast<double>(t - timeDelay);
        const double pi = std::acos(-1.0);
        const double denom = 2.0 * (pi*pi * f0*f0);
        const double gaussian = std::exp((pi*pi * f0*f0) * rt*rt) / denom;
        for (int i = 0; i < 4; ++i) {
            REQUIRE(g1[i] == Catch::Approx(-mt[i] * gaussian));
        }
    }
}

// End of file

