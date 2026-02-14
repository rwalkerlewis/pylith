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

#include "TestPointForce.hh" // USES TestPointForce_Data

namespace pylith {
    class PointForceExplosion2D;
}

class pylith::PointForceExplosion2D {
public:

    // Data factory methods
    static TestPointForce_Data* TriP1(void);

    static TestPointForce_Data* TriP2(void);

    static TestPointForce_Data* QuadQ1(void);

    static TestPointForce_Data* QuadQ2(void);

private:

    PointForceExplosion2D(void); ///< Not implemented
}; // PointForceExplosion2D

// End of file
