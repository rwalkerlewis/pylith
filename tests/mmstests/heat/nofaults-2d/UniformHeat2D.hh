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

#include "TestHeat.hh" // USES TestHeat_Data

namespace pylith {
    class UniformHeat2D;
}

class pylith::UniformHeat2D {
public:

    // Data factory methods
    static TestHeat_Data* TriP1(void);

    static TestHeat_Data* TriP2(void);

    static TestHeat_Data* TriP3(void);

    static TestHeat_Data* QuadQ1(void);

    static TestHeat_Data* QuadQ2(void);

    static TestHeat_Data* QuadQ3(void);

private:

    UniformHeat2D(void); ///< Not implemented
}; // UniformHeat2D

// End of file
