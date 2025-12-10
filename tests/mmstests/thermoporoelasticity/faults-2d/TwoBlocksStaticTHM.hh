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

#include "TestFaultKinThermoporoelasticity.hh" // USES TestFaultKinThermoporoelasticity_Data

namespace pylith {
    class TwoBlocksStaticTHM;
}

class pylith::TwoBlocksStaticTHM {
public:

    // Data factory methods
    static TestFaultKinThermoporoelasticity_Data* TriP2P1P1P1(void);

    static TestFaultKinThermoporoelasticity_Data* TriP3P2P2P2(void);

    static TestFaultKinThermoporoelasticity_Data* QuadQ2Q1Q1Q1(void);

    static TestFaultKinThermoporoelasticity_Data* QuadQ3Q2Q2Q2(void);

private:

    TwoBlocksStaticTHM(void); ///< Not implemented
}; // TwoBlocksStaticTHM

// End of file
