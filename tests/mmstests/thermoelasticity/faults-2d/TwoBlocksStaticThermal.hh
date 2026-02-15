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

#include "TestFaultKinThermoelasticity.hh" // USES TestFaultKinThermoelasticity_Data

namespace pylith {
    class TwoBlocksStaticThermal;
}

class pylith::TwoBlocksStaticThermal {
public:

    // Data factory methods
    static TestFaultKinThermoelasticity_Data* TriP1(void);

    static TestFaultKinThermoelasticity_Data* TriP2(void);

    static TestFaultKinThermoelasticity_Data* TriP3(void);

    static TestFaultKinThermoelasticity_Data* QuadQ1(void);

    static TestFaultKinThermoelasticity_Data* QuadQ2(void);

    static TestFaultKinThermoelasticity_Data* QuadQ3(void);

private:

    TwoBlocksStaticThermal(void); ///< Not implemented
}; // TwoBlocksStaticThermal

// End of file
