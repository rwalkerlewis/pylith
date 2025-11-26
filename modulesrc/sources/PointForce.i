// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
// SWIG interface for PointForce

namespace pylith {
    namespace sources {
        class PointForce;
    } // sources
} // pylith

%{
#include "pylith/sources/PointForce.hh"
%}

%include "pylith/sources/PointForce.hh"

// End of file
