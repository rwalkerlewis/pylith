// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

namespace pylith {
    namespace faults {
        class FaultCohesive;
        class FaultCohesiveKin;
        class FaultCohesiveKinPoro;
        class FaultCohesiveImpulses;
        class AuxiliaryFieldFactory;
        class DiagnosticFieldFactory;
        class DerivedFieldFactory;
        class AuxiliaryFactoryKinematicPoro;

        class KinSrc;
        class KinSrcPoro;
        class KinSrcConstRate;
        class KinSrcStep;
        class KinSrcPoroStep;
        class KinSrcRamp;
        class KinSrcBrune;
        class KinSrcLiuCos;
        class KinSrcTimeHistory;
        class KinSrcAuxiliaryFactory;
        class KinSrcPoroAuxiliaryFactory;

        class TopologyOps;
        class FaultOps;
    } // faults
} // pylith

// End of file
