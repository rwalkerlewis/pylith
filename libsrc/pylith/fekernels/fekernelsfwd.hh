// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/** @file libsrc/fekernels/fekernelsfwd.hh
 *
 * @brief Forward declarations for PyLith fekernels objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_fekernels_fekernelsfwd_hh)
#define pylith_fekernels_fekernelsfwd_hh

namespace pylith {
    namespace fekernels {
        class DispVel;

        class Elasticity;
        class ElasticityPlaneStrain;
        class Elasticity3D;

        class IsotropicLinearElasticity;
        class IsotropicLinearElasticityPlaneStrain;
        class IsotropicLinearElasticity3D;

        class Viscoelasticity;
        class IsotropicLinearMaxwell;
        class IsotropicLinearMaxwellPlaneStrain;
        class IsotropicLinearMaxwell3D;

        class IsotropicLinearGenMaxwell;
        class IsotropicLinearGenMaxwellPlaneStrain;
        class IsotropicLinearGenMaxwell3D;

        class IsotropicPowerLaw;
        class IsotropicPowerLawPlaneStrain;
        class IsotropicPowerLaw3D;
        class IsotropicPowerLawEffectiveStress;

        class IncompressibleElasticity;
        class IsotropicLinearIncompElasticity;
        class IsotropicLinearIncompElasticityPlaneStrain;
        class IsotropicLinearIncompElasticity3D;

        class Poroelasticity;
        class IsotropicLinearPoroelasticity;

        class TimeDependentFn;
        class NeumannTimeDependent;
        class AbsorbingDampers;

        class FaultCohesiveKin;

    } // fekernels
} // pylith

#endif // pylith_fekernels_fekernelsfwd_hh

// End of file
