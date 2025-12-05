// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/IsotropicLinearThermoelasticity.i
 *
 * Python interface to C++ IsotropicLinearThermoelasticity.
 */

namespace pylith {
    namespace materials {
        class IsotropicLinearThermoelasticity : public pylith::materials::RheologyThermoelasticity {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            IsotropicLinearThermoelasticity(void);

            /// Destructor.
            ~IsotropicLinearThermoelasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Include reference stress/strain?
             *
             * @param value Flag indicating to include reference stress/strain.
             */
            void useReferenceState(const bool value);

            /** Include reference stress/strain?
             *
             * @returns True if including reference stress/strain, false otherwise.
             */
            bool useReferenceState(void) const;

        }; // class IsotropicLinearThermoelasticity

    } // materials
} // pylith

// End of file
