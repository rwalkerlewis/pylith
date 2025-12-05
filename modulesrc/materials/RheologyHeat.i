// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/RheologyHeat.i
 *
 * Python interface to C++ RheologyHeat.
 */

namespace pylith {
    namespace materials {
        class RheologyHeat : public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            RheologyHeat(void);

            /// Destructor.
            virtual ~RheologyHeat(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

        }; // class RheologyHeat

    } // materials
} // pylith

// End of file
