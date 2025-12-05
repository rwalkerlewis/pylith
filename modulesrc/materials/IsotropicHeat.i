// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/IsotropicHeat.i
 *
 * Python interface to C++ IsotropicHeat.
 */

namespace pylith {
    namespace materials {
        class IsotropicHeat : public pylith::materials::RheologyHeat {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            IsotropicHeat(void);

            /// Destructor.
            ~IsotropicHeat(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::materials::AuxiliaryFactoryHeat* getAuxiliaryFactory(void);

            /// Add rheology subfields to auxiliary field.
            void addAuxiliarySubfields(void);

            // Note: Kernel methods returning function pointers (PetscPointFn*, PetscPointJacFn*)
            // are not exposed to Python - they are internal implementation details used by
            // the C++ integrator.

        }; // class IsotropicHeat

    } // materials
} // pylith

// End of file
