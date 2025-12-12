// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/RheologyThermoporoelasticity.i
 *
 * Python interface to C++ RheologyThermoporoelasticity.
 */

namespace pylith {
    namespace materials {
        class RheologyThermoporoelasticity : public pylith::utils::PyreComponent {
            // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Destructor.
            virtual ~RheologyThermoporoelasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual pylith::materials::AuxiliaryFactoryThermoporoelasticity* getAuxiliaryFactory(void) = 0;

            /// Add rheology subfields to auxiliary field.
            virtual void addAuxiliarySubfields(void) = 0;

            /** Update kernel constants.
             *
             * @param[inout] kernelConstants Array of constants used in integration kernels.
             * @param[in] dt Current time step.
             */
            virtual void updateKernelConstants(pylith::real_array* kernelConstants,
                                               const PylithReal dt) const;

        }; // class RheologyThermoporoelasticity

    } // materials
} // pylith

// End of file
