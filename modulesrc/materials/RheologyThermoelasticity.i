// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/RheologyThermoelasticity.i
 *
 * Python interface to C++ abstract base class RheologyThermoelasticity.
 */

namespace pylith {
    namespace materials {
        // RheologyThermoelasticity is an abstract base class - no constructor exposed.
        // Concrete implementations (e.g., IsotropicLinearThermoelasticity) are used instead.
        class RheologyThermoelasticity : public pylith::utils::PyreComponent {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Destructor.
            virtual ~RheologyThermoelasticity(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            virtual
            pylith::materials::AuxiliaryFactoryThermoelasticity* getAuxiliaryFactory(void) = 0;

            /// Add rheology subfields to auxiliary field.
            virtual
            void addAuxiliarySubfields(void) = 0;

            /** Update kernel constants.
             *
             * @param[inout] kernelConstants Array of constants used in integration kernels.
             * @param[in] dt Current time step.
             */
            virtual
            void updateKernelConstants(pylith::real_array* kernelConstants,
                                       const PylithReal dt) const;

        }; // class RheologyThermoelasticity

    } // materials
} // pylith

// End of file
