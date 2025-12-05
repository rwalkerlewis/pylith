// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/materials/IsotropicLinearThermoporoelasticity.i
 *
 * Python interface to C++ IsotropicLinearThermoporoelasticity.
 */

namespace pylith {
    namespace materials {
        class IsotropicLinearThermoporoelasticity : public pylith::materials::RheologyThermoporoelasticity {
            // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            IsotropicLinearThermoporoelasticity(void);

            /// Destructor.
            ~IsotropicLinearThermoporoelasticity(void);

            /// Deallocate PETSc and local data structures.
            void deallocate(void);

            /** Use reference stress and strain in computation of stress and strain?
             *
             * @param[in] value Flag indicating to include reference stress and strain.
             */
            void useReferenceState(const bool value);

            /** Use reference stress and strain in computation of stress and strain?
             *
             * @returns True if using reference stress and strain, false otherwise.
             */
            bool useReferenceState(void) const;

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::materials::AuxiliaryFactoryThermoporoelasticity* getAuxiliaryFactory(void);

            /// Add rheology subfields to auxiliary field.
            void addAuxiliarySubfields(void);

        }; // class IsotropicLinearThermoporoelasticity

    } // materials
} // pylith

// End of file
