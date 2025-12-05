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

            /** Get auxiliary factory associated with physics.
             *
             * @return Auxiliary factory for physics object.
             */
            pylith::materials::AuxiliaryFactoryThermoelasticity* getAuxiliaryFactory(void);

            /** Add rheology subfields to auxiliary field.
             */
            void addAuxiliarySubfields(void);

            /** Get stress kernel for LHS residual.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Stress kernel for LHS residual.
             */
            PetscPointFunc getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get stress kernel for RHS residual.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Stress kernel for RHS residual.
             */
            PetscPointFunc getKernelg1u_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get elastic constants kernel for LHS Jacobian.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Elastic constants kernel for LHS Jacobian.
             */
            PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get coupling kernel Jf2_uT for stress-temperature coupling.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Coupling kernel for LHS Jacobian.
             */
            PetscPointJac getKernelJf2uT(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get heat flux kernel for LHS residual (implicit).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Heat flux kernel for LHS residual.
             */
            PetscPointFunc getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get heat flux kernel for RHS residual (explicit).
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Heat flux kernel for RHS residual.
             */
            PetscPointFunc getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get thermal conductivity kernel for LHS Jacobian.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Thermal conductivity kernel for LHS Jacobian.
             */
            PetscPointJac getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get stress kernel for derived field.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Stress kernel for derived field.
             */
            PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get heat flux kernel for derived field.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Heat flux kernel for derived field.
             */
            PetscPointFunc getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const;

        }; // class IsotropicLinearThermoelasticity

    } // materials
} // pylith

// End of file
