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

            /** Get heat flux kernel for LHS residual.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS residual kernel for heat flux (f1).
             */
            PetscPointFn* getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get heat flux kernel for RHS residual.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return RHS residual kernel for heat flux.
             */
            PetscPointFn* getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get thermal conductivity kernel for LHS Jacobian.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return LHS Jacobian kernel for thermal conductivity.
             */
            PetscPointJacFn* getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const;

            /** Get heat flux kernel for derived field.
             *
             * @param[in] coordsys Coordinate system.
             *
             * @return Project kernel for computing heat flux subfield in derived field.
             */
            PetscPointFn* getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const;

        }; // class IsotropicHeat

    } // materials
} // pylith

// End of file
