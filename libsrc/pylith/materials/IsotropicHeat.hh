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

#include "pylith/materials/materialsfwd.hh" // forward declarations

#include "pylith/materials/RheologyHeat.hh" // ISA RheologyHeat

class pylith::materials::IsotropicHeat : public pylith::materials::RheologyHeat {
    friend class TestIsotropicHeat; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
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

    // ============================= LHS Residual ==================================== //

    /** Get heat flux kernel for LHS residual, F(t,s,\dot{s})
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for heat flux (f1).
     */
    PetscPointFn* getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ============================= RHS Residual ==================================== //

    /** Get heat flux kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for heat flux.
     */
    PetscPointFn* getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    // =============================== LHS Jacobian =================================== //

    /** Get thermal conductivity kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel for thermal conductivity.
     */
    PetscPointJacFn* getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ============================ DERIVED FIELDS ========================== //

    /** Get heat flux kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing heat flux subfield in derived field.
     */
    PetscPointFn* getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::materials::AuxiliaryFactoryHeat* _auxiliaryFactory; ///< Factory for auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IsotropicHeat(const IsotropicHeat&); ///< Not implemented.
    const IsotropicHeat& operator=(const IsotropicHeat&); /// Not implemented.

}; // class IsotropicHeat

// End of file
