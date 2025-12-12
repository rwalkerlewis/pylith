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
#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES std::vector
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain::ProjectKernels

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES Coordsys

#include "petscds.h" // USES PetscPointFn*, PetscPointJacFn*

class pylith::materials::RheologyHeat : public pylith::utils::PyreComponent {
    friend class TestIsotropicHeat; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RheologyHeat(void);

    /// Destructor.
    virtual ~RheologyHeat(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::materials::AuxiliaryFactoryHeat* getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual
    void addAuxiliarySubfields(void) = 0;

    // ============================= LHS Residual ==================================== //

    /** Get heat flux kernel for LHS residual, F(t,s,\dot{s})
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS residual kernel for heat flux (f1).
     */
    virtual
    PetscPointFn* getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================= RHS Residual ==================================== //

    /** Get heat flux kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return RHS residual kernel for heat flux.
     */
    virtual
    PetscPointFn* getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // =============================== LHS Jacobian =================================== //

    /** Get thermal conductivity kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return LHS Jacobian kernel for thermal conductivity.
     */
    virtual
    PetscPointJacFn* getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================ DERIVED FIELDS ========================== //

    /** Get heat flux kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Project kernel for computing heat flux subfield in derived field.
     */
    virtual
    PetscPointFn* getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Update kernel constants.
     *
     * @param[inout] kernelConstants Array of constants used in integration kernels.
     * @param[in] dt Current time step.
     */
    virtual
    void updateKernelConstants(pylith::real_array* kernelConstants,
                               const PylithReal dt) const;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    RheologyHeat(const RheologyHeat&); ///< Not implemented.
    const RheologyHeat& operator=(const RheologyHeat&); /// Not implemented.

}; // class RheologyHeat

// End of file
