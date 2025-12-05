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

class pylith::materials::RheologyThermoporoelasticity : public pylith::utils::PyreComponent {
    friend class TestIsotropicLinearThermoporoelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RheologyThermoporoelasticity(void);

    /// Destructor.
    virtual ~RheologyThermoporoelasticity(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::materials::AuxiliaryFactoryThermoporoelasticity* getAuxiliaryFactory(void) = 0;

    /// Add rheology subfields to auxiliary field.
    virtual
    void addAuxiliarySubfields(void) = 0;

    // ============================= RHS (Explicit) ==================================== //

    /** Get stress kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     * @return RHS residual kernel for stress.
     */
    virtual
    PetscPointFn* getKernelg1u_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Darcy flux kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     * @param[in] gravityField Flag for gravity.
     * @return RHS residual kernel for Darcy flux.
     */
    virtual
    PetscPointFn* getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool gravityField) const = 0;

    /** Get heat flux kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     * @return RHS residual kernel for heat flux.
     */
    virtual
    PetscPointFn* getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // =============================== LHS (Implicit) =================================== //

    /** Get stress kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     * @return LHS residual kernel for stress.
     */
    virtual
    PetscPointFn* getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get fluid content kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     * @param[in] useSourceDensity Flag for source density.
     * @return LHS residual kernel for fluid content.
     */
    virtual
    PetscPointFn* getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool useSourceDensity) const = 0;

    /** Get Darcy flux kernel for LHS residual.
     *
     * @param[in] coordsys Coordinate system.
     * @param[in] gravityField Flag for gravity.
     * @return LHS residual kernel for Darcy flux.
     */
    virtual
    PetscPointFn* getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool gravityField) const = 0;

    /** Get heat capacity kernel for LHS residual.
     *
     * @param[in] coordsys Coordinate system.
     * @param[in] useHeatSource Flag for heat source.
     * @return LHS residual kernel for heat capacity.
     */
    virtual
    PetscPointFn* getKernelf0T_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool useHeatSource) const = 0;

    /** Get heat flux kernel for LHS residual.
     *
     * @param[in] coordsys Coordinate system.
     * @return LHS residual kernel for heat flux.
     */
    virtual
    PetscPointFn* getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================= LHS Jacobians ==================================== //

    /** Get elastic stiffness kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf3_uu.
     */
    virtual
    PetscPointJacFn* getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Biot coupling kernel for LHS Jacobian (∂σ/∂p).
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf2_up.
     */
    virtual
    PetscPointJacFn* getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get thermal coupling kernel for LHS Jacobian (∂σ/∂T).
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf2_uT.
     */
    virtual
    PetscPointJacFn* getKernelJf2uT(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get storage coefficient kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf0_pp.
     */
    virtual
    PetscPointJacFn* getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Biot coupling kernel for LHS Jacobian (∂ζ/∂ε_v).
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf0_pe.
     */
    virtual
    PetscPointJacFn* getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get thermal coupling kernel for LHS Jacobian (∂ζ/∂T).
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf0_pT.
     */
    virtual
    PetscPointJacFn* getKernelJf0pT(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get Darcy conductivity kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf3_pp.
     */
    virtual
    PetscPointJacFn* getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get heat capacity kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf0_TT.
     */
    virtual
    PetscPointJacFn* getKernelJf0TT(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get thermal conductivity kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     * @return Jacobian kernel Jf3_TT.
     */
    virtual
    PetscPointJacFn* getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================ DERIVED FIELDS ========================== //

    /** Get Cauchy stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     * @return Project kernel for computing stress subfield in derived field.
     */
    virtual
    PetscPointFn* getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get fluid content kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     * @return Project kernel for computing fluid content subfield in derived field.
     */
    virtual
    PetscPointFn* getKernelFluidContent(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get heat flux kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
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

    RheologyThermoporoelasticity(const RheologyThermoporoelasticity&); ///< Not implemented.
    const RheologyThermoporoelasticity& operator=(const RheologyThermoporoelasticity&); /// Not implemented.

}; // class RheologyThermoporoelasticity

// End of file
