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
#include "pylith/feassemble/feassemblefwd.hh" // USES AuxiliaryFactory
#include "pylith/utils/petscfwd.h" // USES PetscPointFunc, PetscPointJac

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys

class pylith::materials::RheologyThermoelasticity : public pylith::utils::PyreComponent {
    friend class TestRheologyThermoelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    RheologyThermoelasticity(void);

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

    /** Add rheology subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    virtual
    void addAuxiliarySubfields(void) = 0;

    // ============================= Displacement Equation =============================

    /** Get stress kernel for LHS residual, F(t,s,\dot{s}).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Stress kernel for LHS residual.
     */
    virtual
    PetscPointFunc getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get stress kernel for RHS residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Stress kernel for RHS residual.
     */
    virtual
    PetscPointFunc getKernelg1u_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get elastic constants kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Elastic constants kernel for LHS Jacobian.
     */
    virtual
    PetscPointJac getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get coupling kernel for stress-temperature coupling (Jf2_uT).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Coupling kernel for LHS Jacobian.
     */
    virtual
    PetscPointJac getKernelJf2uT(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================= Temperature Equation =============================

    /** Get heat flux kernel for LHS residual (implicit).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Heat flux kernel for LHS residual.
     */
    virtual
    PetscPointFunc getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get heat flux kernel for RHS residual (explicit).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Heat flux kernel for RHS residual.
     */
    virtual
    PetscPointFunc getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get thermal conductivity kernel for LHS Jacobian.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Thermal conductivity kernel for LHS Jacobian.
     */
    virtual
    PetscPointJac getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    // ============================= Two-Way Coupling Kernels =============================

    /** Get thermoelastic heating kernel for temperature equation (f0).
     *
     * For two-way coupling: f0_T includes -T α 3K ∇·v term.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Thermoelastic heating kernel for LHS residual.
     */
    virtual
    PetscPointFunc getKernelf0T_thermoelastic(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get Jacobian kernel for temperature-velocity coupling (Jf3_Tv).
     *
     * For two-way coupling: ∂f0_T/∂(∇v) = -T α 3K I
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Coupling kernel for temperature-velocity Jacobian.
     */
    virtual
    PetscPointJac getKernelJf3Tv(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Get Jacobian kernel for temperature equation with thermoelastic effect (Jf0_TT).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Jacobian kernel for temperature-temperature with thermoelastic effect.
     */
    virtual
    PetscPointJac getKernelJf0TT_thermoelastic(const spatialdata::geocoords::CoordSys* coordsys) const;

    // ============================= Derived Fields =============================

    /** Get stress kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Stress kernel for derived field.
     */
    virtual
    PetscPointFunc getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

    /** Get heat flux kernel for derived field.
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return Heat flux kernel for derived field.
     */
    virtual
    PetscPointFunc getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const = 0;

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

    RheologyThermoelasticity(const RheologyThermoelasticity &); ///< Not implemented.
    const RheologyThermoelasticity& operator=(const RheologyThermoelasticity&); ///< Not implemented

}; // class RheologyThermoelasticity

// End of file
