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
#include "pylith/materials/RheologyThermoporoelasticity.hh" // ISA RheologyThermoporoelasticity

class pylith::materials::IsotropicLinearThermoporoelasticity : public pylith::materials::RheologyThermoporoelasticity {
    friend class TestIsotropicLinearThermoporoelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
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
    pylith::materials::AuxiliaryFactoryThermoporoelasticity* getAuxiliaryFactory(void) override;

    /// Add rheology subfields to auxiliary field.
    void addAuxiliarySubfields(void) override;

    // ============================= RHS (Explicit) ==================================== //

    PetscPointFn* getKernelg1u_explicit(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointFn* getKernelg1p_explicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool gravityField) const override;

    PetscPointFn* getKernelg1T_explicit(const spatialdata::geocoords::CoordSys* coordsys) const override;

    // =============================== LHS (Implicit) =================================== //

    PetscPointFn* getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointFn* getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool useSourceDensity) const override;

    PetscPointFn* getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool gravityField) const override;

    PetscPointFn* getKernelf0T_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                        const bool useHeatSource) const override;

    PetscPointFn* getKernelf1T_implicit(const spatialdata::geocoords::CoordSys* coordsys) const override;

    // ============================= LHS Jacobians ==================================== //

    PetscPointJacFn* getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf2uT(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf0pT(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf0TT(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointJacFn* getKernelJf3TT(const spatialdata::geocoords::CoordSys* coordsys) const override;

    // ============================ DERIVED FIELDS ========================== //

    PetscPointFn* getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointFn* getKernelFluidContent(const spatialdata::geocoords::CoordSys* coordsys) const override;

    PetscPointFn* getKernelHeatFluxVector(const spatialdata::geocoords::CoordSys* coordsys) const override;

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::materials::AuxiliaryFactoryThermoporoelasticity* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    bool _useReferenceState; ///< Flag to use reference stress and strain.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    IsotropicLinearThermoporoelasticity(const IsotropicLinearThermoporoelasticity&); ///< Not implemented.
    const IsotropicLinearThermoporoelasticity& operator=(const IsotropicLinearThermoporoelasticity&); /// Not implemented.

}; // class IsotropicLinearThermoporoelasticity

// End of file
