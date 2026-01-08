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
#include "pylith/materials/AuxiliaryFactoryPoroelasticity.hh" // ISA AuxiliaryFactoryPoroelasticity

/**
 * @brief Factory for auxiliary subfields for black oil poroelasticity rheology.
 *
 * Additional subfields for black oil model:
 * - reference_pressure: Reference pressure for fluid property correlations
 * - fluid_compressibility: Fluid compressibility (1/Pa)
 * - fluid_compressibility_coefficient: Rate of change of compressibility with pressure
 * - viscosity_coefficient: Coefficient for pressure-dependent viscosity
 */
class pylith::materials::AuxiliaryFactoryBlackOilPoroelastic : public pylith::materials::AuxiliaryFactoryPoroelasticity {
    friend class TestAuxiliaryFactoryBlackOilPoroelastic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryBlackOilPoroelastic(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryBlackOilPoroelastic(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addIsotropicPermeability(void);

    /// Add tensor permeability subfield to auxiliary subfields.
    void addTensorPermeability(void);

    /// Add drained Bulk Modulus subfield to auxiliary subfields.
    void addDrainedBulkModulus(void);

    /// Add fluid Biot Coefficient subfield to auxiliary subfields.
    void addBiotCoefficient(void);

    /// Add fluid Biot Modulus subfield to auxiliary subfields.
    void addBiotModulus(void);

    /// Add reference stress subfield to auxiliary fields.
    void addReferenceStress(void);

    /// Add reference strain subfield to auxiliary fields.
    void addReferenceStrain(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    // Black Oil specific fields

    /// Add reference pressure subfield to auxiliary subfields.
    void addReferencePressure(void);

    /// Add fluid compressibility subfield to auxiliary subfields.
    void addFluidCompressibility(void);

    /// Add fluid compressibility coefficient (pressure dependence) subfield.
    void addFluidCompressibilityCoefficient(void);

    /// Add viscosity coefficient (pressure dependence) subfield.
    void addViscosityCoefficient(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryBlackOilPoroelastic(const AuxiliaryFactoryBlackOilPoroelastic &); ///< Not implemented.
    const AuxiliaryFactoryBlackOilPoroelastic& operator=(const AuxiliaryFactoryBlackOilPoroelastic&); ///< Not implemented

}; // class AuxiliaryFactoryBlackOilPoroelastic

// End of file
