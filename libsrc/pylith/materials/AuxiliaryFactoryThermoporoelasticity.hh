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
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::materials::AuxiliaryFactoryThermoporoelasticity : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryThermoporoelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryThermoporoelasticity(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryThermoporoelasticity(void);

    // ============================= Base poroelastic fields =============================

    /// Add solid density subfield to auxiliary subfields.
    void addSolidDensity(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addFluidDensity(void);

    /// Add fluid viscosity subfield to auxiliary subfields.
    void addFluidViscosity(void);

    /// Add porosity subfield to auxiliary subfields.
    void addPorosity(void);

    // ============================= Optional fields =============================

    /// Add body force subfield to auxiliary subfields.
    void addBodyForce(void);

    /// Add gravity subfield to auxiliary subfields.
    void addGravityField(spatialdata::spatialdb::GravityField* gf);

    /// Add source density subfield to auxiliary subfields.
    void addSourceDensity(void);

    /// Add heat source subfield to auxiliary subfields.
    void addHeatSource(void);

    // ============================= Rheology fields =============================

    /// Add Biot coefficient subfield to auxiliary subfields.
    void addBiotCoefficient(void);

    /// Add Biot modulus subfield to auxiliary subfields.
    void addBiotModulus(void);

    /// Add drained bulk modulus subfield to auxiliary subfields.
    void addDrainedBulkModulus(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    /// Add isotropic permeability subfield to auxiliary subfields.
    void addIsotropicPermeability(void);

    /// Add reference temperature subfield to auxiliary subfields.
    void addReferenceTemperature(void);

    /// Add thermal expansion coefficient subfield to auxiliary subfields.
    void addThermalExpansionCoeff(void);

    /// Add fluid thermal expansion subfield to auxiliary subfields.
    void addFluidThermalExpansion(void);

    /// Add thermal conductivity subfield to auxiliary subfields.
    void addThermalConductivity(void);

    /// Add specific heat subfield to auxiliary subfields.
    void addSpecificHeat(void);

    /// Add reference stress subfield to auxiliary subfields.
    void addReferenceStress(void);

    /// Add reference strain subfield to auxiliary subfields.
    void addReferenceStrain(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryThermoporoelasticity(const AuxiliaryFactoryThermoporoelasticity &); ///< Not implemented.
    const AuxiliaryFactoryThermoporoelasticity& operator=(const AuxiliaryFactoryThermoporoelasticity&); ///< Not implemented

}; // class AuxiliaryFactoryThermoporoelasticity

// End of file
