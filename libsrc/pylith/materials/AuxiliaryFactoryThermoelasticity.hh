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

class pylith::materials::AuxiliaryFactoryThermoelasticity : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryThermoelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryThermoelasticity(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryThermoelasticity(void);

    /// Add density subfield to auxiliary subfields.
    void addDensity(void);

    /// Add body force subfield to auxiliary subfields.
    void addBodyForce(void);

    /// Add gravity subfield to auxiliary subfields.
    void addGravityField(spatialdata::spatialdb::GravityField* gf);

    /// Add specific heat subfield to auxiliary subfields.
    void addSpecificHeat(void);

    /// Add thermal conductivity subfield to auxiliary subfields.
    void addThermalConductivity(void);

    /// Add heat source subfield to auxiliary subfields.
    void addHeatSource(void);

    /// Add reference temperature subfield to auxiliary subfields.
    void addReferenceTemperature(void);

    /// Add thermal expansion coefficient subfield to auxiliary subfields.
    void addThermalExpansionCoeff(void);

    /// Add shear modulus subfield to auxiliary subfields.
    void addShearModulus(void);

    /// Add bulk modulus subfield to auxiliary subfields.
    void addBulkModulus(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryThermoelasticity(const AuxiliaryFactoryThermoelasticity &); ///< Not implemented.
    const AuxiliaryFactoryThermoelasticity& operator=(const AuxiliaryFactoryThermoelasticity&); ///< Not implemented

}; // class AuxiliaryFactoryThermoelasticity

// End of file
