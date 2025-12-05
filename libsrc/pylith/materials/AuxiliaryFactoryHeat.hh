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

class pylith::materials::AuxiliaryFactoryHeat : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryHeat; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryHeat(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryHeat(void);

    /// Add density subfield to auxiliary subfields.
    void addDensity(void);

    /// Add specific heat subfield to auxiliary subfields.
    void addSpecificHeat(void);

    /// Add thermal conductivity subfield to auxiliary subfields.
    void addThermalConductivity(void);

    /// Add heat source subfield to auxiliary subfields.
    void addHeatSource(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryHeat(const AuxiliaryFactoryHeat &); ///< Not implemented.
    const AuxiliaryFactoryHeat& operator=(const AuxiliaryFactoryHeat&); ///< Not implemented

}; // class AuxiliaryFactoryHeat

// End of file
