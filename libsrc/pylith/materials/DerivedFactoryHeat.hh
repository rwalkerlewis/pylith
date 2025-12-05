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
#include "pylith/topology/FieldFactory.hh" // ISA FieldFactory

class pylith::materials::DerivedFactoryHeat : public pylith::topology::FieldFactory {
    friend class TestDerivedFactoryHeat; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    DerivedFactoryHeat(void);

    /// Destructor.
    virtual ~DerivedFactoryHeat(void);

    /// Add heat flux subfield to derived field.
    void addHeatFlux(void);

    /// Add subfields using discretizations provided.
    void addSubfields(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    DerivedFactoryHeat(const DerivedFactoryHeat &); ///< Not implemented.
    const DerivedFactoryHeat& operator=(const DerivedFactoryHeat&); ///< Not implemented

}; // class DerivedFactoryHeat

// End of file
