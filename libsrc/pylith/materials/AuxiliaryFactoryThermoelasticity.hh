// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/AuxiliaryFactoryThermoelasticity.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the thermoelastic equations.
 */

#if !defined(pylith_materials_auxiliaryfactorythermoelasticity_hh)
#define pylith_materials_auxiliaryfactorythermoelasticity_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

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

    /** Add gravity subfield to auxiliary subfields.
     *
     * @param[in] gf Gravity field.
     */
    void addGravityField(spatialdata::spatialdb::GravityField* gf);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryThermoelasticity(const AuxiliaryFactoryThermoelasticity &); ///< Not implemented.
    const AuxiliaryFactoryThermoelasticity& operator=(const AuxiliaryFactoryThermoelasticity&); ///< Not implemented

}; // class AuxiliaryFactoryThermoelasticity

#endif // pylith_materials_auxiliaryfactorythermoelasticity_hh

// End of file
