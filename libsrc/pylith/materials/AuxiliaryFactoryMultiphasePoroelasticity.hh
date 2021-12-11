// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/AuxiliaryFactoryMultiphasePoroelasticity.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the multiphaseporoelasticity equation.
 */

#if !defined(pylith_materials_auxiliaryfactorymultiphaseporoelasticity_hh)
#define pylith_materials_auxiliaryfactorymultiphaseporoelasticity_hh

#include "materialsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::materials::AuxiliaryFactoryMultiphasePoroelasticity : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryMultiphasePoroelasticity; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryMultiphasePoroelasticity(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryMultiphasePoroelasticity(void);

    /// Add body force subfield to auxiliary subfields.
    void addBodyForce(void);

    /** Add gravity subfield to auxiliary subfields.
     *
     * @param[in] gf Gravity field.
     */
    void addGravityField(spatialdata::spatialdb::GravityField* gf);

    /// Add porosity subfield to auxiliary subfields.
    void addPorosity(void);

    /// Add solid density subfield to auxiliary subfields.
    void addSolidDensity(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addFluidDensity(void);

    /// Add fluid viscosity subfield to auxiliary subfields.
    void addFluidViscosity(void);

    /// Add reference sourceDensity subfield to auxiliary fields.
    void addSourceDensity(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryMultiphasePoroelasticity(const AuxiliaryFactoryMultiphasePoroelasticity &); ///< Not implemented.
    const AuxiliaryFactoryMultiphasePoroelasticity& operator=(const AuxiliaryFactoryMultiphasePoroelasticity&); ///< Not
                                                                                                                ///< implemented

}; // class AuxiliaryFactoryMultiphasePoroelasticity

#endif // pylith_materials_auxiliaryfactorymultiphasemultiphaseporoelasticity_hh

// End of file
