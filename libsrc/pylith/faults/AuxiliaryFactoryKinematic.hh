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
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/AuxiliaryFactoryKinematic.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for faults.
 */

#if !defined(pylith_faults_auxiliaryfactorykinematic_hh)
#define pylith_faults_auxiliaryfactorykinematic_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::AuxiliaryFactoryKinematic : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryKinematic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryKinematic(void);

    /// Destructor.
    ~AuxiliaryFactoryKinematic(void);

    /// Add fault strike direction subfield to auxiliary field.
    void addStrikeDir(void);

    /// Add fault up-dip direction modulus subfield to auxiliary field.
    void addUpDipDir(void);

    /// Add fault normal direction subfield to auxiliary field.
    void addNormalDir(void);

    /// Add slip subfield to auxiliary field.
    void addSlip(void);

    /// Add slip rate subfield to auxiliary field.
    void addSlipRate(void);

    /// Add slip acceleration subfield to auxiliary field.
    void addSlipAcceleration(void);

    /// Add undrained bulk modulus to auxiliary field.
    void addUndrainedBulkModulus(void);

    /// Add shear modulus to auxiliary field.
    void addShearModulus(void);

    /// Add skempton coefficient to auxiliary field.
    void addSkemptonCoefficient(void);

    /// Add functions for FaultPoroDiffusionCohesivekin

    /// Add layer thickness to auxiliary field.
    void addThickness(void);

    /// Add porosity to auxiliary field.
    void addPorosity(void);

    /// Add betaP to auxiliary field.
    void addBetaP(void);

    /// Add betaSigma to auxiliary field.
    void addBetaSigma(void);

    /// Add fault permeability to auxiliary field.
    void addFaultPermeability(void);

    /// Add fluid viscosity to auxiliary field.
    void addFluidViscosity(void);

    /// Add body force subfield to auxiliary fields.
    void addBodyForce(void);

    /// Add reference source subfield to auxiliary fields.
    void addSource(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryKinematic(const AuxiliaryFactoryKinematic &); ///< Not implemented.
    const AuxiliaryFactoryKinematic& operator=(const AuxiliaryFactoryKinematic&); ///< Not implemented

}; // class AuxiliaryFactoryKinematic

#endif // pylith_faults_auxiliaryfactorykinematic_hh

// End of file
