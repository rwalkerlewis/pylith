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

/** @file libsrc/faults/AuxiliaryFactoryKinematicPoro.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for faults.
 */

#if !defined(pylith_faults_auxiliaryfactorykinematicporo_hh)
#define pylith_faults_auxiliaryfactorykinematicporo_hh

#include "faultsfwd.hh"                          // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::AuxiliaryFactoryKinematicPoro : public pylith::feassemble::AuxiliaryFactory
{
    friend class TestAuxiliaryFactoryKinematic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /// Default constructor.
    AuxiliaryFactoryKinematicPoro(void);

    /// Destructor.
    ~AuxiliaryFactoryKinematicPoro(void);

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

    /// Add skempton coefficient to auxiliary field.
    void addSkemptonCoefficient(void);

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
    AuxiliaryFactoryKinematicPoro(const AuxiliaryFactoryKinematicPoro &);                  ///< Not implemented.
    const AuxiliaryFactoryKinematicPoro &operator=(const AuxiliaryFactoryKinematicPoro &); ///< Not implemented

}; // class AuxiliaryFactoryKinematicPoro

#endif // pylith_faults_auxiliaryfactorykinematicporo_hh

// End of file
