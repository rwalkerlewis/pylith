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

#include <portinfo>

#include "KinSrcPoroStep.hh" // implementation of object methods

#include "pylith/faults/KinSrcPoroAuxiliaryFactory.hh" // USES KinSrcPoroAuxiliaryFactory

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcPoroStep::KinSrcPoroStep(void) {
    pylith::utils::PyreComponent::setName("kinsrcporostep");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcPoroStep::~KinSrcPoroStep(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Thickness time function kernel.
void
pylith::faults::KinSrcPoroStep::thicknessFn(const PylithInt dim,
                                            const PylithInt numS,
                                            const PylithInt numA,
                                            const PylithInt sOff[],
                                            const PylithInt sOff_x[],
                                            const PylithScalar s[],
                                            const PylithScalar s_t[],
                                            const PylithScalar s_x[],
                                            const PylithInt aOff[],
                                            const PylithInt aOff_x[],
                                            const PylithScalar a[],
                                            const PylithScalar a_t[],
                                            const PylithScalar a_x[],
                                            const PylithReal t,
                                            const PylithScalar x[],
                                            const PylithInt numConstants,
                                            const PylithScalar constants[],
                                            PylithScalar thickness[]) {
    // const PylithInt _numA = 2;

    thickness[0] += 0.0;

} // thicknessFn


// ---------------------------------------------------------------------------------------------------------------------
// Porosity time function kernel.
void
pylith::faults::KinSrcPoroStep::porosityFn(const PylithInt dim,
                                           const PylithInt numS,
                                           const PylithInt numA,
                                           const PylithInt sOff[],
                                           const PylithInt sOff_x[],
                                           const PylithScalar s[],
                                           const PylithScalar s_t[],
                                           const PylithScalar s_x[],
                                           const PylithInt aOff[],
                                           const PylithInt aOff_x[],
                                           const PylithScalar a[],
                                           const PylithScalar a_t[],
                                           const PylithScalar a_x[],
                                           const PylithReal t,
                                           const PylithScalar x[],
                                           const PylithInt numConstants,
                                           const PylithScalar constants[],
                                           PylithScalar porosity[]) {
    // const PylithInt _numA = 2;

    porosity[0] += 0.0;

} // porosityFn


// ---------------------------------------------------------------------------------------------------------------------
// Beta_p time function kernel.
void
pylith::faults::KinSrcPoroStep::beta_pFn(const PylithInt dim,
                                         const PylithInt numS,
                                         const PylithInt numA,
                                         const PylithInt sOff[],
                                         const PylithInt sOff_x[],
                                         const PylithScalar s[],
                                         const PylithScalar s_t[],
                                         const PylithScalar s_x[],
                                         const PylithInt aOff[],
                                         const PylithInt aOff_x[],
                                         const PylithScalar a[],
                                         const PylithScalar a_t[],
                                         const PylithScalar a_x[],
                                         const PylithReal t,
                                         const PylithScalar x[],
                                         const PylithInt numConstants,
                                         const PylithScalar constants[],
                                         PylithScalar beta_p[]) {
    // const PylithInt _numA = 2;

    beta_p[0] += 0.0;

} // beta_pFn


// ---------------------------------------------------------------------------------------------------------------------
// Beta_sigma time function kernel.
void
pylith::faults::KinSrcPoroStep::beta_sigmaFn(const PylithInt dim,
                                             const PylithInt numS,
                                             const PylithInt numA,
                                             const PylithInt sOff[],
                                             const PylithInt sOff_x[],
                                             const PylithScalar s[],
                                             const PylithScalar s_t[],
                                             const PylithScalar s_x[],
                                             const PylithInt aOff[],
                                             const PylithInt aOff_x[],
                                             const PylithScalar a[],
                                             const PylithScalar a_t[],
                                             const PylithScalar a_x[],
                                             const PylithReal t,
                                             const PylithScalar x[],
                                             const PylithInt numConstants,
                                             const PylithScalar constants[],
                                             PylithScalar beta_sigma[]) {
    beta_sigma[0] += 0.0;

} // beta_sigmaFn


// ---------------------------------------------------------------------------------------------------------------------
// permeability_normal time function kernel.
void
pylith::faults::KinSrcPoroStep::fault_permeabilityFn(const PylithInt dim,
                                                     const PylithInt numS,
                                                     const PylithInt numA,
                                                     const PylithInt sOff[],
                                                     const PylithInt sOff_x[],
                                                     const PylithScalar s[],
                                                     const PylithScalar s_t[],
                                                     const PylithScalar s_x[],
                                                     const PylithInt aOff[],
                                                     const PylithInt aOff_x[],
                                                     const PylithScalar a[],
                                                     const PylithScalar a_t[],
                                                     const PylithScalar a_x[],
                                                     const PylithReal t,
                                                     const PylithScalar x[],
                                                     const PylithInt numConstants,
                                                     const PylithScalar constants[],
                                                     PylithScalar permeability_normal[]) {
    for (PylithInt i = 0; i < dim * dim; ++i) {
        permeability_normal[i] += 0.0;
    }

} // permeability_normalFn


// ---------------------------------------------------------------------------------------------------------------------
// fluid_viscosity time function kernel.
void
pylith::faults::KinSrcPoroStep::fluid_viscosityFn(const PylithInt dim,
                                                  const PylithInt numS,
                                                  const PylithInt numA,
                                                  const PylithInt sOff[],
                                                  const PylithInt sOff_x[],
                                                  const PylithScalar s[],
                                                  const PylithScalar s_t[],
                                                  const PylithScalar s_x[],
                                                  const PylithInt aOff[],
                                                  const PylithInt aOff_x[],
                                                  const PylithScalar a[],
                                                  const PylithScalar a_t[],
                                                  const PylithScalar a_x[],
                                                  const PylithReal t,
                                                  const PylithScalar x[],
                                                  const PylithInt numConstants,
                                                  const PylithScalar constants[],
                                                  PylithScalar fluid_viscosity[]) {
    fluid_viscosity[0] += 0.0;

} // fluid_viscosityFn


// ---------------------------------------------------------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcPoroStep::slipFn(const PylithInt dim,
                                       const PylithInt numS,
                                       const PylithInt numA,
                                       const PylithInt sOff[],
                                       const PylithInt sOff_x[],
                                       const PylithScalar s[],
                                       const PylithScalar s_t[],
                                       const PylithScalar s_x[],
                                       const PylithInt aOff[],
                                       const PylithInt aOff_x[],
                                       const PylithScalar a[],
                                       const PylithScalar a_t[],
                                       const PylithScalar a_x[],
                                       const PylithReal t,
                                       const PylithScalar x[],
                                       const PylithInt numConstants,
                                       const PylithScalar constants[],
                                       PylithScalar slip[]) {
    const PylithInt _numA = 2;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slip);

    const PylithInt i_initiationTime = 0;
    const PylithInt i_finalSlip = numA - 1;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar *finalSlip = &a[aOff[i_finalSlip]];

    const PylithInt i_originTime = 0;
    // const PylithScalar originTime = constants[i_originTime];
    // const PylithScalar t0 = originTime + initiationTime;
    // Ignoring originTime/setting as zero for sake of expediency
    const PylithScalar originTime = 0.0;
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        for (PylithInt i = 0; i < dim; ++i) {
            slip[i] = finalSlip[i];
        } // for
    } // if
} // slipFn


// ---------------------------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcPoroStep::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional &normalizer,
                                                     const spatialdata::geocoords::CoordSys *cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time function
    // kernel.

    // _auxiliaryFactory->addThickness(); // 0
    // _auxiliaryFactory->addPorosity(); // 1
    // _auxiliaryFactory->addBetaP(); // 2
    // _auxiliaryFactory->addBetaSigma(); // 3
    // _auxiliaryFactory->addPermeabilityTangential(); // 4
    // _auxiliaryFactory->addPermeabilityNormal(); // 5
    // _auxiliaryFactory->addFluidViscosity(); // 6
    _auxiliaryFactory->addInitiationTime(); // 0
    _auxiliaryFactory->addFinalSlip(); // 1

    // Add other kernels
    // _thicknessFnKernel = pylith::faults::KinSrcPoroStep::thicknessFn;                   // 0
    // _porosityFnKernel = pylith::faults::KinSrcPoroStep::porosityFn;                     // 1
    // _beta_pFnKernel = pylith::faults::KinSrcPoroStep::beta_pFn;                         // 2
    // _beta_sigmaFnKernel = pylith::faults::KinSrcPoroStep::beta_sigmaFn;                 // 3
    // _fault_permeabilityFnKernel = pylith::faults::KinSrcPoroStep::fault_permeabilityFn; // 4
    // _fluid_viscosityFnKernel = pylith::faults::KinSrcPoroStep::fluid_viscosityFn;       // 5
    _slipFnKernel = pylith::faults::KinSrcPoroStep::slipFn; // numA - 1
    _slipRateFnKernel = NULL; // Undefined for step function.
    _slipAccFnKernel = NULL; // Undefined for step function.

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// End of file
