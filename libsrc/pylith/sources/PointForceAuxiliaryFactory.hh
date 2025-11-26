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

#include "pylith/sources/sourcesfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::sources::PointForceAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestPointForceAuxiliaryFactory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    PointForceAuxiliaryFactory(void);

    /// Destructor.
    virtual ~PointForceAuxiliaryFactory(void);

    /** Add source location subfield to auxiliary subfields.
     *
     * @param[in] location Source location coordinates.
     * @param[in] size Number of coordinates (spatial dimension).
     */
    void addSourceLocation(const PylithReal* location,
                           const int size);

    /** Add moment tensor subfield to auxiliary subfields.
     *
     * @param[in] momentTensor Moment tensor components.
     * @param[in] size Number of components.
     */
    void addMomentTensor(const PylithReal* momentTensor,
                         const int size);

    /** Add magnitude subfield to auxiliary subfields.
     *
     * @param[in] value Magnitude value.
     */
    void addMagnitude(const PylithReal value);

    /** Add origin time subfield to auxiliary subfields.
     *
     * @param[in] value Origin time.
     */
    void addOriginTime(const PylithReal value);

    /** Add dominant frequency subfield to auxiliary subfields.
     *
     * @param[in] value Dominant frequency.
     */
    void addDominantFrequency(const PylithReal value);

    /** Add time delay subfield to auxiliary subfields.
     *
     * @param[in] value Time delay.
     */
    void addTimeDelay(const PylithReal value);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    PointForceAuxiliaryFactory(const PointForceAuxiliaryFactory&); ///< Not implemented.
    const PointForceAuxiliaryFactory& operator=(const PointForceAuxiliaryFactory&); ///< Not implemented.

}; // class PointForceAuxiliaryFactory

// End of file
