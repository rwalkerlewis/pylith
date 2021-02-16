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

/** @file libsrc/faults/PointSource.hh
 *
 * @brief C++ abstract base class for a moment tensor point source
 *
 */

#if !defined(pylith_faults_pointsource_hh)
#define pylith_faults_pointsource_hh

#include "faultsfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include <string> // HASA std::string

class pylith::faults::PointSource : public pylith::problems::Physics {
    friend class TestPointSource; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    PointSource(void);

    /// Destructor.
    ~PointSource(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set value of label material-id used to identify material cells.
     *
     * @param value Material identifier
     */
    void setPointSourceId(const int value);

    /** Get value of label material-id used to identify material cells.
     *
     * @returns Material identifier
     */
    int getPointSourceId(void) const;

    /** Set descriptive label for material.
     *
     * @param value Label of material.
     */
    void setPointSourceLabel(const char* value);

    /** Get descruptive label of material.
     *
     * @returns Label of material
     */
    const char* getPointSourceLabel(void) const;

    /** Set origin time for source.
     *
     * @param value origin time.
     */
    void setOriginTime(const PylithReal value);

    /** Get origin time for source.
     *
     * @returns origin time of source.
     */
    PylithReal getOriginTime(void) const;

    /** Set user identified full moment tensor.
     *
     * @param vec Reference direction unit vector.
     */
    void setMomentTensor(const PylithReal vec[9]);

    /** Set dominant frequency of Ricker function.
     *
     * @param value dominant frequency of Ricker function.
     */
    void setDominantFrequency(const PylithReal value);

    /** Get dominant frequency of Ricker function.
     *
     * @param value dominant frequency of Ricker function.
     */
    PylithReal getDominantFrequency(const PylithReal value);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ domainMesh Finite-element mesh associated with integration domain.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& domainMesh);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:
    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _momentTensor[9]; ///< Full representation of the moment tensor.
    PylithReal _pointLocation[3]; ///< Cartesian representation of the origin of the point source.
    PylithReal _dominantFrequency[1]; ///< Dominant frequency (Hz) for Ricker source function.
    int _pointSourceId; ///< Identifier for point source cell.
    std::string _sourceLabel; ///< Label identifying point source.
    std::string _momentTensorConvention; ///< Label identifying convention of moment tensor reference frame.
    PylithReal _originTime; ///< Elapsed time to pass for implementation of the point source.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PointSource(const PointSource&); ///< Not implemented
    const PointSource& operator=(const PointSource&); ///< Not implemented

}; // class PointSource

#endif // pylith_faults_pointsource_hh

// End of file
