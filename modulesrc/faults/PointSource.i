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

/** @file modulesrc/faults/PointSource.i
 *
 * @brief Python interface to C++ PointSource object.
 */

namespace pylith {
    namespace faults {
        class PointSource : public pylith::problems::Physics {
            // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////
public:

            /// Default constructor.
            PointSource(void);

            /// Destructor.
            virtual ~PointSource(void);

            /// Deallocate PETSc and local data structures.
            virtual
            void deallocate(void);

            /** Set identifier for point source.
             *
             * @param[in] value Point source identifier
             */
            void setPointSourceId(const int value);

            /** Get identifier for point source.
             *
             * @returns Point source identifier
             */
            int getPointSourceId(void) const;

            /** Set label marking point source.
             *
             * @param[in] value Label of point source (from mesh generator).
             */
            void setPointSourceLabel(const char* value);

            /** Get label marking point source.
             *
             * @returns Label of point source (from mesh generator).
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
            
            // PROTECTED METHODS
            // ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

            /** Update kernel constants.
             *
             * @param[in] dt Current time step.
             */
            void _updateKernelConstants(const PylithReal dt);            

        }; // class PointSource

    } // faults
} // pylith

// End of file
