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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/ElasticIsotropic3D.hh
 *
 * @brief C++ 3-D, isotropic, linear elastic material.
 */

#if !defined(pylith_materials_elasticisotropic3d_hh)
#define pylith_materials_elasticisotropic3d_hh

// Include directives ---------------------------------------------------
#include "ElasticMaterial.hh" // ISA ElasticMaterial

// ElasticIsotropic3D ---------------------------------------------------
/** @brief 3-D, isotropic, linear elastic material. 
 *
 * The physical properties are specified using density, shear-wave
 * speed, and compressional-wave speed. The physical properties are
 * stored internally using density, lambda, and mu, which are directly
 * related to the elasticity constants used in the finite-element
 * integration.
 *
 * $\sigma - \sigma_0 = C (\epsilon - \epsilon_0)
 *
 * This implies that when $\epsilon = \epsilon_0$, $\sigma =
 * \sigma_0$.
 */
class pylith::materials::ElasticIsotropic3D : public ElasticMaterial
{ // class ElasticIsotropic3D
  friend class TestElasticIsotropic3D; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  ElasticIsotropic3D(void);

  /// Destructor
  ~ElasticIsotropic3D(void);

  /** Get stable time step for implicit time integration.
   *
   * Default is MAXDOUBLE (or 1.0e+30 if MAXFLOAT is not defined in math.h).
   *
   * @param mesh Finite-element mesh.
   * @returns Time step
   */
  double stableTimeStepImplicit(const topology::Mesh& mesh);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Compute properties from values in spatial database.
   *
   * Order of values in arrays matches order used in dbValues() and
   * parameterNames().
   *
   * @param propValues Array of property values.
   * @param dbValues Array of database values.
   */
  void _dbToProperties(double* const propValues,
		       const double_array& dbValues);

  /** Nondimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _nondimProperties(double* const values,
			 const int nvalues) const;

  /** Dimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _dimProperties(double* const values,
		      const int nvalues) const;

  /** Compute density from properties.
   *
   * @param density Array for density.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   */
  void _calcDensity(double* const density,
		    const double* properties,
		    const int numProperties,
		    const double* stateVars,
		    const int numStateVars);

  /** Compute stress tensor from properties and state variables. If
   * the state variables are from the previous time step, then the
   * computeStateVars flag should be set to true so that the state
   * variables are updated (but not stored) when computing the
   * stresses.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   * @param computeStateVars Flag indicating to compute updated state variables.
   */
  void _calcStress(double* const stress,
		   const int stressSize,
		   const double* properties,
		   const int numProperties,
		   const double* stateVars,
		   const int numStateVars,
		   const double* totalStrain,
		   const int strainSize,
		   const double* initialStress,
		   const int initialStressSize,
		   const double* initialStrain,
		   const int initialStrainSize,
		   const bool computeStateVars);

  /** Compute derivatives of elasticity matrix from properties.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _calcElasticConsts(double* const elasticConsts,
			  const int numElasticConsts,
			  const double* properties,
			  const int numProperties,
			  const double* stateVars,
			  const int numStateVars,
			  const double* totalStrain,
			  const int strainSize,
			  const double* initialStress,
			  const int initialStressSize,
			  const double* initialStrain,
			  const int initialStrainSize);

  /** Get stable time step for implicit time integration.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   *
   * @returns Time step
   */
  double _stableTimeStepImplicit(const double* properties,
				 const int numProperties,
				 const double* stateVars,
				 const int numStateVars) const;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  static const int p_density;
  static const int p_mu;
  static const int p_lambda;
  static const int db_density;
  static const int db_vs;
  static const int db_vp;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  ElasticIsotropic3D(const ElasticIsotropic3D&);

  /// Not implemented
  const ElasticIsotropic3D& operator=(const ElasticIsotropic3D&);

}; // class ElasticIsotropic3D

#include "ElasticIsotropic3D.icc" // inline methods

#endif // pylith_materials_elasticisotropic3d_hh


// End of file 