// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application elasticitylgdeformapp.

#include "ElasticityImplicitLgDeformGravData3DQuadratic.hh"

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_spaceDim = 3;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_cellDim = 3;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_numVertices = 10;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_numCells = 1;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_numBasis = 10;

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_numQuadPts = 4;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_matType = "ElasticIsotropic3D";

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_matDBFilename = "data/elasticisotropic3d.spatialdb";

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_matId = 0;

const char* pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_matLabel = "elastic isotropic 3-D";

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_dt =   1.00000000e-02;

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_gravityVec[] = {
  0.00000000e+00,  0.00000000e+00, -1.00000000e+08,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_vertices[] = {
 -5.00000000e-01, -2.00000000e+00, -1.00000000e+00,
  2.00000000e+00, -2.00000000e+00, -5.00000000e-01,
  1.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  2.00000000e-01,  5.00000000e-01,  2.00000000e+00,
  1.50000000e+00, -5.00000000e-01, -2.50000000e-01,
  2.50000000e-01, -5.00000000e-01, -5.00000000e-01,
  7.50000000e-01, -2.00000000e+00, -7.50000000e-01,
 -1.50000000e-01, -7.50000000e-01,  5.00000000e-01,
  1.10000000e+00, -7.50000000e-01,  7.50000000e-01,
  6.00000000e-01,  7.50000000e-01,  1.00000000e+00,
};

const int pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_cells[] = {
0,1,2,3,4,5,6,7,8,9,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_verticesRef[] = {
 -1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
  1.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  1.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -1.00000000e+00,
 -1.00000000e+00,  0.00000000e+00, -1.00000000e+00,
  0.00000000e+00, -1.00000000e+00, -1.00000000e+00,
 -1.00000000e+00, -1.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -1.00000000e+00,  0.00000000e+00,
 -1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_quadPts[] = {
 -8.00000000e-01, -8.00000000e-01, -8.00000000e-01,
  5.00000000e-01, -8.00000000e-01, -8.00000000e-01,
 -8.00000000e-01,  5.00000000e-01, -8.00000000e-01,
 -8.00000000e-01, -8.00000000e-01,  5.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_quadWts[] = {
  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,  3.33333333e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_basis[] = {
  2.80000000e-01, -8.00000000e-02, -8.00000000e-02,
 -8.00000000e-02,  4.00000000e-02,  2.80000000e-01,
  2.80000000e-01,  2.80000000e-01,  4.00000000e-02,
  4.00000000e-02, -4.50000000e-02,  3.75000000e-01,
 -8.00000000e-02, -8.00000000e-02,  3.00000000e-01,
  2.00000000e-02,  1.50000000e-01,  2.00000000e-02,
  3.00000000e-01,  4.00000000e-02, -4.50000000e-02,
 -8.00000000e-02,  3.75000000e-01, -8.00000000e-02,
  3.00000000e-01,  1.50000000e-01,  2.00000000e-02,
  2.00000000e-02,  4.00000000e-02,  3.00000000e-01,
 -4.50000000e-02, -8.00000000e-02, -8.00000000e-02,
  3.75000000e-01,  4.00000000e-02,  2.00000000e-02,
  2.00000000e-02,  1.50000000e-01,  3.00000000e-01,
  3.00000000e-01,};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_basisDerivRef[] = {
 -9.00000000e-01, -9.00000000e-01, -9.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  2.00000000e-01,  2.00000000e-01,  0.00000000e+00,
 -2.00000000e-01,  1.20000000e+00, -2.00000000e-01,
  1.20000000e+00, -2.00000000e-01, -2.00000000e-01,
 -2.00000000e-01, -2.00000000e-01,  1.20000000e+00,
  2.00000000e-01,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  2.00000000e-01,  2.00000000e-01,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
  1.00000000e+00,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  2.00000000e-01,  1.50000000e+00,  0.00000000e+00,
 -2.00000000e-01, -1.00000000e-01, -2.00000000e-01,
 -1.40000000e+00, -1.50000000e+00, -1.50000000e+00,
 -2.00000000e-01, -2.00000000e-01, -1.00000000e-01,
  2.00000000e-01,  0.00000000e+00,  1.50000000e+00,
  0.00000000e+00,  2.00000000e-01,  2.00000000e-01,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  1.00000000e+00,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00, -3.00000000e-01,
  1.50000000e+00,  2.00000000e-01,  0.00000000e+00,
 -1.50000000e+00, -1.40000000e+00, -1.50000000e+00,
 -1.00000000e-01, -2.00000000e-01, -2.00000000e-01,
 -2.00000000e-01, -2.00000000e-01, -1.00000000e-01,
  2.00000000e-01,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  2.00000000e-01,  1.50000000e+00,
  4.00000000e-01,  4.00000000e-01,  4.00000000e-01,
 -3.00000000e-01,  0.00000000e+00,  0.00000000e+00,
  0.00000000e+00, -3.00000000e-01,  0.00000000e+00,
  0.00000000e+00,  0.00000000e+00,  1.00000000e+00,
  2.00000000e-01,  2.00000000e-01,  0.00000000e+00,
 -2.00000000e-01, -1.00000000e-01, -2.00000000e-01,
 -1.00000000e-01, -2.00000000e-01, -2.00000000e-01,
 -1.50000000e+00, -1.50000000e+00, -1.40000000e+00,
  1.50000000e+00,  0.00000000e+00,  2.00000000e-01,
  0.00000000e+00,  1.50000000e+00,  2.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_fieldTIncr[] = {
  3.00000000e-01, -4.00000000e-01, -4.00000000e-01,
 -6.00000000e-01,  8.00000000e-01,  2.00000000e-01,
  5.00000000e-01,  5.00000000e-01,  7.00000000e-01,
 -7.00000000e-01, -5.00000000e-01, -7.00000000e-01,
 -6.00000000e-01, -3.00000000e-01,  8.00000000e-01,
 -4.00000000e-01, -8.00000000e-01, -5.00000000e-01,
  7.00000000e-01,  8.00000000e-01, -5.00000000e-01,
 -5.00000000e-01, -5.00000000e-01, -7.00000000e-01,
 -3.00000000e-01, -9.00000000e-01,  8.00000000e-01,
 -1.00000000e-01,  5.00000000e-01, -9.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_fieldT[] = {
  1.00000000e-01, -2.00000000e-01, -6.00000000e-01,
 -3.00000000e-01,  4.00000000e-01,  9.00000000e-01,
  6.00000000e-01,  8.00000000e-01,  5.00000000e-01,
 -8.00000000e-01, -6.00000000e-01, -8.00000000e-01,
 -0.00000000e+00, -2.00000000e-01,  6.00000000e-01,
 -4.00000000e-01, -7.00000000e-01, -2.00000000e-01,
  7.00000000e-01,  6.00000000e-01, -1.00000000e-01,
 -4.00000000e-01, -3.00000000e-01, -3.00000000e-01,
 -7.00000000e-01, -6.00000000e-01,  1.00000000e-01,
 -9.00000000e-01,  3.00000000e-01, -8.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_fieldTmdt[] = {
  2.00000000e-01, -3.00000000e-01, -1.00000000e-01,
 -4.00000000e-01,  2.00000000e-01,  3.00000000e-01,
 -5.00000000e-01,  2.00000000e-01,  2.00000000e-01,
 -3.00000000e-01, -8.00000000e-01, -3.00000000e-01,
 -5.00000000e-01, -9.00000000e-01,  4.00000000e-01,
 -3.00000000e-01, -6.00000000e-01, -8.00000000e-01,
  9.00000000e-01,  5.00000000e-01, -2.00000000e-01,
 -7.00000000e-01, -2.00000000e-01, -9.00000000e-01,
 -5.00000000e-01, -8.00000000e-01,  4.00000000e-01,
 -4.00000000e-01,  5.00000000e-01, -7.00000000e-01,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_valsResidual[] = {
  2.83531193e+11,  6.85360097e+11, -7.60839301e+11,
  8.17609945e+11, -2.46059924e+11,  5.52982643e+11,
 -1.61369461e+12, -3.33728123e+11, -2.00813431e+12,
 -3.86612970e+11,  1.55209617e+12, -1.59708576e+12,
 -8.25945978e+11, -2.06081757e+11, -9.32629697e+11,
  2.39916408e+12,  1.36233180e+12,  2.35919362e+12,
 -1.83049591e+12, -1.27565741e+12,  2.69052281e+11,
  2.60556997e+11, -1.53217232e+12,  1.31884398e+12,
 -8.95741060e+11,  3.53530714e+12, -4.30653354e+12,
  1.79162831e+12, -3.54139567e+12,  4.39369174e+12,
};

const PylithScalar pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::_valsJacobian[] = {
  4.65661925e+11,  1.45207670e+11,  8.48382894e+10,
  1.24772472e+11,  7.41894782e+10, -4.42698676e+10,
  7.18213188e+10,  5.09242305e+09,  4.05530373e+10,
  9.25950209e+10,  3.01383081e+09,  4.41940539e+10,
  1.40659877e+11,  7.81646933e+09,  8.55770103e+10,
  7.55814313e+10,  1.00158408e+11, -5.39884553e+09,
 -7.04119227e+11, -2.55717104e+11, -3.23748072e+10,
 -3.27406977e+11, -3.58266405e+10, -1.14846242e+11,
  1.90099476e+11, -7.46675761e+10,  1.18407263e+10,
 -1.29665315e+11,  3.07330418e+10, -7.01133548e+10,
  1.45207670e+11,  4.17879082e+11,  2.69285551e+10,
  7.41894782e+10,  7.12501045e+10,  1.69326571e+10,
  5.09242305e+09,  7.79843128e+10,  1.24078373e+09,
  3.01383081e+09,  1.26067202e+11, -2.62782048e+10,
  7.81646933e+09,  1.18102474e+11,  2.96832254e+10,
  1.00158408e+11,  2.69890946e+10,  6.52474095e+10,
 -2.55717104e+11, -5.58070068e+11, -2.09279922e+10,
 -3.58266405e+10, -4.64117178e+11,  6.84864638e+10,
 -7.46675761e+10,  3.28253126e+11, -2.44433188e+11,
  3.07330418e+10, -1.44338150e+11,  8.31202909e+10,
  8.48382894e+10,  2.69285551e+10,  4.15301734e+11,
 -4.42698676e+10,  1.69326571e+10,  7.33047953e+10,
  4.05530373e+10,  1.24078373e+09,  1.27489728e+11,
  4.41940539e+10, -2.62782048e+10,  1.51395370e+11,
  8.55770103e+10,  2.96832254e+10,  2.11321369e+11,
 -5.39884553e+09,  6.52474095e+10, -1.46463462e+11,
 -3.23748072e+10, -2.09279922e+10, -5.59853491e+11,
 -1.14846242e+11,  6.84864638e+10, -5.09095160e+11,
  1.18407263e+10, -2.44433188e+11,  4.97279724e+11,
 -7.01133548e+10,  8.31202909e+10, -2.60680607e+11,
  1.24772472e+11,  7.41894782e+10, -4.42698676e+10,
  3.18994605e+11,  6.08636475e+10, -5.28381253e+10,
 -5.95541375e+10,  1.26084960e+10, -4.24114682e+10,
 -8.37120894e+10,  2.17101919e+10, -9.24436085e+08,
 -2.46254169e+11, -4.29845346e+10, -4.41618213e+10,
  3.69282312e+11,  1.07296224e+11,  1.48157570e+11,
 -5.47553983e+11, -1.78849586e+11,  9.36194643e+10,
  7.95730371e+10, -7.67652780e+10,  4.79024142e+10,
 -2.15276680e+11,  8.00666055e+10, -1.77270676e+11,
  2.59728633e+11, -5.81352446e+10,  7.21969464e+10,
  7.41894782e+10,  7.12501045e+10,  1.69326571e+10,
  6.08636475e+10,  3.25267756e+11, -1.22309091e+11,
  1.26084960e+10, -5.31861531e+10,  1.29709419e+10,
  2.17101919e+10, -1.52532444e+11,  1.00603438e+11,
 -4.29845346e+10, -2.62321656e+11,  6.79528717e+10,
  1.07296224e+11,  3.27288566e+11,  5.30178659e+10,
 -1.78849586e+11, -4.32757567e+11,  7.65998486e+10,
 -7.67652780e+10,  1.29553761e+11, -1.26591592e+11,
  8.00666055e+10, -3.89057716e+11,  1.34140735e+11,
 -5.81352446e+10,  4.36495348e+11, -2.13317675e+11,
 -4.42698676e+10,  1.69326571e+10,  7.33047953e+10,
 -5.28381253e+10, -1.22309091e+11,  4.88759132e+11,
 -4.24114682e+10,  1.29709419e+10, -6.74464625e+10,
 -9.24436085e+08,  1.00603438e+11, -1.98419387e+11,
 -4.41618213e+10,  6.79528717e+10, -3.98133382e+11,
  1.48157570e+11,  5.30178659e+10,  3.01368983e+11,
  9.36194643e+10,  7.65998486e+10, -5.77601119e+11,
  4.79024142e+10, -1.26591592e+11,  1.94445962e+11,
 -1.77270676e+11,  1.34140735e+11, -2.74232321e+11,
  7.21969464e+10, -2.13317675e+11,  4.57953800e+11,
  7.18213188e+10,  5.09242305e+09,  4.05530373e+10,
 -5.95541375e+10,  1.26084960e+10, -4.24114682e+10,
  1.10135125e+12,  3.25720553e+11,  3.62390508e+11,
  2.94017649e+11,  1.23894657e+10,  8.94735240e+10,
  1.34408493e+10, -4.81163047e+10,  1.60257602e+11,
 -6.81749121e+11, -2.46273575e+11, -2.38131042e+11,
  1.52813392e+11,  6.38089974e+10,  2.92728877e+10,
 -1.19071391e+10,  6.22561142e+10, -3.62572653e+10,
  2.38104621e+11, -5.73710641e+10, -4.42540169e+09,
 -1.11833868e+12, -1.30115106e+11, -3.60722383e+11,
  5.09242305e+09,  7.79843128e+10,  1.24078373e+09,
  1.26084960e+10, -5.31861531e+10,  1.29709419e+10,
  3.25720553e+11,  1.09654928e+12,  1.25558927e+11,
  1.23894657e+10,  4.02203557e+11, -1.25992638e+11,
 -4.81163047e+10, -2.03145182e+11,  1.24652000e+11,
 -2.46273575e+11, -6.32348620e+11, -1.80580051e+11,
  6.38089974e+10,  9.79869768e+10,  3.14800157e+10,
  6.22561142e+10, -2.38293966e+10,  8.52080054e+10,
 -5.73710641e+10,  4.92506284e+11, -2.94720961e+11,
 -1.30115106e+11, -1.25472106e+12,  2.20182976e+11,
  4.05530373e+10,  1.24078373e+09,  1.27489728e+11,
 -4.24114682e+10,  1.29709419e+10, -6.74464625e+10,
  3.62390508e+11,  1.25558927e+11,  1.03669910e+12,
  8.94735240e+10, -1.25992638e+11,  4.11450418e+11,
  1.60257602e+11,  1.24652000e+11, -3.17810955e+09,
 -2.38131042e+11, -1.80580051e+11, -6.28819588e+11,
  2.92728877e+10,  3.14800157e+10,  3.14543861e+10,
 -3.62572653e+10,  8.52080054e+10, -1.29228969e+11,
 -4.42540169e+09, -2.94720961e+11,  5.35392012e+11,
 -3.60722383e+11,  2.20182976e+11, -1.31381251e+12,
  9.25950209e+10,  3.01383081e+09,  4.41940539e+10,
 -8.37120894e+10,  2.17101919e+10, -9.24436085e+08,
  2.94017649e+11,  1.23894657e+10,  8.94735240e+10,
  4.97528833e+11, -1.14626728e+10,  6.59882894e+10,
  6.65491224e+10,  8.47757035e+09,  1.31257632e+10,
  7.66331934e+10,  8.55959553e+10, -3.73685357e+07,
 -3.56709943e+10, -2.76429617e+10, -6.19572817e+10,
 -3.22925718e+11, -5.74029189e+09, -8.83327493e+10,
  3.17445492e+11, -1.22169413e+11,  1.00263616e+11,
 -9.02460510e+11,  3.58283248e+10, -1.61793411e+11,
  3.01383081e+09,  1.26067202e+11, -2.62782048e+10,
  2.17101919e+10, -1.52532444e+11,  1.00603438e+11,
  1.23894657e+10,  4.02203557e+11, -1.25992638e+11,
 -1.14626728e+10,  8.65931363e+11, -3.82661325e+11,
  8.47757035e+09, -1.28217470e+10,  1.35369129e+10,
  8.55959553e+10,  3.26503144e+10,  1.22944799e+11,
 -2.76429617e+10, -3.20097931e+10, -5.88793091e+10,
 -5.74029189e+09, -4.10305623e+11,  1.11266648e+11,
 -1.22169413e+11,  7.48622075e+11, -4.24857358e+11,
  3.58283248e+10, -1.56780490e+12,  6.70317037e+11,
  4.41940539e+10, -2.62782048e+10,  1.51395370e+11,
 -9.24436085e+08,  1.00603438e+11, -1.98419387e+11,
  8.94735240e+10, -1.25992638e+11,  4.11450418e+11,
  6.59882894e+10, -3.82661325e+11,  8.70711287e+11,
  1.31257632e+10,  1.35369129e+10,  9.07016473e+10,
 -3.73685357e+07,  1.22944799e+11, -5.25010034e+10,
 -6.19572817e+10, -5.88793091e+10,  1.08121328e+10,
 -8.83327493e+10,  1.11266648e+11, -5.16970075e+11,
  1.00263616e+11, -4.24857358e+11,  7.55013813e+11,
 -1.61793411e+11,  6.70317037e+11, -1.52219420e+12,
  1.40659877e+11,  7.81646933e+09,  8.55770103e+10,
 -2.46254169e+11, -4.29845346e+10, -4.41618213e+10,
  1.34408493e+10, -4.81163047e+10,  1.60257602e+11,
  6.65491224e+10,  8.47757035e+09,  1.31257632e+10,
  1.16270901e+12,  1.84156516e+11,  1.94618725e+11,
 -5.72995180e+11, -7.97684286e+10, -3.20409939e+11,
  2.64947735e+10,  1.68329992e+10, -9.87964535e+09,
 -2.49889410e+11, -3.10663177e+10, -7.33658238e+10,
 -1.71088244e+11, -6.30056465e+10,  1.74314405e+11,
 -1.69626626e+11,  4.76576773e+10, -1.80076276e+11,
  7.81646933e+09,  1.18102474e+11,  2.96832254e+10,
 -4.29845346e+10, -2.62321656e+11,  6.79528717e+10,
 -4.81163047e+10, -2.03145182e+11,  1.24652000e+11,
  8.47757035e+09, -1.28217470e+10,  1.35369129e+10,
  1.84156516e+11,  1.16609381e+12, -2.28629589e+11,
 -7.97684286e+10, -4.65673799e+11, -9.51163063e+10,
  1.68329992e+10,  7.30677747e+10, -6.02206307e+10,
 -3.10663177e+10, -2.30487506e+11,  7.81937180e+09,
 -6.30056465e+10, -2.28430845e+11,  1.98763507e+11,
  4.76576773e+10,  4.56166756e+10, -5.84413625e+10,
  8.55770103e+10,  2.96832254e+10,  2.11321369e+11,
 -4.41618213e+10,  6.79528717e+10, -3.98133382e+11,
  1.60257602e+11,  1.24652000e+11, -3.17810955e+09,
  1.31257632e+10,  1.35369129e+10,  9.07016473e+10,
  1.94618725e+11, -2.28629589e+11,  1.70825144e+12,
 -3.20409939e+11, -9.51163063e+10, -7.97801243e+11,
 -9.87964535e+09, -6.02206307e+10,  7.70360435e+10,
 -7.33658238e+10,  7.81937180e+09, -2.93644299e+11,
  1.74314405e+11,  1.98763507e+11, -4.10328847e+11,
 -1.80076276e+11, -5.84413625e+10, -1.84224621e+11,
  7.55814313e+10,  1.00158408e+11, -5.39884553e+09,
  3.69282312e+11,  1.07296224e+11,  1.48157570e+11,
 -6.81749121e+11, -2.46273575e+11, -2.38131042e+11,
  7.66331934e+10,  8.55959553e+10, -3.73685357e+07,
 -5.72995180e+11, -7.97684286e+10, -3.20409939e+11,
  2.45592464e+12,  8.72069906e+11,  7.29084849e+11,
 -9.75215619e+11, -4.34128069e+11, -2.56910788e+11,
 -5.51910015e+11, -3.31945132e+11, -1.04698326e+11,
 -5.96253782e+11, -1.03426849e+11, -1.13624994e+11,
  4.00702141e+11,  3.04215606e+10,  1.61968884e+11,
  1.00158408e+11,  2.69890946e+10,  6.52474095e+10,
  1.07296224e+11,  3.27288566e+11,  5.30178659e+10,
 -2.46273575e+11, -6.32348620e+11, -1.80580051e+11,
  8.55959553e+10,  3.26503144e+10,  1.22944799e+11,
 -7.97684286e+10, -4.65673799e+11, -9.51163063e+10,
  8.72069906e+11,  2.27130607e+12,  5.76476353e+11,
 -4.34128069e+11, -7.36984695e+11, -2.23203143e+11,
 -3.31945132e+11, -4.91827419e+11, -3.22783678e+11,
 -1.03426849e+11, -7.11586437e+11,  1.31379586e+11,
  3.04215606e+10,  3.80186919e+11, -1.27382836e+11,
 -5.39884553e+09,  6.52474095e+10, -1.46463462e+11,
  1.48157570e+11,  5.30178659e+10,  3.01368983e+11,
 -2.38131042e+11, -1.80580051e+11, -6.28819588e+11,
 -3.73685357e+07,  1.22944799e+11, -5.25010034e+10,
 -3.20409939e+11, -9.51163063e+10, -7.97801243e+11,
  7.29084849e+11,  5.76476353e+11,  2.14357815e+12,
 -2.56910788e+11, -2.23203143e+11, -4.02976629e+11,
 -1.04698326e+11, -3.22783678e+11, -2.22817271e+11,
 -1.13624994e+11,  1.31379586e+11, -7.44361231e+11,
  1.61968884e+11, -1.27382836e+11,  5.50793300e+11,
 -7.04119227e+11, -2.55717104e+11, -3.23748072e+10,
 -5.47553983e+11, -1.78849586e+11,  9.36194643e+10,
  1.52813392e+11,  6.38089974e+10,  2.92728877e+10,
 -3.56709943e+10, -2.76429617e+10, -6.19572817e+10,
  2.64947735e+10,  1.68329992e+10, -9.87964535e+09,
 -9.75215619e+11, -4.34128069e+11, -2.56910788e+11,
  1.89896017e+12,  6.56740815e+11, -8.16178375e+10,
  3.38393444e+11,  1.64966431e+11,  1.12967809e+11,
 -3.12949675e+10,  2.44157356e+10,  1.92718403e+11,
 -1.22806992e+11, -3.04272562e+10,  1.41617951e+10,
 -2.55717104e+11, -5.58070068e+11, -2.09279922e+10,
 -1.78849586e+11, -4.32757567e+11,  7.65998486e+10,
  6.38089974e+10,  9.79869768e+10,  3.14800157e+10,
 -2.76429617e+10, -3.20097931e+10, -5.88793091e+10,
  1.68329992e+10,  7.30677747e+10, -6.02206307e+10,
 -4.34128069e+11, -7.36984695e+11, -2.23203143e+11,
  6.56740815e+11,  1.40195492e+12, -6.50991546e+10,
  1.64966431e+11,  4.34421272e+11,  6.47670884e+10,
  2.44157356e+10, -1.07992730e+11,  2.24564280e+11,
 -3.04272562e+10, -1.39616091e+11,  3.09189976e+10,
 -3.23748072e+10, -2.09279922e+10, -5.59853491e+11,
  9.36194643e+10,  7.65998486e+10, -5.77601119e+11,
  2.92728877e+10,  3.14800157e+10,  3.14543861e+10,
 -6.19572817e+10, -5.88793091e+10,  1.08121328e+10,
 -9.87964535e+09, -6.02206307e+10,  7.70360435e+10,
 -2.56910788e+11, -2.23203143e+11, -4.02976629e+11,
 -8.16178375e+10, -6.50991546e+10,  1.51957659e+12,
  1.12967809e+11,  6.47670884e+10,  3.90219017e+11,
  1.92718403e+11,  2.24564280e+11, -4.04297886e+11,
  1.41617951e+10,  3.09189976e+10, -8.43690419e+10,
 -3.27406977e+11, -3.58266405e+10, -1.14846242e+11,
  7.95730371e+10, -7.67652780e+10,  4.79024142e+10,
 -1.19071391e+10,  6.22561142e+10, -3.62572653e+10,
 -3.22925718e+11, -5.74029189e+09, -8.83327493e+10,
 -2.49889410e+11, -3.10663177e+10, -7.33658238e+10,
 -5.51910015e+11, -3.31945132e+11, -1.04698326e+11,
  3.38393444e+11,  1.64966431e+11,  1.12967809e+11,
  1.31141204e+12,  1.23926458e+11,  3.32914096e+11,
 -6.26226600e+11,  2.47475905e+11, -2.25050603e+11,
  3.60887342e+11, -1.17281247e+11,  1.48766690e+11,
 -3.58266405e+10, -4.64117178e+11,  6.84864638e+10,
 -7.67652780e+10,  1.29553761e+11, -1.26591592e+11,
  6.22561142e+10, -2.38293966e+10,  8.52080054e+10,
 -5.74029189e+09, -4.10305623e+11,  1.11266648e+11,
 -3.10663177e+10, -2.30487506e+11,  7.81937180e+09,
 -3.31945132e+11, -4.91827419e+11, -3.22783678e+11,
  1.64966431e+11,  4.34421272e+11,  6.47670884e+10,
  1.23926458e+11,  1.60701662e+12, -2.64550313e+11,
  2.47475905e+11, -9.92458260e+11,  6.53136619e+11,
 -1.17281247e+11,  4.42033730e+11, -2.76758613e+11,
 -1.14846242e+11,  6.84864638e+10, -5.09095160e+11,
  4.79024142e+10, -1.26591592e+11,  1.94445962e+11,
 -3.62572653e+10,  8.52080054e+10, -1.29228969e+11,
 -8.83327493e+10,  1.11266648e+11, -5.16970075e+11,
 -7.33658238e+10,  7.81937180e+09, -2.93644299e+11,
 -1.04698326e+11, -3.22783678e+11, -2.22817271e+11,
  1.12967809e+11,  6.47670884e+10,  3.90219017e+11,
  3.32914096e+11, -2.64550313e+11,  1.69650893e+12,
 -2.25050603e+11,  6.53136619e+11, -1.33623672e+12,
  1.48766690e+11, -2.76758613e+11,  7.26818588e+11,
  1.90099476e+11, -7.46675761e+10,  1.18407263e+10,
 -2.15276680e+11,  8.00666055e+10, -1.77270676e+11,
  2.38104621e+11, -5.73710641e+10, -4.42540169e+09,
  3.17445492e+11, -1.22169413e+11,  1.00263616e+11,
 -1.71088244e+11, -6.30056465e+10,  1.74314405e+11,
 -5.96253782e+11, -1.03426849e+11, -1.13624994e+11,
 -3.12949675e+10,  2.44157356e+10,  1.92718403e+11,
 -6.26226600e+11,  2.47475905e+11, -2.25050603e+11,
  1.97967373e+12, -2.24128696e+11,  2.63978666e+11,
 -1.08518305e+12,  2.92810999e+11, -2.22744142e+11,
 -7.46675761e+10,  3.28253126e+11, -2.44433188e+11,
  8.00666055e+10, -3.89057716e+11,  1.34140735e+11,
 -5.73710641e+10,  4.92506284e+11, -2.94720961e+11,
 -1.22169413e+11,  7.48622075e+11, -4.24857358e+11,
 -6.30056465e+10, -2.28430845e+11,  1.98763507e+11,
 -1.03426849e+11, -7.11586437e+11,  1.31379586e+11,
  2.44157356e+10, -1.07992730e+11,  2.24564280e+11,
  2.47475905e+11, -9.92458260e+11,  6.53136619e+11,
 -2.24128696e+11,  3.09010212e+12, -1.54522896e+12,
  2.92810999e+11, -2.22995762e+12,  1.16725574e+12,
  1.18407263e+10, -2.44433188e+11,  4.97279724e+11,
 -1.77270676e+11,  1.34140735e+11, -2.74232321e+11,
 -4.42540169e+09, -2.94720961e+11,  5.35392012e+11,
  1.00263616e+11, -4.24857358e+11,  7.55013813e+11,
  1.74314405e+11,  1.98763507e+11, -4.10328847e+11,
 -1.13624994e+11,  1.31379586e+11, -7.44361231e+11,
  1.92718403e+11,  2.24564280e+11, -4.04297886e+11,
 -2.25050603e+11,  6.53136619e+11, -1.33623672e+12,
  2.63978666e+11, -1.54522896e+12,  3.55709452e+12,
 -2.22744142e+11,  1.16725574e+12, -2.17532306e+12,
 -1.29665315e+11,  3.07330418e+10, -7.01133548e+10,
  2.59728633e+11, -5.81352446e+10,  7.21969464e+10,
 -1.11833868e+12, -1.30115106e+11, -3.60722383e+11,
 -9.02460510e+11,  3.58283248e+10, -1.61793411e+11,
 -1.69626626e+11,  4.76576773e+10, -1.80076276e+11,
  4.00702141e+11,  3.04215606e+10,  1.61968884e+11,
 -1.22806992e+11, -3.04272562e+10,  1.41617951e+10,
  3.60887342e+11, -1.17281247e+11,  1.48766690e+11,
 -1.08518305e+12,  2.92810999e+11, -2.22744142e+11,
  2.50676306e+12, -1.01492749e+11,  5.98355251e+11,
  3.07330418e+10, -1.44338150e+11,  8.31202909e+10,
 -5.81352446e+10,  4.36495348e+11, -2.13317675e+11,
 -1.30115106e+11, -1.25472106e+12,  2.20182976e+11,
  3.58283248e+10, -1.56780490e+12,  6.70317037e+11,
  4.76576773e+10,  4.56166756e+10, -5.84413625e+10,
  3.04215606e+10,  3.80186919e+11, -1.27382836e+11,
 -3.04272562e+10, -1.39616091e+11,  3.09189976e+10,
 -1.17281247e+11,  4.42033730e+11, -2.76758613e+11,
  2.92810999e+11, -2.22995762e+12,  1.16725574e+12,
 -1.01492749e+11,  4.03210515e+12, -1.49589456e+12,
 -7.01133548e+10,  8.31202909e+10, -2.60680607e+11,
  7.21969464e+10, -2.13317675e+11,  4.57953800e+11,
 -3.60722383e+11,  2.20182976e+11, -1.31381251e+12,
 -1.61793411e+11,  6.70317037e+11, -1.52219420e+12,
 -1.80076276e+11, -5.84413625e+10, -1.84224621e+11,
  1.61968884e+11, -1.27382836e+11,  5.50793300e+11,
  1.41617951e+10,  3.09189976e+10, -8.43690419e+10,
  1.48766690e+11, -2.76758613e+11,  7.26818588e+11,
 -2.22744142e+11,  1.16725574e+12, -2.17532306e+12,
  5.98355251e+11, -1.49589456e+12,  3.80503835e+12,
};

pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::ElasticityImplicitLgDeformGravData3DQuadratic(void)
{ // constructor
  spaceDim = _spaceDim;
  cellDim = _cellDim;
  numVertices = _numVertices;
  numCells = _numCells;
  numBasis = _numBasis;
  numQuadPts = _numQuadPts;
  matType = const_cast<char*>(_matType);
  matDBFilename = const_cast<char*>(_matDBFilename);
  matId = _matId;
  matLabel = const_cast<char*>(_matLabel);
  dt = _dt;
  gravityVec = const_cast<PylithScalar*>(_gravityVec);
  vertices = const_cast<PylithScalar*>(_vertices);
  cells = const_cast<int*>(_cells);
  verticesRef = const_cast<PylithScalar*>(_verticesRef);
  quadPts = const_cast<PylithScalar*>(_quadPts);
  quadWts = const_cast<PylithScalar*>(_quadWts);
  basis = const_cast<PylithScalar*>(_basis);
  basisDerivRef = const_cast<PylithScalar*>(_basisDerivRef);
  fieldTIncr = const_cast<PylithScalar*>(_fieldTIncr);
  fieldT = const_cast<PylithScalar*>(_fieldT);
  fieldTmdt = const_cast<PylithScalar*>(_fieldTmdt);
  valsResidual = const_cast<PylithScalar*>(_valsResidual);
  valsJacobian = const_cast<PylithScalar*>(_valsJacobian);
} // constructor

pylith::feassemble::ElasticityImplicitLgDeformGravData3DQuadratic::~ElasticityImplicitLgDeformGravData3DQuadratic(void)
{}


// End of file
