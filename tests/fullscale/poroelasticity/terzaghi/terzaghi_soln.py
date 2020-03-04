# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2018 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/linearporoelasticity/nofaults-2d/terzaghi_soln.py
#
# @brief Analytical solution to Terzaghi's problem.
#  
#        -1000 Pa
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         Uy=0
#
# Dirichlet boundary conditions
#   Ux(+-1,0) = 0
#   Uy(x,-1) = 0
# Neumann boundary conditions
#   \tau_normal(x,+1) = -1000*Pa
#

import numpy


# Physical properties
p_porosity               =     0.19 # frac,    phi
p_density                =   2500.0 # kg/m**3, rhos
p_fluid_density          =     1000 # kg/m**3, rhof
p_fluid_viscosity        =    0.001 # Pa*s,    mu_f
p_bulk_modulus           =   3.6e10 # Pa,      K_sg
p_shear_modulus          =      6e9 # Pa,      mu or G
p_biot_coefficient       =    0.778 # -,       alpha
p_isotropic_permeability = 1.90e-13 # m**2,    k
p_fluid_bulk_modulus     =      2e9 # Pa,      K_fl
p_external_force         =    -1000 # Pa,      F
p_iterations             =      200 # Number of 'infinite' iterations, N

p_bulk_density = (1.0 - p_porosity)*p_density + p_porosity*p_fluid_density
p_drained_bulk_modulus = -1.0*p_bulk_modulus * (p_biot_coefficient - 1.0)


p_mu = p_shear_modulus
p_lambda = p_drained_bulk_modulus - 2.0/3.0 * p_shear_modulus

gacc = 9.80665  # m/s
ymax = +1.0
ymin = -1.0  # m

# Height of column, m
H = ymax - ymin
L = H

# Drained Matrix Bulk Modulus, Pa
K_dr = (-1*p_bulk_modulus)*(p_biot_coefficient - 1)

# Confined compressibility of the Porous Medium
m_v = 1 / (K_dr + 4/3 * p_shear_modulus)

# Drained storage coefficient
#S = (p_biot_coefficient - p_porosity)/p_bulk_modulus + p_porosity/p_fluid_bulk_modulus
S = p_porosity/p_fluid_bulk_modulus + 3/(3*p_drained_bulk_modulus + 2*p_shear_modulus)
# Consolidation coefficient
#c_v = p_isotropic_permeability / (p_fluid_viscosity * (p_biot_coefficient**2 * m_v + S))
c = (p_isotropic_permeability/p_fluid_viscosity) / S

# Undrained Fluid Pressure Response
p_0 = (p_biot_coefficient*m_v)/(p_biot_coefficient**2*m_v+S) * p_external_force

# Undrained response of consolidation
u_z0 = (m_v - (p_biot_coefficient**2*m_v**2)/(p_biot_coefficient**2*m_v + S))*p_external_force*H

# Steady State Consolidation
u_z_inf = m_v*p_external_force*H 

# Drained Poisson Ratio
nu = (3*p_drained_bulk_modulus - 2*p_shear_modulus) / (2*(3*p_drained_bulk_modulus + p_shear_modulus))

# Undrained Bulk Modulus
p_undrained_bulk_modulus = p_drained_bulk_modulus + (p_fluid_bulk_modulus*(p_bulk_modulus - p_drained_bulk_modulus)**2) / \
                           (p_fluid_bulk_modulus*(p_bulk_modulus - p_drained_bulk_modulus) + p_porosity*p_bulk_modulus*(p_bulk_modulus - p_fluid_bulk_modulus))

# Undrained Poisson Ratio
nu_u = (3*p_undrained_bulk_modulus - 2*p_shear_modulus) / (2*(3*p_undrained_bulk_modulus + p_shear_modulus))

# Poroelastic Stress Coefficient
eta = (p_biot_coefficient*(1 - 2*nu))/(2*(1-nu))

# Different terms for displacement calculation
H_tm = K_dr + 4.0*p_shear_modulus/3.0
M    = p_bulk_modulus / (p_biot_coefficient - p_porosity*(1.0 - p_bulk_modulus/p_fluid_bulk_modulus))
K_ud = p_biot_coefficient**2.0 * M + K_dr
a    = 1.0 / H_tm
ai   = 1.0 / (K_ud + 4.0*p_shear_modulus/3.0)
Cf   = (p_isotropic_permeability/p_fluid_viscosity) / (a * p_biot_coefficient**2.0 + 1.0/M)

# Time steps
tsteps = numpy.arange(1.0, 5.01, 1.0) + 1 # sec

# Time interval
tR = 1.0 # sec

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to gravitational body forces (no initial stress).
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 200

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.trace_strain,
            "porosity": self.porosity,
            "density": self.density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "bulk_modulus": self.bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "isotropic_permeability": self.isotropic_permeability,
            "fluid_bulk_modulus": self.fluid_bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": self.zero_vector,
        }
        self.key = None
        return

    def getField(self, name, pts):
        if self.key is None:
            field = self.fields[name](pts)
        else:
            field = self.fields[name][self.key](pts)
        return field

    # Material Tests
    def porosity(self, locs):
        """
        Compute porosity field at locations.
        """
        (npts, dim) = locs.shape
        porosity = p_porosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return porosity

    def density(self, locs):
        """
        Compute density field at locations.
        """
        (npts, dim) = locs.shape
        density = p_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def fluid_density(self, locs):
        """
        Compute fluid density field at locations.
        """
        (npts, dim) = locs.shape
        fluid_density = p_fluid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_density

    def fluid_viscosity(self, locs):
        """
        Compute fluid viscosity field at locations.
        """
        (npts, dim) = locs.shape
        fluid_viscosity = p_fluid_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_viscosity

    def shear_modulus(self, locs):
        """
        Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_shear_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def bulk_modulus(self, locs):
        """
        Compute bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        bulk_modulus = p_bulk_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bulk_modulus

    def biot_coefficient(self, locs):
        """
        Compute biot coefficient field at locations.
        """
        (npts, dim) = locs.shape
        biot_coefficient = p_biot_coefficient * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coefficient  

    def isotropic_permeability(self, locs):
        """
        Compute isotropic permeability field at locations.
        """
        (npts, dim) = locs.shape
        isotropic_permeability = p_isotropic_permeability * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return isotropic_permeability

    def fluid_bulk_modulus(self, locs):
        """
        Compute bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        fluid_bulk_modulus = p_fluid_bulk_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return  fluid_bulk_modulus

    # Solution Tests    

    def pressure(self,locs):
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts), dtype=numpy.float64)
        z = locs[:,1] - ymax
        t_track = 0
#        for t in tsteps:
#            for i in numpy.arange(1,p_iterations+1):
#                pressure[t_track,:] += 4/numpy.pi * p_0 * 1/(2*i -1) * numpy.sin( ((2*i - 1)*numpy.pi) / (2*H) * z) * numpy.exp( -1*( ((2*i-1)*numpy.pi)/(2*H) )**2 * c_v * t) 
#            t_track += 1
        z_star = z/L
        t_star = (c*t) / (4*L**2)
        pressure[t_track,:] += (p_external_force * eta)/(p_shear_modulus*S) * self.F1(z_star, t_star)
        t_track += 1
        return pressure
    
    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        #strain = self.strain(locs)
        disp = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
        z = locs[:,1]  - ymax
        t_track = 0
        for t in tsteps:
#            term = lambda m, z,t : 1.0/(2.0*m+1.0)**2 \
#                   * numpy.exp( - ((2.0*m+1.0)*numpy.pi/(2.0*H))**2 * Cf * t ) \
#                   * numpy.cos( (2.0*m+1.0)*numpy.pi/(2.0*H) * z )
#            disp[t_track,:,1] = a*p_external_force*(H-z) - 8.0/(numpy.pi**2) * (a-ai)*(H)*p_external_force\
#                * numpy.sum( [term(m, z,t) for m in numpy.arange(0.0,30.0) ], axis=0)
#            t_track += 1
            z_star = z/L
            t_star = (c*t) / (4*L**2)
            disp[t_track,:,1] = ( (p_external_force * L * (1 - 2*nu_u))/(2*p_shear_modulus*(1 - nu_u)) )*(1 - z_star) + \
                                ( (p_external_force * L * (nu_u - nu)) / (2*p_shear_modulus*(1 - nu_u)*(1 - nu)) ) * \
                                self.F2(z_star, t_star)
            t_track += 1

        return disp

    def trace_strain(self,locs):
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        return trace_strain

    # Additional Functions

    def F1(self, z_star, t_star):
        F1 = 0
        for m in numpy.arange(1,self.ITERATIONS*2+1,2):
            F1 += 4/(m*numpy.pi) * numpy.sin((m*numpy.pi*z_star)/2)*numpy.exp(-m**2*numpy.pi**2*t_star)
        return F1

    def F2(self, z_star, t_star):
        F2 = 0
        for m in numpy.arange(1,self.ITERATIONS*2+1,2):
            F2 += 8/(m**2 * numpy.pi**2)  * numpy.cos( (m*numpy.pi*z_star)/2 ) * (1 - numpy.exp(-m**2 * numpy.pi**2 * t_star) )    
        return F2


    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def strain(self, locs):
        """
        Compute strain field at locations.
        """
        eyy = p_density * gacc * (locs[:, 1] - ymax) / (p_lambda + 2 * p_mu)
        exx = 0
        ezz = 0
        exy = 0

        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, :, 0] = exx
        strain[:, :, 1] = eyy
        strain[:, :, 2] = ezz
        strain[:, :, 3] = exy
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        syy = p_density * gacc * (locs[:, 1] - ymax)
        sxx = p_lambda / (p_lambda + 2 * p_mu) * syy
        sxy = 0.0
        szz = p_lambda / (2 * p_lambda + 2 * p_mu) * (sxx + syy)

        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0, :, 0] = sxx
        stress[0, :, 1] = syy
        stress[0, :, 2] = szz
        stress[0, :, 3] = sxy
        return stress

    def gacc(self, locs):
        """Compute gravitational acceleration at locations.
        """
        (npts, dim) = locs.shape
        gravacc = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        gravacc[0, :, 1] = -gacc
        return gravacc


# End of file
