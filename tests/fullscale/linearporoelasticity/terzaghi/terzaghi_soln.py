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
# @file tests/fullscale/linearporoelasticity/terzaghi/terzaghi_soln.py
#
# @brief Analytical solution to Terzaghi's problem.
#
#          Uy=0
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         +1 Pa
#
# Dirichlet boundary conditions
#   Ux(+-1,0) = 0
#   Uy(x,+1) = 0
# Neumann boundary conditions
#   \tau_normal(x,+1) = 1*Pa

import numpy

# Physical properties
p_solid_density = 2500 # kg / m**3
p_fluid_density = 1000 # kg / m**3
p_fluid_dynamic_viscosity = 1.0 # Pa*s
p_G = 3.0 # Pa
p_K_u = 9.76 # Pa
p_alpha = 0.6 # -
p_M = 16.0 # Pa
p_isotropic_permeability = 1.5 # m**2

p_zmax = 10.0 # m
p_zmin = 0.0 # m
p_ymax = 10.0 # m
p_ymin = 0.0 # m
p_xmax = 10.0 # m
p_xmin = 0.0 # m
p_vertical_stress = 1.0 # Pa

# Height of column, m
H = p_zmax - p_zmin
L = H

p_K_d = p_K_u - p_alpha*p_alpha*p_M # Pa,      Cheng (B.5)
p_nu = (3.0*p_K_d - 2.0*p_G) / (2.0*(3.0*p_K_d + p_G)) # -,       Cheng (B.8)
p_nu_u = (3.0*p_K_u - 2.0*p_G) / (2.0*(3.0*p_K_u + p_G)) # -,       Cheng (B.9)
p_eta = (3.0*p_alpha*p_G) /(3.0*p_K_d + 4.0*p_G) #  -,       Cheng (B.11)
p_S = (3.0*p_K_u + 4.0*p_G) / (p_M*(3.0*p_K_d + 4.0*p_G)) # Pa^{-1}, Cheng (B.14)
p_c = (p_isotropic_permeability / p_fluid_dynamic_viscosity) / p_S # m^2 / s, Cheng (B.16)

# Time steps
tsteps = numpy.arange(0.0, 0.0057333334, 0.0028666667) # sec

# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """
    Analytical solution to Terzaghi's problem
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 500

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.trace_strain,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "undrained_bulk_modulus": self.undrained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "initial_amplitude": {
                "x_neg": self.zero_vector,
                "x_pos": self.zero_vector,
                "y_neg_neu": self.y_neg_neu,
                "y_neg_dir": self.zero_scalar,
                "y_pos": self.zero_vector,
            }
        }
        self.key = None
        return

    def getField(self, name, pts):
        if self.key is None:
            field = self.fields[name](pts)
        else:
            field = self.fields[name][self.key](pts)
        return field

    def zero_scalar(self, locs):
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]        
        return numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        return numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)

    def solid_density(self, locs):
        """
        Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        solid_density = p_solid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return solid_density

    def fluid_density(self, locs):
        """
        Compute fluid density field at locations.
        """
        (npts, dim) = locs.shape
        fluid_density = p_fluid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_density

    def shear_modulus(self, locs):
        """
        Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_G * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def fluid_viscosity(self, locs):
        """
        Compute fluid_viscosity field at locations.
        """
        (npts, dim) = locs.shape
        fluid_viscosity = p_fluid_dynamic_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_viscosity

    def undrained_bulk_modulus(self, locs):
        """
        Compute undrained bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        undrained_bulk_modulus = p_K_u * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return undrained_bulk_modulus

    def biot_coefficient(self, locs):
        """
        Compute biot coefficient field at locations.
        """
        (npts, dim) = locs.shape
        biot_coefficient = p_alpha * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coefficient

    def biot_modulus(self, locs):
        """
        Compute biot modulus field at locations.
        """
        (npts, dim) = locs.shape
        biot_modulus = p_M * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_modulus

    def isotropic_permeability(self, locs):
        """
        Compute isotropic permeability field at locations.
        """
        (npts, dim) = locs.shape
        isotropic_permeability = p_isotropic_permeability * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return isotropic_permeability

    def displacement(self, locs):
        """
        Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
        z = locs[:,1] 
        t_track = 0
        z_star = z/L

        for t in tsteps:
            t_star = (p_c*t) / ( (2*L)**2 )
            displacement[t_track,:,1] = ( (p_vertical_stress*L*(1.0 - 2.0*p_nu_u)) / (2.0*p_G*(1.0 - p_nu_u)) ) \
                                        * (1.0 - z_star) + ((p_vertical_stress*L*(p_nu_u - p_nu)) / (2.0*p_G*(1.0 - p_nu_u)*(1.0 - p_nu))) \
                                        * self.F2(z_star, t_star)
            t_track += 1

        return displacement

    def pressure(self, locs):
        """
        Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts), dtype=numpy.float64)
        z = locs[:,1] - p_ymax
        t_track = 0

        for t in tsteps:
            z_star = z/L
            t_star = (c*t) / (4*L**2)
            pressure[t_track,:] = (p_vertical_stress * p_eta) / (p_G * p_S) * self.F1(z_star, t_star)
            t_track += 1

        return pressure

    def trace_strain(self, locs):
        """
        Compute trace strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts), dtype=numpy.float64)
        z = locs[:,1] - p_ymax
        t_track = 0

        for t in tsteps:
            z_star = z/L
            t_star = (c*t) / (4*L**2)
            trace_strain[t_track,:] = -((p_vertical_stress*L*(1.0 - 2.0*p_nu_u)) / (2.0*p_G*(1.0 - p_nu_u)*L)) \
                                      + ((p_vertical_stress*L*(p_nu_u - p_nu)) / (2.0*p_G*(1.0 - p_nu_u)*(1.0 - p_nu)))*self.F3(z_star, t_star)
            t_track += 1

        return trace_strain

    # Series functions

    def F1(self, z_star, t_star):
        F1 = 0
        for m in numpy.arange(1,2*self.ITERATIONS+1,2):
            F1 += 4./(m*numpy.pi) * numpy.sin(0.5*m*numpy.pi*z_star)*numpy.exp( -(m*numpy.pi)**2*t_star)
        return F1

    def F2(self, z_star, t_star):
        F2 = 0
        for m in numpy.arange(1,2*self.ITERATIONS+1,2):
            F2 += ( 8. / (m*numpy.pi)**2 )  * numpy.cos(0.5*m*numpy.pi*z_star) * (1. - numpy.exp( -(m * numpy.pi)**2 * t_star) )
        return F2

    def F3(self, z_star, t_star):
        F3 = 0
        for m in numpy.arange(1,2*self.ITERATIONS+1,2):
            F3 += (-4.0 / (m*numpy.pi*L)) * numpy.sin(0.5*m*numpy.pi*z_star) * (1.0 - numpy.exp( -(m*numpy.pi)*t_star))
        return F3

    def strain(self, locs):
        """
        Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_zz = 0.0
        e_xy = 0.0

        strain = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, :, 0] = exx
        strain[:, :, 1] = eyy
        strain[:, :, 2] = ezz
        strain[:, :, 3] = exy
        return strain

    def stress(self, locs):
        """
        Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        p_poisson_ratio = (3*p_K_d - 2*p_G) / (2*(3*p_K_d + p_G))
        trace_strain = self.trace_strain(locs)
        pressure = self.pressure(locs)
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_xy = 0.0

        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:, :, 0] = ( (2*p_G*p_poisson_ratio) / (1 - 2*p_poisson_ratio) )*trace_strain + 2*p_G*e_xx - p_alpha*pressure
        stress[:, :, 1] = ( (2*p_G*p_poisson_ratio) / (1 - 2*p_poisson_ratio) )*trace_strain + 2*p_G*e_yy - p_alpha*pressure
        stress[:, :, 2] = ( (2*p_G*p_poisson_ratio) / (1 - 2*p_poisson_ratio) )*trace_strain - p_alpha*pressure
        stress[:, :, 3] = 2*p_G*e_xy
        return stress

    def y_neg_neu(self, locs):
        """Compute initial traction at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        traction = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:, :, 0] = 0.0
        traction[:, :, 1] = p_vertical_stress
        return traction
        
    def initial_displacement(self, locs):
        """
        Compute initial displacement at locations
        """
        (npts, dim) = locs.shape
        displacement = numpy.zeros((1, npts, dim), dtype=numpy.float64)
        z = locs[:,1] - p_ymax 
        z_star = z/L        
        
        displacement[0,:,1] = ( (p_vertical_stress*L*(1.0 - 2.0*p_nu_u) ) / (2.0*p_G*(1.0 - p_nu_u)) ) * (1.0 - z_star)
        return displacement

    def initial_pressure(self, locs):
        """
        Compute initial pressure at locations
        """
        (npts, dim) = locs.shape
        pressure = numpy.zeros((ntpts, npts), dtype=numpy.float64)
        z = locs[:,1] - p_ymax 
        z_star = z/L        

        pressure[0,:] = (p_vertical_stress * p_eta) / (p_G * p_S)
        
        return pressure
        
    def initial_trace_strain(self, locs):
        """
        Compute initial trace strain field at locations.
        """
        (npts, dim) = locs.shape

        trace_strain = numpy.zeros((1, npts), dtype=numpy.float64)
        z = locs[:,1] - p_ymax
        z_star = z/L

        trace_strain[0,:] = -((p_vertical_stress(1.0 - 2.0*p_nu_u)) / (2.0*p_G*(1.0 - p_nu_u)))

        return trace_strain        
        

# End of file
