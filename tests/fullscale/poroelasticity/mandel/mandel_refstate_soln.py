# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2022 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#
# @file tests/fullscale/poroelasticity/mandel/mandel_refstate_soln.py
#
# @brief Analytical solution to Mandel's problem.
# Owing to the symmetry of the problem, we only need consider the quarter
# domain case.
#
#           -F
#        ----------
#        |        |
#  Ux=0  |        | P=0
#        |        |
#        |        |
#        ----------
#          Uy=0
#
# Dirichlet boundary conditions
#   Ux(0,y) = 0
#   Uy(x,0) = 0
# Neumann boundary conditions
#   \tau_normal(x,ymax) = -1*Pa

import numpy

# Physical properties
p_solid_density = 2500  # kg / m**3
p_fluid_density = 1000  # kg / m**3
p_fluid_viscosity = 1.0  # Pa*s
p_shear_modulus = 3.0  # Pa
p_solid_bulk_modulus = 10.0  # Pa
p_fluid_bulk_modulus = 8.0  # Pa
p_drained_bulk_modulus = 4.0  # Pa
# p_undrained_bulk_modulus = 2.6941176470588233 # Pa
p_biot_coefficient = 0.6  # -
p_porosity = 0.1
# p_biot_modulus = 4.705882352941176 # Pa
p_isotropic_permeability = 1.5  # m**2

p_biot_modulus = 1.0 / ( p_porosity / p_fluid_bulk_modulus + (p_biot_coefficient - p_porosity) / p_solid_bulk_modulus )
p_bulk_density = p_solid_density*(1.0 - p_porosity) + p_fluid_density * p_porosity

p_mu = p_shear_modulus
p_undrained_bulk_modulus = p_drained_bulk_modulus + p_biot_coefficient * p_biot_coefficient * p_biot_modulus  # Pa,      Cheng (B.5)
p_drained_lambda = p_drained_bulk_modulus - (2.0 * p_shear_modulus) / 3.0
p_undrained_lambda = p_undrained_bulk_modulus - (2.0 * p_shear_modulus) / 3.0
p_skempton_coefficient = (p_biot_coefficient * p_biot_modulus) / p_undrained_bulk_modulus
p_constant_stress_storage = p_undrained_bulk_modulus / (p_biot_modulus * (p_undrained_bulk_modulus - p_biot_coefficient**2 * p_biot_modulus))
p_drained_poisson_ratio = p_drained_lambda / (2.0 * (p_drained_lambda + p_mu))
p_undrained_poisson_ratio = p_undrained_lambda / (2.0 * (p_undrained_lambda + p_mu))
p_drained_youngs_modulus = 2.0 * p_shear_modulus * (1.0 * p_drained_poisson_ratio)
p_undrained_youngs_modulus = 2.0 * p_shear_modulus * (1.0 * p_undrained_poisson_ratio)

xmin = 0.0  # m
xmax = 10.0  # m
ymin = 0.0  # m
ymax = 1.0  # m


# Height of column, m
a = (xmax - xmin)
b = (ymax - ymin)


vertical_stress = 1.0  # Pa
F = vertical_stress*a

# p_biot_modulus = 1.0 / (p_porosity / p_fluid_bulk_modulus + (p_biot_coefficient - p_porosity) / p_solid_bulk_modulus)  # Pa
# p_undrained_bulk_modulus = p_drained_bulk_modulus + p_biot_coefficient * p_biot_coefficient * p_biot_modulus  # Pa,      Cheng (B.5)
# # p_drained_bulk_modulus = p_undrained_bulk_modulus - p_biot_coefficient*p_biot_coefficient*p_biot_modulus # Pa,      Cheng (B.5)
# p_drained_poisson_ratio = (3.0 * p_drained_bulk_modulus - 2.0 * p_shear_modulus) / (2.0 * (3.0 * p_drained_bulk_modulus + p_shear_modulus))  # -,       Cheng (B.8)
# p_undrained_poisson_ratio = (3.0 * p_undrained_bulk_modulus - 2.0 * p_shear_modulus) / (2.0 * (3.0 * p_undrained_bulk_modulus + p_shear_modulus))  # -,       Cheng (B.9)
eta = (3.0 * p_biot_coefficient * p_shear_modulus) / (3.0 * p_drained_bulk_modulus + 4.0 * p_shear_modulus)  # -,       Cheng (B.11)
S = (3.0 * p_undrained_bulk_modulus + 4.0 * p_shear_modulus) / (p_biot_modulus * (3.0 * p_drained_bulk_modulus + 4.0 * p_shear_modulus))  # Pa^{-1}, Cheng (B.14)
c = (p_isotropic_permeability / p_fluid_viscosity) / S  # m^2 / s, Cheng (B.16)
B = (3. * (p_undrained_poisson_ratio - p_drained_poisson_ratio)) / (p_biot_coefficient * (1. - 2. * p_drained_poisson_ratio) * (1. + p_undrained_poisson_ratio))

M_xx = p_undrained_lambda + 2.0 * p_shear_modulus
M_xy = p_undrained_lambda
M_xz = p_undrained_lambda
M_yx = p_undrained_lambda
M_yy = p_undrained_lambda + 2.0 * p_shear_modulus
M_yz = p_undrained_lambda
M_zx = p_undrained_lambda
M_zy = p_undrained_lambda
M_zz = p_undrained_lambda + 2.0 * p_shear_modulus

A1 = ((p_biot_coefficient**2) * M_yy - 2.0 * p_biot_coefficient*p_biot_coefficient*M_xy + p_biot_coefficient**2 * M_xx) / (p_biot_coefficient*M_xx - p_biot_coefficient*M_xy) \
     + (M_xx*M_yy - M_xy**2) / (p_biot_modulus*(p_biot_coefficient*M_xx - p_biot_coefficient*M_xy))
A2 = (p_biot_coefficient*M_xx - p_biot_coefficient*M_xy) / M_xx

# Time steps
ts = 0.0028666667  # sec
nts = 1
tsteps = numpy.arange(0.0, ts * nts, ts) + ts  # sec


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to Mandel's problem with Reference Stress and Strain
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 300
    EPS = 1e-25

    def __init__(self):
        self.fields = {
            "displacement": self.zero_vector_time,
            "pressure": self.zero_scalar,
            "trace_strain": self.trace_strain,
            "porosity": self.porosity,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "shear_modulus": self.shear_modulus,
            "undrained_bulk_modulus": self.undrained_bulk_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "reference_stress": self.stress,
            "reference_strain": self.strain,
            "cauchy_stress": self.zero_tensor,
            "cauchy_strain": self.zero_tensor,                               
            "initial_amplitude": {
                "x_neg": self.zero_vector,
                "x_pos": self.zero_scalar,
                "y_neg": self.zero_vector,
                "y_pos": self.initial_displacement,
            }
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name in "initial_amplitude":
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def zero_scalar(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, 1), dtype=numpy.float64)

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def zero_tensor(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)

    def zero_vector_time(self, locs):
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        return numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)

    def ones_scalar(self, locs):
        (npts, dim) = locs.shape
        return numpy.ones((1, npts, 1), dtype=numpy.float64)

    def ones_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.ones((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def solid_density(self, locs):
        """Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        solid_density = p_solid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return solid_density

    def fluid_density(self, locs):
        """Compute fluid density field at locations.
        """
        (npts, dim) = locs.shape
        fluid_density = p_fluid_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_density

    def porosity(self, locs):
        """Compute solid_density field at locations.
        """
        (npts, dim) = locs.shape
        porosity = p_porosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return porosity

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_shear_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def fluid_viscosity(self, locs):
        """Compute fluid_viscosity field at locations.
        """
        (npts, dim) = locs.shape
        fluid_viscosity = p_fluid_viscosity * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_viscosity

    def undrained_bulk_modulus(self, locs):
        """Compute undrained bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        undrained_bulk_modulus = p_undrained_bulk_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return undrained_bulk_modulus

    def drained_bulk_modulus(self, locs):
        """Compute undrained bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        drained_bulk_modulus = p_drained_bulk_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return drained_bulk_modulus

    def biot_coefficient(self, locs):
        """Compute biot coefficient field at locations.
        """
        (npts, dim) = locs.shape
        biot_coefficient = p_biot_coefficient * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coefficient

    def biot_modulus(self, locs):
        """Compute biot modulus field at locations.
        """
        (npts, dim) = locs.shape
        biot_modulus = p_biot_modulus * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_modulus

    def isotropic_permeability(self, locs):
        """Compute isotropic permeability field at locations.
        """
        (npts, dim) = locs.shape
        isotropic_permeability = p_isotropic_permeability * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return isotropic_permeability

    def y_pos(self, locs):
        """Compute initial traction at locations.

        :TODO: If this is the initial traction, then it should be a single time point (0).
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        t_track = 0

        displacement = numpy.zeros((ntpts, npts, self.SPACE_DIM), dtype=numpy.float64)
        displacement[:, :, 0] = 0.0
        for t in tsteps:
            displacement[t_track, :, 1] = F
            t_track += 1
        return displacement

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:
            A_x = 0.0
            B_x = 0.0

            for n in numpy.arange(1, self.ITERATIONS + 1, 1):
                a_n = zeroArray[n - 1]
                A_x += (numpy.sin(a_n) * numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                    numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))
                B_x += (numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                    numpy.sin((a_n * x) / a) * numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))

            displacement[t_track, :, 0] = ((F * p_drained_poisson_ratio) / (2.0 * p_shear_modulus * a) - (F * p_undrained_poisson_ratio) / (p_shear_modulus * a) * A_x) * x + F / p_shear_modulus * B_x
            displacement[t_track, :, 1] = (-1 * (F * (1.0 - p_drained_poisson_ratio)) / (2 * p_shear_modulus * a) + (F * (1 - p_undrained_poisson_ratio)) / (p_shear_modulus * a) * A_x) * z
            t_track += 1

        return displacement

    # def scaled_displacement(self, locs):
    #     """ Generate scaled y pos y displacement"""
    #     ts = 0.0028666667  # sec
    #     nts_scaled = 200
    #     tsteps_scaled = numpy.arange(0.0, ts * nts_scaled, ts)  # sec

    #     (npts, dim) = locs.shape
    #     ntpts = tsteps_scaled.shape[0]
    #     displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
    #     x = locs[:, 0]
    #     z = locs[:, 1]
    #     t_track = 0
    #     zeroArray = self.mandelZeros()

    #     for t in tsteps_scaled:
    #         A_x = 0.0
    #         B_x = 0.0

    #         for n in numpy.arange(1, self.ITERATIONS + 1, 1):
    #             a_n = zeroArray[n - 1]
    #             A_x += (numpy.sin(a_n) * numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
    #                 numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))
    #             B_x += (numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
    #                 numpy.sin((a_n * x) / a) * numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))

    #         displacement[t_track, :, 0] = ((F * p_drained_poisson_ratio) / (2.0 * p_shear_modulus * a) - (F * p_undrained_poisson_ratio) / (p_shear_modulus * a) * A_x) * x + F / p_shear_modulus * B_x
    #         displacement[t_track, :, 1] = (-1 * (F * (1.0 - p_drained_poisson_ratio)) / (2 * p_shear_modulus * a) + (F * (1 - p_undrained_poisson_ratio)) / (p_shear_modulus * a) * A_x) * z
    #         t_track += 1

    #     return displacement


    def pressure(self, locs):
        """Compute pressure field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:

            if t == 0.0:
                pressure[t_track, :] = (1. / (3. * a)) * (B * (1. + p_undrained_poisson_ratio)) * F
            else:
                p = 0.0
                for n in numpy.arange(1, self.ITERATIONS + 1, 1):
                    x_n = zeroArray[n - 1]
                    p += (numpy.sin(x_n) / (x_n - numpy.sin(x_n) * numpy.cos(x_n))) * \
                        (numpy.cos((x_n * x) / a) - numpy.cos(x_n)) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))
                pressure[t_track, :, 0] = ((2.0 * (F * B * (1.0 + p_undrained_poisson_ratio))) / (3.0 * a)) * p
            t_track += 1

        return pressure

    def trace_strain(self, locs):
        """Compute trace strain field at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:

            eps_A = 0.0
            eps_B = 0.0
            eps_C = 0.0

            for i in numpy.arange(1, self.ITERATIONS+1,1):
                x_n = zeroArray[i-1]
                eps_A += (x_n * numpy.exp( (-1.0*x_n*x_n*c*t)/(a*a)) * numpy.cos(x_n)*numpy.cos( (x_n*x)/a)) / (a * (x_n - numpy.sin(x_n)*numpy.cos(x_n)))
                eps_B += ( numpy.exp( (-1.0*x_n*x_n*c*t)/(a*a)) * numpy.sin(x_n)*numpy.cos(x_n)) / (x_n - numpy.sin(x_n)*numpy.cos(x_n))
                eps_C += ( numpy.exp( (-1.0*x_n*x_n*c*t)/(x_n*x_n)) * numpy.sin(x_n)*numpy.cos(x_n)) / (x_n - numpy.sin(x_n)*numpy.cos(x_n))

            trace_strain[t_track,:,0] = (F/p_shear_modulus)*eps_A + ( (F*p_drained_poisson_ratio)/(2.0*p_shear_modulus*a)) - eps_B/(p_shear_modulus*a) - (F*(1.0-p_drained_poisson_ratio))/(2/0*p_shear_modulus*a) + eps_C/(p_shear_modulus*a)
            t_track += 1

        return trace_strain

    # Series functions

    def mandelZeros(self):
        """Compute roots for analytical mandel problem solutions
        """
        zeroArray = numpy.zeros(self.ITERATIONS)
        x0 = 0

        for i in numpy.arange(1, self.ITERATIONS + 1, 1):
            a1 = x0 + numpy.pi/4
            a2 = x0 + numpy.pi/2 - 10000*2.2204e-16
            am = a1
            for j in numpy.arange(0, self.ITERATIONS, 1):
                y1 = numpy.tan(a1) - ((1.0 - p_drained_poisson_ratio) / (p_undrained_poisson_ratio - p_drained_poisson_ratio)) * a1
                y2 = numpy.tan(a2) - ((1.0 - p_drained_poisson_ratio) / (p_undrained_poisson_ratio - p_drained_poisson_ratio)) * a2
                am = (a1 + a2) / 2.0
                ym = numpy.tan(am) - (1 - p_drained_poisson_ratio) / (p_undrained_poisson_ratio - p_drained_poisson_ratio) * am
                if ((ym * y1) > 0):
                    a1 = am
                else:
                    a2 = am
                if (numpy.abs(y2) < self.EPS):
                    am = a2
            zeroArray[i - 1] = am
            x0 += numpy.pi
        return zeroArray

    def strain(self, locs):
        """Compute strain field at locations. 
           Uses reference stress, so strain should be zero.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        # pressure = self.pressure(locs)[:,:,0]
        # stress = self.stress(locs)

        # sxx = stress[:,:,0]
        # syy = stress[:,:,1]

        # exx = (M_yy*sxx - M_xy*syy + p_biot_coefficient*M_yy*pressure - p_biot_coefficient*M_xy*pressure) / (M_xx*M_yy - M_xy**2)
        # eyy = (M_xx*syy - M_xy*sxx + p_biot_coefficient*M_xx*pressure - p_biot_coefficient*M_xy*pressure) / (M_xx*M_yy - M_xy**2)        
        # ezz = 0.0
        # exy = 0.0

        strain = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        # strain[:, :, 0] = exx
        # strain[:, :, 1] = eyy
        # strain[:, :, 2] = ezz
        # strain[:, :, 3] = exy
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]

        ntpts = tsteps.shape[0]
        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        pressure = self.pressure(locs)[:,:,0]

        syy = self.sigma_zz(locs)
        sxx = numpy.zeros((ntpts, npts))
        sxy = numpy.zeros((ntpts, npts))
        # szz = p_undrained_lambda / (2 * p_undrained_lambda + 2 * p_mu) * (sxx + syy)
        szz = 0.5 * p_drained_lambda / (p_drained_lambda + p_shear_modulus) * (sxx + syy)
        # szz = p_drained_poisson_ratio * (sxx + syy) - p_biot_coefficient*(1.0 - 2.0*p_drained_poisson_ratio) * pressure

        stress[:,:,0] = sxx
        stress[:,:,1] = syy
        stress[:,:,2] = szz
        stress[:,:,3] = sxy

        return stress

    def initial_traction(self, locs):
        """Compute traction at locations.

        :TODO: If this is the initial traction, then it should be a single time point (0).
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        traction = numpy.zeros((npts), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        t = 0.0

        sigma_zz_A = 0.0
        sigma_zz_B = 0.0

        for i in numpy.arange(1, self.ITERATIONS + 1, 1):
            x_n = zeroArray[i - 1]
            sigma_zz_A += (numpy.sin(x_n) / (x_n - numpy.sin(x_n) * numpy.cos(x_n))) * \
                numpy.cos((x_n * x) / a) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))
            sigma_zz_B += ((numpy.sin(x_n) * numpy.cos(x_n)) / (x_n - numpy.sin(x_n) *
                                                                numpy.cos(x_n))) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))

        traction[:] = -(F / a) - ((2.0 * F) / a * A2 / A1) * sigma_zz_A + ((2.0 * F) / a) * sigma_zz_B

        return traction

    def initial_stress(self, locs):
        """Compute stress at locations.

        :TODO: If this is the initial traction, then it should be a single time point (0).
        """
        (npts, dim) = locs.shape

        syy = self.initial_traction(locs)

        stress = numpy.zeros((npts, self.TENSOR_SIZE), dtype=numpy.float64)
        sxx = numpy.zeros(npts)
        sxy = numpy.zeros(npts)
        # szz = p_undrained_lambda / (2 * p_undrained_lambda + 2 * p_mu) * (sxx + syy)
        szz = 0.5 * p_drained_lambda / (p_drained_lambda + p_shear_modulus) * (sxx + syy)
        # szz = p_drained_poisson_ratio * (sxx + syy) - p_biot_coefficient*(1.0 - 2.0*p_drained_poisson_ratio) * pressure

        stress[:,0] = sxx
        stress[:,1] = syy
        stress[:,2] = szz
        stress[:,3] = sxy

        return stress

    def initial_strain(self, locs):
        """Compute initial strain field at locations.
        """
        (npts, dim) = locs.shape

        exx = 0.0
        eyy = 0.0
        ezz = 0.0
        exy = 0.0

        strain = numpy.zeros((npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, 0] = exx
        strain[:, 1] = eyy
        strain[:, 2] = ezz
        strain[:, 3] = exy
        return strain

    def initial_displacement(self, locs):
        """Compute initial displacement at locations
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((1, npts, dim), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        t = 0.0
        A_x = 0.0
        B_x = 0.0

        for n in numpy.arange(1, self.ITERATIONS + 1, 1):
            a_n = zeroArray[n - 1]
            A_x += (numpy.sin(a_n) * numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))
            B_x += (numpy.cos(a_n) / (a_n - numpy.sin(a_n) * numpy.cos(a_n))) * \
                numpy.sin((a_n * x) / a) * numpy.exp(-1.0 * (a_n * a_n * c * t) / (a * a))

        displacement[t_track, :, 0] = ((F * p_drained_poisson_ratio) / (2.0 * p_shear_modulus * a) - (F * p_undrained_poisson_ratio) / (p_shear_modulus * a) * A_x) * x + F / p_shear_modulus * B_x
        displacement[t_track, :, 1] = (-1 * (F * (1.0 - p_drained_poisson_ratio)) / (2 * p_shear_modulus * a) + (F * (1 - p_undrained_poisson_ratio)) / (p_shear_modulus * a) * A_x) * z
        t_track += 1
        return displacement

    def initial_pressure(self, locs):
        """Compute initial pressure at locations
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((1, npts, 1), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        t = 0.0

        if t == 0.0:
            pressure[t_track, :] = (1. / (3. * a)) * (B * (1. + p_undrained_poisson_ratio)) * F
        else:
            p = 0.0
            for n in numpy.arange(1, self.ITERATIONS + 1, 1):
                x_n = zeroArray[n - 1]
                p += (numpy.sin(x_n) / (x_n - numpy.sin(x_n) * numpy.cos(x_n))) * \
                    (numpy.cos((x_n * x) / a) - numpy.cos(x_n)) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))
            pressure[t_track, :, 0] = ((2.0 * (F * B * (1.0 + p_undrained_poisson_ratio))) / (3.0 * a)) * p
        t_track += 1


        return pressure

    def initial_trace_strain(self, locs):
        """Compute initial trace strain field at locations.
        """
        (npts, dim) = locs.shape
        zeroArray = self.mandelZeros()
        trace_strain = numpy.zeros((1, npts), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t = 0.0

        trace_strain[0, :] = 0.0

        return trace_strain

    def sigma_zz(self, locs):
        """Compute traction at locations.
        """
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        traction = numpy.zeros((ntpts, npts), dtype=numpy.float64)
        x = locs[:, 0]
        z = locs[:, 1]
        t_track = 0
        zeroArray = self.mandelZeros()

        for t in tsteps:

            sigma_zz_A = 0.0
            sigma_zz_B = 0.0

            for i in numpy.arange(1, self.ITERATIONS + 1, 1):
                x_n = zeroArray[i - 1]
                sigma_zz_A += (numpy.sin(x_n) / (x_n - numpy.sin(x_n) * numpy.cos(x_n))) * \
                    numpy.cos((x_n * x) / a) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))
                sigma_zz_B += ((numpy.sin(x_n) * numpy.cos(x_n)) / (x_n - numpy.sin(x_n) *
                                                                    numpy.cos(x_n))) * numpy.exp(-1.0 * (x_n * x_n * c * t) / (a * a))

            traction[t_track, :] = -(F / a) - ((2.0 * F) / a * A2 / A1) * sigma_zz_A + ((2.0 * F) / a) * sigma_zz_B
            t_track += 1

        return traction


# End of file
