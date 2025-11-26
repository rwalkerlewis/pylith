# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""
Point moment tensor source for elastodynamics simulations.

A moment tensor point source represents a seismic source at a single point
in the domain. The source is characterized by:
- Location (x, y, z coordinates)
- Moment tensor components (Mxx, Myy, Mzz, Mxy, Mxz, Myz)
- Magnitude (scalar multiplier)
- Source time function (temporal evolution)

The equivalent body force from a moment tensor point source is:
    f_i(x, t) = -M_{ij,j} * delta(x - x_s) * S(t)

where M_{ij} is the moment tensor, delta is the Dirac delta function,
x_s is the source location, and S(t) is the source time function.
"""

from pylith.utils.PetscComponent import PetscComponent
import numpy as np


def validateLocation(value):
    """Validate source location (must be a list of 2 or 3 floats)."""
    msg = "Location must be a list of 2 or 3 coordinates."
    if not isinstance(value, list):
        raise ValueError(msg)
    if len(value) not in (2, 3):
        raise ValueError(msg)
    try:
        coords = list(map(float, value))
    except (TypeError, ValueError):
        raise ValueError(msg)
    return coords


def validateMomentTensor(value):
    """Validate moment tensor components (must be a list of 3, 4, or 6 floats)."""
    msg = "Moment tensor must be a list of 3 (2D), 4 (2D with shear), or 6 (3D) components."
    if not isinstance(value, list):
        raise ValueError(msg)
    if len(value) not in (3, 4, 6):
        raise ValueError(msg)
    try:
        components = list(map(float, value))
    except (TypeError, ValueError):
        raise ValueError(msg)
    return components


class PointMomentTensor(PetscComponent):
    """
    Point moment tensor source for elastodynamics simulations.

    This class implements a point source characterized by a moment tensor
    at a specified location with a specified time function.

    The moment tensor is specified in Voigt notation:
    - 2D: [Mxx, Myy, Mxy] (plane strain)
    - 3D: [Mxx, Myy, Mzz, Mxy, Mxz, Myz]

    Common moment tensor configurations:
    - Explosion: [1, 1, 1, 0, 0, 0] (isotropic)
    - Vertical dipole: [0, 0, 1, 0, 0, 0]
    - Double couple (strike-slip): [1, -1, 0, 0, 0, 0]
    """

    DOC_CONFIG = {
        "cfg": """
            # Point source with explosion mechanism
            [pylithapp.problem.sources.explosion]
            location = [0.0*km, -5.0*km, 0.0*km]
            moment_tensor = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
            magnitude = 1.0e+18*N*m

            time_function = pylith.sources.SourceTimeRicker
            time_function.origin_time = 0.0*s
            time_function.peak_frequency = 2.0*Hz
            time_function.delay = 0.5*s

            # Point source with double-couple (strike-slip) mechanism
            [pylithapp.problem.sources.strike_slip]
            location = [10.0*km, -10.0*km, 0.0*km]
            moment_tensor = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0]
            magnitude = 5.0e+17*N*m

            time_function = pylith.sources.SourceTimeGaussian
            time_function.origin_time = 0.5*s
            time_function.sigma = 0.2*s
        """
    }

    import pythia.pyre.inventory
    from pythia.pyre.units.length import meter, kilometer
    from pythia.pyre.units.SI import newton

    location = pythia.pyre.inventory.list("location", default=[0.0, 0.0, 0.0], validator=validateLocation)
    location.meta["tip"] = "Source location coordinates [x, y, z] (use units, e.g., [0.0*km, -5.0*km, 0.0*km])."

    momentTensor = pythia.pyre.inventory.list(
        "moment_tensor", default=[1.0, 1.0, 1.0, 0.0, 0.0, 0.0], validator=validateMomentTensor
    )
    momentTensor.meta["tip"] = (
        "Moment tensor components in Voigt notation. "
        "2D: [Mxx, Myy, Mxy], 3D: [Mxx, Myy, Mzz, Mxy, Mxz, Myz]. "
        "Values are normalized by magnitude."
    )

    magnitude = pythia.pyre.inventory.dimensional("magnitude", default=1.0e+18 * newton * meter)
    magnitude.meta["tip"] = "Seismic moment magnitude (in N*m)."

    from .SourceTimeFunction import SourceTimeRicker

    timeFunction = pythia.pyre.inventory.facility(
        "time_function", family="source_time_function", factory=SourceTimeRicker
    )
    timeFunction.meta["tip"] = "Source time function."

    # Optional: characteristic length for delta function regularization
    characteristicLength = pythia.pyre.inventory.dimensional(
        "characteristic_length", default=0.0 * meter
    )
    characteristicLength.meta["tip"] = (
        "Characteristic length for regularizing the delta function. "
        "If 0 or not set, uses local element size."
    )

    def __init__(self, name="pointmomenttensor"):
        """Constructor."""
        PetscComponent.__init__(self, name, facility="point_source")
        self._locationNd = None
        self._momentTensorNd = None
        self._magnitudeNd = None

    def preinitialize(self, problem):
        """Do pre-initialization setup."""
        from pylith.mpi.Communicator import mpi_is_root

        if mpi_is_root():
            self._info.log(
                f"Pre-initializing point moment tensor source '{self.aliases[-1]}'."
            )

        # Get scales for nondimensionalization
        scales = problem.scales
        lengthScale = scales.getLengthScale()
        rigidityScale = scales.getRigidityScale()
        timeScale = scales.getTimeScale()

        # Nondimensionalize location
        self._locationNd = np.array([float(x) for x in self.location])
        self._locationNd /= lengthScale

        # Nondimensionalize moment tensor
        # Moment tensor has units of [Force * Length] = [Rigidity * Length^3]
        momentScale = rigidityScale * lengthScale ** 3
        self._magnitudeNd = self.magnitude.value / momentScale

        # Normalize moment tensor components (dimensionless direction)
        self._momentTensorNd = np.array([float(m) for m in self.momentTensor])
        norm = np.sqrt(np.sum(self._momentTensorNd ** 2))
        if norm > 0:
            self._momentTensorNd /= norm

        # Initialize time function
        self.timeFunction.preinitialize(problem)

    def verifyConfiguration(self, mesh):
        """Verify configuration is valid for the given mesh."""
        from pylith.mpi.Communicator import mpi_is_root

        spaceDim = mesh.getCoordSys().getSpaceDim()

        # Verify location dimensions
        if len(self.location) != spaceDim:
            raise ValueError(
                f"Source location has {len(self.location)} components, "
                f"but mesh has {spaceDim} dimensions."
            )

        # Verify moment tensor dimensions
        mtLen = len(self.momentTensor)
        if spaceDim == 2:
            if mtLen not in (3, 4):
                raise ValueError(
                    f"For 2D problems, moment tensor must have 3 or 4 components, got {mtLen}."
                )
        elif spaceDim == 3:
            if mtLen != 6:
                raise ValueError(
                    f"For 3D problems, moment tensor must have 6 components, got {mtLen}."
                )

        if mpi_is_root():
            self._info.log(
                f"Point source '{self.aliases[-1]}' configuration verified for {spaceDim}D mesh."
            )

    def getLocation(self):
        """Get the source location (nondimensionalized).
        
        Returns:
            numpy.ndarray: Source location coordinates.
        """
        return self._locationNd

    def getMomentTensor(self):
        """Get the normalized moment tensor components.
        
        Returns:
            numpy.ndarray: Moment tensor components in Voigt notation.
        """
        return self._momentTensorNd

    def getMagnitude(self):
        """Get the source magnitude (nondimensionalized).
        
        Returns:
            float: Nondimensionalized magnitude.
        """
        return self._magnitudeNd

    def getMomentTensorFull(self, spaceDim):
        """Get the full moment tensor matrix.
        
        Args:
            spaceDim: Spatial dimension (2 or 3).
            
        Returns:
            numpy.ndarray: Full moment tensor as a matrix.
        """
        mt = self._momentTensorNd
        
        if spaceDim == 2:
            # 2D: [Mxx, Myy, Mxy] or [Mxx, Myy, Mzz, Mxy]
            if len(mt) == 3:
                return np.array([
                    [mt[0], mt[2]],
                    [mt[2], mt[1]]
                ])
            else:  # len == 4, ignore Mzz
                return np.array([
                    [mt[0], mt[3]],
                    [mt[3], mt[1]]
                ])
        else:
            # 3D: [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
            return np.array([
                [mt[0], mt[3], mt[4]],
                [mt[3], mt[1], mt[5]],
                [mt[4], mt[5], mt[2]]
            ])

    def evaluateTimeFunction(self, t):
        """Evaluate the source time function at time t.
        
        Args:
            t: Time value.
            
        Returns:
            float: Source time function value.
        """
        return self.timeFunction.evaluate(t)

    def computeEquivalentForce(self, x, t, spaceDim, charLength=None):
        """
        Compute the equivalent body force at point x and time t.
        
        For a moment tensor point source, the equivalent body force is:
            f_i = -M_{ij} * dphi/dx_j * S(t)
        
        where dphi/dx_j is the spatial derivative of a regularized delta function.
        
        Args:
            x: Coordinates of the evaluation point.
            t: Time value.
            spaceDim: Spatial dimension.
            charLength: Characteristic length for delta function regularization.
            
        Returns:
            numpy.ndarray: Equivalent force vector at point x.
        """
        # Get regularization length
        if charLength is None or charLength <= 0:
            charLength = 1.0  # Default, should be set based on element size
        
        # Compute distance from source
        xs = self._locationNd[:spaceDim]
        r = x[:spaceDim] - xs
        dist = np.sqrt(np.sum(r ** 2))
        
        # Regularized delta function derivative (Gaussian approximation)
        # phi(r) = (1 / (sqrt(2*pi) * sigma)^dim) * exp(-r^2 / (2*sigma^2))
        # dphi/dr = -r / sigma^2 * phi(r)
        sigma = charLength
        if spaceDim == 2:
            norm = 1.0 / (2.0 * np.pi * sigma ** 2)
        else:
            norm = 1.0 / ((2.0 * np.pi) ** 1.5 * sigma ** 3)
        
        phi = norm * np.exp(-dist ** 2 / (2.0 * sigma ** 2))
        
        # Spatial derivative of delta function
        if dist > 1e-10 * sigma:
            dphi_dr = -r / (sigma ** 2) * phi
        else:
            dphi_dr = np.zeros(spaceDim)
        
        # Get moment tensor matrix
        M = self.getMomentTensorFull(spaceDim) * self._magnitudeNd
        
        # Evaluate time function
        stf = self.evaluateTimeFunction(t)
        
        # Compute equivalent force: f_i = -M_{ij} * dphi/dx_j * S(t)
        force = -np.dot(M, dphi_dr) * stf
        
        return force

    def getBodyForceFieldFunction(self, mesh, t):
        """
        Get a function that computes the body force field at any point.
        
        This can be used to set up the body force contribution in the
        elasticity equation.
        
        Args:
            mesh: The finite element mesh.
            t: Current time.
            
        Returns:
            callable: Function that takes coordinates and returns force vector.
        """
        spaceDim = mesh.getCoordSys().getSpaceDim()
        
        # Estimate characteristic length from mesh (simplified)
        charLength = self.characteristicLength.value
        if charLength <= 0:
            # Use a default based on expected element size
            # In practice, this should be computed from the mesh
            charLength = 0.1  # nondimensionalized
        
        def forceFunction(x):
            return self.computeEquivalentForce(x, t, spaceDim, charLength)
        
        return forceFunction

    def _configure(self):
        """Set members based on inventory."""
        PetscComponent._configure(self)


def point_source():
    """Factory associated with PointMomentTensor."""
    return PointMomentTensor()


# End of file
