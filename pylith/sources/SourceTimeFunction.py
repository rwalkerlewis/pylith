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
Source time functions for point sources in elastodynamics.

Provides various source time functions (STF) including:
- Step function
- Ramp function
- Gaussian function
- Ricker wavelet
- User-defined time history
"""

from pylith.utils.PetscComponent import PetscComponent
import numpy as np


class SourceTimeFunction(PetscComponent):
    """
    Abstract base class for source time functions.
    
    A source time function S(t) modulates the amplitude of a point source over time.
    The total source contribution is: f(x, t) = M * delta(x - x_s) * S(t)
    where M is the source magnitude/moment tensor and x_s is the source location.
    """

    import pythia.pyre.inventory
    from pythia.pyre.units.time import second

    originTime = pythia.pyre.inventory.dimensional("origin_time", default=0.0 * second)
    originTime.meta["tip"] = "Origin time for the source time function."

    def __init__(self, name="sourcetimefunction"):
        """Constructor."""
        PetscComponent.__init__(self, name, facility="source_time_function")

    def preinitialize(self, problem):
        """Do pre-initialization setup."""
        self._originTimeValue = self.originTime.value

    def evaluate(self, t):
        """
        Evaluate the source time function at time t.
        
        Args:
            t: Time value (in seconds, dimensional).
            
        Returns:
            float: Value of the source time function at time t.
        """
        raise NotImplementedError("Implement evaluate() in derived class.")

    def evaluateDerivative(self, t):
        """
        Evaluate the time derivative of the source time function at time t.
        
        Args:
            t: Time value (in seconds, dimensional).
            
        Returns:
            float: Value of the time derivative at time t.
        """
        raise NotImplementedError("Implement evaluateDerivative() in derived class.")


class SourceTimeStep(SourceTimeFunction):
    """
    Step source time function: S(t) = H(t - t0) * magnitude
    
    where H is the Heaviside step function, t0 is the origin time,
    and magnitude is the amplitude of the step.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.sources.my_source.time_function]
            origin_time = 0.0*s
            magnitude = 1.0
        """
    }

    import pythia.pyre.inventory

    magnitude = pythia.pyre.inventory.float("magnitude", default=1.0)
    magnitude.meta["tip"] = "Magnitude of the step function."

    def __init__(self, name="sourcetimestep"):
        """Constructor."""
        SourceTimeFunction.__init__(self, name)

    def evaluate(self, t):
        """Evaluate the step function at time t."""
        if t >= self._originTimeValue:
            return self.magnitude
        return 0.0

    def evaluateDerivative(self, t):
        """Evaluate the derivative of the step function (delta function approximation)."""
        # The derivative of a step function is a delta function
        # For numerical purposes, return 0 (the derivative is infinite at t0)
        return 0.0


def source_time_function_step():
    """Factory associated with SourceTimeStep."""
    return SourceTimeStep()


class SourceTimeRamp(SourceTimeFunction):
    """
    Ramp source time function: linear ramp over a specified rise time.
    
    S(t) = 0                                       for t < t0
    S(t) = magnitude * (t - t0) / rise_time        for t0 <= t < t0 + rise_time
    S(t) = magnitude                               for t >= t0 + rise_time
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.sources.my_source.time_function]
            origin_time = 0.0*s
            rise_time = 1.0*s
            magnitude = 1.0
        """
    }

    import pythia.pyre.inventory
    from pythia.pyre.units.time import second

    riseTime = pythia.pyre.inventory.dimensional(
        "rise_time", default=1.0 * second, validator=pythia.pyre.inventory.greater(0.0 * second)
    )
    riseTime.meta["tip"] = "Rise time for the ramp function."

    magnitude = pythia.pyre.inventory.float("magnitude", default=1.0)
    magnitude.meta["tip"] = "Final magnitude of the ramp function."

    def __init__(self, name="sourcetimeramp"):
        """Constructor."""
        SourceTimeFunction.__init__(self, name)

    def preinitialize(self, problem):
        """Do pre-initialization setup."""
        SourceTimeFunction.preinitialize(self, problem)
        self._riseTimeValue = self.riseTime.value

    def evaluate(self, t):
        """Evaluate the ramp function at time t."""
        t0 = self._originTimeValue
        tr = self._riseTimeValue
        
        if t < t0:
            return 0.0
        elif t < t0 + tr:
            return self.magnitude * (t - t0) / tr
        else:
            return self.magnitude

    def evaluateDerivative(self, t):
        """Evaluate the derivative of the ramp function at time t."""
        t0 = self._originTimeValue
        tr = self._riseTimeValue
        
        if t0 <= t < t0 + tr:
            return self.magnitude / tr
        return 0.0


def source_time_function_ramp():
    """Factory associated with SourceTimeRamp."""
    return SourceTimeRamp()


class SourceTimeGaussian(SourceTimeFunction):
    """
    Gaussian source time function.
    
    S(t) = magnitude * exp(-((t - t0 - center_time) / sigma)^2 / 2)
    
    This is useful for smooth impulsive sources.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.sources.my_source.time_function]
            origin_time = 0.0*s
            center_time = 0.5*s
            sigma = 0.1*s
            magnitude = 1.0
        """
    }

    import pythia.pyre.inventory
    from pythia.pyre.units.time import second

    centerTime = pythia.pyre.inventory.dimensional("center_time", default=0.5 * second)
    centerTime.meta["tip"] = "Center time of the Gaussian (time of peak amplitude)."

    sigma = pythia.pyre.inventory.dimensional(
        "sigma", default=0.1 * second, validator=pythia.pyre.inventory.greater(0.0 * second)
    )
    sigma.meta["tip"] = "Standard deviation (width) of the Gaussian."

    magnitude = pythia.pyre.inventory.float("magnitude", default=1.0)
    magnitude.meta["tip"] = "Peak magnitude of the Gaussian."

    def __init__(self, name="sourcetimegaussian"):
        """Constructor."""
        SourceTimeFunction.__init__(self, name)

    def preinitialize(self, problem):
        """Do pre-initialization setup."""
        SourceTimeFunction.preinitialize(self, problem)
        self._centerTimeValue = self.centerTime.value
        self._sigmaValue = self.sigma.value

    def evaluate(self, t):
        """Evaluate the Gaussian function at time t."""
        t0 = self._originTimeValue
        tc = self._centerTimeValue
        sig = self._sigmaValue
        
        if t < t0:
            return 0.0
        
        tau = (t - t0 - tc) / sig
        return self.magnitude * np.exp(-tau * tau / 2.0)

    def evaluateDerivative(self, t):
        """Evaluate the derivative of the Gaussian function at time t."""
        t0 = self._originTimeValue
        tc = self._centerTimeValue
        sig = self._sigmaValue
        
        if t < t0:
            return 0.0
        
        tau = (t - t0 - tc) / sig
        return -self.magnitude * tau / sig * np.exp(-tau * tau / 2.0)


def source_time_function_gaussian():
    """Factory associated with SourceTimeGaussian."""
    return SourceTimeGaussian()


class SourceTimeRicker(SourceTimeFunction):
    """
    Ricker wavelet source time function (Mexican hat wavelet).
    
    S(t) = magnitude * (1 - 2 * pi^2 * f^2 * (t - t0 - delay)^2) * 
           exp(-pi^2 * f^2 * (t - t0 - delay)^2)
    
    The Ricker wavelet is the second derivative of a Gaussian and is commonly
    used in seismic wave propagation simulations because it has a well-defined
    peak frequency.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.sources.my_source.time_function]
            origin_time = 0.0*s
            peak_frequency = 1.0*Hz
            delay = 1.0*s
            magnitude = 1.0
        """
    }

    import pythia.pyre.inventory
    from pythia.pyre.units.time import second
    from pythia.pyre.units.SI import hertz

    peakFrequency = pythia.pyre.inventory.dimensional(
        "peak_frequency", default=1.0 * hertz, validator=pythia.pyre.inventory.greater(0.0 * hertz)
    )
    peakFrequency.meta["tip"] = "Peak (central) frequency of the Ricker wavelet."

    delay = pythia.pyre.inventory.dimensional("delay", default=1.0 * second)
    delay.meta["tip"] = "Time delay to shift the wavelet (to avoid starting with non-zero amplitude)."

    magnitude = pythia.pyre.inventory.float("magnitude", default=1.0)
    magnitude.meta["tip"] = "Peak magnitude of the Ricker wavelet."

    def __init__(self, name="sourcetimericker"):
        """Constructor."""
        SourceTimeFunction.__init__(self, name)

    def preinitialize(self, problem):
        """Do pre-initialization setup."""
        SourceTimeFunction.preinitialize(self, problem)
        self._peakFrequencyValue = self.peakFrequency.value
        self._delayValue = self.delay.value

    def evaluate(self, t):
        """Evaluate the Ricker wavelet at time t."""
        t0 = self._originTimeValue
        f = self._peakFrequencyValue
        delay = self._delayValue
        
        if t < t0:
            return 0.0
        
        tau = t - t0 - delay
        pi2_f2_tau2 = (np.pi * f * tau) ** 2
        return self.magnitude * (1.0 - 2.0 * pi2_f2_tau2) * np.exp(-pi2_f2_tau2)

    def evaluateDerivative(self, t):
        """Evaluate the derivative of the Ricker wavelet at time t."""
        t0 = self._originTimeValue
        f = self._peakFrequencyValue
        delay = self._delayValue
        
        if t < t0:
            return 0.0
        
        tau = t - t0 - delay
        pi2_f2 = (np.pi * f) ** 2
        pi2_f2_tau2 = pi2_f2 * tau * tau
        
        # Derivative of Ricker wavelet
        return self.magnitude * tau * pi2_f2 * (2.0 * pi2_f2_tau2 - 3.0) * np.exp(-pi2_f2_tau2)


def source_time_function_ricker():
    """Factory associated with SourceTimeRicker."""
    return SourceTimeRicker()


class SourceTimeHistory(SourceTimeFunction):
    """
    User-defined source time function from a time history file.
    
    The time history file should have two columns: time and amplitude.
    Linear interpolation is used between data points.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.sources.my_source.time_function]
            origin_time = 0.0*s
            time_history = my_source_history.timedb
        """
    }

    import pythia.pyre.inventory

    timeHistoryFile = pythia.pyre.inventory.str("time_history", default="")
    timeHistoryFile.meta["tip"] = "Filename for time history data."

    def __init__(self, name="sourcetimehistory"):
        """Constructor."""
        SourceTimeFunction.__init__(self, name)
        self._times = None
        self._values = None

    def preinitialize(self, problem):
        """Do pre-initialization setup."""
        SourceTimeFunction.preinitialize(self, problem)
        
        if not self.timeHistoryFile:
            raise ValueError("Time history filename is required for SourceTimeHistory.")
        
        self._loadTimeHistory()

    def _loadTimeHistory(self):
        """Load the time history from file."""
        import os
        
        if not os.path.exists(self.timeHistoryFile):
            raise IOError(f"Time history file not found: {self.timeHistoryFile}")
        
        data = np.loadtxt(self.timeHistoryFile, comments=["#", "//"])
        if data.ndim == 1:
            data = data.reshape(-1, 2)
        
        self._times = data[:, 0]
        self._values = data[:, 1]

    def evaluate(self, t):
        """Evaluate the time history at time t using linear interpolation."""
        t0 = self._originTimeValue
        
        if t < t0:
            return 0.0
        
        t_rel = t - t0
        return np.interp(t_rel, self._times, self._values)

    def evaluateDerivative(self, t):
        """Evaluate the derivative using finite differences."""
        dt = 1e-6  # Small time step for numerical derivative
        return (self.evaluate(t + dt) - self.evaluate(t - dt)) / (2.0 * dt)


def source_time_function_history():
    """Factory associated with SourceTimeHistory."""
    return SourceTimeHistory()


# End of file
