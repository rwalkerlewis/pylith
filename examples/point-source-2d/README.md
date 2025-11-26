# Point Source Example (2D)

This example demonstrates how to use moment tensor point sources in
elastodynamic simulations with PyLith. Point sources are useful for
simulating earthquake sources without explicitly modeling a fault.

## Features

- Point moment tensor sources with user-defined:
  - Location (x, y coordinates)
  - Moment tensor components
  - Magnitude (seismic moment)
  - Source time function (Ricker wavelet, Gaussian, step, ramp, or custom)

## Simulations

### Step 1: Explosion Source

A simple isotropic (explosion) source demonstrating:
- Location at the center of the domain
- Isotropic moment tensor [Mxx=1, Myy=1, Mzz=1, Mxy=0]
- Ricker wavelet source time function

### Step 2: Strike-slip Source

A double-couple source typical of strike-slip earthquakes:
- Location offset from center
- Strike-slip moment tensor [Mxx=1, Myy=-1, Mxy=0]
- Gaussian source time function

### Step 3: Multiple Sources

Multiple point sources in the same simulation:
- Two sources at different locations
- Different source mechanisms
- Different source time functions

## Running the Simulations

```bash
# Step 1: Explosion source
pylith step01_explosion.cfg

# Step 2: Strike-slip source  
pylith step02_strikeslip.cfg

# Step 3: Multiple sources
pylith step03_multiple.cfg
```

## Configuration

The point sources are configured in the `[pylithapp.problem.sources]` section.
Key parameters include:

- `location`: Source coordinates [x, y] or [x, y, z]
- `moment_tensor`: Moment tensor components in Voigt notation
- `magnitude`: Seismic moment (N*m)
- `time_function`: Source time function type and parameters

See the individual step configuration files for detailed examples.

## Moment Tensor Convention

The moment tensor is specified in Voigt notation:
- 2D: [Mxx, Myy, Mxy]
- 3D: [Mxx, Myy, Mzz, Mxy, Mxz, Myz]

Common source mechanisms:
- Explosion: [1, 1, 1, 0, 0, 0] (isotropic)
- Strike-slip: [1, -1, 0, 0, 0, 0]
- Vertical dipole: [0, 0, 1, 0, 0, 0]

## Source Time Functions

Available source time functions:
- `SourceTimeStep`: Step (Heaviside) function
- `SourceTimeRamp`: Linear ramp function
- `SourceTimeGaussian`: Gaussian pulse
- `SourceTimeRicker`: Ricker wavelet (Mexican hat)
- `SourceTimeHistory`: User-defined from file
