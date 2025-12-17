#!/usr/bin/env python3
"""
Example demonstrating PyLith dynamic poroelasticity scales.

This example shows how to use the DynamicPoroelasticity scales
to nondimensionalize dynamic coupled mechanical-hydraulic problems.
"""

import numpy as np


def main():
    """Demonstrate dynamic poroelasticity scales."""
    
    # Import PyLith scales modules
    from pylith.scales.DynamicPoroelasticity import DynamicPoroelasticity
    from pylith.scales.ElasticityScales import ElasticityScales
    
    print("=" * 70)
    print("PyLith Dynamic Poroelasticity Scales Example")
    print("=" * 70)
    
    # =========================================================================
    # Example 1: Using the DynamicPoroelasticity convenience class
    # =========================================================================
    print("\nExample 1: Using DynamicPoroelasticity convenience class")
    print("-" * 70)
    
    # Create normalizer instance
    normalizer = DynamicPoroelasticity()
    normalizer._configure()
    
    # Display default scales
    print(f"Length scale:            {normalizer.getLengthScale():.3e} m")
    print(f"Displacement scale:      {normalizer.getDisplacementScale():.3e} m")
    print(f"Rigidity scale:          {normalizer.getRigidityScale():.3e} Pa")
    print(f"Time scale:              {normalizer.getTimeScale():.3e} s")
    print(f"                         ({normalizer.getTimeScale():.3f} s)")
    
    # Compute derived scales
    print("\nDerived scales:")
    print(f"Stress scale:            {ElasticityScales.getStressScale(normalizer):.3e} Pa")
    print(f"Strain scale:            {ElasticityScales.getStrainScale(normalizer):.3e}")
    print(f"Fluid pressure scale:    {ElasticityScales.getFluidPressureScale(normalizer):.3e} Pa")
    print(f"Body force scale:        {ElasticityScales.getBodyForceScale(normalizer):.3e} Pa/m")
    print(f"Density scale:           {ElasticityScales.getDensityScale(normalizer):.3e} kg/m³")
    print(f"Velocity scale:          {ElasticityScales.getVelocityScale(normalizer):.3e} m/s")
    print(f"Acceleration scale:      {ElasticityScales.getAccelerationScale(normalizer):.3e} m/s²")
    print(f"Viscosity scale:         {ElasticityScales.getViscosityScale(normalizer):.3e} Pa*s")
    print(f"Permeability scale:      {ElasticityScales.getPermeabilityScale(normalizer):.3e} m²")
    
    # =========================================================================
    # Example 2: Custom dynamic poroelasticity problem
    # =========================================================================
    print("\n\nExample 2: Seismic wave propagation in saturated medium")
    print("-" * 70)
    
    from pylith.scales.Scales import Scales
    
    # Create custom scales for seismic wave propagation
    # Problem: Earthquake waves propagating through saturated sedimentary basin
    custom_scales = Scales()
    
    # Physical parameters
    length_scale = 10.0e3  # 10 km (basin size)
    shear_wave_speed = 1.5e3  # 1.5 km/s (soft sediments)
    permeability = 1.0e-13  # m² (sedimentary rocks)
    viscosity = 1.0e-3  # Pa*s (water)
    rigidity = 5.0e9  # Pa (5 GPa for soft sediments)
    
    # Set scales
    ElasticityScales.setDynamicPoroelasticity(
        custom_scales,
        lengthScale=length_scale,
        velocityScale=shear_wave_speed,
        permeability=permeability,
        viscosity=viscosity,
        rigidity=rigidity
    )
    
    # Compute time scale (wave travel time)
    time_scale = length_scale / shear_wave_speed
    
    # Compute fluid diffusion time scale for comparison
    time_fluid = ElasticityScales.computePoroelasticityTimeScale(
        viscosity, permeability, length_scale, rigidity
    )
    
    print(f"Basin parameters:")
    print(f"  Basin size:             {length_scale/1e3:.1f} km")
    print(f"  Shear wave speed:       {shear_wave_speed/1e3:.2f} km/s")
    print(f"  Shear modulus:          {rigidity/1e9:.1f} GPa")
    print(f"  Permeability:           {permeability:.2e} m²")
    print(f"  Fluid viscosity:        {viscosity:.3f} Pa*s")
    
    print(f"\nTime scales:")
    print(f"  Wave propagation:       {time_scale:.3e} s ({time_scale:.2f} s)")
    print(f"  Fluid diffusion:        {time_fluid:.3e} s ({time_fluid/31557600:.1f} years)")
    print(f"  Ratio (diffusion/wave): {time_fluid/time_scale:.2e}")
    
    print(f"\n→ Fluid diffusion is {time_fluid/time_scale:.1e}x slower than wave propagation")
    print("  Undrained response dominates during seismic wave passage")
    
    # =========================================================================
    # Example 3: Biot slow wave analysis
    # =========================================================================
    print("\n\nExample 3: Biot slow wave propagation")
    print("-" * 70)
    
    # Parameters for Biot slow wave (diffusive wave)
    L = 100.0  # 100 m observation scale
    phi = 0.3  # Porosity
    k_perm = 1.0e-12  # m² permeability
    mu_fluid = 1.0e-3  # Pa*s
    K_s = 40.0e9  # Pa (solid grain bulk modulus)
    K_f = 2.2e9  # Pa (fluid bulk modulus)
    G = 10.0e9  # Pa (shear modulus)
    rho_s = 2650.0  # kg/m³ (solid density)
    rho_f = 1000.0  # kg/m³ (fluid density)
    
    # Compute composite properties
    rho = (1 - phi) * rho_s + phi * rho_f  # Bulk density
    
    # Biot-Willis coefficient
    alpha = 1.0 - (K_s / (3.0 * G))  # Approximation for this case
    
    # Skempton's coefficient (for typical values)
    B = 0.8
    
    # Fast P-wave speed (drained)
    K_drained = (2.0/3.0) * G  # Approximation
    c_p_fast = np.sqrt((K_drained + (4.0/3.0) * G) / rho)
    
    # Slow wave diffusion coefficient
    D_slow = (k_perm * G) / (mu_fluid * phi)
    
    # Characteristic frequencies
    omega_c = D_slow / (L**2)  # Critical frequency
    f_c = omega_c / (2 * np.pi)
    
    print(f"Poroelastic medium properties:")
    print(f"  Porosity:                  {phi:.2f}")
    print(f"  Permeability:              {k_perm:.2e} m²")
    print(f"  Shear modulus:             {G/1e9:.1f} GPa")
    print(f"  Biot-Willis coefficient:   {alpha:.2f}")
    print(f"  Skempton coefficient:      {B:.2f}")
    
    print(f"\nWave speeds:")
    print(f"  Fast P-wave:               {c_p_fast:.0f} m/s")
    print(f"  Slow wave diffusivity:     {D_slow:.3e} m²/s")
    
    print(f"\nCharacteristic scales:")
    print(f"  Length scale:              {L:.1f} m")
    print(f"  Critical frequency:        {f_c:.2f} Hz")
    print(f"  Critical period:           {1.0/f_c:.2f} s")
    
    print(f"\n→ For frequencies >> {f_c:.2f} Hz: undrained (elastic) response")
    print(f"→ For frequencies << {f_c:.2f} Hz: drained response")
    
    # =========================================================================
    # Example 4: Comparing static and dynamic poroelasticity
    # =========================================================================
    print("\n\nExample 4: Quasi-static vs Dynamic poroelasticity comparison")
    print("-" * 70)
    
    from pylith.scales.QuasistaticPoroelasticity import QuasistaticPoroelasticity
    from pylith.scales.DynamicElasticity import DynamicElasticity
    
    L = 5.0e3  # 5 km
    
    # Quasi-static poroelasticity (diffusion-controlled)
    quasi_poro = QuasistaticPoroelasticity()
    quasi_poro._configure()
    t_quasi = quasi_poro.getTimeScale()
    
    # Dynamic elasticity (wave propagation)
    dynamic_elast = DynamicElasticity()
    dynamic_elast._configure()
    t_dynamic = dynamic_elast.getTimeScale()
    
    # Dynamic poroelasticity (wave propagation with fluid coupling)
    dynamic_poro = DynamicPoroelasticity()
    dynamic_poro._configure()
    t_dynamic_poro = dynamic_poro.getTimeScale()
    
    print(f"For length scale L = {L/1e3:.0f} km:")
    print(f"\nTime scales:")
    print(f"  Quasi-static poroelasticity:  {t_quasi:.3e} s ({t_quasi/31557600:.1f} yr)")
    print(f"  Dynamic elasticity:           {t_dynamic:.3e} s ({t_dynamic:.2f} s)")
    print(f"  Dynamic poroelasticity:       {t_dynamic_poro:.3e} s ({t_dynamic_poro:.2f} s)")
    
    print(f"\nTime scale ratios:")
    print(f"  Quasi-static / Dynamic:       {t_quasi/t_dynamic:.2e}")
    print(f"  Dynamic poro / Dynamic elast: {t_dynamic_poro/t_dynamic:.2f}")
    
    print("\n→ Dynamic poroelasticity uses wave propagation time scale")
    print("→ Quasi-static poroelasticity uses diffusion time scale")
    print(f"→ Diffusion is {t_quasi/t_dynamic:.1e}x slower than wave propagation")
    
    # =========================================================================
    # Example 5: Liquefaction time scale analysis
    # =========================================================================
    print("\n\nExample 5: Soil liquefaction during earthquake shaking")
    print("-" * 70)
    
    # Loose sand parameters
    liquefaction_scales = Scales()
    
    L_soil = 10.0  # 10 m soil layer thickness
    k_sand = 1.0e-11  # m² (loose sand)
    mu_water = 1.0e-3  # Pa*s
    G_sand = 50.0e6  # Pa (50 MPa - very loose sand)
    c_s_sand = 150.0  # m/s (shear wave speed in loose sand)
    
    # Earthquake parameters
    f_earthquake = 2.0  # Hz (typical earthquake frequency)
    duration_shaking = 20.0  # seconds
    
    # Compute time scales
    t_wave = L_soil / c_s_sand  # Wave travel time
    t_diffusion = (mu_water * L_soil**2) / (k_sand * G_sand)  # Fluid diffusion
    period_earthquake = 1.0 / f_earthquake
    
    print(f"Loose sand layer properties:")
    print(f"  Layer thickness:           {L_soil:.1f} m")
    print(f"  Permeability:              {k_sand:.2e} m²")
    print(f"  Shear modulus:             {G_sand/1e6:.1f} MPa")
    print(f"  Shear wave speed:          {c_s_sand:.0f} m/s")
    
    print(f"\nEarthquake characteristics:")
    print(f"  Dominant frequency:        {f_earthquake:.1f} Hz")
    print(f"  Period:                    {period_earthquake:.2f} s")
    print(f"  Shaking duration:          {duration_shaking:.0f} s")
    
    print(f"\nTime scales:")
    print(f"  Wave travel time:          {t_wave:.3f} s")
    print(f"  Fluid diffusion time:      {t_diffusion:.1f} s")
    print(f"  Earthquake period:         {period_earthquake:.2f} s")
    
    print(f"\nLiquefaction analysis:")
    if t_diffusion > duration_shaking:
        print(f"→ Diffusion time ({t_diffusion:.1f} s) > Shaking duration ({duration_shaking:.0f} s)")
        print("  LIQUEFACTION RISK: Undrained conditions during shaking")
        print("  Pore pressure cannot dissipate quickly enough")
        print("  Effective stress can drop to near zero")
    else:
        print(f"→ Diffusion time ({t_diffusion:.1f} s) < Shaking duration ({duration_shaking:.0f} s)")
        print("  Lower liquefaction risk: Partial drainage during shaking")
    
    if t_wave << period_earthquake:
        print(f"→ Wave travel time ({t_wave:.3f} s) << Earthquake period ({period_earthquake:.2f} s)")
        print("  Quasi-static loading approximation valid within layer")
    
    # Dimensionless group
    Pi_liquefaction = t_diffusion / duration_shaking
    print(f"\nDimensionless liquefaction number: {Pi_liquefaction:.2f}")
    if Pi_liquefaction > 1.0:
        print("  π > 1: Undrained response → HIGH LIQUEFACTION POTENTIAL")
    else:
        print("  π < 1: Partially drained → REDUCED LIQUEFACTION POTENTIAL")
    
    print("\n" + "=" * 70)
    print("Example completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
