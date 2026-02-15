#!/usr/bin/env python3
"""
Example demonstrating PyLith thermoporoelasticity scales.

This example shows how to use the QuasistaticThermoporoelasticity scales
to nondimensionalize a coupled thermal-mechanical-hydraulic problem.
"""

import numpy as np


def main():
    """Demonstrate thermoporoelasticity scales."""
    
    # Import PyLith scales modules
    from pylith.scales.QuasistaticThermoporoelasticity import QuasistaticThermoporoelasticity
    from pylith.scales.ElasticityScales import ElasticityScales
    
    print("=" * 70)
    print("PyLith Thermoporoelasticity Scales Example")
    print("=" * 70)
    
    # =========================================================================
    # Example 1: Using the QuasistaticThermoporoelasticity convenience class
    # =========================================================================
    print("\nExample 1: Using QuasistaticThermoporoelasticity convenience class")
    print("-" * 70)
    
    # Create normalizer instance
    normalizer = QuasistaticThermoporoelasticity()
    normalizer._configure()
    
    # Display default scales
    print(f"Length scale:            {normalizer.getLengthScale():.3e} m")
    print(f"Displacement scale:      {normalizer.getDisplacementScale():.3e} m")
    print(f"Rigidity scale:          {normalizer.getRigidityScale():.3e} Pa")
    print(f"Time scale:              {normalizer.getTimeScale():.3e} s")
    print(f"                         ({normalizer.getTimeScale()/31557600:.3e} years)")
    print(f"Temperature scale:       {normalizer.getTemperatureScale():.3e} K")
    
    # Compute derived scales
    print("\nDerived scales:")
    print(f"Stress scale:            {ElasticityScales.getStressScale(normalizer):.3e} Pa")
    print(f"Strain scale:            {ElasticityScales.getStrainScale(normalizer):.3e}")
    print(f"Fluid pressure scale:    {ElasticityScales.getFluidPressureScale(normalizer):.3e} Pa")
    print(f"Body force scale:        {ElasticityScales.getBodyForceScale(normalizer):.3e} Pa/m")
    print(f"Density scale:           {ElasticityScales.getDensityScale(normalizer):.3e} kg/m³")
    print(f"Velocity scale:          {ElasticityScales.getVelocityScale(normalizer):.3e} m/s")
    print(f"Viscosity scale:         {ElasticityScales.getViscosityScale(normalizer):.3e} Pa*s")
    print(f"Permeability scale:      {ElasticityScales.getPermeabilityScale(normalizer):.3e} m²")
    print(f"Heat flux scale:         {ElasticityScales.getHeatFluxScale(normalizer):.3e} W/m²")
    
    # =========================================================================
    # Example 2: Custom thermoporoelasticity scales
    # =========================================================================
    print("\n\nExample 2: Custom thermoporoelasticity problem")
    print("-" * 70)
    
    from pylith.scales.Scales import Scales
    
    # Create custom scales for a specific problem
    # Problem: Thermal pressurization on a fault zone at 10 km depth
    custom_scales = Scales()
    
    # Physical parameters
    length_scale = 10.0e3  # 10 km (depth to fault)
    permeability = 1.0e-17  # m^2 (very low permeability fault gouge)
    viscosity = 1.0e-3  # Pa*s (water viscosity)
    rigidity = 30.0e9  # Pa (30 GPa shear modulus)
    thermal_conductivity = 2.0  # W/(m*K) for fault zone
    density = 2600.0  # kg/m³
    specific_heat = 900.0  # J/(kg*K)
    
    # Set scales
    ElasticityScales.setQuasistaticThermoporoelasticity(
        custom_scales,
        lengthScale=length_scale,
        permeability=permeability,
        viscosity=viscosity,
        rigidity=rigidity,
        thermalConductivity=thermal_conductivity,
        density=density,
        specificHeat=specific_heat
    )
    
    # Compute individual time scales
    time_poro = ElasticityScales.computePoroelasticityTimeScale(
        viscosity, permeability, length_scale, rigidity
    )
    time_thermal = ElasticityScales.computeThermoelasticityTimeScale(
        length_scale, thermal_conductivity, density, specific_heat
    )
    time_combined = ElasticityScales.computeThermoporoelasticityTimeScale(
        length_scale, permeability, viscosity, rigidity,
        thermal_conductivity, density, specific_heat
    )
    
    print(f"Problem parameters:")
    print(f"  Depth to fault zone:    {length_scale/1e3:.1f} km")
    print(f"  Permeability:           {permeability:.2e} m²")
    print(f"  Fluid viscosity:        {viscosity:.3f} Pa*s")
    print(f"  Shear modulus:          {rigidity/1e9:.1f} GPa")
    print(f"  Thermal conductivity:   {thermal_conductivity:.1f} W/(m*K)")
    print(f"  Density:                {density:.0f} kg/m³")
    print(f"  Specific heat:          {specific_heat:.0f} J/(kg*K)")
    
    print(f"\nTime scales:")
    print(f"  Poroelastic:            {time_poro:.3e} s ({time_poro/31557600:.1f} years)")
    print(f"  Thermal:                {time_thermal:.3e} s ({time_thermal/31557600:.1f} years)")
    print(f"  Combined (minimum):     {time_combined:.3e} s ({time_combined/31557600:.1f} years)")
    
    if time_poro < time_thermal:
        print("\n→ Fluid diffusion controls the time scale (faster than heat diffusion)")
    else:
        print("\n→ Heat diffusion controls the time scale (faster than fluid diffusion)")
    
    # =========================================================================
    # Example 3: Nondimensionalization analysis for thermal pressurization
    # =========================================================================
    print("\n\nExample 3: Thermal pressurization analysis")
    print("-" * 70)
    
    # Dimensional problem parameters for co-seismic thermal pressurization
    L = 0.1  # 0.1 m (fault zone thickness)
    mu_fluid = 1.0e-3  # Pa*s
    kappa = 1.0e-18  # m² (very low permeability)
    G = 30.0e9  # Pa
    k_thermal = 1.5  # W/(m*K)
    rho = 2700.0  # kg/m³
    c = 800.0  # J/(kg*K)
    alpha_T = 3.0e-5  # 1/K thermal expansion
    Lambda = 0.6  # Skempton's coefficient
    
    # Compute time scales
    t_poro = (mu_fluid * L**2) / (kappa * G)
    t_thermal = (rho * c * L**2) / k_thermal
    
    # Compute pressure change from temperature rise
    Delta_T = 100.0  # K temperature rise from frictional heating
    Delta_p = Lambda * 3.0 * alpha_T * G * Delta_T  # Pressure change
    
    print(f"Fault zone parameters:")
    print(f"  Thickness:                 {L*1e3:.1f} mm")
    print(f"  Permeability:              {kappa:.2e} m²")
    print(f"  Skempton coefficient:      {Lambda:.2f}")
    print(f"  Temperature rise:          {Delta_T:.0f} K")
    
    print(f"\nTime scales:")
    print(f"  Poroelastic diffusion:     {t_poro:.3e} s ({t_poro:.2f} s)")
    print(f"  Thermal diffusion:         {t_thermal:.3e} s ({t_thermal:.2f} s)")
    print(f"  Ratio (thermal/poro):      {t_thermal/t_poro:.2f}")
    
    print(f"\nPressure response:")
    print(f"  Thermal pressure rise:     {Delta_p/1e6:.1f} MPa")
    
    # Compute nondimensional groups
    Pi_poro = (mu_fluid * L**2) / (kappa * G * t_poro)
    Pi_thermal = (rho * c * L**2) / (k_thermal * t_thermal)
    
    print(f"\nNondimensional groups:")
    print(f"  Poroelastic Biot number:   {Pi_poro:.2f}")
    print(f"  Thermal Biot number:       {Pi_thermal:.2f}")
    
    if t_thermal < t_poro:
        print("\n→ Heat diffuses faster: undrained heating regime")
        print("   Pore pressure rises quickly, then slowly dissipates")
    else:
        print("\n→ Fluid diffuses faster: drained heating regime")
        print("   Pore pressure dissipates during heating")
    
    # =========================================================================
    # Example 4: Comparing all time scales
    # =========================================================================
    print("\n\nExample 4: Multi-physics time scale comparison")
    print("-" * 70)
    
    from pylith.scales.QuasistaticElasticity import QuasistaticElasticity
    from pylith.scales.QuasistaticPoroelasticity import QuasistaticPoroelasticity
    from pylith.scales.QuasistaticThermoelasticity import QuasistaticThermoelasticity
    
    L = 50.0e3  # 50 km
    
    # Elasticity time scale (prescribed)
    elast_scales = QuasistaticElasticity()
    elast_scales._configure()
    t_elastic = elast_scales.getTimeScale()
    
    # Poroelasticity time scale (fluid diffusion)
    poro_scales = QuasistaticPoroelasticity()
    poro_scales._configure()
    t_poro = poro_scales.getTimeScale()
    
    # Thermoelasticity time scale (thermal diffusion)
    thermo_scales = QuasistaticThermoelasticity()
    thermo_scales._configure()
    t_thermo = thermo_scales.getTimeScale()
    
    # Thermoporoelasticity time scale (minimum of thermal and fluid diffusion)
    thermoporo_scales = QuasistaticThermoporoelasticity()
    thermoporo_scales._configure()
    t_thermoporo = thermoporo_scales.getTimeScale()
    
    print(f"For length scale L = {L/1e3:.0f} km:")
    print(f"\nTime scales:")
    print(f"  Elasticity (prescribed):          {t_elastic:.3e} s ({t_elastic/31557600:.1f} yr)")
    print(f"  Poroelasticity (fluid diffusion): {t_poro:.3e} s ({t_poro/31557600:.1f} yr)")
    print(f"  Thermoelasticity (heat diffusion):{t_thermo:.3e} s ({t_thermo/31557600:.1f} yr)")
    print(f"  Thermoporoelasticity (coupled):   {t_thermoporo:.3e} s ({t_thermoporo/31557600:.1f} yr)")
    
    print(f"\nTime scale ratios:")
    print(f"  Thermal / Poroelastic:            {t_thermo/t_poro:.2f}")
    print(f"  Thermoporo / Poroelastic:         {t_thermoporo/t_poro:.2f}")
    print(f"  Thermoporo / Thermal:             {t_thermoporo/t_thermo:.2f}")
    
    print("\n→ Thermoporoelastic time scale is controlled by the fastest")
    print("  diffusion process (minimum of fluid and thermal diffusion)")
    
    # =========================================================================
    # Example 5: Fault zone thermal pressurization scenario
    # =========================================================================
    print("\n\nExample 5: Realistic fault zone thermal pressurization")
    print("-" * 70)
    
    # Realistic parameters for a mature fault zone
    fault_scales = Scales()
    
    L_fault = 0.05  # 5 cm fault core thickness
    k_fault = 5.0e-19  # m² (very low permeability)
    mu_water = 1.0e-3  # Pa*s
    G_fault = 20.0e9  # Pa (damaged zone)
    k_therm = 1.2  # W/(m*K) (clay-rich gouge)
    rho_fault = 2400.0  # kg/m³
    c_fault = 850.0  # J/(kg*K)
    
    # Compute time scales
    t_poro_fault = (mu_water * L_fault**2) / (k_fault * G_fault)
    t_thermal_fault = (rho_fault * c_fault * L_fault**2) / k_therm
    
    # Seismic parameters
    V_slip = 1.0  # m/s slip velocity
    tau_friction = 50.0e6  # Pa (50 MPa friction stress)
    duration = 5.0  # seconds of slip
    
    # Energy released
    Q_per_area = tau_friction * V_slip * duration  # J/m²
    temperature_rise = Q_per_area / (rho_fault * c_fault * L_fault)
    
    print(f"Fault zone properties:")
    print(f"  Core thickness:              {L_fault*100:.1f} cm")
    print(f"  Permeability:                {k_fault:.2e} m²")
    print(f"  Thermal conductivity:        {k_therm:.1f} W/(m*K)")
    
    print(f"\nSeismic slip parameters:")
    print(f"  Slip velocity:               {V_slip:.1f} m/s")
    print(f"  Friction stress:             {tau_friction/1e6:.0f} MPa")
    print(f"  Slip duration:               {duration:.1f} s")
    print(f"  Total slip:                  {V_slip*duration:.1f} m")
    
    print(f"\nThermal response:")
    print(f"  Heat generated:              {Q_per_area/1e6:.1f} MJ/m²")
    print(f"  Temperature rise:            {temperature_rise:.0f} K")
    
    print(f"\nDiffusion time scales:")
    print(f"  Poroelastic (fluid):         {t_poro_fault:.3e} s ({t_poro_fault:.1f} s)")
    print(f"  Thermal (heat):              {t_thermal_fault:.3e} s ({t_thermal_fault:.1f} s)")
    
    if t_thermal_fault < duration:
        print("\n→ Heat diffuses during slip (drained thermal regime)")
    else:
        print("\n→ Heat is trapped during slip (undrained thermal regime)")
        
    if t_poro_fault < duration:
        print("→ Fluid pressure equilibrates during slip (drained fluid regime)")
    else:
        print("→ Fluid pressure trapped during slip (undrained fluid regime)")
        print("   THERMAL PRESSURIZATION EXPECTED!")
    
    print("\n" + "=" * 70)
    print("Example completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
