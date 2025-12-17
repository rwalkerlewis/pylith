// =================================================================================================
// This code is part of SpatialData, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/spatialdata).
//
// Copyright (c) 2010-2025, University of California, Davis and the SpatialData Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file modulesrc/units/Scales.i
 *
 * @brief SWIG interface for C++ Scales object.
 */

namespace pylith {
    namespace scales {
        class ElasticityScales {
public:

            // PUBLIC METHODS /////////////////////////////////////////////////

            /** Set defaults scales for quasi-static elasticity.
             *
             * @param[inout] Scales for nondimensionalization.
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] timeScale Default time scale in seconds.
             *
             */
            static
            void setQuasistaticElasticity(pylith::scales::Scales* scales,
                                          const double lengthScale=100.0e+3,
                                          const double timeScale=31557600.0e+2);

            /** Set defaults scales for dynamic elasticity.
             *
             * @param[inout] Scales for nondimensionalization.
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] velocityScale Default velocity scale in m/s.
             *
             */
            static
            void setDynamicElasticity(pylith::scales::Scales* scales,
                                      const double lengthScale=100.0e+3,
                                      const double velocityScale=3.0e+3);

            /** Set defaults scales for dynamic poroelasticity.
             *
             * @param[inout] Scales for nondimensionalization.
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] velocityScale Default velocity scale in m/s.
             * @param[in] permeability Default permeability scale in m^2.
             * @param[in] viscosity Default viscosity scale in Pa*s.
             * @param[in] rigidity Default rigidity scale in Pa.
             */
            static
            void setDynamicPoroelasticity(pylith::scales::Scales* scales,
                                          const double lengthScale=100.0e+3,
                                          const double velocityScale=3.0e+3,
                                          const double permeability=1.0e-12,
                                          const double viscosity=1.0e-3,
                                          const double rigidity=25.0e+9);

            /** Set defaults scales for quasi-static poroelasticity.
             *
             * @param[inout] Scales for nondimensionalization.
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] permeability Default permeability scale in m^2.
             * @param[in] viscosity Default viscosity scale in Pa*s.
             * @param[in] rigidity Default rigidity scale in Pa.
             */
            static
            void setQuasistaticPoroelasticity(pylith::scales::Scales* scales,
                                              const double lengthScale=100.0e+3,
                                              const double permeability=1.0e-12,
                                              const double viscosity=1.0e-3,
                                              const double rigidity=25.0e+9);

            /** Set time scale for poroelasticity.
             *
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] permeability Default permeability scale in m^2.
             * @param[in] viscosity Default viscosity scale in Pa*s.
             * @param[in] rigidity Default rigidity scale in Pa.
             *
             * @returns Time scale in seconds.
             */
            static
            double computePoroelasticityTimeScale(const double viscosity,
                                                  const double permeability,
                                                  const double length,
                                                  const double rigidity);

            /** Set defaults scales for quasi-static thermoelasticity.
             *
             * @param[inout] Scales for nondimensionalization.
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] thermalConductivity Default thermal conductivity in W/(m*K).
             * @param[in] density Default density in kg/m^3.
             * @param[in] specificHeat Default specific heat capacity in J/(kg*K).
             */
            static
            void setQuasistaticThermoelasticity(pylith::scales::Scales* scales,
                                                const double lengthScale=100.0e+3,
                                                const double thermalConductivity=2.5,
                                                const double density=2500.0,
                                                const double specificHeat=1000.0);

            /** Compute time scale for thermoelasticity.
             *
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] thermalConductivity Default thermal conductivity in W/(m*K).
             * @param[in] density Default density in kg/m^3.
             * @param[in] specificHeat Default specific heat capacity in J/(kg*K).
             *
             * @returns Time scale in seconds.
             */
            static
            double computeThermoelasticityTimeScale(const double lengthScale,
                                                    const double thermalConductivity,
                                                    const double density,
                                                    const double specificHeat);

            /** Set defaults scales for quasi-static thermoporoelasticity.
             *
             * @param[inout] Scales for nondimensionalization.
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] permeability Default permeability scale in m^2.
             * @param[in] viscosity Default viscosity scale in Pa*s.
             * @param[in] rigidity Default rigidity scale in Pa.
             * @param[in] thermalConductivity Default thermal conductivity in W/(m*K).
             * @param[in] density Default density in kg/m^3.
             * @param[in] specificHeat Default specific heat capacity in J/(kg*K).
             */
            static
            void setQuasistaticThermoporoelasticity(pylith::scales::Scales* scales,
                                                    const double lengthScale=100.0e+3,
                                                    const double permeability=1.0e-12,
                                                    const double viscosity=1.0e-3,
                                                    const double rigidity=25.0e+9,
                                                    const double thermalConductivity=2.5,
                                                    const double density=2500.0,
                                                    const double specificHeat=1000.0);

            /** Set defaults scales for heat transfer.
             *
             * @param[inout] Scales for nondimensionalization.
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] thermalConductivity Default thermal conductivity in W/(m*K).
             * @param[in] density Default density in kg/m^3.
             * @param[in] specificHeat Default specific heat capacity in J/(kg*K).
             */
            static
            void setHeat(pylith::scales::Scales* scales,
                         const double lengthScale=100.0e+3,
                         const double thermalConductivity=2.5,
                         const double density=2500.0,
                         const double specificHeat=1000.0);

            /** Compute time scale for thermoporoelasticity.
             *
             * Chooses minimum of poroelastic and thermal diffusion time scales.
             *
             * @param[in] lengthScale Default length scale in meters.
             * @param[in] permeability Default permeability scale in m^2.
             * @param[in] viscosity Default viscosity scale in Pa*s.
             * @param[in] rigidity Default rigidity scale in Pa.
             * @param[in] thermalConductivity Default thermal conductivity in W/(m*K).
             * @param[in] density Default density in kg/m^3.
             * @param[in] specificHeat Default specific heat capacity in J/(kg*K).
             *
             * @returns Time scale in seconds.
             */
            static
            double computeThermoporoelasticityTimeScale(const double lengthScale,
                                                        const double permeability,
                                                        const double viscosity,
                                                        const double rigidity,
                                                        const double thermalConductivity,
                                                        const double density,
                                                        const double specificHeat);

            /** Get value to nondimensionalize stress.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Stress scale in Pa (SI units).
             */
            static
            double getStressScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize pressure.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Fluid pressure scale in Pascals (SI units).
             */
            static
            double getFluidPressureScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize strain.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Strain scale.
             */
            static
            double getStrainScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize body force.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Body force scale in Pa/m (SI units).
             */
            static
            double getBodyForceScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize density.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Density scale in kg/m^3 (SI units).
             */
            static
            double getDensityScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize velocity.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Velocity scale in m/s (SI units).
             */
            static
            double getVelocityScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize acceleration.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Acceleration scale in m/s^2 (SI units).
             */
            static
            double getAccelerationScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize fluid viscosity.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Fluid viscosity scale in Pa*s (SI units).
             */
            static
            double getViscosityScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize permeability.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Permeability scale in m/s^2 (SI units).
             */
            static
            double getPermeabilityScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize temperature.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Temperature scale in K (SI units).
             */
            static
            double getTemperatureScale(const pylith::scales::Scales& scales);

            /** Get value to nondimensionalize heat flux.
             *
             * @param[in] Scales for nondimensionalization.
             * @returns Heat flux scale in W/m^2 (SI units).
             */
            static
            double getHeatFluxScale(const pylith::scales::Scales& scales);

        }; // class ElasticityScales

    } // units
} // spatialdata

// End of file
