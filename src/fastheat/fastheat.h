//
// Created for FastHeat
//

#ifndef FASTHEAT_H
#define FASTHEAT_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "../../algorithms/rasteralgorithms/shadowmaps.h"
#include "../../common/threading/parallel.h"
#include "../../geospatial/map/map.h"

namespace FastHeatModel
{

    // -------------------------------------------------------------------------
    // Psychrometrics: Vapor pressures and Relative Humidity
    // -------------------------------------------------------------------------

    // Saturation vapor pressure [hPa] using Magnus formula
    inline double saturationVaporPressureHPa(double tempC)
    {
        return 6.112 * std::exp((17.67 * tempC) / (tempC + 243.5));
    }

    // Actual vapor pressure [hPa] from temperature [C] and RH [%]
    inline double actualVaporPressureHPa(double tempC, double rhPercent)
    {
        double es = saturationVaporPressureHPa(tempC);
        return std::max(0.1, (std::clamp(rhPercent, 1.0, 100.0) / 100.0) * es);
    }

    // Natural wet bulb temperature [C] using Stull (2011) formula
    inline double wetBulbTemperatureStull(double tempC, double rhPercent)
    {
        double T = tempC;
        double RH = std::clamp(rhPercent, 1.0, 100.0);

        double Tw = T * std::atan(0.151977 * std::sqrt(RH + 8.313659))
                  + std::atan(T + RH)
                  - std::atan(RH - 1.676331)
                  + 0.00391838 * std::pow(RH, 1.5) * std::atan(0.023101 * RH)
                  - 4.686035;

        return Tw;
    }

    // Conditioned summer afternoon relative humidity based on latitude and longitude
    // Centered around standard 32% (ea ~ 20 hPa) for conditioned summer afternoons in mid-latitudes,
    // adjusted for tropical (higher) and arid/desert belts.
    inline double conditionedSummerRH(double latitudeDeg, double /*longitudeDeg*/)
    {
        double absLat = std::abs(latitudeDeg);
        double rh = 32.0;

        if (absLat < 15.0)
        {
            // Humid tropics (equatorial zone)
            rh = 65.0 - (absLat / 15.0) * 33.0; // 65% down to 32%
        }
        else if (absLat < 35.0)
        {
            // Subtropical / desert / Mediterranean summer dry belt
            rh = 30.0 + ((absLat - 15.0) / 20.0) * 2.0; // 30% to 32%
        }
        else if (absLat < 55.0)
        {
            // Temperate mid-latitudes
            rh = 32.0 + ((absLat - 35.0) / 20.0) * 16.0; // 32% to 48%
        }
        else
        {
            // Subpolar / high latitudes
            rh = 48.0 + std::min(20.0, (absLat - 55.0) * 0.8); // 48% to 65%
        }

        return std::clamp(rh, 15.0, 85.0);
    }

    // GEV-based placeholder temperature for return period heat hazard
    inline double placeholderGEVTemperature(double latitudeDeg, double returnPeriodYears)
    {
        double absLat = std::abs(latitudeDeg);
        double rp = std::clamp(returnPeriodYears, 2.0, 1000.0);

        // Baseline summer max temperature at sea level
        double tBase = 38.0 - 0.22 * absLat;
        if (absLat >= 15.0 && absLat <= 40.0)
        {
            tBase += 4.0; // Continental subtropical summer heating
        }

        // GEV-like return period increase
        double deltaTrp = 1.8 * std::log(std::max(1.5, rp / 2.0));
        double temp = tBase + deltaTrp;

        return std::clamp(temp, 20.0, 52.0);
    }

    // -------------------------------------------------------------------------
    // Universal Thermal Climate Index (UTCI) - Official 6th Order Polynomial
    // COST Action 730 / Broede et al. (2012)
    //
    // Inputs:
    //   Ta     : Air temperature [deg C]
    //   Tmrt   : Mean radiant temperature [deg C]
    //   v10    : Wind speed at 10m [m/s]
    //   eaHPa  : Water vapor pressure [hPa]
    // -------------------------------------------------------------------------
    inline double calculateUTCI(double Ta, double Tmrt, double v10, double eaHPa)
    {
        double d_tr = Tmrt - Ta;
        double vel = std::clamp(v10, 0.5, 17.0);
        double pa_pr = std::clamp(eaHPa * 0.1, 0.05, 5.0); // convert hPa to kPa
        double ta = std::clamp(Ta, -50.0, 50.0);

        // Pre-calculate powers so we can re-use them
        double ta2 = ta * ta;
        double ta3 = ta2 * ta;
        double ta4 = ta3 * ta;
        double ta5 = ta4 * ta;
        double ta6 = ta5 * ta;

        double vel2 = vel * vel;
        double vel3 = vel2 * vel;
        double vel4 = vel3 * vel;
        double vel5 = vel4 * vel;
        double vel6 = vel5 * vel;

        double d_tr2 = d_tr * d_tr;
        double d_tr3 = d_tr2 * d_tr;
        double d_tr4 = d_tr3 * d_tr;
        double d_tr5 = d_tr4 * d_tr;
        double d_tr6 = d_tr5 * d_tr;

        double pa_pr2 = pa_pr * pa_pr;
        double pa_pr3 = pa_pr2 * pa_pr;
        double pa_pr4 = pa_pr3 * pa_pr;
        double pa_pr5 = pa_pr4 * pa_pr;
        double pa_pr6 = pa_pr5 * pa_pr;

        double utci = ta +
            0.607562052 +
            -0.0227712343 * ta +
            8.06470249e-4 * ta2 +
            -1.54271372e-4 * ta3 +
            -3.24651735e-6 * ta4 +
            7.32602852e-8 * ta5 +
            1.35959073e-9 * ta6 +
            -2.25836520 * vel +
            0.0880326035 * ta * vel +
            0.00216844454 * ta2 * vel +
            -1.53347087e-5 * ta3 * vel +
            -5.72983704e-7 * ta4 * vel +
            -2.55090145e-9 * ta5 * vel +
            -0.751269505 * vel2 +
            -0.00408350271 * ta * vel2 +
            -5.21670675e-5 * ta2 * vel2 +
            1.94544667e-6 * ta3 * vel2 +
            1.14099531e-8 * ta4 * vel2 +
            0.158137256 * vel3 +
            -6.57263143e-5 * ta * vel3 +
            2.22697524e-7 * ta2 * vel3 +
            -4.16117031e-8 * ta3 * vel3 +
            -0.0127762753 * vel4 +
            9.66891875e-6 * ta * vel4 +
            2.52785852e-9 * ta2 * vel4 +
            4.56306672e-4 * vel5 +
            -1.74202546e-7 * ta * vel5 +
            -5.91491269e-6 * vel6 +
            0.398374029 * d_tr +
            1.83945314e-4 * ta * d_tr +
            -1.73754510e-4 * ta2 * d_tr +
            -7.60781159e-7 * ta3 * d_tr +
            3.77830287e-8 * ta4 * d_tr +
            5.43079673e-10 * ta5 * d_tr +
            -0.0200518269 * vel * d_tr +
            8.92859837e-4 * ta * vel * d_tr +
            3.45433048e-6 * ta2 * vel * d_tr +
            -3.77925774e-7 * ta3 * vel * d_tr +
            -1.69699377e-9 * ta4 * vel * d_tr +
            1.69992415e-4 * vel2 * d_tr +
            -4.99204314e-5 * ta * vel2 * d_tr +
            2.47417178e-7 * ta2 * vel2 * d_tr +
            1.07596466e-8 * ta3 * vel2 * d_tr +
            8.49242932e-5 * vel3 * d_tr +
            1.35191328e-6 * ta * vel3 * d_tr +
            -6.21531254e-9 * ta2 * vel3 * d_tr +
            -4.99410301e-6 * vel4 * d_tr +
            -1.89489258e-8 * ta * vel4 * d_tr +
            8.15300114e-8 * vel5 * d_tr +
            7.55043090e-4 * d_tr2 +
            -5.65095215e-5 * ta * d_tr2 +
            -4.52166564e-7 * ta2 * d_tr2 +
            2.46688878e-8 * ta3 * d_tr2 +
            2.42674348e-10 * ta4 * d_tr2 +
            1.54547250e-4 * vel * d_tr2 +
            5.24110970e-6 * ta * vel * d_tr2 +
            -8.75874982e-8 * ta2 * vel * d_tr2 +
            -1.50743064e-9 * ta3 * vel * d_tr2 +
            -1.56236307e-5 * vel2 * d_tr2 +
            -1.33895614e-7 * ta * vel2 * d_tr2 +
            2.49709824e-9 * ta2 * vel2 * d_tr2 +
            6.51711721e-7 * vel3 * d_tr2 +
            1.94960053e-9 * ta * vel3 * d_tr2 +
            -1.00361113e-8 * vel4 * d_tr2 +
            -1.21206673e-5 * d_tr3 +
            -2.18203660e-7 * ta * d_tr3 +
            7.51269482e-9 * ta2 * d_tr3 +
            9.79063848e-11 * ta3 * d_tr3 +
            1.25006734e-6 * vel * d_tr3 +
            -1.81584736e-9 * ta * vel * d_tr3 +
            -3.52197671e-10 * ta2 * vel * d_tr3 +
            -3.36514630e-8 * vel2 * d_tr3 +
            1.35908359e-10 * ta * vel2 * d_tr3 +
            4.17032620e-10 * vel3 * d_tr3 +
            -1.30369025e-9 * d_tr4 +
            4.13908461e-10 * ta * d_tr4 +
            9.22652254e-12 * ta2 * d_tr4 +
            -5.08220384e-9 * vel * d_tr4 +
            -2.24730961e-11 * ta * vel * d_tr4 +
            1.17139133e-10 * vel2 * d_tr4 +
            6.62154879e-10 * d_tr5 +
            4.03863260e-13 * ta * d_tr5 +
            1.95087203e-12 * vel * d_tr5 +
            -4.73602469e-12 * d_tr6 +
            5.12733497 * pa_pr +
            -0.312788561 * ta * pa_pr +
            -0.0196701861 * ta2 * pa_pr +
            9.99690870e-4 * ta3 * pa_pr +
            9.51738512e-6 * ta4 * pa_pr +
            -4.66426341e-7 * ta5 * pa_pr +
            0.548050612 * vel * pa_pr +
            -0.00330552823 * ta * vel * pa_pr +
            -0.00164119440 * ta2 * vel * pa_pr +
            -5.16670694e-6 * ta3 * vel * pa_pr +
            9.52692432e-7 * ta4 * vel * pa_pr +
            -0.0429223622 * vel2 * pa_pr +
            0.00500845667 * ta * vel2 * pa_pr +
            1.00601257e-6 * ta2 * vel2 * pa_pr +
            -1.81748644e-6 * ta3 * vel2 * pa_pr +
            -1.25813502e-3 * vel3 * pa_pr +
            -1.79330391e-4 * ta * vel3 * pa_pr +
            2.34994441e-6 * ta2 * vel3 * pa_pr +
            1.29735808e-4 * vel4 * pa_pr +
            1.29064870e-6 * ta * vel4 * pa_pr +
            -2.28558686e-6 * vel5 * pa_pr +
            -0.0369476348 * d_tr * pa_pr +
            0.00162325322 * ta * d_tr * pa_pr +
            -3.14279680e-5 * ta2 * d_tr * pa_pr +
            2.59835559e-6 * ta3 * d_tr * pa_pr +
            -4.77136523e-8 * ta4 * d_tr * pa_pr +
            8.64203390e-3 * vel * d_tr * pa_pr +
            -6.87405181e-4 * ta * vel * d_tr * pa_pr +
            -9.13863872e-6 * ta2 * vel * d_tr * pa_pr +
            5.15916806e-7 * ta3 * vel * d_tr * pa_pr +
            -3.59217476e-5 * vel2 * d_tr * pa_pr +
            3.28696511e-5 * ta * vel2 * d_tr * pa_pr +
            -7.10542454e-7 * ta2 * vel2 * d_tr * pa_pr +
            -1.24382300e-5 * vel3 * d_tr * pa_pr +
            -7.38584400e-9 * ta * vel3 * d_tr * pa_pr +
            2.20609296e-7 * vel4 * d_tr * pa_pr +
            -7.32469180e-4 * d_tr2 * pa_pr +
            -1.87381964e-5 * ta * d_tr2 * pa_pr +
            4.80925239e-6 * ta2 * d_tr2 * pa_pr +
            -8.75492040e-8 * ta3 * d_tr2 * pa_pr +
            2.77862930e-5 * vel * d_tr2 * pa_pr +
            -5.06004592e-6 * ta * vel * d_tr2 * pa_pr +
            1.14325367e-7 * ta2 * vel * d_tr2 * pa_pr +
            2.53016723e-6 * vel2 * d_tr2 * pa_pr +
            -1.72857035e-8 * ta * vel2 * d_tr2 * pa_pr +
            -3.95079398e-8 * vel3 * d_tr2 * pa_pr +
            -3.59413173e-7 * d_tr3 * pa_pr +
            7.04388046e-7 * ta * d_tr3 * pa_pr +
            -1.89309167e-8 * ta2 * d_tr3 * pa_pr +
            -4.79768731e-7 * vel * d_tr3 * pa_pr +
            7.96079978e-9 * ta * vel * d_tr3 * pa_pr +
            1.62897058e-9 * vel2 * d_tr3 * pa_pr +
            3.94367674e-8 * d_tr4 * pa_pr +
            -1.18566247e-9 * ta * d_tr4 * pa_pr +
            3.34678041e-10 * vel * d_tr4 * pa_pr +
            -1.15606447e-10 * d_tr5 * pa_pr +
            -2.80626406 * pa_pr2 +
            0.548712484 * ta * pa_pr2 +
            -0.00399428410 * ta2 * pa_pr2 +
            -9.54009191e-4 * ta3 * pa_pr2 +
            1.93090978e-5 * ta4 * pa_pr2 +
            -0.308806365 * vel * pa_pr2 +
            0.0116952364 * ta * vel * pa_pr2 +
            4.95271903e-4 * ta2 * vel * pa_pr2 +
            -1.90710882e-5 * ta3 * vel * pa_pr2 +
            0.00210787756 * vel2 * pa_pr2 +
            -6.98445738e-4 * ta * vel2 * pa_pr2 +
            2.30109073e-5 * ta2 * vel2 * pa_pr2 +
            4.17856590e-4 * vel3 * pa_pr2 +
            -1.27043871e-5 * ta * vel3 * pa_pr2 +
            -3.04620472e-6 * vel4 * pa_pr2 +
            0.0514507424 * d_tr * pa_pr2 +
            -0.00432510997 * ta * d_tr * pa_pr2 +
            8.99281156e-5 * ta2 * d_tr * pa_pr2 +
            -7.14663943e-7 * ta3 * d_tr * pa_pr2 +
            -2.66016305e-4 * vel * d_tr * pa_pr2 +
            2.63789586e-4 * ta * vel * d_tr * pa_pr2 +
            -7.01199003e-6 * ta2 * vel * d_tr * pa_pr2 +
            -1.06823306e-4 * vel2 * d_tr * pa_pr2 +
            3.61341136e-6 * ta * vel2 * d_tr * pa_pr2 +
            2.29748967e-7 * vel3 * d_tr * pa_pr2 +
            3.04788893e-4 * d_tr2 * pa_pr2 +
            -6.42070836e-5 * ta * d_tr2 * pa_pr2 +
            1.16257971e-6 * ta2 * d_tr2 * pa_pr2 +
            7.68023384e-6 * vel * d_tr2 * pa_pr2 +
            -5.47446896e-7 * ta * vel * d_tr2 * pa_pr2 +
            -3.59937910e-8 * vel2 * d_tr2 * pa_pr2 +
            -4.36497725e-6 * d_tr3 * pa_pr2 +
            1.68737969e-7 * ta * d_tr3 * pa_pr2 +
            2.67489271e-8 * vel * d_tr3 * pa_pr2 +
            3.23926897e-9 * d_tr4 * pa_pr2 +
            -0.0353874123 * pa_pr3 +
            -0.221201190 * ta * pa_pr3 +
            0.0155126038 * ta2 * pa_pr3 +
            -2.63917279e-4 * ta3 * pa_pr3 +
            0.0453433455 * vel * pa_pr3 +
            -0.00432943862 * ta * vel * pa_pr3 +
            1.45389826e-4 * ta2 * vel * pa_pr3 +
            2.17508610e-4 * vel2 * pa_pr3 +
            -6.66724702e-5 * ta * vel2 * pa_pr3 +
            3.33217140e-5 * vel3 * pa_pr3 +
            -0.00226921615 * d_tr * pa_pr3 +
            3.80261982e-4 * ta * d_tr * pa_pr3 +
            -5.45314314e-9 * ta2 * d_tr * pa_pr3 +
            -7.96355448e-4 * vel * d_tr * pa_pr3 +
            2.53458034e-5 * ta * vel * d_tr * pa_pr3 +
            -6.31223658e-6 * vel2 * d_tr * pa_pr3 +
            3.02122035e-4 * d_tr2 * pa_pr3 +
            -4.77403547e-6 * ta * d_tr2 * pa_pr3 +
            1.73825715e-6 * vel * d_tr2 * pa_pr3 +
            -4.09087898e-7 * d_tr3 * pa_pr3 +
            0.614155345 * pa_pr4 +
            -0.0616755931 * ta * pa_pr4 +
            0.00133374846 * ta2 * pa_pr4 +
            0.00355375387 * vel * pa_pr4 +
            -5.13027851e-4 * ta * vel * pa_pr4 +
            1.02449757e-4 * vel2 * pa_pr4 +
            -0.00148526421 * d_tr * pa_pr4 +
            -4.11469183e-5 * ta * d_tr * pa_pr4 +
            -6.80434415e-6 * vel * d_tr * pa_pr4 +
            -9.77675906e-6 * d_tr2 * pa_pr4 +
            0.0882773108 * pa_pr5 +
            -0.00301859306 * ta * pa_pr5 +
            0.00104452989 * vel * pa_pr5 +
            2.47090539e-4 * d_tr * pa_pr5 +
            0.00148348065 * pa_pr6;

        return utci;
    }

    // -------------------------------------------------------------------------
    // SOLWEIG Mean Radiant Temperature (Tmrt) calculation
    // Computes shortwave and longwave radiation fluxes on a human body
    // -------------------------------------------------------------------------
    inline double calculateTmrtSOLWEIG(
        double Ta,
        double RH,
        double svf,
        float beamFactor,
        double dni,
        double dhi,
        double solarElevationDeg,
        int luClass = 50,
        float canopyHeight = 0.0f,
        double albedoOverride = -1.0)
    {
        const double sigma = 5.670374e-8; // Stefan-Boltzmann constant [W/m^2/K^4]
        const double ak = 0.70;           // Shortwave absorption of human body
        const double al = 0.97;           // Longwave emissivity/absorption of human body

        double TaK = Ta + 273.15;
        double ea = actualVaporPressureHPa(Ta, RH);

        // 1. Shortwave radiation fluxes on standing human cylinder
        double elevRad = Solar::degToRad(std::max(0.0, solarElevationDeg));
        double elevDeg = std::max(0.0, solarElevationDeg);

        // Projected area factor for a standing human (Fanger 1970 / SOLWEIG)
        // fp(theta) = 0.308 * cos(theta_deg * (0.998 - theta_deg^2 / 22600.0))
        double fAngRad = Solar::degToRad(elevDeg * (0.998 - (elevDeg * elevDeg) / 22600.0));
        double fp = std::clamp(0.308 * std::cos(fAngRad), 0.08, 0.308);

        // Direct beam shortwave flux on human
        double psi = std::clamp(static_cast<double>(beamFactor), 0.0, 1.0);
        double Kdir = psi * dni * fp;

        // Diffuse sky shortwave flux on human (cylinder view factor to sky ~ 0.5 * SVF)
        double Ksky = 0.5 * svf * dhi;

        // Ground-reflected shortwave flux
        double Igh_ground = psi * dni * std::sin(elevRad) + svf * dhi;

        // Ground albedo and excess ground surface heating depending on surface type / land use class
        // Supported classes:
        // 80: Water
        // 30, 40, 90: Grass, Herbaceous, Crops
        // 20: Shrubs
        // 10, 95: Trees, Forest, Mangroves
        // 50: Built-up, Concrete
        // 60, default: Asphalt, Bare ground
        double alphaGround = 0.18;
        double deltaTgBase = 20.0;
        double epsGround = 0.95;

        if (luClass == 80) // Water
        {
            alphaGround = 0.06;
            deltaTgBase = 0.0; // Water doesn't overheat like dry urban pavement
            epsGround = 0.98;
        }
        else if (luClass == 30 || luClass == 40 || luClass == 90) // Grass / lawn / meadow / crops
        {
            alphaGround = 0.22;
            deltaTgBase = 6.0; // Evaporative cooling / transpiration
            epsGround = 0.97;
        }
        else if (luClass == 20) // Shrubs / brush
        {
            alphaGround = 0.20;
            deltaTgBase = 8.0;
            epsGround = 0.96;
        }
        else if (luClass == 10 || luClass == 95) // Trees / Forest
        {
            alphaGround = 0.16;
            deltaTgBase = 3.0;
            epsGround = 0.97;
        }
        else if (luClass == 50) // Built-up / concrete
        {
            alphaGround = 0.28;
            deltaTgBase = 18.0;
            epsGround = 0.93;
        }
        else // Asphalt / bare ground / default
        {
            alphaGround = 0.15;
            deltaTgBase = 26.0;
            epsGround = 0.95;
        }

        // Differentiate effects of canopy height:
        // - Mature trees (canopyHeight >= 4.0m): overhead shade & heavy understory cooling (deltaTgBase -> 3.0K)
        // - Shrubs (1.5m <= canopyHeight < 4.0m): low bushes with intermediate cooling (deltaTgBase -> 8.0K)
        // - Grass / ground (< 1.5m): underlying land use governs surface heating (e.g. grass = 6.0K)
        if (luClass != 80 && canopyHeight >= 1.5f)
        {
            if (canopyHeight >= 4.0f)
            {
                double fveg = std::clamp(static_cast<double>(canopyHeight) / 12.0, 0.4, 1.0);
                deltaTgBase = (1.0 - fveg) * deltaTgBase + fveg * 3.0;
                alphaGround = (1.0 - fveg) * alphaGround + fveg * 0.16;
                epsGround = (1.0 - fveg) * epsGround + fveg * 0.97;
            }
            else
            {
                // Shrubs / bushes: intermediate ground cooling and albedo
                double fshrub = std::clamp((static_cast<double>(canopyHeight) - 1.5) / 2.5, 0.0, 1.0);
                deltaTgBase = (1.0 - fshrub) * std::min(deltaTgBase, 12.0) + fshrub * 8.0;
                alphaGround = (1.0 - fshrub) * alphaGround + fshrub * 0.20;
                epsGround = (1.0 - fshrub) * epsGround + fshrub * 0.96;
            }
        }

        if (albedoOverride >= 0.0)
        {
            alphaGround = std::clamp(albedoOverride, 0.02, 0.95);
            // High albedo / cool materials reflect solar radiation, significantly reducing excess ground surface heating
            deltaTgBase = std::max(0.0, 26.0 * (1.0 - alphaGround) / 0.85);
        }

        double Kground = 0.5 * alphaGround * Igh_ground;

        // Wall-reflected shortwave flux from surrounding buildings/trees
        // In shadow, surroundings also do not receive direct beam
        double alphaWall = 0.20;
        double Kwalls = 0.5 * (1.0 - svf) * alphaWall * (dhi + 0.5 * psi * dni * std::sin(elevRad));

        double S_short = ak * (Kdir + Ksky + Kground + Kwalls);

        // 2. Longwave radiation fluxes on human
        // Clear-sky atmospheric emissivity (Prata 1996)
        double xi = 46.5 * (ea / TaK);
        double epsSky = 1.0 - (1.0 + xi) * std::exp(-std::sqrt(1.2 + 3.0 * xi));
        epsSky = std::clamp(epsSky, 0.65, 0.98);

        // Downwelling longwave from sky
        double Lsky = svf * epsSky * sigma * std::pow(TaK, 4.0);

        // Ground surface temperature and upwelling longwave
        // Unshaded ground heats up in intense summer afternoon sun, modulated by land surface type
        double deltaTg = deltaTgBase * psi * std::clamp(dni / 850.0, 0.0, 1.3);
        double TgK = TaK + deltaTg;
        double Lground = 0.5 * epsGround * sigma * std::pow(TgK, 4.0);

        // Surrounding walls/canopy longwave
        double deltaTw = 6.0 * psi;
        double TwK = TaK + deltaTw;
        double epsWall = 0.90;
        double Lwalls = 0.5 * (1.0 - svf) * epsWall * sigma * std::pow(TwK, 4.0);

        double S_long = al * (0.5 * Lsky + Lground + Lwalls);

        // Total radiant absorption density Sstr [W/m^2]
        double Sstr = S_short + S_long;

        // Stefan-Boltzmann inversion for Mean Radiant Temperature
        double TmrtK = std::pow(std::max(1.0, Sstr / (al * sigma)), 0.25);
        double TmrtC = TmrtK - 273.15;

        return TmrtC;
    }

    // -------------------------------------------------------------------------
    // Wet Bulb Globe Temperature (WBGT)
    // ISO 7243: WBGT = 0.7 * Tnw + 0.2 * Tg + 0.1 * Ta
    // -------------------------------------------------------------------------
    inline double calculateWBGT(double Ta, double RH, double Tmrt, double v10)
    {
        double Tnw = wetBulbTemperatureStull(Ta, RH);

        // Black globe temperature Tg from thermal equilibrium with Tmrt, Ta, and v10
        // Diameter D = 0.05m black globe
        double D = 0.05;
        double v = std::max(0.2, v10);
        double hc = 6.3 * std::pow(v, 0.6) / std::pow(D, 0.4); // Convective coefficient [W/m^2/K]

        const double sigma = 5.670374e-8;
        const double epsG = 0.95;

        double TaK = Ta + 273.15;
        double TmrtK = Tmrt + 273.15;
        double TbarK = 0.5 * (TaK + TmrtK);
        double radCoeff = 4.0 * epsG * sigma * std::pow(TbarK, 3.0);

        // Steady-state globe temperature
        double TgK = (hc * TaK + radCoeff * TmrtK) / (hc + radCoeff);
        double Tg = TgK - 273.15;

        double wbgt = 0.7 * Tnw + 0.2 * Tg + 0.1 * Ta;
        return wbgt;
    }

    // -------------------------------------------------------------------------
    // Tygron / DPRA Heat Stress Module Methodology
    // Based on DPRA Heat Stress Report & Tygron Heat Stress Module Theory
    // -------------------------------------------------------------------------

    // 1. Surface Properties: Vegetation Fraction and Bowen Ratio
    // Per DPRA report:
    // - Mature Trees (canopyHeight >= 4.0m): F_veg = 1.0, Bowen ratio = 0.4
    // - Shrubs / Bushes (1.5m <= canopyHeight < 4.0m): F_veg = 0.8, Bowen ratio = 0.4
    // - Ground / grass / crops (canopyHeight < 1.5m): governed by land use class:
    //   * Nature / Trees / Forest (LU 10, 95): F_veg = 1.0, Bowen ratio = 0.4
    //   * Shrubs / Brush (LU 20): F_veg = 0.8, Bowen ratio = 0.4
    //   * Grass / Lawn / Meadow / Agriculture (LU 30, 40, 90): F_veg = 1.0, Bowen ratio = 0.4
    //   * Water (LU 80): Daytime evaporative cooling behaves as vegetation (F_veg = 1.0, Bowen ratio = 0.4)
    //   * Built-up / Paved / Asphalt / Concrete (LU 50, 60, default): F_veg = 0.0, Bowen ratio = 3.0
    inline void getSurfaceProperties(int luClass, float canopyHeight, double& fveg, double& bowenRatio)
    {
        // 1. Overhead arboreal tree canopy
        if (canopyHeight >= 4.0f)
        {
            fveg = 1.0;
            bowenRatio = 0.4;
            return;
        }
        // 2. Low woody shrubs and hedges (1.5m - 4.0m)
        else if (canopyHeight >= 1.5f)
        {
            fveg = 0.8;
            bowenRatio = 0.4;
            return;
        }

        // 3. Ground surface cover (canopyHeight < 1.5m: grass plots, water, pavement)
        if (luClass == 80) // Water (evaporative cooling during day)
        {
            fveg = 1.0;
            bowenRatio = 0.4;
        }
        else if (luClass == 10 || luClass == 95) // Forest / Trees land use
        {
            fveg = 1.0;
            bowenRatio = 0.4;
        }
        else if (luClass == 20) // Shrubs land use
        {
            fveg = 0.8;
            bowenRatio = 0.4;
        }
        else if (luClass == 30 || luClass == 40 || luClass == 90) // Grass, herbaceous, crops
        {
            fveg = 1.0;
            bowenRatio = 0.4;
        }
        else // Built-up (50), asphalt/pavement (60), bare ground, default
        {
            fveg = 0.0;
            bowenRatio = 3.0;
        }
    }

    // 2. Fast 2D Integral Image (Prefix Sum) Neighborhood Averaging
    // Runs in O(rows * cols) time independent of window size.
    // Handles missing values (MV) and domain boundary padding cleanly.
    inline std::shared_ptr<cTMapf> calculateNeighborhoodAverage(
        const std::shared_ptr<cTMapf>& inputMap,
        double windowMeters,
        double dx)
    {
        if (!inputMap) return nullptr;
        const int rows = inputMap->nrRows();
        const int cols = inputMap->nrCols();
        if (rows <= 0 || cols <= 0) return nullptr;

        auto result = Solar::makeLike(inputMap);
        int rad = std::max(1, static_cast<int>(std::round(0.5 * windowMeters / std::max(0.1, dx))));

        // Allocate prefix sum and count tables of size (rows + 1) * (cols + 1)
        const size_t stride = static_cast<size_t>(cols + 1);
        std::vector<double> sumTab(static_cast<size_t>(rows + 1) * stride, 0.0);
        std::vector<int> cntTab(static_cast<size_t>(rows + 1) * stride, 0);

        for (int r = 0; r < rows; ++r)
        {
            double rowSum = 0.0;
            int rowCnt = 0;
            const size_t currRowOffset = static_cast<size_t>(r + 1) * stride;
            const size_t prevRowOffset = static_cast<size_t>(r) * stride;

            for (int c = 0; c < cols; ++c)
            {
                float v = inputMap->data[r][c];
                if (!pcr::isMV(v) && std::isfinite(v))
                {
                    rowSum += static_cast<double>(v);
                    rowCnt += 1;
                }
                sumTab[currRowOffset + (c + 1)] = sumTab[prevRowOffset + (c + 1)] + rowSum;
                cntTab[currRowOffset + (c + 1)] = cntTab[prevRowOffset + (c + 1)] + rowCnt;
            }
        }

        // Parallel evaluation of window average for each cell
        GetTPool().parallel_for(0, rows, [&](int r)
        {
            int r0 = std::max(0, r - rad);
            int r1 = std::min(rows - 1, r + rad);

            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(inputMap->data[r][c]))
                {
                    pcr::setMV(result->data[r][c]);
                    continue;
                }

                int c0 = std::max(0, c - rad);
                int c1 = std::min(cols - 1, c + rad);

                double s = sumTab[static_cast<size_t>(r1 + 1) * stride + (c1 + 1)]
                         - sumTab[static_cast<size_t>(r0) * stride + (c1 + 1)]
                         - sumTab[static_cast<size_t>(r1 + 1) * stride + c0]
                         + sumTab[static_cast<size_t>(r0) * stride + c0];

                int n = cntTab[static_cast<size_t>(r1 + 1) * stride + (c1 + 1)]
                      - cntTab[static_cast<size_t>(r0) * stride + (c1 + 1)]
                      - cntTab[static_cast<size_t>(r1 + 1) * stride + c0]
                      + cntTab[static_cast<size_t>(r0) * stride + c0];

                if (n > 0)
                {
                    result->data[r][c] = static_cast<float>(s / static_cast<double>(n));
                }
                else
                {
                    result->data[r][c] = inputMap->data[r][c];
                }
            }
        });
        GetTPool().waitext();

        return result;
    }

    // 3. DPRA Sun Daily Motion Factor f_sun,t
    // Sourced from DPRA Heat Stress Report (Section B1.2.7) / Tygron Sun Motion Table (Table 5)
    // Interpolates UHI_max diurnal fraction based on local solar hour.
    inline double sunDailyMotionFactor(double solarHour)
    {
        double h = std::fmod(solarHour, 24.0);
        if (h < 0.0) h += 24.0;

        static const double dpraTable[25] = {
            1.000, // 00:00
            0.888, // 01:00
            0.728, // 02:00
            0.609, // 03:00
            0.490, // 04:00
            0.256, // 05:00
            0.079, // 06:00
            0.007, // 07:00
            -0.020, // 08:00
            -0.006, // 09:00
            0.010, // 10:00
            0.033, // 11:00
            0.056, // 12:00
            0.082, // 13:00
            0.128, // 14:00
            0.184, // 15:00
            0.270, // 16:00
            0.366, // 17:00
            0.506, // 18:00
            0.651, // 19:00
            0.803, // 20:00
            0.901, // 21:00
            0.958, // 22:00
            0.983, // 23:00
            1.000  // 24:00
        };

        int idx = static_cast<int>(std::floor(h));
        double frac = h - static_cast<double>(idx);
        int nextIdx = (idx + 1) % 24;
        double val = (1.0 - frac) * dpraTable[idx] + frac * dpraTable[nextIdx];
        return std::clamp(val, 0.0, 1.0);
    }

    // 4. Urban Heat Island (UHI_max) Calculation (Theeuwes et al. 2016 / DPRA / Tygron)
    //
    // Formula:
    //   UHI_max = (2 - S_vf_avg - F_veg_avg) * ( (S * (T_max - T_min)^3) / U )^(1/4)
    // where:
    //   S = Q_gl_avg / (C_air * rho_air)   [K m/s]
    //   rho_air = p / (R_specific * (T_station + 273.15))   [kg/m^3]
    //   R_specific = 287.058 J/(kg K)
    //   C_air = 1007.0 J/(kg K)
    //   p = atmospheric pressure [Pa] (barometrically adjusted for elevation)
    //   T_max - T_min = daily diurnal temperature range [K] (default: 15 K)
    //   U = 10m reference wind speed [m/s]
    inline double calculateUHIMax(
        double svfAvg,
        double fvegAvg,
        double tStationC,
        double elevMeters,
        double qglAvg,
        double deltaT,
        double wind10m)
    {
        double p = 101325.0 * std::exp(-elevMeters / 8434.0);
        double tK = tStationC + 273.15;
        double rhoAir = p / (287.058 * tK);
        double s = qglAvg / (1007.0 * rhoAir);

        double u = std::max(0.5, wind10m);
        double dt = std::max(2.0, deltaT);

        double tempEffect = std::pow((s * dt * dt * dt) / u, 0.25);
        double factor = std::clamp(2.0 - svfAvg - fvegAvg, 0.0, 2.0);

        return factor * tempEffect;
    }

    // 5. Physiological Equivalent Temperature (PET) Calculation
    // Based on DPRA Heat Stress Report & Tygron PET Formulas
    //
    // In Sun (direct radiation present, beamFactor > 0):
    //   PET_sun = -13.26 + 1.25 * Ta + 0.011 * Q_gl - 3.37 * ln(u) + 0.078 * Tw
    //           + 0.0055 * Q_gl * ln(u) + 5.56 * sin(phi) - 0.0103 * Q_gl * ln(u) * sin(phi)
    //           + 0.546 * Bb + 1.94 * Svf
    //
    // In Shade / Night (beamFactor == 0 or sun below horizon):
    //   PET_shade = -12.14 + 1.25 * Ta - 1.47 * ln(u) + 0.060 * Tw
    //             + 0.015 * Svf * Q_d + 0.0060 * (1 - Svf) * sigma * (Ta + 273.15)^4
    //
    // Partial shade is blended via:
    //   PET = (1 - psi) * PET_shade + psi * PET_sun
    inline double calculatePET(
        double Ta,               // Atmospheric temperature [deg C] (including UHI)
        double RH,               // Relative humidity [%]
        double svf,              // Sky View Factor [0 to 1]
        float beamFactor,        // Direct beam shadow factor [0 to 1]
        double ghi,              // Global solar radiation Q_gl [W/m^2]
        double solarElevDeg,     // Solar altitude angle phi in degrees
        double windSpeed10m,     // Reference 10m wind speed [m/s]
        double bowenRatio,       // Bowen ratio Bb (0.4 for veg/water, 3.0 for urban)
        double dhi = -1.0)       // Diffuse radiation Q_d [W/m^2] (optional)
    {
        double elevRad = Solar::degToRad(std::max(0.0, solarElevDeg));
        double sinElev = std::sin(elevRad);

        // Wind speed at 1.2m pedestrian level
        double u12 = std::max(0.25, windSpeed10m * 0.50);
        double lnU = std::log(u12);

        // Natural wet-bulb temperature Tw
        double Tw = wetBulbTemperatureStull(Ta, RH);

        // Diffuse radiation Qd
        double Qd = 0.0;
        if (solarElevDeg > 0.0 && ghi > 0.0)
        {
            if (dhi >= 0.0)
            {
                Qd = dhi;
            }
            else
            {
                // Tygron atmospheric transmissivity & transmission factor
                double tauA = ghi / std::max(10.0, 1367.0 * sinElev);
                double fTrans = 0.2;
                if (tauA < 0.3)
                {
                    fTrans = 1.0;
                }
                else if (tauA < 0.7)
                {
                    fTrans = 1.6 - 2.0 * tauA;
                }
                Qd = ghi * fTrans;
            }
        }

        // PET Shade and Night formula
        const double sigma = 5.670374e-8;
        double TaK = Ta + 273.15;
        double taK4 = TaK * TaK * TaK * TaK;
        double petShade = -12.14
                        + 1.25 * Ta
                        - 1.47 * lnU
                        + 0.060 * Tw
                        + 0.015 * svf * Qd
                        + 0.0060 * (1.0 - svf) * sigma * taK4;

        if (solarElevDeg <= 0.0 || ghi <= 0.0)
        {
            return petShade;
        }

        // PET Sun formula
        double petSun = -13.26
                      + 1.25 * Ta
                      + 0.011 * ghi
                      - 3.37 * lnU
                      + 0.078 * Tw
                      + 0.0055 * ghi * lnU
                      + 5.56 * sinElev
                      - 0.0103 * ghi * lnU * sinElev
                      + 0.546 * bowenRatio
                      + 1.94 * svf;

        float psi = std::clamp(beamFactor, 0.0f, 1.0f);
        double pet = (1.0 - static_cast<double>(psi)) * petShade + static_cast<double>(psi) * petSun;
        return pet;
    }

    // -------------------------------------------------------------------------
    // FastHeat Core Model Simulation
    // -------------------------------------------------------------------------
    struct HeatModelOutputs
    {
        std::shared_ptr<cTMapf> UTCI;       // Peak UTCI [deg C]
        std::shared_ptr<cTMapf> WBGT;       // Peak WBGT [deg C]
        std::shared_ptr<cTMapf> SUNHOURS;   // Non-shaded daily sun hours [hours]
        std::shared_ptr<cTMapf> WBGT_H38;   // Number of hours in day with heat stress >= 38 C
        std::shared_ptr<cTMapf> SVF;        // Sky View Factor [-]
        std::shared_ptr<cTMapf> TMRT;       // Peak Mean Radiant Temperature [deg C]
        std::shared_ptr<cTMapf> TMRT_AVG;   // Daylight average Mean Radiant Temperature [deg C]
        std::shared_ptr<cTMapf> UTCI_AVG;   // Daylight average UTCI [deg C]
        std::shared_ptr<cTMapf> PET;        // Peak Physiological Equivalent Temperature (PET) [deg C]
        std::shared_ptr<cTMapf> PET_AVG;    // Daylight average PET [deg C]
        std::shared_ptr<cTMapf> UHI;        // Urban Heat Island effect UHI_max [deg C]
        std::shared_ptr<cTMapf> TA;         // Peak Atmospheric Temperature (T_station + UHI*f_sun) [deg C]
    };

    inline int RunFastHeat(
        const std::shared_ptr<cTMapf>& DEM,
        const std::shared_ptr<cTMapf>& BH,
        const std::shared_ptr<cTMapf>& CH,
        const std::shared_ptr<cTMapf>& TEMP,
        const std::shared_ptr<cTMapf>& RH,
        double windSpeedV10,
        float canopyTransmissivity,
        double latitudeDeg,
        double longitudeDeg,
        int year,
        int month,
        int day,
        double peakHourUTC,
        HeatModelOutputs& outputs,
        const std::shared_ptr<cTMapf>& LU = nullptr,
        const std::shared_ptr<cTMapf>& ALBEDO = nullptr)
    {
        if (!DEM)
        {
            std::cerr << "FastHeat Error: DEM raster is null." << std::endl;
            return -1;
        }

        const int rows = DEM->data.nr_rows();
        const int cols = DEM->data.nr_cols();
        if (rows <= 0 || cols <= 0)
        {
            std::cerr << "FastHeat Error: DEM has invalid dimensions." << std::endl;
            return -1;
        }

        std::cout << "[FastHeat] Starting simulation for " << rows << "x" << cols << " domain..." << std::endl;
        std::cout << "[FastHeat] Domain center lat/lng: (" << latitudeDeg << ", " << longitudeDeg << ")" << std::endl;
        std::cout << "[FastHeat] Date: " << year << "-" << month << "-" << day << " at UTC hour: " << peakHourUTC << std::endl;

        // 1. Calculate Sky View Factor (SVF)
        std::cout << "[FastHeat] Calculating Sky View Factor (SVF) with 16 directional horizon sweeps..." << std::endl;
        auto svf = Solar::skyViewFactor(DEM, BH, CH, canopyTransmissivity, 16);
        outputs.SVF = svf;

        // 2. Calculate daily non-shaded sun hours
        std::cout << "[FastHeat] Calculating daily non-shaded sun hours (canopy transmissivity = " << canopyTransmissivity << ")..." << std::endl;
        auto sunhours = Solar::dailySunHoursWithCanopy(DEM, BH, CH, canopyTransmissivity, latitudeDeg, longitudeDeg, year, month, day, 16);
        outputs.SUNHOURS = sunhours;

        // 3. Precompute surface properties: Vegetation fraction and Bowen ratio
        std::cout << "[FastHeat] Deriving surface properties (Vegetation fraction and Bowen ratio) per DPRA guidelines..." << std::endl;
        auto fvegMap = Solar::makeLike(DEM);
        auto bowenMap = Solar::makeLike(DEM);

        GetTPool().parallel_for(0, rows, [&](int r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(fvegMap->data[r][c]);
                    pcr::setMV(bowenMap->data[r][c]);
                    continue;
                }

                int luVal = (LU && !pcr::isMV(LU->data[r][c])) ? static_cast<int>(std::round(LU->data[r][c])) : 50;
                float chVal = (CH && !pcr::isMV(CH->data[r][c])) ? std::max(0.0f, CH->data[r][c]) : 0.0f;

                double fv = 0.0;
                double bw = 3.0;
                getSurfaceProperties(luVal, chVal, fv, bw);
                fvegMap->data[r][c] = static_cast<float>(fv);
                bowenMap->data[r][c] = static_cast<float>(bw);
            }
        });
        GetTPool().waitext();

        // 4. Neighborhood averaging for DPRA UHI calculation (700m window)
        double dx = std::abs(DEM->cellSizeX());
        if (dx <= 0.0) dx = std::abs(DEM->cellSize());
        if (dx <= 0.0) dx = 2.5;

        std::cout << "[FastHeat] Calculating 700m neighborhood averages for DPRA UHI (SVF and Vegetation Fraction)..." << std::endl;
        auto svfAvg = calculateNeighborhoodAverage(svf, 700.0, dx);
        auto fvegAvg = calculateNeighborhoodAverage(fvegMap, 700.0, dx);

        // 5. Clear-sky solar radiation at peak afternoon time
        double avgElev = 0.0;
        int validCount = 0;
        for (int r = 0; r < rows; r += 10)
        {
            for (int c = 0; c < cols; c += 10)
            {
                if (!pcr::isMV(DEM->data[r][c]))
                {
                    avgElev += DEM->data[r][c];
                    validCount++;
                }
            }
        }
        if (validCount > 0) avgElev /= validCount;

        Solar::ClearSkyRadiation peakRad = Solar::clearSkySolarRadiation(latitudeDeg, longitudeDeg, year, month, day, peakHourUTC, avgElev);
        std::cout << "[FastHeat] Peak solar radiation: DNI (Idir) = " << peakRad.dni << " W/m^2, DHI (Idiff) = "
                  << peakRad.dhi << " W/m^2, GHI = " << peakRad.ghi << " W/m^2 (elev = " << peakRad.elevationDeg << " deg)" << std::endl;

        // Peak shadow sweep with canopy transmissivity
        std::cout << "[FastHeat] Performing dual-surface shadow sweep at peak condition..." << std::endl;
        auto peakBeamFactor = Solar::shadowSweepCanopy(DEM, BH, CH, peakRad.azimuthDeg, peakRad.elevationDeg, canopyTransmissivity);

        // Allocate output rasters
        outputs.UTCI = Solar::makeLike(DEM);
        outputs.WBGT = Solar::makeLike(DEM);
        outputs.TMRT = Solar::makeLike(DEM);
        outputs.WBGT_H38 = Solar::makeLike(DEM);
        outputs.TMRT_AVG = Solar::makeLike(DEM);
        outputs.UTCI_AVG = Solar::makeLike(DEM);
        outputs.PET = Solar::makeLike(DEM);
        outputs.PET_AVG = Solar::makeLike(DEM);
        outputs.UHI = Solar::makeLike(DEM);
        outputs.TA = Solar::makeLike(DEM);

        // Precompute daylight steps for hours above 38 C and average daylight metrics
        const double daylight = Solar::daylightHours(latitudeDeg, year, month, day);
        const int nDiurnalSteps = 14;
        const double dtDiurnal = daylight > 0.0 ? daylight / static_cast<double>(nDiurnalSteps) : 0.0;
        const double sunrise = 12.0 - daylight / 2.0;

        struct DiurnalStep
        {
            double hourUTC;
            double localHour;
            Solar::ClearSkyRadiation rad;
            std::shared_ptr<cTMapf> beamFactor;
            double tempFactor; // Normalized diurnal temperature curve factor [0.0, 1.0]
        };

        std::vector<DiurnalStep> diurnalSteps;
        double sumDaylightGHI = 0.0;
        if (daylight > 0.0)
        {
            diurnalSteps.reserve(nDiurnalSteps);
            for (int s = 0; s < nDiurnalSteps; ++s)
            {
                double solarHr = sunrise + (static_cast<double>(s) + 0.5) * dtDiurnal;
                double utcHr = solarHr - (longitudeDeg / 15.0);
                while (utcHr < 0.0) utcHr += 24.0;
                while (utcHr >= 24.0) utcHr -= 24.0;

                Solar::ClearSkyRadiation sRad = Solar::clearSkySolarRadiation(latitudeDeg, longitudeDeg, year, month, day, utcHr, avgElev);
                if (sRad.elevationDeg > 0.0)
                {
                    DiurnalStep step;
                    step.hourUTC = utcHr;
                    step.localHour = solarHr;
                    step.rad = sRad;
                    step.beamFactor = Solar::shadowSweepCanopy(DEM, BH, CH, sRad.azimuthDeg, sRad.elevationDeg, canopyTransmissivity);

                    // Diurnal sinusoidal temperature factor peaking at 15:00 local time
                    double timeFromPeak = std::abs(step.localHour - 15.0);
                    step.tempFactor = std::max(0.0, 1.0 - (timeFromPeak / 8.0) * (timeFromPeak / 8.0));
                    sumDaylightGHI += sRad.ghi;
                    diurnalSteps.push_back(step);
                }
            }
        }

        // Daily average global radiation Q_gl_avg during daylight hours (default: 340 W/m^2 per DPRA report)
        double qglAvg = (!diurnalSteps.empty()) ? (sumDaylightGHI / static_cast<double>(diurnalSteps.size())) : 340.0;

        // Peak local solar time and sun daily motion factor f_sun
        double peakLocalHour = peakHourUTC + (longitudeDeg / 15.0);
        while (peakLocalHour < 0.0) peakLocalHour += 24.0;
        while (peakLocalHour >= 24.0) peakLocalHour -= 24.0;
        double fsunPeak = sunDailyMotionFactor(peakLocalHour);

        std::cout << "[FastHeat] Computing peak biometeorology (UHI, Ta, PET, UTCI, WBGT, Tmrt)..." << std::endl;

        // Parallel processing across rows
        GetTPool().parallel_for(0, rows, [&](int r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(outputs.UTCI->data[r][c]);
                    pcr::setMV(outputs.WBGT->data[r][c]);
                    pcr::setMV(outputs.TMRT->data[r][c]);
                    pcr::setMV(outputs.WBGT_H38->data[r][c]);
                    pcr::setMV(outputs.TMRT_AVG->data[r][c]);
                    pcr::setMV(outputs.UTCI_AVG->data[r][c]);
                    pcr::setMV(outputs.PET->data[r][c]);
                    pcr::setMV(outputs.PET_AVG->data[r][c]);
                    pcr::setMV(outputs.UHI->data[r][c]);
                    pcr::setMV(outputs.TA->data[r][c]);
                    continue;
                }

                // Reference station temperature and RH at pixel
                double tStation = (TEMP && !pcr::isMV(TEMP->data[r][c])) ? static_cast<double>(TEMP->data[r][c]) : 35.0;
                double rh = (RH && !pcr::isMV(RH->data[r][c])) ? static_cast<double>(RH->data[r][c]) : 32.0;

                double svfVal = outputs.SVF ? static_cast<double>(outputs.SVF->data[r][c]) : 1.0;
                float beam = peakBeamFactor ? peakBeamFactor->data[r][c] : 1.0f;

                double svfAvgVal = (svfAvg && !pcr::isMV(svfAvg->data[r][c])) ? static_cast<double>(svfAvg->data[r][c]) : svfVal;
                double fvegAvgVal = (fvegAvg && !pcr::isMV(fvegAvg->data[r][c])) ? static_cast<double>(fvegAvg->data[r][c]) : 0.0;

                // Land use / surface type class and Bowen ratio for pixel
                int luVal = (LU && !pcr::isMV(LU->data[r][c])) ? static_cast<int>(std::round(LU->data[r][c])) : 50;
                float chVal = (CH && !pcr::isMV(CH->data[r][c])) ? std::max(0.0f, CH->data[r][c]) : 0.0f;
                double bowenVal = (bowenMap && !pcr::isMV(bowenMap->data[r][c])) ? static_cast<double>(bowenMap->data[r][c]) : 3.0;
                double albVal = (ALBEDO && !pcr::isMV(ALBEDO->data[r][c])) ? static_cast<double>(ALBEDO->data[r][c]) : -1.0;

                // 1. Urban Heat Island (UHI_max) calculation per Theeuwes et al. (2016) / DPRA
                double elev = static_cast<double>(DEM->data[r][c]);
                double deltaT = 15.0; // Standard diurnal range (DPRA default: 34 - 19 C)
                double uhiMax = calculateUHIMax(svfAvgVal, fvegAvgVal, tStation, elev, qglAvg, deltaT, windSpeedV10);
                outputs.UHI->data[r][c] = static_cast<float>(uhiMax);

                // 2. Peak Atmospheric Temperature Ta incorporating UHI effect:
                // Ta,t = Tstation,t + UHI_max * f_sun,t
                double taPeak = tStation + uhiMax * fsunPeak;
                outputs.TA->data[r][c] = static_cast<float>(taPeak);

                // 3. Peak Tmrt via SOLWEIG using localized atmospheric temperature Ta
                double tmrt = calculateTmrtSOLWEIG(taPeak, rh, svfVal, beam, peakRad.dni, peakRad.dhi, peakRad.elevationDeg, luVal, chVal, albVal);
                outputs.TMRT->data[r][c] = static_cast<float>(tmrt);

                // 4. Peak UTCI using Ta
                double ea = actualVaporPressureHPa(taPeak, rh);
                double utci = calculateUTCI(taPeak, tmrt, windSpeedV10, ea);
                outputs.UTCI->data[r][c] = static_cast<float>(utci);

                // 5. Peak WBGT using Ta
                double wbgt = calculateWBGT(taPeak, rh, tmrt, windSpeedV10);
                outputs.WBGT->data[r][c] = static_cast<float>(wbgt);

                // 6. Peak PET using Tygron / DPRA formulas
                double pet = calculatePET(taPeak, rh, svfVal, beam, peakRad.ghi, peakRad.elevationDeg, windSpeedV10, bowenVal, peakRad.dhi);
                outputs.PET->data[r][c] = static_cast<float>(pet);

                // 7. Diurnal cycle integration for daylight averages and heat stress hours
                double sumTmrt = 0.0;
                double sumUtci = 0.0;
                double sumPet = 0.0;
                float hoursAbove38 = 0.0f;
                for (const auto& step : diurnalSteps)
                {
                    double fsunStep = sunDailyMotionFactor(step.localHour);
                    double tStationStep = tStation - 10.0 * (1.0 - step.tempFactor);
                    double taStep = tStationStep + uhiMax * fsunStep;

                    float stepBeam = step.beamFactor ? step.beamFactor->data[r][c] : 1.0f;
                    double tmrtStep = calculateTmrtSOLWEIG(taStep, rh, svfVal, stepBeam, step.rad.dni, step.rad.dhi, step.rad.elevationDeg, luVal, chVal, albVal);
                    double eaStep = actualVaporPressureHPa(taStep, rh);
                    double utciStep = calculateUTCI(taStep, tmrtStep, windSpeedV10, eaStep);
                    double petStep = calculatePET(taStep, rh, svfVal, stepBeam, step.rad.ghi, step.rad.elevationDeg, windSpeedV10, bowenVal, step.rad.dhi);

                    sumTmrt += tmrtStep;
                    sumUtci += utciStep;
                    sumPet += petStep;

                    if (utciStep >= 38.0)
                    {
                        hoursAbove38 += static_cast<float>(dtDiurnal);
                    }
                }

                outputs.WBGT_H38->data[r][c] = std::clamp(hoursAbove38, 0.0f, static_cast<float>(daylight));

                if (!diurnalSteps.empty())
                {
                    double nSteps = static_cast<double>(diurnalSteps.size());
                    outputs.TMRT_AVG->data[r][c] = static_cast<float>(sumTmrt / nSteps);
                    outputs.UTCI_AVG->data[r][c] = static_cast<float>(sumUtci / nSteps);
                    outputs.PET_AVG->data[r][c] = static_cast<float>(sumPet / nSteps);
                }
                else
                {
                    outputs.TMRT_AVG->data[r][c] = static_cast<float>(tmrt);
                    outputs.UTCI_AVG->data[r][c] = static_cast<float>(utci);
                    outputs.PET_AVG->data[r][c] = static_cast<float>(pet);
                }
            }
        });

        GetTPool().waitext();

        std::cout << "[FastHeat] Simulation completed successfully." << std::endl;
        return 0;
    }

} // namespace FastHeatModel

#endif // FASTHEAT_H
