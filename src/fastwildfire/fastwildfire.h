//
// FastWildFire (FastFire) Model
// Directional GEV Extreme Fire Weather, Thermodynamic RH, Dual-Scenario Fosberg & McArthur,
// and Fast Sweeping Method (FSM) Eikonal Landscape Risk Solver
//

#ifndef FASTWILDFIRE_H
#define FASTWILDFIRE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include "../../common/threading/parallel.h"
#include "../../geospatial/map/map.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace FastWildfireModel
{

    // -------------------------------------------------------------------------
    // Constants and Directional Sector Definitions
    // -------------------------------------------------------------------------
    constexpr int NUM_SECTORS = 12;

    // Azimuth of sector k (k in 0..11): theta_k = k * 30 degrees (0 = North, 90 = East, 180 = South, 270 = West)
    inline double getSectorAzimuthDeg(int sector_idx)
    {
        return static_cast<double>((sector_idx % NUM_SECTORS) * 30);
    }

    // Downwind azimuth of sector k: fire spreads in downwind direction (theta_k + 180)
    inline double getDownwindAzimuthDeg(int sector_idx)
    {
        double dw = getSectorAzimuthDeg(sector_idx) + 180.0;
        while (dw >= 360.0) dw -= 360.0;
        return dw;
    }

    // -------------------------------------------------------------------------
    // GEV Extreme Value Distribution Quantile
    // -------------------------------------------------------------------------
    // Coles (2001) / standard GEV inverse CDF:
    // Q(R; mu, sigma, xi) = mu + (sigma / xi) * (y_R^(-xi) - 1)  if xi != 0
    //                     = mu - sigma * ln(y_R)                 if xi == 0
    // where y_R = -ln(1 - 1/R)
    inline double calculateGEVQuantile(double R_years, double mu, double sigma, double xi)
    {
        double R = std::max(1.0001, R_years);
        double p = 1.0 - (1.0 / R);
        double y = -std::log(p);

        if (std::abs(xi) > 1e-5)
        {
            return mu + (sigma / xi) * (std::pow(y, -xi) - 1.0);
        }
        else
        {
            return mu - sigma * std::log(y);
        }
    }

    // -------------------------------------------------------------------------
    // Thermodynamic Relative Humidity (RH) via Magnus-Tetens Formulation
    // -------------------------------------------------------------------------
    // Computes physically consistent RH (%) from Temperature (deg C) and Dew Point (deg C),
    // bounded to extreme dry floor RH_min = 3.0% and upper limit 100.0%.
    // Physical constraint: dew point cannot exceed ambient dry bulb temperature.
    inline double calculateRelativeHumidity(double tempC, double dewPointC)
    {
        double t = std::clamp(tempC, -30.0, 65.0);
        double td = std::clamp(dewPointC, -45.0, t - 0.1);

        double num = (17.625 * td) / (243.04 + td) - (17.625 * t) / (243.04 + t);
        double rh = 100.0 * std::exp(num);

        return std::clamp(rh, 3.0, 100.0);
    }

    // Climatological warm-season dew point estimation with latitude, elevation,
    // and soil moisture boundary-layer coupling (used when explicit dew point rasters are absent)
    inline double estimateWarmSeasonDewPoint(double latitudeDeg, double elevM = 0.0, double sm_volumetric = -1.0)
    {
        double absLat = std::abs(latitudeDeg);
        double td_base = 14.0;

        if (absLat < 20.0)
        {
            // Tropical maritime / humid equatorial zone: Warm ocean SSTs (26-29 C) sustain
            // high maritime boundary layer water vapor (typically 21.0 C to 24.5 C at sea level).
            double t_trop = 1.0 - (absLat / 20.0);
            td_base = 18.5 + 5.0 * t_trop; // 23.5 C at equator down to 18.5 C at 20 deg
        }
        else if (absLat >= 20.0 && absLat <= 35.0)
        {
            // Subtropical belt (Mediterranean / arid desert transition): 14.0 C down to 11.0 C
            td_base = 14.0 - ((absLat - 20.0) / 15.0) * 3.0;
        }
        else if (absLat > 35.0 && absLat <= 60.0)
        {
            // Temperate zone: 11.0 C down to 5.0 C
            td_base = 11.0 - ((absLat - 35.0) / 25.0) * 6.0;
        }
        else
        {
            // Boreal / polar
            td_base = std::max(0.0, 5.0 - ((absLat - 60.0) / 30.0) * 5.0);
        }

        // Soil moisture land-atmosphere feedback:
        // Moist soils (sm > 0.25) sustain high boundary-layer humidity via evapotranspiration.
        // Arid soils (sm < 0.10) depress dew point due to strong sensible heating.
        if (sm_volumetric > 0.0)
        {
            double sm = std::clamp(sm_volumetric, 0.04, 0.45);
            double delta_td = (sm - 0.20) * 14.0; // [-2.24 C to +3.5 C]
            td_base += delta_td;
        }

        // Environmental lapse rate correction for elevation (~5.5 C / 1000m)
        td_base -= (std::max(0.0, elevM) / 1000.0) * 5.5;

        return std::clamp(td_base, -20.0, 26.5);
    }

    // -------------------------------------------------------------------------
    // Effective Sustained Fire-Weather Wind Speed
    // -------------------------------------------------------------------------
    // Translates ambient 10m wind speed [km/h] (which for 5-50 year return periods
    // can reach 20 to 50 m/s = 72 to 180 km/h in storm/cyclone regimes) into the
    // effective sustained wind speed at fuel/flame scale.
    // In wildland fire science (McArthur 1967, Noble 1980, Rothermel 1972, Cheney 1998, Sullivan 2009):
    // - Below 45 km/h (~12.5 m/s): linear/unattenuated open wind response.
    // - Above 45 km/h: physical flame detachment, fuel bed convective cooling/blowout, and canopy
    //   boundary layer decoupling cause the effective combustion-driving wind to transition
    //   sublinearly, asymptotically approaching the physical flame-scale limit of extreme firestorms (~60 km/h).
    inline double calculateEffectiveFireWindKmH(double windKmH)
    {
        double V = std::max(0.0, windKmH);
        if (V <= 45.0)
        {
            return V;
        }
        // Smooth tanh saturation above 45 km/h towards physical flame-level limit (~60 km/h):
        return 45.0 + 15.0 * std::tanh((V - 45.0) / 15.0);
    }

    // -------------------------------------------------------------------------
    // Fosberg Fire Weather Index (FFWI)
    // -------------------------------------------------------------------------
    // Inputs:
    //   tempC   : Air Temperature [deg C]
    //   rhPct   : Relative Humidity [%]
    //   windKmH : 10m Wind Speed [km/h]
    // Output:
    //   FFWI dimensionless index [0, ~100+]
    inline double calculateFFWI(double tempC, double rhPct, double windKmH)
    {
        double TF = 1.8 * tempC + 32.0;
        double V_eff = calculateEffectiveFireWindKmH(windKmH);
        double Umph = 0.621371 * V_eff;
        double rh = std::clamp(rhPct, 1.0, 100.0);

        // Equilibrium Moisture Content m
        double m = 0.0;
        if (rh < 10.0)
        {
            m = 0.03229 + 0.281073 * rh - 0.000578 * rh * TF;
        }
        else if (rh <= 50.0)
        {
            m = 2.22749 + 0.160107 * rh - 0.01478 * TF;
        }
        else
        {
            m = 21.0606 + 0.005565 * rh * rh - 0.00035 * rh * TF - 0.483199 * rh;
        }
        m = std::max(0.01, m);

        // Moisture Damping factor eta (extinction at m >= 30%)
        double xm = m / 30.0;
        double eta = 0.0;
        if (xm < 1.0)
        {
            eta = 1.0 - 2.0 * xm + 1.5 * xm * xm - 0.5 * xm * xm * xm;
        }
        eta = std::clamp(eta, 0.0, 1.0);

        double ffwi = (eta * std::sqrt(1.0 + Umph * Umph)) / 0.3002;
        return std::max(0.0, ffwi);
    }

    // -------------------------------------------------------------------------
    // McArthur Forest Fire Danger Index (FFDI Mark 5)
    // -------------------------------------------------------------------------
    // Inputs:
    //   tempC         : Air Temperature [deg C]
    //   rhPct         : Relative Humidity [%]
    //   windKmH       : 10m Wind Speed [km/h]
    //   droughtFactor : Drought Factor [1.0 - 10.0]
    // Output:
    //   FFDI Mark 5 index [0, ~100+]
    inline double calculateFFDI(double tempC, double rhPct, double windKmH, double droughtFactor = 10.0)
    {
        double df = std::clamp(droughtFactor, 1.0, 10.0);
        double rh = std::clamp(rhPct, 3.0, 100.0);
        double V_eff = calculateEffectiveFireWindKmH(windKmH);
        double lnDF = std::log(df);

        double exponent = -0.450 + 0.987 * lnDF - 0.0345 * rh + 0.0338 * tempC + 0.0234 * V_eff;
        double ffdi = 2.0 * std::exp(exponent);
        return std::max(0.0, ffdi);
    }

    // -------------------------------------------------------------------------
    // McArthur Drought Factor DF from Volumetric Soil Moisture & Annual Stats
    // -------------------------------------------------------------------------
    // In Australian & international fire science (Griffiths 1998, Finkele et al. 2006):
    // The Drought Factor (1.0 to 10.0) represents the proportion of fine fuel available
    // for combustion, driven by soil water depletion.
    // In perennially moist/tropical regions (e.g. Caribbean with SM ~ 0.35-0.40 m3/m3),
    // DF drops to 1.0 - 2.0 (fuel availability ~ 0).
    // In dry summer / Mediterranean regimes (SM <= 0.08 m3/m3), DF reaches 9.5 - 10.0.
    inline double calculateDroughtFactorFromSM(double sm_volumetric, double sm_min = -1.0, double sm_var = -1.0)
    {
        double sm = std::clamp(sm_volumetric, 0.02, 0.60);
        // Effective saturation ratio: residual ~0.06, field capacity ~0.36
        double s_rel = std::clamp((sm - 0.06) / (0.36 - 0.06), 0.0, 1.0);
        double df_base = 1.0 + 9.0 * std::pow(1.0 - s_rel, 1.5);

        // Modulate with annual minimum soil moisture:
        // In the Caribbean, the driest month still has high soil moisture (sm_min > 0.22).
        // If the minimum soil moisture throughout the entire annual cycle remains high,
        // fuels never cure completely, establishing a ceiling on fuel availability:
        if (sm_min > 0.18)
        {
            double min_rel = std::clamp((sm_min - 0.06) / (0.36 - 0.06), 0.0, 1.0);
            double df_cap = 1.0 + 9.0 * std::pow(1.0 - min_rel, 2.0);
            df_base = std::min(df_base, df_cap);
        }

        return std::clamp(df_base, 1.0, 10.0);
    }

    // -------------------------------------------------------------------------
    // Soil Moisture Live & Dead Fuel Damping Factor (eta_SM)
    // -------------------------------------------------------------------------
    // Based on Rothermel (1972) / NFDRS / Krueger et al. (2015 "Soil moisture as
    // a predictor of wildfire danger", Int. J. Wildland Fire).
    // When soil moisture is high (> 0.25 m3/m3), live fuel moisture (LFMC) remains > 120%
    // and duff/litter layers stay humid, suppressing combustion rate of spread:
    // eta_SM in [0.02, 1.0].
    inline double calculateSoilMoistureDampingFactor(double sm_volumetric, double sm_min = -1.0, double sm_var = -1.0)
    {
        double sm = std::clamp(sm_volumetric, 0.02, 0.60);
        // Sigmoid damping centered at ~0.22 m3/m3 with transition scale 0.045
        double eta = 1.0 / (1.0 + std::exp((sm - 0.22) / 0.045));

        // Incorporate seasonal variation and annual minimum:
        // Perennially humid regions with sm_min > 0.20 and low seasonal amplitude
        // maintain lush understory foliage year-round, dampening fire spread even further:
        if (sm_min > 0.20)
        {
            double eta_min_cap = 1.0 / (1.0 + std::exp((sm_min - 0.20) / 0.040));
            eta = std::min(eta, eta_min_cap);
        }

        return std::clamp(eta, 0.02, 1.0);
    }

    // -------------------------------------------------------------------------
    // Normalized Hybrid Hazard Parameter (H_hyb) in [0, 1]
    // -------------------------------------------------------------------------
    // Maps standard meteorological danger indices (FFWI and FFDI on standard 0-100+ scale)
    // to a strictly bounded [0, 1] hazard rating with smooth, non-saturating asymptotic
    // growth that preserves full spatial contrast across the entire landscape.
    inline double calculateHybridHazard(double ffwi, double ffdi, double omega = 0.5)
    {
        // Smooth non-clipping mapping:
        // - idx = 0  -> 0.00
        // - idx = 25 -> 0.25 (High)
        // - idx = 50 -> 0.49 (Very High)
        // - idx = 75 -> 0.66 (Severe/Extreme)
        // - idx = 100 -> 0.78 (Catastrophic benchmark)
        // - idx = 150 -> 0.92
        // - idx = 250 -> 0.99
        auto normIndex = [](double idx) -> double {
            if (idx <= 0.0) return 0.0;
            double z = idx / 70.0;
            return std::clamp(1.0 - std::exp(-std::pow(z, 1.2)), 0.0, 1.0);
        };
        double norm_ffwi = normIndex(ffwi);
        double norm_ffdi = normIndex(ffdi);
        double w = std::clamp(omega, 0.0, 1.0);
        double h_hyb = w * norm_ffwi + (1.0 - w) * norm_ffdi;
        return std::clamp(h_hyb, 0.0, 1.0);
    }

    // -------------------------------------------------------------------------
    // Directional Factors: Wind Alignment & Topographic Slope
    // -------------------------------------------------------------------------
    // Wind Alignment & Lateral Flanking Factor A_W:
    // Uses the industry-standard Anderson (1983) / Alexander (1985) fire ellipse length-to-breadth ratio:
    // LB = 1.0 + 0.0776 * V_kmh
    // eccentricity e = sqrt(1.0 - 1.0 / (LB^2))
    // Ellipse polar radius relative to head spread:
    // Lambda_wind(phi) = (1.0 - e) / (1.0 - e * cos(phi))
    inline double calculateWindAlignmentFactor(double stepAzimuthDeg, double downwindAzimuthDeg, double windKmH, double cw = 0.0776)
    {
        double diff = std::abs(stepAzimuthDeg - downwindAzimuthDeg);
        while (diff > 180.0) diff = std::abs(diff - 360.0);
        double phiRad = diff * (M_PI / 180.0);

        // Sublinear wind saturation above 60 km/h bounds the ellipse length-to-breadth ratio (LB <= 8.0)
        // for extreme storm/cyclone winds (10 to 50 m/s = 36 to 180 km/h), preventing the ellipse from
        // collapsing into an unphysical 1-pixel needle.
        double V_eff = calculateEffectiveFireWindKmH(windKmH);
        double lb = 1.0 + cw * V_eff;
        double e_sq = std::max(0.0, 1.0 - 1.0 / (lb * lb));
        double e = std::sqrt(e_sq);

        double num = 1.0 - e;
        double den = std::max(1e-4, 1.0 - e * std::cos(phiRad));
        double aw = num / den;

        // Minimum lateral / backing spread floor (turbulent diffusion and radiation)
        return std::clamp(aw, 0.06, 1.0);
    }

    // -------------------------------------------------------------------------
    // Topographic Slope & Cliff Barrier Factor A_Z
    // -------------------------------------------------------------------------
    // Inputs:
    //   elevFromM : Elevation at origin cell [m]
    //   elevToM   : Elevation at destination cell [m]
    //   stepDistM : Horizontal distance between cells [m]
    //   cs        : McArthur / Noble slope exponent (default 3.533)
    // 
    // Physical Mechanics (Rothermel 1972, Cheney 1981, USFS Andrews 2012):
    // 1. Moderate slopes (0 to 30 deg): Flames preheat uphill fuels, accelerating
    //    forward rate of spread according to McArthur/Noble: exp(cs * tanBeta).
    // 2. Downhill travel (tanBeta < 0): Buoyant convection and flame tilt point away
    //    from the fuel bed, drastically retarding spread down to 0.15 - 0.35x.
    // 3. Steep terrain and cliffs (> 35 deg, and especially 45-75+ deg):
    //    - Plume detachment: Convective heat lifts vertically away from the face.
    //    - Rockiness / scree: Fuel continuity breaks down on vertical cliffs/ravines.
    //    - Topographic barrier penalty: Cliffs act as prominent firebreaks and 
    //      impassable barriers to contiguous flame propagation.
    inline double calculateTopographicSlopeFactor(double elevFromM, double elevToM, double stepDistM, double cs = 3.533)
    {
        if (stepDistM <= 1e-4) return 1.0;
        double tanBeta = (elevToM - elevFromM) / stepDistM;
        double absSlope = std::abs(tanBeta);
        double slopeDeg = std::atan(absSlope) * (180.0 / M_PI);

        // Uphill acceleration / downhill deceleration:
        // Uphill: capped at ~32 deg (tanBeta ~ 0.62) where convective flame preheating saturates.
        // Downhill: smooth retarding factor down to ~0.20x.
        double dirSlope = std::clamp(tanBeta, -0.65, 0.60);
        double az_dir = std::exp(cs * dirSlope);

        // Cliff and ruggedness barrier penalty P_cliff:
        // Sharp barrier penalty begins at 30 deg and accelerates up to cliffs (45-75+ deg).
        // Slopes > 50 deg reduce spread potential by 75% to 95%, acting as natural firebreaks.
        double p_cliff = 1.0;
        if (slopeDeg > 30.0)
        {
            double excess = (slopeDeg - 30.0) / 10.0;
            p_cliff = 1.0 / (1.0 + excess * excess);
        }

        return az_dir * p_cliff;
    }

    // -------------------------------------------------------------------------
    // Fuel Flammability from Land Use (ESA WorldCover)
    // -------------------------------------------------------------------------
    inline float defaultFuelFromLandUse(int luClass)
    {
        switch (luClass)
        {
            case 10: return 0.80f; // Tree cover (dense wildland fuels, continuous canopy & duff)
            case 95: return 0.30f; // Mangroves (wetland woody)
            case 20: return 0.70f; // Shrubland (heathland/chaparral/scrub, volatile oils)
            case 30: return 0.45f; // Grassland (rapid flashy surface fuels)
            case 40: return 0.20f; // Cropland (seasonal agriculture, low curing/sparse)
            case 90: return 0.05f; // Herbaceous wetland (humid boundary layer)
            case 50: return 0.02f; // Built-up / urban / streets (pavement, masonry, non-burnable roads)
            case 60: return 0.01f; // Bare / sparse vegetation (minimal fuel)
            case 70: return 0.0f;  // Snow and ice (non-flammable)
            case 80: return 0.0f;  // Permanent water bodies (non-flammable)
            default: return 0.15f;
        }
    }

    // -------------------------------------------------------------------------
    // Structures: Parameters and Model Outputs
    // -------------------------------------------------------------------------
    struct WildfireParameters
    {
        double return_period_ext = 50.0;    // Target extreme return period [years]
        double return_period_comp = 2.0;    // Companion baseline return period [years]
        double drought_factor = 10.0;       // McArthur Drought Factor (1 - 10)
        double tau_0 = 3600.0;              // Characteristic spatial decay scale [seconds / resistance units]
        double omega = 0.5;                 // Fosberg vs McArthur blend weight [0.0 - 1.0]
        double delta_T_cc = 0.0;            // Climate change additive temperature offset [deg C]
        double gamma_U = 1.0;               // Climate change multiplicative wind factor
        double delta_Td_cc = 0.0;           // Climate change additive dew point offset [deg C]
        double cw = 0.0776;                 // Anderson fire ellipse length-to-breadth wind coefficient [h/km]
        double cs = 3.533;                  // Topographic slope exponent
        int num_sweeps = 6;                 // Fast Sweeping iterations
        int single_sector = -1;             // -1 for all 12 sectors, or 0..11
        bool wildland_as_source = true;     // Initialize upwind fuel beds as sources
        float fuel_source_threshold = 0.4f; // Minimum fuel flammability to act as wildland ignition bed
        bool use_soil_moisture = true;      // Incorporate soil moisture and annual variation damping
        bool log_verbose = true;
    };

    struct SectorResults
    {
        int sector_index = 0;
        double azimuth_deg = 0.0;
        double downwind_azimuth_deg = 0.0;

        std::shared_ptr<cTMapf> Wind_ext = nullptr;      // 10m wind speed extreme [km/h]
        std::shared_ptr<cTMapf> Wind_comp = nullptr;     // 10m wind speed companion [km/h]

        std::shared_ptr<cTMapf> FFWI_H = nullptr;        // Fosberg under Heat-dominant regime
        std::shared_ptr<cTMapf> FFWI_W = nullptr;        // Fosberg under Wind-dominant regime
        std::shared_ptr<cTMapf> FFWI_max = nullptr;      // Envelope max Fosberg

        std::shared_ptr<cTMapf> FFDI_H = nullptr;        // McArthur under Heat-dominant regime
        std::shared_ptr<cTMapf> FFDI_W = nullptr;        // McArthur under Wind-dominant regime
        std::shared_ptr<cTMapf> FFDI_max = nullptr;      // Envelope max McArthur

        std::shared_ptr<cTMapf> H_hyb_H = nullptr;       // Hybrid Hazard Heat regime
        std::shared_ptr<cTMapf> H_hyb_W = nullptr;       // Hybrid Hazard Wind regime
        std::shared_ptr<cTMapf> H_hyb_max = nullptr;     // Hybrid Hazard envelope max

        std::shared_ptr<cTMapf> Regime_k = nullptr;      // Dominant regime: 0 = Heat (H), 1 = Wind (W)
        std::shared_ptr<cTMapf> S_max = nullptr;         // Local spread potential S_max
        std::shared_ptr<cTMapf> Resistance_T = nullptr;  // Cumulative resistance surface T from FSM
        std::shared_ptr<cTMapf> ConnectedRisk_J = nullptr;// Connected Landscape Risk J in [0, 1]
    };

    struct WildfireModelOutputs
    {
        // Thermodynamic state
        std::shared_ptr<cTMapf> Temp_ext = nullptr;          // Extreme temperature [deg C]
        std::shared_ptr<cTMapf> Temp_comp = nullptr;         // Companion temperature [deg C]
        std::shared_ptr<cTMapf> DewPoint = nullptr;          // Dew point [deg C]
        std::shared_ptr<cTMapf> RH_ext = nullptr;            // Relative humidity during extreme heat [%]
        std::shared_ptr<cTMapf> RH_comp = nullptr;           // Relative humidity during companion heat [%]

        // Soil moisture & fuel drying state
        std::shared_ptr<cTMapf> Soil_Moisture = nullptr;     // Volumetric soil moisture [m3/m3]
        std::shared_ptr<cTMapf> Drought_Factor = nullptr;    // Realized McArthur Drought Factor DF [1.0 - 10.0]
        std::shared_ptr<cTMapf> SM_Damping = nullptr;        // Soil moisture spread damping factor eta_SM [0.02 - 1.0]

        // Landscape fuel
        std::shared_ptr<cTMapf> Fuel_Flammability = nullptr; // Fuel flammability F(x) in [0.01, 1.0]

        // Directional sector suite (12 directions)
        std::vector<SectorResults> sectors;

        // Omnidirectional Envelopes
        std::shared_ptr<cTMapf> FFWI_env = nullptr;          // Omnidirectional envelope realized FFWI (spatially modulated)
        std::shared_ptr<cTMapf> FFDI_env = nullptr;          // Omnidirectional envelope realized FFDI (spatially modulated)
        std::shared_ptr<cTMapf> FFWI_met_env = nullptr;      // Omnidirectional envelope pure meteorological FFWI
        std::shared_ptr<cTMapf> FFDI_met_env = nullptr;      // Omnidirectional envelope pure meteorological FFDI
        std::shared_ptr<cTMapf> H_hyb_env = nullptr;         // Omnidirectional envelope max Hybrid Hazard
        std::shared_ptr<cTMapf> S_env = nullptr;             // Omnidirectional envelope max Spread Potential
        std::shared_ptr<cTMapf> ConnectedRisk_env = nullptr; // Omnidirectional Connected Risk J_env [0, 1]
        std::shared_ptr<cTMapf> Dominant_Regime_env = nullptr; // Dominant regime of max sector (0=Heat, 1=Wind)

        // Adaptation Benefit
        std::shared_ptr<cTMapf> Adaptation_Benefit_dJ = nullptr; // Delta J = J_baseline - J_adapted
    };

    // -------------------------------------------------------------------------
    // Fast Sweeping Method (FSM) Anisotropic Eikonal Solver
    // -------------------------------------------------------------------------
    // Fast Sweeping Method (FSM) 8-Stencil Anisotropic Eikonal Solver
    // -------------------------------------------------------------------------
    // Solves |nabla T| = R on a 2D grid with anisotropic directional resistance
    // across all 8 cardinal and lateral/diagonal directions (N, NE, E, SE, S, SW, W, NW)
    // using the Anderson (1983) fire ellipse and 8 triangular simplex updates
    // for seamless rotational symmetry without Cartesian diamond artifacts.
    inline std::shared_ptr<cTMapf> SolveFastSweepingEikonal(
        const std::shared_ptr<cTMapf>& DEM,
        const std::shared_ptr<cTMapf>& S_map,
        const std::shared_ptr<cTMapf>& Wind_map,
        double downwindAzimuthDeg,
        const std::shared_ptr<cTMapf>& IgnitionMask,
        const WildfireParameters& params)
    {
        if (!DEM || !S_map) return nullptr;

        const int rows = DEM->nrRows();
        const int cols = DEM->nrCols();
        const double dx = std::max(0.1, std::abs(DEM->cellSizeX()));
        const double dy = std::max(0.1, std::abs(DEM->cellSizeY()));
        const double diagDist = std::sqrt(dx * dx + dy * dy);

        // Precompute diagonal step compass azimuths (from neighbor to current cell)
        double az_NW = std::atan2(dx, -dy) * (180.0 / M_PI); // 135 deg
        if (az_NW < 0.0) az_NW += 360.0;
        double az_NE = std::atan2(-dx, -dy) * (180.0 / M_PI); // 225 deg
        if (az_NE < 0.0) az_NE += 360.0;
        double az_SW = std::atan2(dx, dy) * (180.0 / M_PI); // 45 deg
        if (az_SW < 0.0) az_SW += 360.0;
        double az_SE = std::atan2(-dx, dy) * (180.0 / M_PI); // 315 deg
        if (az_SE < 0.0) az_SE += 360.0;

        std::shared_ptr<cTMapf> T_map = DEM->GetCopySP0();
        T_map->setAllMV();

        const double cw = params.cw;
        const double cs = params.cs;
        const int numSweeps = std::max(4, params.num_sweeps);

        // 1. Initialize T: 0.0 at source cells, MV elsewhere
        // Compute wind direction components to determine upwind physical neighbors
        double windAzimuthDeg = downwindAzimuthDeg - 180.0;
        while (windAzimuthDeg < 0.0) windAzimuthDeg += 360.0;
        while (windAzimuthDeg >= 360.0) windAzimuthDeg -= 360.0;

        double radWind = windAzimuthDeg * (M_PI / 180.0);
        double sinW = std::sin(radWind);
        double cosW = std::cos(radWind);

        // In raster coordinates, row decreases Northwards (-r), col increases Eastwards (+c):
        int dr_up = (std::abs(cosW) > 0.35) ? -static_cast<int>(std::round(cosW / std::abs(cosW))) : 0;
        int dc_up = (std::abs(sinW) > 0.35) ? static_cast<int>(std::round(sinW / std::abs(sinW))) : 0;

        auto isBarrierCell = [&](int row_idx, int col_idx) -> bool
        {
            if (row_idx < 0 || row_idx >= rows || col_idx < 0 || col_idx >= cols) return false;
            if (pcr::isMV(DEM->data[row_idx][col_idx]) || DEM->data[row_idx][col_idx] <= 0.0f) return true;
            if (pcr::isMV(S_map->data[row_idx][col_idx]) || S_map->data[row_idx][col_idx] <= 0.03f) return true;
            return false;
        };

        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]) || DEM->data[r][c] <= 0.0f)
                {
                    continue; // Water / ocean is NEVER an ignition or fuel source
                }

                float sVal = S_map->data[r][c];
                if (pcr::isMV(sVal) || sVal < 0.02f)
                {
                    continue; // Non-flammable / bare cells cannot ignite or propagate
                }

                bool is_source = false;

                // Explicit ignition source mask (targeted fire spread scenario)
                if (IgnitionMask)
                {
                    if (!pcr::isMV(IgnitionMask->data[r][c]) && IgnitionMask->data[r][c] > 0.0f)
                    {
                        is_source = true;
                    }
                }
                else
                {
                    // Autohazard barrier-shielding mode:
                    // Only cells immediately downwind of a real physical barrier (water, rock, fuelbreak)
                    // are initialized as barrier-shadow sources (T = 0.0).
                    // Continuous forest across domain borders does NOT act as a barrier, preventing artificial edge halos!
                    bool upwind_has_barrier = false;
                    if (isBarrierCell(r + dr_up, c + dc_up)) upwind_has_barrier = true;
                    if (dr_up != 0 && isBarrierCell(r + dr_up, c)) upwind_has_barrier = true;
                    if (dc_up != 0 && isBarrierCell(r, c + dc_up)) upwind_has_barrier = true;

                    if (upwind_has_barrier)
                    {
                        is_source = true;
                    }
                }

                if (is_source)
                {
                    T_map->data[r][c] = 0.0f;
                }
            }
        }

        // Helper lambda to compute directional resistance between neighbor (rn, cn) -> current (r, c)
        auto computeResistance = [&](int rn, int cn, int rc, int cc, double stepDistM, double stepAzimuthDeg) -> double
        {
            if (DEM->data[rn][cn] <= 0.0f || DEM->data[rc][cc] <= 0.0f) return 1e9; // Water is impassable

            float sFrom = S_map->data[rn][cn];
            if (pcr::isMV(sFrom) || sFrom <= 1e-4f) return 1e9;

            float zFrom = DEM->data[rn][cn];
            float zTo = DEM->data[rc][cc];
            if (pcr::isMV(zFrom) || pcr::isMV(zTo)) return 1e9;

            double wKmH = 20.0;
            if (Wind_map && !pcr::isMV(Wind_map->data[rn][cn]))
            {
                wKmH = static_cast<double>(Wind_map->data[rn][cn]);
            }

            double aw = calculateWindAlignmentFactor(stepAzimuthDeg, downwindAzimuthDeg, wKmH, cw);
            double az = calculateTopographicSlopeFactor(zFrom, zTo, stepDistM, cs);

            double effectiveSpread = static_cast<double>(sFrom) * aw * az;
            return 1.0 / std::max(1e-4, effectiveSpread);
        };

        // 2. Perform alternating 8-stencil Gauss-Seidel sweeps
        for (int sweep = 0; sweep < numSweeps; ++sweep)
        {
            int r_start = (sweep % 4 == 0 || sweep % 4 == 1) ? 0 : rows - 1;
            int r_end   = (sweep % 4 == 0 || sweep % 4 == 1) ? rows : -1;
            int r_step  = (sweep % 4 == 0 || sweep % 4 == 1) ? 1 : -1;

            int c_start = (sweep % 4 == 0 || sweep % 4 == 3) ? 0 : cols - 1;
            int c_end   = (sweep % 4 == 0 || sweep % 4 == 3) ? cols : -1;
            int c_step  = (sweep % 4 == 0 || sweep % 4 == 3) ? 1 : -1;

            for (int r = r_start; r != r_end; r += r_step)
            {
                for (int c = c_start; c != c_end; c += c_step)
                {
                    if (pcr::isMV(DEM->data[r][c]) || DEM->data[r][c] <= 0.0f)
                    {
                        continue; // Skip water cells completely
                    }
                    if (!pcr::isMV(T_map->data[r][c]) && T_map->data[r][c] == 0.0f)
                    {
                        continue;
                    }

                    double T_cand_best = -1.0;

                    // Fetch raw arrival times at 8 adjacent neighbors:
                    double val_N = (r > 0 && !pcr::isMV(T_map->data[r - 1][c])) ? static_cast<double>(T_map->data[r - 1][c]) : -1.0;
                    double val_S = (r + 1 < rows && !pcr::isMV(T_map->data[r + 1][c])) ? static_cast<double>(T_map->data[r + 1][c]) : -1.0;
                    double val_W = (c > 0 && !pcr::isMV(T_map->data[r][c - 1])) ? static_cast<double>(T_map->data[r][c - 1]) : -1.0;
                    double val_E = (c + 1 < cols && !pcr::isMV(T_map->data[r][c + 1])) ? static_cast<double>(T_map->data[r][c + 1]) : -1.0;

                    double val_NW = (r > 0 && c > 0 && !pcr::isMV(T_map->data[r - 1][c - 1])) ? static_cast<double>(T_map->data[r - 1][c - 1]) : -1.0;
                    double val_NE = (r > 0 && c + 1 < cols && !pcr::isMV(T_map->data[r - 1][c + 1])) ? static_cast<double>(T_map->data[r - 1][c + 1]) : -1.0;
                    double val_SW = (r + 1 < rows && c > 0 && !pcr::isMV(T_map->data[r + 1][c - 1])) ? static_cast<double>(T_map->data[r + 1][c - 1]) : -1.0;
                    double val_SE = (r + 1 < rows && c + 1 < cols && !pcr::isMV(T_map->data[r + 1][c + 1])) ? static_cast<double>(T_map->data[r + 1][c + 1]) : -1.0;

                    // 1. 1D Discrete Ray Arrivals along 8 compass directions
                    if (val_N >= 0.0)
                    {
                        double T_1D = val_N + dy * computeResistance(r - 1, c, r, c, dy, 180.0);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }
                    if (val_S >= 0.0)
                    {
                        double T_1D = val_S + dy * computeResistance(r + 1, c, r, c, dy, 0.0);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }
                    if (val_W >= 0.0)
                    {
                        double T_1D = val_W + dx * computeResistance(r, c - 1, r, c, dx, 90.0);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }
                    if (val_E >= 0.0)
                    {
                        double T_1D = val_E + dx * computeResistance(r, c + 1, r, c, dx, 270.0);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }
                    if (val_NW >= 0.0)
                    {
                        double T_1D = val_NW + diagDist * computeResistance(r - 1, c - 1, r, c, diagDist, az_NW);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }
                    if (val_NE >= 0.0)
                    {
                        double T_1D = val_NE + diagDist * computeResistance(r - 1, c + 1, r, c, diagDist, az_NE);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }
                    if (val_SW >= 0.0)
                    {
                        double T_1D = val_SW + diagDist * computeResistance(r + 1, c - 1, r, c, diagDist, az_SW);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }
                    if (val_SE >= 0.0)
                    {
                        double T_1D = val_SE + diagDist * computeResistance(r + 1, c + 1, r, c, diagDist, az_SE);
                        if (T_cand_best < 0.0 || T_1D < T_cand_best) T_cand_best = T_1D;
                    }

                    // 2. 2D Continuous Planar Wavefront Solver for the 4 Cartesian Quadrants:
                    // Solves ((T - u)/dx)^2 + ((T - v)/dy)^2 = R^2 for incoming wavefronts between u and v
                    auto solveQuadrant2D = [&](double u, double v, double R_quad)
                    {
                        if (u < 0.0 || v < 0.0 || R_quad <= 0.0) return;
                        double inv_dx2 = 1.0 / (dx * dx);
                        double inv_dy2 = 1.0 / (dy * dy);
                        double A = inv_dx2 + inv_dy2;
                        double B = -2.0 * (u * inv_dx2 + v * inv_dy2);
                        double C = (u * u) * inv_dx2 + (v * v) * inv_dy2 - (R_quad * R_quad);
                        double discr = B * B - 4.0 * A * C;
                        if (discr >= 0.0)
                        {
                            double T_quad = (-B + std::sqrt(discr)) / (2.0 * A);
                            // Causality condition: arrival must be later than both incoming wave points
                            if (T_quad > u && T_quad > v)
                            {
                                if (T_cand_best < 0.0 || T_quad < T_cand_best) T_cand_best = T_quad;
                            }
                        }
                    };

                    // Northwest Quadrant (W, N)
                    if (val_W >= 0.0 && val_N >= 0.0)
                    {
                        double R_NW = computeResistance(r - 1, c - 1, r, c, diagDist, az_NW);
                        solveQuadrant2D(val_W, val_N, R_NW);
                    }
                    // Northeast Quadrant (E, N)
                    if (val_E >= 0.0 && val_N >= 0.0)
                    {
                        double R_NE = computeResistance(r - 1, c + 1, r, c, diagDist, az_NE);
                        solveQuadrant2D(val_E, val_N, R_NE);
                    }
                    // Southwest Quadrant (W, S)
                    if (val_W >= 0.0 && val_S >= 0.0)
                    {
                        double R_SW = computeResistance(r + 1, c - 1, r, c, diagDist, az_SW);
                        solveQuadrant2D(val_W, val_S, R_SW);
                    }
                    // Southeast Quadrant (E, S)
                    if (val_E >= 0.0 && val_S >= 0.0)
                    {
                        double R_SE = computeResistance(r + 1, c + 1, r, c, diagDist, az_SE);
                        solveQuadrant2D(val_E, val_S, R_SE);
                    }

                    // 3. 8 Triangular Simplex Updates between adjacent cardinal and diagonal neighbors:
                    // Guarantees C2-smooth isotropic / rotated elliptical wavefronts in all 360 degrees
                    // eliminating Cartesian rectangular grid artifacts.
                    auto solveSimplexHoriz = [&](double Tc, double Td, double R_tri)
                    {
                        if (Tc < 0.0 || Td < 0.0 || R_tri <= 0.0) return;
                        double p = (Tc - Td) / dx;
                        if (std::abs(p) <= R_tri)
                        {
                            double q = std::sqrt(std::max(0.0, R_tri * R_tri - p * p));
                            double T_simp = Tc + dy * q;
                            if (T_simp >= Tc && T_simp >= Td)
                            {
                                if (T_cand_best < 0.0 || T_simp < T_cand_best) T_cand_best = T_simp;
                            }
                        }
                    };

                    auto solveSimplexVert = [&](double Tc, double Td, double R_tri)
                    {
                        if (Tc < 0.0 || Td < 0.0 || R_tri <= 0.0) return;
                        double q = (Tc - Td) / dy;
                        if (std::abs(q) <= R_tri)
                        {
                            double p = std::sqrt(std::max(0.0, R_tri * R_tri - q * q));
                            double T_simp = Tc + dx * p;
                            if (T_simp >= Tc && T_simp >= Td)
                            {
                                if (T_cand_best < 0.0 || T_simp < T_cand_best) T_cand_best = T_simp;
                            }
                        }
                    };

                    if (val_N >= 0.0 && val_NE >= 0.0) solveSimplexHoriz(val_N, val_NE, computeResistance(r - 1, c + 1, r, c, diagDist, az_NE));
                    if (val_N >= 0.0 && val_NW >= 0.0) solveSimplexHoriz(val_N, val_NW, computeResistance(r - 1, c - 1, r, c, diagDist, az_NW));
                    if (val_S >= 0.0 && val_SE >= 0.0) solveSimplexHoriz(val_S, val_SE, computeResistance(r + 1, c + 1, r, c, diagDist, az_SE));
                    if (val_S >= 0.0 && val_SW >= 0.0) solveSimplexHoriz(val_S, val_SW, computeResistance(r + 1, c - 1, r, c, diagDist, az_SW));

                    if (val_E >= 0.0 && val_NE >= 0.0) solveSimplexVert(val_E, val_NE, computeResistance(r - 1, c + 1, r, c, diagDist, az_NE));
                    if (val_E >= 0.0 && val_SE >= 0.0) solveSimplexVert(val_E, val_SE, computeResistance(r + 1, c + 1, r, c, diagDist, az_SE));
                    if (val_W >= 0.0 && val_NW >= 0.0) solveSimplexVert(val_W, val_NW, computeResistance(r - 1, c - 1, r, c, diagDist, az_NW));
                    if (val_W >= 0.0 && val_SW >= 0.0) solveSimplexVert(val_W, val_SW, computeResistance(r + 1, c - 1, r, c, diagDist, az_SW));

                    if (T_cand_best >= 0.0)
                    {
                        if (pcr::isMV(T_map->data[r][c]) || T_cand_best < static_cast<double>(T_map->data[r][c]))
                        {
                            T_map->data[r][c] = static_cast<float>(T_cand_best);
                        }
                    }
                }
            }
        }

        return T_map;
    }

    // -------------------------------------------------------------------------
    // Main RunFastWildfire Simulation Pipeline
    // -------------------------------------------------------------------------
    inline int RunFastWildfire(
        const std::shared_ptr<cTMapf>& DEM,
        const std::shared_ptr<cTMapf>& Fuel_input,
        const std::shared_ptr<cTMapf>& LU,
        const std::shared_ptr<cTMapf>& Temp_ext_in,
        const std::shared_ptr<cTMapf>& Temp_comp_in,
        const std::shared_ptr<cTMapf>& DewPoint_in,
        const std::vector<std::shared_ptr<cTMapf>>& Wind_ext_dirs,
        const std::vector<std::shared_ptr<cTMapf>>& Wind_comp_dirs,
        const std::vector<std::shared_ptr<cTMapf>>& FastWind_dirs,
        const std::shared_ptr<cTMapf>& IgnitionMask,
        const WildfireParameters& params,
        WildfireModelOutputs& outputs,
        const std::shared_ptr<cTMapf>& SoilMoisture_in = nullptr,
        const std::shared_ptr<cTMapf>& SoilMoisture_min_in = nullptr,
        const std::shared_ptr<cTMapf>& SoilMoisture_var_in = nullptr)
    {
        if (!DEM)
        {
            std::cerr << "[FastWildfire] Error: DEM raster is required." << std::endl;
            return -1;
        }

        const int rows = DEM->nrRows();
        const int cols = DEM->nrCols();

        if (params.log_verbose)
        {
            std::cout << "[FastWildfire] Initializing simulation on domain (" << rows << "x" << cols << ")..." << std::endl;
            std::cout << "               Return Period Ext:  " << params.return_period_ext << " years" << std::endl;
            std::cout << "               Return Period Comp: " << params.return_period_comp << " years" << std::endl;
            std::cout << "               Drought Factor DF:  " << params.drought_factor << std::endl;
            std::cout << "               Hybrid Weight omega: " << params.omega << std::endl;
            std::cout << "               FSM Decay Scale tau: " << params.tau_0 << " s" << std::endl;
        }

        // Compute reference elevation for terrain lapse rate downscaling
        double elev_sum = 0.0;
        int elev_cnt = 0;
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (!pcr::isMV(DEM->data[r][c]) && DEM->data[r][c] > 0.0f)
                {
                    elev_sum += DEM->data[r][c];
                    elev_cnt++;
                }
            }
        }
        double ref_elev = (elev_cnt > 0) ? (elev_sum / elev_cnt) : 0.0;

        // 1. Prepare Fuel Flammability Layer F(x) in [0.01, 1.0] (0.0 for water)
        outputs.Fuel_Flammability = DEM->GetCopySP0();
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(outputs.Fuel_Flammability->data[r][c]);
                    continue;
                }

                // Water / ocean detection
                bool is_water = (DEM->data[r][c] <= 0.0f);
                if (LU && !pcr::isMV(LU->data[r][c]))
                {
                    int luClass = static_cast<int>(std::round(LU->data[r][c]));
                    if (luClass == 80) is_water = true;
                }

                if (is_water)
                {
                    outputs.Fuel_Flammability->data[r][c] = 0.0f;
                    continue;
                }

                float fVal = 0.30f;
                if (Fuel_input && !pcr::isMV(Fuel_input->data[r][c]))
                {
                    fVal = Fuel_input->data[r][c];
                }
                else if (LU && !pcr::isMV(LU->data[r][c]))
                {
                    int luClass = static_cast<int>(std::round(LU->data[r][c]));
                    fVal = defaultFuelFromLandUse(luClass);
                }

                outputs.Fuel_Flammability->data[r][c] = std::clamp(fVal, 0.01f, 1.0f);
            }
        }

        // 2. Prepare Temperature Quantiles: T_ext and T_comp with topographic lapse rate (-6.5 C / 1000m)
        outputs.Temp_ext = DEM->GetCopySP0();
        outputs.Temp_comp = DEM->GetCopySP0();

        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(outputs.Temp_ext->data[r][c]);
                    pcr::setMV(outputs.Temp_comp->data[r][c]);
                    continue;
                }

                float tExt = 36.0f;
                float tComp = 30.0f;

                if (Temp_ext_in && !pcr::isMV(Temp_ext_in->data[r][c]))
                {
                    tExt = Temp_ext_in->data[r][c];
                }
                if (Temp_comp_in && !pcr::isMV(Temp_comp_in->data[r][c]))
                {
                    tComp = Temp_comp_in->data[r][c];
                }
                else if (Temp_ext_in && !pcr::isMV(Temp_ext_in->data[r][c]))
                {
                    // Default companion baseline is ~6 C cooler than target extreme
                    tComp = tExt - 6.0f;
                }

                // Apply climate change delta
                tExt += static_cast<float>(params.delta_T_cc);
                tComp += static_cast<float>(params.delta_T_cc);

                // Topographic environmental lapse rate (-6.5 C per 1000m)
                // Downscale temperature according to local elevation relative to domain reference:
                double z = std::max(0.0, static_cast<double>(DEM->data[r][c]));
                double lapse_dT = -0.0065 * (z - ref_elev);
                tExt += static_cast<float>(lapse_dT);
                tComp += static_cast<float>(lapse_dT);

                outputs.Temp_ext->data[r][c] = tExt;
                outputs.Temp_comp->data[r][c] = tComp;
            }
        }

        // 3. Prepare Soil Moisture, Dynamic Drought Factor DF, and Live/Dead Fuel Damping eta_SM
        // Computed prior to Dew Point so that surface moisture can inform boundary-layer humidity.
        outputs.Soil_Moisture = DEM->GetCopySP0();
        outputs.Drought_Factor = DEM->GetCopySP0();
        outputs.SM_Damping = DEM->GetCopySP0();

        LSMVector2 center = DEM->GetPixelLngLat(rows / 2, cols / 2);
        double domainLat = std::isfinite(center.y) ? center.y : 35.0;
        double absLat = std::abs(domainLat);

        // Climatic default when soil moisture layer is absent:
        // - Deep tropics (|lat| < 20 deg): default DF ~ 3.5, eta_SM ~ 0.40, SM ~ 0.30
        // - Subtropical/arid (|lat| 20-35 deg): default DF ~ 8.5, eta_SM ~ 0.85, SM ~ 0.12
        // - Temperate (|lat| > 35 deg): default DF ~ 5.5, eta_SM ~ 0.60, SM ~ 0.22
        float default_sm = 0.20f;
        float default_df = static_cast<float>(params.drought_factor);
        float default_eta = 0.70f;
        if ((!SoilMoisture_in || !params.use_soil_moisture) && params.drought_factor >= 9.99)
        {
            if (absLat < 20.0)
            {
                default_sm = 0.30f;
                default_df = 3.5f;
                default_eta = 0.35f;
            }
            else if (absLat >= 20.0 && absLat <= 35.0)
            {
                default_sm = 0.12f;
                default_df = 8.5f;
                default_eta = 0.85f;
            }
            else
            {
                default_sm = 0.22f;
                default_df = 5.5f;
                default_eta = 0.60f;
            }
        }

        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]) || DEM->data[r][c] <= 0.0f)
                {
                    pcr::setMV(outputs.Soil_Moisture->data[r][c]);
                    pcr::setMV(outputs.Drought_Factor->data[r][c]);
                    pcr::setMV(outputs.SM_Damping->data[r][c]);
                    continue;
                }

                if (params.use_soil_moisture && SoilMoisture_in && !pcr::isMV(SoilMoisture_in->data[r][c]))
                {
                    float sm_val = SoilMoisture_in->data[r][c];
                    float sm_min_val = (SoilMoisture_min_in && !pcr::isMV(SoilMoisture_min_in->data[r][c]))
                                       ? SoilMoisture_min_in->data[r][c] : -1.0f;
                    float sm_var_val = (SoilMoisture_var_in && !pcr::isMV(SoilMoisture_var_in->data[r][c]))
                                       ? SoilMoisture_var_in->data[r][c] : -1.0f;

                    double df_val = calculateDroughtFactorFromSM(sm_val, sm_min_val, sm_var_val);
                    double eta_val = calculateSoilMoistureDampingFactor(sm_val, sm_min_val, sm_var_val);

                    outputs.Soil_Moisture->data[r][c] = sm_val;
                    outputs.Drought_Factor->data[r][c] = static_cast<float>(df_val);
                    outputs.SM_Damping->data[r][c] = static_cast<float>(eta_val);
                }
                else
                {
                    outputs.Soil_Moisture->data[r][c] = default_sm;
                    outputs.Drought_Factor->data[r][c] = default_df;
                    outputs.SM_Damping->data[r][c] = default_eta;
                }
            }
        }

        // 4. Prepare Dew Point T_d and Compute Thermodynamic RH_ext and RH_comp
        outputs.DewPoint = DEM->GetCopySP0();
        outputs.RH_ext = DEM->GetCopySP0();
        outputs.RH_comp = DEM->GetCopySP0();

        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(outputs.DewPoint->data[r][c]);
                    pcr::setMV(outputs.RH_ext->data[r][c]);
                    pcr::setMV(outputs.RH_comp->data[r][c]);
                    continue;
                }

                double z = std::max(0.0, static_cast<double>(DEM->data[r][c]));
                float td = 10.0f;
                if (DewPoint_in && !pcr::isMV(DewPoint_in->data[r][c]))
                {
                    td = DewPoint_in->data[r][c] - static_cast<float>(0.0020 * (z - ref_elev));
                }
                else
                {
                    double local_sm = (!outputs.Soil_Moisture || pcr::isMV(outputs.Soil_Moisture->data[r][c]))
                                      ? static_cast<double>(default_sm) : static_cast<double>(outputs.Soil_Moisture->data[r][c]);
                    td = static_cast<float>(estimateWarmSeasonDewPoint(domainLat, DEM->data[r][c], local_sm));
                }
                td += static_cast<float>(params.delta_Td_cc);
                outputs.DewPoint->data[r][c] = td;

                float tExt = outputs.Temp_ext->data[r][c];
                float tComp = outputs.Temp_comp->data[r][c];

                double rhExt = calculateRelativeHumidity(tExt, td);
                double rhComp = calculateRelativeHumidity(tComp, td);

                outputs.RH_ext->data[r][c] = static_cast<float>(rhExt);
                outputs.RH_comp->data[r][c] = static_cast<float>(rhComp);
            }
        }

        // 5. Initialize Omnidirectional Envelope Rasters
        outputs.FFWI_env = DEM->GetCopySP0();
        outputs.FFDI_env = DEM->GetCopySP0();
        outputs.FFWI_met_env = DEM->GetCopySP0();
        outputs.FFDI_met_env = DEM->GetCopySP0();
        outputs.H_hyb_env = DEM->GetCopySP0();
        outputs.S_env = DEM->GetCopySP0();
        outputs.ConnectedRisk_env = DEM->GetCopySP0();
        outputs.Dominant_Regime_env = DEM->GetCopySP0();

        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(outputs.FFWI_env->data[r][c]);
                    pcr::setMV(outputs.FFDI_env->data[r][c]);
                    pcr::setMV(outputs.FFWI_met_env->data[r][c]);
                    pcr::setMV(outputs.FFDI_met_env->data[r][c]);
                    pcr::setMV(outputs.H_hyb_env->data[r][c]);
                    pcr::setMV(outputs.S_env->data[r][c]);
                    pcr::setMV(outputs.ConnectedRisk_env->data[r][c]);
                    pcr::setMV(outputs.Dominant_Regime_env->data[r][c]);
                }
                else
                {
                    outputs.FFWI_env->data[r][c] = 0.0f;
                    outputs.FFDI_env->data[r][c] = 0.0f;
                    outputs.FFWI_met_env->data[r][c] = 0.0f;
                    outputs.FFDI_met_env->data[r][c] = 0.0f;
                    outputs.H_hyb_env->data[r][c] = 0.0f;
                    outputs.S_env->data[r][c] = 0.0f;
                    outputs.ConnectedRisk_env->data[r][c] = 0.0f;
                    outputs.Dominant_Regime_env->data[r][c] = 0.0f;
                }
            }
        }

        outputs.sectors.resize(NUM_SECTORS);

        int startSec = (params.single_sector >= 0 && params.single_sector < NUM_SECTORS) ? params.single_sector : 0;
        int endSec = (params.single_sector >= 0 && params.single_sector < NUM_SECTORS) ? params.single_sector + 1 : NUM_SECTORS;

        if (params.log_verbose)
        {
            std::cout << "[FastWildfire] Simulating " << (endSec - startSec) << " directional sector(s)..." << std::endl;
        }

        const double dx = std::max(0.1, std::abs(DEM->cellSizeX()));
        const double dy = std::max(0.1, std::abs(DEM->cellSizeY()));

        // 5. Evaluate Directional Sectors (k in 0..11)
        for (int k = startSec; k < endSec; ++k)
        {
            SectorResults& sec = outputs.sectors[k];
            sec.sector_index = k;
            sec.azimuth_deg = getSectorAzimuthDeg(k);
            sec.downwind_azimuth_deg = getDownwindAzimuthDeg(k);

            sec.Wind_ext = DEM->GetCopySP0();
            sec.Wind_comp = DEM->GetCopySP0();

            sec.FFWI_H = DEM->GetCopySP0();
            sec.FFWI_W = DEM->GetCopySP0();
            sec.FFWI_max = DEM->GetCopySP0();

            sec.FFDI_H = DEM->GetCopySP0();
            sec.FFDI_W = DEM->GetCopySP0();
            sec.FFDI_max = DEM->GetCopySP0();

            sec.H_hyb_H = DEM->GetCopySP0();
            sec.H_hyb_W = DEM->GetCopySP0();
            sec.H_hyb_max = DEM->GetCopySP0();

            sec.Regime_k = DEM->GetCopySP0();
            sec.S_max = DEM->GetCopySP0();
            sec.ConnectedRisk_J = DEM->GetCopySP0();

            const std::shared_ptr<cTMapf>& wExtDir = (k < static_cast<int>(Wind_ext_dirs.size())) ? Wind_ext_dirs[k] : nullptr;
            const std::shared_ptr<cTMapf>& wCompDir = (k < static_cast<int>(Wind_comp_dirs.size())) ? Wind_comp_dirs[k] : nullptr;
            const std::shared_ptr<cTMapf>& fwDir = (k < static_cast<int>(FastWind_dirs.size())) ? FastWind_dirs[k] : nullptr;

            for (int r = 0; r < rows; ++r)
            {
                for (int c = 0; c < cols; ++c)
                {
                    if (pcr::isMV(DEM->data[r][c]))
                    {
                        pcr::setMV(sec.Wind_ext->data[r][c]);
                        pcr::setMV(sec.Wind_comp->data[r][c]);
                        pcr::setMV(sec.FFWI_H->data[r][c]);
                        pcr::setMV(sec.FFWI_W->data[r][c]);
                        pcr::setMV(sec.FFWI_max->data[r][c]);
                        pcr::setMV(sec.FFDI_H->data[r][c]);
                        pcr::setMV(sec.FFDI_W->data[r][c]);
                        pcr::setMV(sec.FFDI_max->data[r][c]);
                        pcr::setMV(sec.H_hyb_H->data[r][c]);
                        pcr::setMV(sec.H_hyb_W->data[r][c]);
                        pcr::setMV(sec.H_hyb_max->data[r][c]);
                        pcr::setMV(sec.Regime_k->data[r][c]);
                        pcr::setMV(sec.S_max->data[r][c]);
                        pcr::setMV(sec.ConnectedRisk_J->data[r][c]);
                        continue;
                    }

                    // 10m wind speeds for sector k [km/h]
                    double vExt = 45.0;
                    double vComp = 18.0;

                    if (fwDir && !pcr::isMV(fwDir->data[r][c]))
                    {
                        // Use terrain-resolved FastWind field
                        vExt = static_cast<double>(fwDir->data[r][c]);
                        vComp = vExt * 0.55;
                    }
                    else
                    {
                        if (wExtDir && !pcr::isMV(wExtDir->data[r][c]))
                        {
                            vExt = static_cast<double>(wExtDir->data[r][c]);
                        }
                        if (wCompDir && !pcr::isMV(wCompDir->data[r][c]))
                        {
                            vComp = static_cast<double>(wCompDir->data[r][c]);
                        }
                        else
                        {
                            vComp = vExt * 0.55;
                        }
                    }

                    // Climate change wind scaling factor gamma_U
                    vExt *= params.gamma_U;
                    vComp *= params.gamma_U;

                    // Annual maximum wind quantiles from all-season ERA5-Land reflect mid-latitude winter storm gales (Nov-Feb).
                    // Condition sustained fire wind to warm-season convective/synoptic extremes in temperate latitudes:
                    if (absLat > 35.0)
                    {
                        double lat_scale = std::clamp((absLat - 35.0) / 25.0, 0.0, 1.0);
                        double season_wind_factor = 1.0 - 0.28 * lat_scale;
                        vExt *= season_wind_factor;
                        vComp *= season_wind_factor;
                    }

                    sec.Wind_ext->data[r][c] = static_cast<float>(vExt);
                    sec.Wind_comp->data[r][c] = static_cast<float>(vComp);

                    bool is_water = (DEM->data[r][c] <= 0.0f);
                    if (LU && !pcr::isMV(LU->data[r][c]))
                    {
                        int luClass = static_cast<int>(std::round(LU->data[r][c]));
                        if (luClass == 80) is_water = true;
                    }
                    if (outputs.Fuel_Flammability->data[r][c] <= 0.005f) is_water = true;

                    if (is_water)
                    {
                        sec.FFWI_H->data[r][c] = 0.0f;
                        sec.FFWI_W->data[r][c] = 0.0f;
                        sec.FFWI_max->data[r][c] = 0.0f;

                        sec.FFDI_H->data[r][c] = 0.0f;
                        sec.FFDI_W->data[r][c] = 0.0f;
                        sec.FFDI_max->data[r][c] = 0.0f;

                        sec.H_hyb_H->data[r][c] = 0.0f;
                        sec.H_hyb_W->data[r][c] = 0.0f;
                        sec.H_hyb_max->data[r][c] = 0.0f;

                        sec.Regime_k->data[r][c] = 0.0f;
                        sec.S_max->data[r][c] = 0.0f;
                        continue;
                    }

                    double tExt = outputs.Temp_ext->data[r][c];
                    double tComp = outputs.Temp_comp->data[r][c];
                    double rhExt = outputs.RH_ext->data[r][c];
                    double rhComp = outputs.RH_comp->data[r][c];

                    // Determine local Drought Factor DF and Fuel Damping eta_SM from real soil moisture
                    double df_local = params.drought_factor;
                    double eta_sm = 1.0;
                    if (outputs.Soil_Moisture && !pcr::isMV(outputs.Soil_Moisture->data[r][c]))
                    {
                        df_local = outputs.Drought_Factor->data[r][c];
                        eta_sm   = outputs.SM_Damping->data[r][c];
                    }

                    // Scenario A — Heat-Dominant Regime (H):
                    // Extreme T + lowest RH + companion wind
                    double ffwi_H = calculateFFWI(tExt, rhExt, vComp);
                    double ffdi_H = calculateFFDI(tExt, rhExt, vComp, df_local);
                    double hyb_H  = calculateHybridHazard(ffwi_H, ffdi_H, params.omega);

                    // Scenario B — Wind-Dominant Regime (W):
                    // Companion T + companion RH + extreme wind
                    double ffwi_W = calculateFFWI(tComp, rhComp, vExt);
                    double ffdi_W = calculateFFDI(tComp, rhComp, vExt, df_local);
                    double hyb_W  = calculateHybridHazard(ffwi_W, ffdi_W, params.omega);

                    // Dominant Regime: compare unclipped raw hazard
                    double raw_H = params.omega * ffwi_H + (1.0 - params.omega) * ffdi_H;
                    double raw_W = params.omega * ffwi_W + (1.0 - params.omega) * ffdi_W;
                    float regime = (raw_W > raw_H) ? 1.0f : 0.0f;

                    float fuelVal = outputs.Fuel_Flammability->data[r][c];
                    double ffwi_max = std::max(ffwi_H, ffwi_W);
                    double ffdi_max = std::max(ffdi_H, ffdi_W);
                    double hyb_max  = std::max(hyb_H, hyb_W);

                    // Standard Pure Meteorological Fire Danger Indices:
                    sec.FFWI_H->data[r][c] = static_cast<float>(ffwi_H);
                    sec.FFWI_W->data[r][c] = static_cast<float>(ffwi_W);
                    sec.FFWI_max->data[r][c] = static_cast<float>(ffwi_max);

                    sec.FFDI_H->data[r][c] = static_cast<float>(ffdi_H);
                    sec.FFDI_W->data[r][c] = static_cast<float>(ffdi_W);
                    sec.FFDI_max->data[r][c] = static_cast<float>(ffdi_max);

                    if (outputs.FFWI_met_env && sec.FFWI_max->data[r][c] > outputs.FFWI_met_env->data[r][c])
                    {
                        outputs.FFWI_met_env->data[r][c] = sec.FFWI_max->data[r][c];
                    }
                    if (outputs.FFDI_met_env && sec.FFDI_max->data[r][c] > outputs.FFDI_met_env->data[r][c])
                    {
                        outputs.FFDI_met_env->data[r][c] = sec.FFDI_max->data[r][c];
                    }

                    // Normalized Hybrid Hazard in [0, 1] modulated by live/dead fuel moisture factor eta_sm:
                    sec.H_hyb_H->data[r][c] = static_cast<float>(hyb_H * eta_sm);
                    sec.H_hyb_W->data[r][c] = static_cast<float>(hyb_W * eta_sm);
                    sec.H_hyb_max->data[r][c] = static_cast<float>(hyb_max * eta_sm);

                    sec.Regime_k->data[r][c] = regime;

                    // Base Forward Spread Potential S_base (without isotropic slope):
                    // Passed into the Fast Sweeping Eikonal solver, which computes directional slope
                    // Az and wind alignment Aw dynamically along propagation vectors without double multiplication.
                    float s_base = static_cast<float>(hyb_max * fuelVal * eta_sm);
                    sec.S_max->data[r][c] = std::clamp(s_base, 0.0001f, 1.0f);
                }
            }

            // Spatial Anisotropic Fast Sweeping Method (FSM) Solver
            sec.Resistance_T = SolveFastSweepingEikonal(
                DEM,
                sec.S_max,
                sec.Wind_ext,
                sec.downwind_azimuth_deg,
                IgnitionMask,
                params);

            // Convert Cumulative Resistance T to Normalized Connected Risk J in [0, 1]
            for (int r = 0; r < rows; ++r)
            {
                for (int c = 0; c < cols; ++c)
                {
                    if (pcr::isMV(DEM->data[r][c]))
                    {
                        pcr::setMV(sec.ConnectedRisk_J->data[r][c]);
                        continue;
                    }

                    bool is_water = (DEM->data[r][c] <= 0.0f);
                    if (LU && !pcr::isMV(LU->data[r][c]))
                    {
                        int luClass = static_cast<int>(std::round(LU->data[r][c]));
                        if (luClass == 80) is_water = true;
                    }
                    if (outputs.Fuel_Flammability->data[r][c] <= 0.005f) is_water = true;

                    if (is_water)
                    {
                        sec.ConnectedRisk_J->data[r][c] = 0.0f;
                        outputs.FFWI_env->data[r][c] = 0.0f;
                        outputs.FFDI_env->data[r][c] = 0.0f;
                        outputs.H_hyb_env->data[r][c] = 0.0f;
                        outputs.S_env->data[r][c] = 0.0f;
                        outputs.ConnectedRisk_env->data[r][c] = 0.0f;
                        outputs.Dominant_Regime_env->data[r][c] = 0.0f;
                        continue;
                    }

                    // Local isotropic steady-state combustion equilibrium hazard S_local
                    // Realized local spread potential incorporates fuel flammability, slope acceleration, and soil moisture damping
                    double hyb_max_cell = static_cast<double>(sec.H_hyb_max->data[r][c]);
                    float fuelVal_cell = outputs.Fuel_Flammability->data[r][c];
                    double eta_sm_cell = outputs.SM_Damping ? static_cast<double>(outputs.SM_Damping->data[r][c]) : 1.0;

                    // Slope factor McArthur / Noble for local reporting:
                    double dz_dx = (DEM->data[r][std::min(cols - 1, c + 1)] - DEM->data[r][std::max(0, c - 1)]) / (2.0 * dx);
                    double dz_dy = (DEM->data[std::min(rows - 1, r + 1)][c] - DEM->data[std::max(0, r - 1)][c]) / (2.0 * dy);
                    double slopeGrad = std::sqrt(dz_dx * dz_dx + dz_dy * dz_dy);
                    double slopeDeg = std::atan(slopeGrad) * (180.0 / M_PI);
                    double slopeMult = std::clamp(std::exp(0.045 * slopeDeg), 1.0, 2.5);
                    if (slopeDeg > 30.0)
                    {
                        double excess = (slopeDeg - 30.0) / 10.0;
                        slopeMult *= (1.0 / (1.0 + excess * excess));
                    }

                    // Local combustion carrying capacity (hyb_max_cell already incorporates soil moisture damping eta_sm):
                    float s_local = std::clamp(static_cast<float>(hyb_max_cell * fuelVal_cell * slopeMult), 0.0f, 1.0f);
                    sec.S_max->data[r][c] = s_local; // Realized sector spread potential

                    // Transmitted exposure from upwind fuels via spatial decay kernel
                    float J_trans = 1.0f; // Default full landscape connectivity
                    float J_val = 0.0f;

                    if (IgnitionMask)
                    {
                        // Specific ignition footprint: T is distance from active ignition point
                        if (sec.Resistance_T && !pcr::isMV(sec.Resistance_T->data[r][c]))
                        {
                            double T_val = static_cast<double>(sec.Resistance_T->data[r][c]);
                            if (T_val >= 0.0)
                            {
                                J_trans = static_cast<float>(std::exp(-T_val / params.tau_0));
                            }
                            else
                            {
                                J_trans = 0.0f;
                            }
                        }
                        else
                        {
                            J_trans = 0.0f; // Unreached by the fire front
                        }
                        J_val = J_trans * s_local;
                    }
                    else
                    {
                        // Autohazard barrier-shielding mode:
                        // T is distance from the nearest upwind physical barrier (water, bare rock, fuelbreak).
                        // In continuous wildland fuel (no upwind barrier), T is MV -> J_trans = 1.0 (unbroken connectivity).
                        // Behind physical barriers, J_shield decays with distance -> J_trans = 1 - J_shield.
                        if (sec.Resistance_T && !pcr::isMV(sec.Resistance_T->data[r][c]))
                        {
                            double T_val = static_cast<double>(sec.Resistance_T->data[r][c]);
                            if (T_val >= 0.0)
                            {
                                float J_shield = static_cast<float>(std::exp(-T_val / params.tau_0));
                                J_trans = std::clamp(1.0f - J_shield, 0.0f, 1.0f);
                            }
                        }
                        // Landscape connected wildfire risk:
                        // Baseline in-situ ignition/spread potential (0.25) plus upwind transmitted continuity (0.75 * J_trans)
                        float transmission_factor = 0.25f + 0.75f * J_trans;
                        J_val = s_local * transmission_factor;
                    }

                    sec.ConnectedRisk_J->data[r][c] = std::clamp(J_val, 0.0f, 1.0f);

                    // Realized Landscape Danger Modulation for FFDI and FFWI:
                    // Modulates the atmospheric fire danger scale across space according to:
                    // 1. Land use fuel flammability F(x) (streets ~0.02, crops ~0.20, grass ~0.45, forest ~0.80+)
                    // 2. Landscape connectivity / transmitted continuity (shielding behind barriers/concrete)
                    // 3. Live & dead fuel moisture damping eta_SM from soil moisture
                    // 4. Topographic slope acceleration slopeMult
                    float trans_factor = (IgnitionMask) ? J_trans : (0.25f + 0.75f * J_trans);
                    float landscape_factor = 0.0f;
                    if (!is_water && fuelVal_cell > 0.005f)
                    {
                        landscape_factor = std::clamp(
                            (fuelVal_cell / 0.80f) * trans_factor * (static_cast<float>(eta_sm_cell) / 0.75f) * static_cast<float>(slopeMult),
                            0.0f,
                            1.5f
                        );
                    }

                    sec.FFWI_H->data[r][c] *= landscape_factor;
                    sec.FFWI_W->data[r][c] *= landscape_factor;
                    sec.FFWI_max->data[r][c] *= landscape_factor;

                    sec.FFDI_H->data[r][c] *= landscape_factor;
                    sec.FFDI_W->data[r][c] *= landscape_factor;
                    sec.FFDI_max->data[r][c] *= landscape_factor;

                    // Update Omnidirectional Envelopes
                    if (sec.FFWI_max->data[r][c] > outputs.FFWI_env->data[r][c])
                    {
                        outputs.FFWI_env->data[r][c] = sec.FFWI_max->data[r][c];
                    }
                    if (sec.FFDI_max->data[r][c] > outputs.FFDI_env->data[r][c])
                    {
                        outputs.FFDI_env->data[r][c] = sec.FFDI_max->data[r][c];
                    }
                    if (sec.H_hyb_max->data[r][c] > outputs.H_hyb_env->data[r][c])
                    {
                        outputs.H_hyb_env->data[r][c] = sec.H_hyb_max->data[r][c];
                        outputs.Dominant_Regime_env->data[r][c] = sec.Regime_k->data[r][c];
                    }
                    if (s_local > outputs.S_env->data[r][c])
                    {
                        outputs.S_env->data[r][c] = s_local;
                    }
                    if (sec.ConnectedRisk_J->data[r][c] > outputs.ConnectedRisk_env->data[r][c])
                    {
                        outputs.ConnectedRisk_env->data[r][c] = sec.ConnectedRisk_J->data[r][c];
                    }

                }
            }
        }

        // Post-processing sweep: guarantee water cells are cleanly set to 0.0 across all outputs
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (!pcr::isMV(DEM->data[r][c]))
                {
                    bool is_water = (DEM->data[r][c] <= 0.0f);
                    if (LU && !pcr::isMV(LU->data[r][c]) && static_cast<int>(std::round(LU->data[r][c])) == 80) is_water = true;
                    if (outputs.Fuel_Flammability->data[r][c] <= 0.005f) is_water = true;

                    if (is_water)
                    {
                        outputs.Fuel_Flammability->data[r][c] = 0.0f;
                        outputs.FFWI_env->data[r][c] = 0.0f;
                        outputs.FFDI_env->data[r][c] = 0.0f;
                        if (outputs.FFWI_met_env) outputs.FFWI_met_env->data[r][c] = 0.0f;
                        if (outputs.FFDI_met_env) outputs.FFDI_met_env->data[r][c] = 0.0f;
                        outputs.H_hyb_env->data[r][c] = 0.0f;
                        outputs.S_env->data[r][c] = 0.0f;
                        outputs.ConnectedRisk_env->data[r][c] = 0.0f;
                        outputs.Dominant_Regime_env->data[r][c] = 0.0f;

                        if (outputs.Soil_Moisture) outputs.Soil_Moisture->data[r][c] = 0.0f;
                        if (outputs.Drought_Factor) outputs.Drought_Factor->data[r][c] = 0.0f;
                        if (outputs.SM_Damping) outputs.SM_Damping->data[r][c] = 0.0f;
                    }
                }
            }
        }

        if (params.log_verbose)
        {
            std::cout << "[FastWildfire] Simulation completed successfully." << std::endl;
        }

        return 0;
    }

} // namespace FastWildfireModel

#endif // FASTWILDFIRE_H
