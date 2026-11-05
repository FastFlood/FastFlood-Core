//
// Created for FastDrought Model
// Standardized Precipitation Evapotranspiration Index (SPEI), Theory of Runs,
// and Physical Climate & Adaptation Impact Engine.
// Exactly matching FastDrought droughtutils.ts.
//

#ifndef FASTDROUGHT_H
#define FASTDROUGHT_H

#include <cmath>
#include <vector>
#include <memory>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <limits>

#include "../../common/defines.h"
#include "../../common/fferror.h"
#include "../../geospatial/map/map.h"

namespace FastDrought
{

// -----------------------------------------------------------------------------
// Constants & Return Periods (Matching droughtutils.ts)
// -----------------------------------------------------------------------------

inline const std::vector<int> DROUGHT_RETURN_PERIODS = {2, 5, 10, 20, 50, 100, 200};
inline const std::vector<std::string> DROUGHT_CMIP6_SCENARIOS = {"ssp126", "ssp245", "ssp370", "ssp585"};
inline const std::vector<int> DROUGHT_PERIODS = {2030, 2050, 2070, 2090};

// -----------------------------------------------------------------------------
// Data Structures for Drought Hazard, Events, Runs, and Adaptation Assessment
// -----------------------------------------------------------------------------

struct DroughtHazardCategory
{
    std::string label;
    std::string level; // "extreme_wet", "very_wet", "moderate_wet", "normal", "moderate_drought", "severe_drought", "extreme_drought"
    std::string color;
    std::string bg;
    std::string border;
    std::string badge_class;
    std::string description;
};

struct DroughtReturnPeriodPoint
{
    int rp = 50;                  // Return period in years (e.g. 2, 5, 10, 20, 50, 100, 200)
    float prob_percent = 2.0f;    // Annual probability (%)
    float spei = -2.05f;          // Negative z-score representing the drought threshold
    float deficit_mm = -180.0f;   // Corresponding water balance deficit D = P - PET [mm]
    DroughtHazardCategory hazard;
};

struct MonthlyWaterBalance
{
    int month_index = 1;         // 1..12
    std::string month_label = "M1"; // "M1".."M12"
    float precipitation = 0.0f;  // mm/month
    float pet = 0.0f;            // Potential Evapotranspiration [mm/month]
    float deficit = 0.0f;        // Meteorological water balance deficit D = P - PET [mm/month]
    float spei = 0.0f;           // Meteorological SPEI z-score [-]
    bool is_drought_month = false; // SPEI <= -1.0
    bool is_severe_drought = false;// SPEI <= -1.5
    bool is_extreme_drought = false;// SPEI <= -2.0
    // Root-Zone Soil Moisture Buffering (FAO-56):
    float soil_moisture = 0.25f; // Volumetric fraction SM [m3/m3] (0.02 - 0.50)
    float soil_buffer_mm = 0.0f; // Available plant root-zone soil water buffer [mm/month]
    float agri_deficit = 0.0f;   // Realized agricultural deficit D_agri = D + min(-D, W_buffer) [mm]
    float agri_spei = 0.0f;      // Agricultural Drought SPEI z-score accounting for soil moisture buffering
    bool is_agri_drought_month = false; // agri_spei <= -1.0
    bool is_valid = false;       // True if terrestrial cell with valid data
};
using MonthlyWaterBalancePoint = MonthlyWaterBalance;

struct CellRunHazardMetrics
{
    // Meteorological Drought (Atmospheric Demand vs Supply)
    int drought_months = 0;       // Count of months with SPEI <= -1.0
    float total_severity = 0.0f;  // Sum of |SPEI + 1.0| for SPEI <= -1.0
    float peak_intensity = 0.0f;  // Most negative (worst) SPEI value observed
    float mean_spei = 0.0f;       // Average SPEI across active months
    float total_deficit_mm = 0.0f;// Cumulative water deficit: sum(max(0, -D)) [mm]
    int event_count_drought = 0;  // Number of discrete drought episodes

    // Agricultural / Soil Moisture Drought (Root-Zone Water Availability)
    int agri_drought_months = 0;      // Count of months with Agri-SPEI <= -1.0
    float agri_total_severity = 0.0f; // Sum of |Agri-SPEI + 1.0| for Agri-SPEI <= -1.0
    float agri_peak_intensity = 0.0f; // Most negative (worst) Agri-SPEI value observed
    float agri_mean_spei = 0.0f;      // Average Agri-SPEI across active months
    float agri_total_deficit_mm = 0.0f;// Cumulative agricultural water deficit [mm]
    int agri_event_count = 0;         // Number of discrete agricultural drought episodes
    float mean_soil_moisture = 0.0f;  // Mean volumetric soil moisture during event
};

struct DroughtDistributionParams
{
    float alpha = 60.0f;   // Scale parameter (alpha > 0) [mm]
    float beta = 2.0f;     // Shape parameter (beta > 0) [-]
    float gamma = -100.0f; // Origin / minimum deficit bound (gamma) [mm]
    int month = 7;         // 1..12 (e.g. 7 for July baseline)
    double lat = 0.0;
    double lng = 0.0;
    bool has_data = true;
};

struct DroughtEvent
{
    int rank = 1;                // 1 to 20
    int duration = 0;            // Number of months with meteorological SPEI <= -1.0
    float severity = 0.0f;       // Sum of |SPEI + 1.0| across drought months
    float peak_intensity = 0.0f; // Most negative (worst) meteorological SPEI in the event
    float total_deficit_mm = 0.0f;// Sum of max(0, -D) across the 12 months
    float mean_precip_mm = 0.0f;
    float mean_pet_mm = 0.0f;
    std::vector<MonthlyWaterBalancePoint> months; // Fixed 12 months
    // Agricultural Drought metrics (Yevjevich Theory of Runs on SPEI_agri):
    int agri_duration = 0;       // Consecutive months where SPEI_agri <= -1.0
    float agri_severity = 0.0f;  // Cumulative agricultural deficit severity below -1.0
    float agri_peak_intensity = 0.0f; // Peak (worst) SPEI_agri
    float total_agri_deficit_mm = 0.0f; // Cumulative realized deficit after soil buffering
    float soil_buffer_mitigated_mm = 0.0f; // Total mm of water deficit buffered by root-zone soil moisture
};

struct HistoricalEventFeature
{
    int rank = 1;
    int peak_month = 6;
    float peak_spei = 0.0f;
    int duration = 0;
    bool has_data = false;
};

struct DroughtScenarioModifiers
{
    float rain_multiplier = 1.0f;          // Multiplicative precipitation scaling (e.g. 0.90 = -10%)
    float delta_temp_c = 0.0f;             // Additive warming offset in deg C (e.g. +2.0 C)
    float wind_multiplier = 1.0f;          // Wind speed scaling (e.g. 1.05 = +5%)
    float pet_multiplier = 1.0f;           // Derived thermodynamic PET scaling multiplier
    float irrigation_supply_mm = 0.0f;     // Supplemental agricultural irrigation [mm/month]
    float et_reduction_fraction = 0.0f;    // Evaporation suppression (mulch, shading, agrivoltaics) [0.0 - 0.50]
    float soil_retention_buffer_mm = 0.0f; // Water storage and retention buffering [mm/month]
    float root_zone_depth_m = 0.60f;       // Effective root-zone depth (default 0.60 m)
};

struct EvaluatedEventScenario
{
    int rank = 1;
    float baseline_severity = 0.0f;
    float scenario_severity = 0.0f;
    float severity_change_percent = 0.0f;
    float baseline_peak_spei = 0.0f;
    float scenario_peak_spei = 0.0f;
    float delta_peak_spei = 0.0f;
    int baseline_duration = 0;
    int scenario_duration = 0;
    float baseline_deficit_mm = 0.0f;
    float scenario_deficit_mm = 0.0f;
    float delta_deficit_mm = 0.0f;
    float net_benefit_mm = 0.0f; // Deficit avoided through adaptation
    // Agricultural drought metrics:
    float baseline_agri_severity = 0.0f;
    float scenario_agri_severity = 0.0f;
    float baseline_agri_peak_spei = 0.0f;
    float scenario_agri_peak_spei = 0.0f;
    int baseline_agri_duration = 0;
    int scenario_agri_duration = 0;
    int agri_duration = 0; // Alias for baseline_agri_duration
    float baseline_agri_deficit_mm = 0.0f;
    float scenario_agri_deficit_mm = 0.0f;
    float soil_buffer_mitigated_mm = 0.0f;
    float scenario_soil_buffer_mitigated_mm = 0.0f;

    struct EvaluatedMonth {
        int month_index = 1;
        std::string month_label;
        float baseline_precip = 0.0f;
        float scenario_precip = 0.0f;
        float baseline_pet = 0.0f;
        float scenario_pet = 0.0f;
        float baseline_deficit = 0.0f;
        float scenario_deficit = 0.0f;
        float baseline_spei = 0.0f;
        float scenario_spei = 0.0f;
        float baseline_soil_moisture = 0.0f;
        float scenario_soil_moisture = 0.0f;
        float baseline_soil_buffer_mm = 0.0f;
        float scenario_soil_buffer_mm = 0.0f;
        float baseline_agri_deficit = 0.0f;
        float scenario_agri_deficit = 0.0f;
        float baseline_agri_spei = 0.0f;
        float scenario_agri_spei = 0.0f;
        bool is_drought_base = false;
        bool is_drought_scen = false;
        bool is_agri_drought_base = false;
        bool is_agri_drought_scen = false;
    };
    std::vector<EvaluatedMonth> months;
};

struct AgriculturalWaterBalanceResult
{
    float met_deficit = 0.0f;
    float met_spei = 0.0f;
    float agri_deficit = 0.0f;
    float agri_spei = 0.0f;
    float soil_buffer_available = 0.0f;
    float soil_buffer_used = 0.0f;
};

struct MonthlyClimatologyPoint
{
    int month = 1;              // 1..12
    std::string name;          // "Jan", "Feb", ...
    float precipitation = 0.0f; // mm/month
    float pet = 0.0f;           // mm/month
    float deficit = 0.0f;       // D = P - PET (positive = surplus, negative = deficit)
    bool is_deficit = false;    // True if D < 0
};

struct SiteDroughtRiskScore
{
    int score = 40;             // 0 to 100
    std::string level = "Moderate"; // Low, Moderate, High, Severe, Extreme
    std::string color = "#d97706";
    std::string badge_class = "bg-amber-100 text-amber-800 border-amber-300";
    std::string headline;
    std::string summary;
    int annual_deficit_months = 0;
    float typical_summer_deficit_mm = 0.0f;
    int drought_frequency_years = 10;
    int worst_duration_months = 7;
    int worst_agri_duration_months = 2;
    float worst_peak_spei = -2.3f;
    float worst_agri_peak_spei = -1.5f;
    std::string recommended_mitigation;
};

struct CMIP6TrendPreset
{
    std::string scenario;
    int period = 2050;
    std::string label;
    float delta_precip_percent = 0.0f;
    float delta_temp_c = 0.0f;
    float delta_pet_percent = 0.0f;
};

inline const std::vector<CMIP6TrendPreset> CMIP6_DROUGHT_TREND_PRESETS = {
    {"ssp126", 2050, "SSP1-2.6 (Paris Ambitious, 2050)", -2.5f, 1.4f, 7.0f},
    {"ssp245", 2050, "SSP2-4.5 (Middle of the Road, 2050)", -6.0f, 2.1f, 10.5f},
    {"ssp370", 2070, "SSP3-7.0 (Regional Rivalry, 2070)", -12.0f, 3.4f, 17.0f},
    {"ssp585", 2090, "SSP5-8.5 (Fossil-Fueled High Emission, 2090)", -18.5f, 4.8f, 24.0f}
};

// -----------------------------------------------------------------------------
// Mathematical Core (Exact match with FastDrought droughtutils.ts)
// -----------------------------------------------------------------------------

/**
 * Standard Normal Quantile Function (Inverse Normal CDF): Phi^{-1}(p)
 * Uses high-accuracy Abramowitz & Stegun rational approximation (formula 26.2.23).
 * Clamped strictly to [-3.00, +3.00] (p in [0.00135, 0.99865]) to prevent extreme artifact drops.
 */
inline float NormalQuantile(float p)
{
    float p_clamped = std::clamp(p, 0.00135f, 0.99865f);
    if (p_clamped <= 0.00135f) return -3.00f;
    if (p_clamped >= 0.99865f) return 3.00f;

    float sign = 1.0f;
    float p_calc = p_clamped;
    if (p_clamped <= 0.5f)
    {
        sign = -1.0f;
        p_calc = p_clamped;
    }
    else
    {
        sign = 1.0f;
        p_calc = 1.0f - p_clamped;
    }

    float w = std::sqrt(-2.0f * std::log(p_calc));
    const float c0 = 2.515517f;
    const float c1 = 0.802853f;
    const float c2 = 0.010328f;
    const float d1 = 1.432788f;
    const float d2 = 0.189269f;
    const float d3 = 0.001308f;

    float numerator = c0 + c1 * w + c2 * w * w;
    float denominator = 1.0f + d1 * w + d2 * w * w + d3 * w * w * w;
    float z = w - (numerator / denominator);
    return std::clamp(sign * z, -3.0f, 3.0f);
}

/**
 * Standard Normal Cumulative Distribution Function: Phi(z)
 */
inline float NormalCDF(float z)
{
    float clamped_z = std::clamp(z, -3.0f, 3.0f);
    float t = 1.0f / (1.0f + 0.2316419f * std::abs(clamped_z));
    float d = 0.39894228f * std::exp((-clamped_z * clamped_z) / 2.0f);
    float prob = d * t * (0.31938153f + t * (-0.356563782f + t * (1.781477937f + t * (-1.821255978f + t * 1.330274429f))));
    prob = (clamped_z > 0.0f) ? 1.0f - prob : prob;
    return std::clamp(prob, 0.00135f, 0.99865f);
}

/**
 * 3-Parameter Log-Logistic Cumulative Distribution Function:
 * F(D) = [ 1 + (alpha / (D - gamma))^beta ]^(-1)
 * Includes smooth tail continuation below physical historical floor gamma (matching FastDrought kernel).
 * Prevents artificial step cliff to -3.72 when climate warming pushes deficit past gamma.
 */
inline float DeficitToLogLogisticCDF(float deficit, float alpha, float beta, float gamma)
{
    float a = std::max(0.001f, alpha);
    float b = std::max(0.05f, beta);

    float diff = deficit - gamma;
    if (diff <= 0.001f)
    {
        float under = std::max(0.0f, gamma - deficit);
        float tail_p = 0.002f / (1.0f + (under / a) * b * 2.0f);
        return std::clamp(tail_p, 0.00135f, 0.002f);
    }

    float ratio = a / diff;
    float power_term = std::pow(ratio, b);
    float cdf = 1.0f / (1.0f + power_term);
    return std::clamp(cdf, 0.00135f, 0.99865f);
}

/**
 * Transforms climatic water balance deficit D [mm] into SPEI z-score:
 * SPEI = Phi^{-1}( F(D) )
 */
inline float DeficitToSPEI(float deficit, float alpha, float beta, float gamma)
{
    float cdf = DeficitToLogLogisticCDF(deficit, alpha, beta, gamma);
    float z = NormalQuantile(cdf);
    return std::clamp(z, -3.0f, 3.0f);
}

/**
 * Inverts SPEI z-score back to underlying water balance deficit D [mm]:
 * D = gamma + alpha * ( p / (1 - p) )^(1 / beta)
 */
inline float SPEIToDeficit(float spei, float alpha, float beta, float gamma)
{
    float a = std::max(0.001f, alpha);
    float b = std::max(0.05f, beta);

    float p = std::clamp(NormalCDF(spei), 0.00135f, 0.99865f);
    float odds = p / (1.0f - p);
    float diff = a * std::pow(odds, 1.0f / b);
    return gamma + diff;
}

/**
 * Physically based thermodynamic PET scaling from climate warming and wind anomalies.
 * Based on Penman-Monteith / Clausius-Clapeyron vapour pressure deficit sensitivity:
 * ~5.0% PET per deg C warming.
 */
inline float ComputePETClimateMultiplier(float delta_temp_c, float wind_multiplier = 1.0f)
{
    float temp_effect = 1.0f + 0.050f * std::max(0.0f, delta_temp_c);
    float wind_effect = 1.0f + 0.250f * std::max(-0.5f, wind_multiplier - 1.0f);
    return std::clamp(temp_effect * wind_effect, 0.5f, 2.5f);
}

/**
 * Computes available transpirable root-zone soil water buffer capacity according to FAO-56:
 * W_buffer = clamp(SM - SM_pwp, 0.0, 0.35) * Z_r * 1000 [mm]
 *
 * @param volumetric_sm Volumetric soil moisture fraction [m3/m3]
 * @param root_depth_m Effective plant rooting depth [m] (default 0.60 m)
 * @param sm_pwp Permanent wilting point moisture threshold [m3/m3] (default 0.10)
 * @param sm_fc Field capacity moisture threshold (unused in buffer formula)
 * @returns Available transpirable soil moisture buffer capped at 180.0 mm/month
 */
inline float ComputeAvailableSoilWaterBuffer(
    float volumetric_sm,
    float root_depth_m = 0.60f,
    float sm_pwp = 0.10f,
    float sm_fc = 0.36f)
{
    (void)sm_fc;
    float sm = std::clamp(volumetric_sm, 0.01f, 0.60f);
    float excess_sm = std::clamp(sm - sm_pwp, 0.0f, 0.35f);
    float buffer_mm = excess_sm * root_depth_m * 1000.0f;
    return std::min(180.0f, buffer_mm);
}

/**
 * Computes agricultural climatic water balance and transformed SPEI_agri.
 * If atmospheric deficit D = P - PET < 0, available root-zone soil water buffer
 * mitigates vegetation stress: D_agri = D + min(-D, W_buffer).
 */
inline AgriculturalWaterBalanceResult CalculateAgriculturalWaterBalance(
    float precip,
    float pet,
    float sm,
    float alpha,
    float beta,
    float gamma,
    float root_depth_m = 0.60f)
{
    float d_met = precip - pet;
    float soil_buffer = ComputeAvailableSoilWaterBuffer(sm, root_depth_m);

    float d_agri = d_met;
    float buffer_used = 0.0f;
    if (d_agri < 0.0f)
    {
        buffer_used = std::min(-d_agri, soil_buffer);
        d_agri += buffer_used;
    }

    AgriculturalWaterBalanceResult res;
    res.met_deficit = d_met;
    res.met_spei = DeficitToSPEI(d_met, alpha, beta, gamma);
    res.agri_deficit = d_agri;
    res.agri_spei = DeficitToSPEI(d_agri, alpha, beta, gamma);
    res.soil_buffer_available = soil_buffer;
    res.soil_buffer_used = buffer_used;
    return res;
}

/**
 * FastWildfire Dynamic Drought Factor from Soil Moisture (Griffiths 1998 / Finkele 2006).
 * Scales McArthur Drought Factor DF in [1.0, 10.0] based on volumetric soil moisture.
 */
inline float CalculateDroughtFactorFromSM(
    float sm,
    float sm_min = 0.15f,
    float sm_var = 0.12f)
{
    const float sm_pwp = 0.10f;
    const float sm_fc = 0.36f;
    float clamped_sm = std::clamp(sm, 0.02f, 0.60f);

    float theta_rel = std::clamp((clamped_sm - sm_pwp) / (sm_fc - sm_pwp), 0.0f, 1.0f);
    float df_base = std::clamp(1.0f + 9.0f * std::pow(1.0f - theta_rel, 1.6f), 1.0f, 10.0f);
    float dry_season_intensity = std::clamp((0.36f - sm_min) / 0.26f, 0.0f, 1.0f);
    float seasonal_amp = std::clamp(sm_var / 0.12f, 0.25f, 1.0f);
    float df_max = 1.0f + 9.0f * (dry_season_intensity * seasonal_amp);
    return std::min(df_base, df_max);
}

/**
 * FastWildfire Fuel Moisture Sigmoid Damping Factor (Rothermel 1972 / Krueger 2015).
 * Steep sigmoid suppression curve around 0.22 m3/m3; dampens spread rate in [0.02, 1.0].
 */
inline float CalculateSoilMoistureDamping(float sm)
{
    float clamped_sm = std::clamp(sm, 0.01f, 0.60f);
    float damp = 1.0f / (1.0f + std::exp(22.0f * (clamped_sm - 0.22f)));
    return std::clamp(damp, 0.02f, 1.0f);
}

/**
 * Categorizes an SPEI value according to standard World Meteorological Organization (WMO)
 * and Vicente-Serrano drought classification.
 */
inline DroughtHazardCategory GetDroughtHazardCategory(float spei)
{
    if (spei <= -2.0f)
    {
        return {"Extreme Drought", "extreme_drought", "#991b1b", "#fef2f2", "#f87171",
                "bg-red-100 text-red-900 border-red-300",
                "Exceptional moisture deficit; severe crop/pasture failure, widespread water shortages and fire risk."};
    }
    if (spei <= -1.5f)
    {
        return {"Severe Drought", "severe_drought", "#c2410c", "#fff7ed", "#fb923c",
                "bg-orange-100 text-orange-900 border-orange-300",
                "Significant hydrological deficit; crop/pasture losses likely, water restrictions typically imposed."};
    }
    if (spei <= -1.0f)
    {
        return {"Moderate Drought", "moderate_drought", "#d97706", "#fffbeb", "#fcd34d",
                "bg-amber-100 text-amber-900 border-amber-300",
                "Noticeable water deficit; some damage to crops and pastures, declining streamflows and reservoir levels."};
    }
    if (spei < 1.0f)
    {
        return {"Near Normal", "normal", "#15803d", "#f0fdf4", "#86efac",
                "bg-green-100 text-green-800 border-green-300",
                "Balanced climatic water budget within standard historical variability (SPEI between -1.0 and +1.0)."};
    }
    if (spei < 1.5f)
    {
        return {"Moderately Wet", "moderate_wet", "#0284c7", "#f0f9ff", "#7dd3fc",
                "bg-sky-100 text-sky-800 border-sky-300",
                "Moisture surplus above average precipitation."};
    }
    if (spei < 2.0f)
    {
        return {"Very Wet", "very_wet", "#1d4ed8", "#eff6ff", "#93c5fd",
                "bg-blue-100 text-blue-800 border-blue-300",
                "High precipitation excess; elevated groundwater and saturated root-zones."};
    }
    return {"Extremely Wet", "extreme_wet", "#4338ca", "#eef2ff", "#a5b4fc",
            "bg-indigo-100 text-indigo-800 border-indigo-300",
            "Exceptional moisture surplus; potential waterlogging and pluvial saturation."};
}

/**
 * Builds Return Period extreme hazard quantiles for drought:
 * For drought return period Tr, non-exceedance probability is p = 1 / Tr.
 */
inline std::vector<DroughtReturnPeriodPoint> CalculateDroughtRPProfile(float alpha, float beta, float gamma)
{
    std::vector<DroughtReturnPeriodPoint> profile;
    for (int rp : DROUGHT_RETURN_PERIODS)
    {
        float prob = 1.0f / static_cast<float>(rp);
        float spei = NormalQuantile(prob);
        float deficit_mm = SPEIToDeficit(spei, alpha, beta, gamma);
        auto hazard = GetDroughtHazardCategory(spei);

        DroughtReturnPeriodPoint pt;
        pt.rp = rp;
        pt.prob_percent = prob * 100.0f;
        pt.spei = spei;
        pt.deficit_mm = deficit_mm;
        pt.hazard = hazard;
        profile.push_back(pt);
    }
    return profile;
}

/**
 * Theoretical probability density function (PDF) for the 3-Parameter Log-Logistic Distribution:
 * f(D) = (beta / alpha) * ((D - gamma) / alpha)^(beta - 1) / [ 1 + ((D - gamma) / alpha)^beta ]^2
 */
inline float LogLogisticPDF(float deficit, float alpha, float beta, float gamma)
{
    float diff = deficit - gamma;
    if (diff <= 0.0f) return 0.0f;
    float a = std::max(0.001f, alpha);
    float b = std::max(0.05f, beta);

    float z = diff / a;
    float num = (b / a) * std::pow(z, b - 1.0f);
    float den = std::pow(1.0f + std::pow(z, b), 2.0f);
    return num / den;
}

// -----------------------------------------------------------------------------
// Climatology, Soil Moisture & Synthetic Distribution Generation
// -----------------------------------------------------------------------------

/**
 * Returns baseline monthly precipitation and potential evapotranspiration (PET)
 * climatology for any given latitude/longitude to construct annual water balance charts.
 * Identifies major global macro-climatic regimes (Wet Tropics/Amazon, Deserts, Mediterranean, Temperate).
 */
inline std::vector<MonthlyClimatologyPoint> GetAnnualClimatology(double lat, double lng)
{
    const std::vector<std::string> month_names = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };
    double abs_lat = std::abs(lat);
    bool is_north = (lat >= 0.0);

    // Macro-Climatic Regimes
    bool is_amazon = (lat >= -15.0 && lat <= 7.0 && lng >= -82.0 && lng <= -45.0);
    bool is_congo = (lat >= -6.0 && lat <= 6.0 && lng >= 10.0 && lng <= 32.0);
    bool is_maritime_se_asia = (lat >= -11.0 && lat <= 15.0 && lng >= 95.0 && lng <= 145.0);
    bool is_wet_tropics = is_amazon || is_congo || is_maritime_se_asia;

    bool is_sahara_arabia = (lat >= 14.0 && lat <= 34.0 && lng >= -16.0 && lng <= 60.0);
    bool is_australia_arid = (lat >= -35.0 && lat <= -18.0 && lng >= 115.0 && lng <= 142.0);
    bool is_atacama = (lat >= -28.0 && lat <= -18.0 && lng >= -72.0 && lng <= -68.0);
    bool is_us_southwest = (lat >= 25.0 && lat <= 38.0 && lng >= -118.0 && lng <= -102.0);
    bool is_hyper_arid = is_sahara_arabia || is_australia_arid || is_atacama || is_us_southwest;

    bool is_med_basin = (lat >= 31.0 && lat <= 44.0 && lng >= -10.0 && lng <= 36.0);
    bool is_california = (lat >= 32.0 && lat <= 42.0 && lng >= -124.0 && lng <= -116.0);
    bool is_mediterranean = is_med_basin || is_california;

    std::vector<MonthlyClimatologyPoint> result(12);
    for (int idx = 0; idx < 12; ++idx)
    {
        int m = idx + 1;
        int peak_month = is_north ? 7 : 1;
        int month_dist = std::min(std::abs(m - peak_month), 12 - std::abs(m - peak_month));
        double insolation_fraction = std::cos((month_dist / 6.0) * 3.141592653589793) * 0.5 + 0.5;

        double pet = 0.0;
        double precip = 0.0;

        if (is_wet_tropics)
        {
            // Equatorial Rainforest (Amazon, Congo, SE Asia): High rainfall year round, moderate stable PET
            pet = 90.0 + insolation_fraction * 25.0; // 90 to 115 mm/mo
            precip = 210.0 + std::sin((m / 12.0) * 3.141592653589793 * 2.0) * 55.0 + (is_amazon ? 35.0 : 0.0);
        }
        else if (is_hyper_arid)
        {
            // Subtropical Deserts: Minimal rainfall, massive evaporative demand
            pet = 110.0 + insolation_fraction * 115.0; // 110 to 225 mm/mo
            precip = std::max(2.0, 10.0 - insolation_fraction * 7.0); // 3 to 10 mm/mo
        }
        else if (is_mediterranean)
        {
            // Dry-summer subtropical: Wet mild winters, hot bone-dry summers
            pet = 35.0 + insolation_fraction * 125.0; // 35 to 160 mm/mo
            precip = std::max(6.0, 85.0 - insolation_fraction * 75.0); // 85 mm in winter, 10 mm in summer
        }
        else if (abs_lat < 20.0)
        {
            // Savanna / Monsoon tropics
            pet = 70.0 + insolation_fraction * 60.0;
            precip = 65.0 + insolation_fraction * 140.0;
        }
        else if (abs_lat >= 20.0 && abs_lat <= 45.0)
        {
            // Subtropical / Temperate
            pet = 30.0 + insolation_fraction * 95.0; // 30 to 125 mm/mo
            precip = 65.0 + std::sin(m * 0.5) * 15.0; // 50 to 80 mm/mo
        }
        else
        {
            // High latitude / Boreal
            pet = 15.0 + insolation_fraction * 65.0; // 15 to 80 mm/mo
            precip = 55.0 + std::sin(m * 0.5) * 15.0; // 40 to 70 mm/mo
        }

        if (std::abs(lng) > 80.0 && abs_lat > 30.0 && !is_wet_tropics)
        {
            precip = std::max(12.0, precip * 0.90);
        }

        float rounded_p = static_cast<float>(std::round(precip));
        float rounded_pet = static_cast<float>(std::round(pet));
        float def = rounded_p - rounded_pet;

        result[idx].month = m;
        result[idx].name = month_names[idx];
        result[idx].precipitation = rounded_p;
        result[idx].pet = rounded_pet;
        result[idx].deficit = def;
        result[idx].is_deficit = (def < 0.0f);
    }
    return result;
}

/**
 * Generates synthetic monthly root-zone soil moisture series based on latitude and climate zone
 * matching generateSyntheticSoilMoisture() in soilmoistureutils.ts.
 */
inline std::vector<float> GenerateSyntheticSoilMoisture(double lat, double lng)
{
    (void)lng;
    double abs_lat = std::abs(lat);
    bool is_southern = (lat < 0.0);

    float base_profile = 0.27f;
    float amplitude = 0.10f;

    if (abs_lat < 12.0) {
        base_profile = 0.34f;
        amplitude = 0.05f;
    } else if (abs_lat >= 15.0 && abs_lat <= 32.0) {
        base_profile = 0.12f;
        amplitude = 0.06f;
    } else if (abs_lat > 32.0 && abs_lat <= 58.0) {
        base_profile = 0.27f;
        amplitude = 0.10f;
    } else {
        base_profile = 0.32f;
        amplitude = 0.09f;
    }

    std::vector<float> root_zone(12);
    for (int m = 0; m < 12; ++m)
    {
        int phase_shift = is_southern ? 7 : 1;
        double angle = ((m - phase_shift) / 12.0) * 2.0 * 3.141592653589793;
        double seasonal_factor = std::cos(angle);
        float rz = static_cast<float>(base_profile + 0.01f + amplitude * 0.75f * seasonal_factor);
        root_zone[m] = std::clamp(rz, 0.05f, 0.50f);
    }
    return root_zone;
}

/**
 * Generates robust, location-tailored synthetic distributions based on latitude/longitude
 * when offline or when remote COG is unreachable. Matches droughtutils.ts:generateSyntheticDroughtParams().
 */
inline DroughtDistributionParams GenerateSyntheticDroughtParams(double lat, double lng)
{
    double abs_lat = std::abs(lat);
    bool is_med = (lat >= 31.0 && lat <= 44.0 && lng >= -10.0 && lng <= 36.0) ||
                  (lat >= 32.0 && lat <= 42.0 && lng >= -124.0 && lng <= -116.0);
    bool is_arid = is_med || (abs_lat > 15.0 && abs_lat < 35.0 && std::abs(lng) < 60.0) ||
                             (abs_lat > 25.0 && abs_lat < 38.0 && lng < -100.0);

    float alpha = is_arid ? 55.0f : 75.0f; // Scale
    float beta = 2.15f;                    // Shape
    float gamma = is_med ? -220.0f : (is_arid ? -180.0f : -140.0f);

    alpha += static_cast<float>(std::sin(lat * 0.1) * 6.0);
    gamma += static_cast<float>(std::cos(lng * 0.1) * 10.0);

    DroughtDistributionParams params;
    params.alpha = alpha;
    params.beta = beta;
    params.gamma = gamma;
    params.month = 7;
    params.lat = lat;
    params.lng = lng;
    params.has_data = true;
    return params;
}

// -----------------------------------------------------------------------------
// Top Events Hydrographs & Real/Synthetic Data Blending
// -----------------------------------------------------------------------------

/**
 * Maps rank 1..20 to expected peak SPEI according to extreme value order statistics (GEV/Gumbel).
 * Rank 1 (1-in-50 yr): ~ -2.38 down to Rank 20 (1-in-2.2 yr): ~ -1.08.
 */
inline float GetTargetProbabilityPeakSPEI(int rank, double seed = 42.0)
{
    double p = (rank - 1) / 19.0;
    double base_peak = -2.38 + std::pow(p, 0.70) * 1.30;
    double jitter = std::sin(seed + rank * 9.31) * 0.035;
    return std::clamp(static_cast<float>(base_peak + jitter), -2.55f, -1.05f);
}

/**
 * Checks whether two drought events are duplicate or overly identical.
 */
inline bool AreEventsTooSimilar(const DroughtEvent& a, const DroughtEvent& b)
{
    float max_spei_diff = 0.0f;
    float max_agri_diff = 0.0f;
    for (size_t m = 0; m < 12; ++m)
    {
        if (m < a.months.size() && m < b.months.size())
        {
            max_spei_diff = std::max(max_spei_diff, std::abs(a.months[m].spei - b.months[m].spei));
            max_agri_diff = std::max(max_agri_diff, std::abs(a.months[m].agri_spei - b.months[m].agri_spei));
        }
    }
    if (max_spei_diff < 0.08f && max_agri_diff < 0.08f) return true;

    bool dur_match = (a.duration == b.duration);
    float peak_diff = std::abs(a.peak_intensity - b.peak_intensity);
    float sev_diff = std::abs(a.severity - b.severity);
    if (dur_match && peak_diff < 0.05f && sev_diff < 0.18f) return true;

    return false;
}

/**
 * Checks whether an event sequence has sufficient discrete event diversity.
 */
inline bool HasSignificantEventOverlap(const std::vector<DroughtEvent>& events)
{
    if (events.size() < 5) return false;
    int dupe_count = 0;
    for (size_t i = 0; i < events.size(); ++i)
    {
        for (size_t j = i + 1; j < events.size(); ++j)
        {
            if (AreEventsTooSimilar(events[i], events[j]))
            {
                dupe_count++;
            }
        }
    }
    return dupe_count > 3;
}

/**
 * Rescales an authentic event whose peak was artificially clamped to -3.00 or -3.72 in the raw COG raster,
 * bringing it into a physically realistic extreme range (-2.35 to -2.45).
 */
inline DroughtEvent TemperExtremeClampedEvent(
    const DroughtEvent& ev,
    float target_peak,
    const DroughtDistributionParams& params)
{
    if (ev.peak_intensity > -2.65f) return ev;

    float alpha = params.alpha;
    float beta = params.beta;
    float gamma = params.gamma;

    float current_peak_excursion = std::abs(ev.peak_intensity + 1.0f);
    float target_peak_excursion = std::abs(target_peak + 1.0f);
    float scale = std::min(1.0f, target_peak_excursion / std::max(0.01f, current_peak_excursion));

    int duration = 0;
    float severity = 0.0f;
    float min_spei = 0.0f;
    float total_deficit = 0.0f;
    int agri_duration = 0;
    float agri_severity = 0.0f;
    float min_agri_spei = 0.0f;
    float total_agri_deficit = 0.0f;
    float total_buffer_mitigated = 0.0f;

    std::vector<MonthlyWaterBalancePoint> months = ev.months;
    for (auto& m : months)
    {
        float spei = m.spei;
        if (spei <= -1.0f)
        {
            float exc = std::abs(spei + 1.0f);
            spei = -1.0f - exc * scale;
        }
        spei = std::clamp(spei, -2.60f, 3.0f);
        float d = SPEIToDeficit(spei, alpha, beta, gamma);
        float pet = m.pet;
        float p = std::max(0.0f, pet + d);

        float sm = m.soil_moisture > 0.0f ? m.soil_moisture : 0.25f;
        float raw_soil_buffer = ComputeAvailableSoilWaterBuffer(sm, 0.60f);
        float buffer_mitigation_fraction = std::clamp(raw_soil_buffer / 160.0f, 0.15f, 0.25f);
        float agri_spei = spei;
        float buf_used = 0.0f;
        float agri_d = d;
        if (spei < 0.0f)
        {
            agri_spei = spei * (1.0f - buffer_mitigation_fraction);
            buf_used = std::min(std::abs(d) * buffer_mitigation_fraction, raw_soil_buffer);
            agri_d = d + buf_used;
        }

        if (spei <= -1.0f)
        {
            duration++;
            severity += std::abs(spei + 1.0f);
        }
        if (spei < min_spei) min_spei = spei;
        if (d < 0.0f) total_deficit += std::abs(d);

        if (agri_spei <= -1.0f)
        {
            agri_duration++;
            agri_severity += std::abs(agri_spei + 1.0f);
        }
        if (agri_spei < min_agri_spei) min_agri_spei = agri_spei;
        if (agri_d < 0.0f) total_agri_deficit += std::abs(agri_d);
        total_buffer_mitigated += buf_used;

        m.precipitation = p;
        m.deficit = d;
        m.spei = spei;
        m.agri_deficit = agri_d;
        m.agri_spei = agri_spei;
        m.is_drought_month = (spei <= -1.0f);
        m.is_severe_drought = (spei <= -1.5f);
        m.is_extreme_drought = (spei <= -2.0f);
        m.is_agri_drought_month = (agri_spei <= -1.0f);
    }

    DroughtEvent res = ev;
    res.duration = std::max(1, duration);
    res.severity = severity;
    res.peak_intensity = min_spei;
    res.total_deficit_mm = total_deficit;
    res.months = months;
    res.agri_duration = std::max(0, agri_duration);
    res.agri_severity = agri_severity;
    res.agri_peak_intensity = min_agri_spei;
    res.total_agri_deficit_mm = total_agri_deficit;
    res.soil_buffer_mitigated_mm = total_buffer_mitigated;
    return res;
}

/**
 * Builds a distinct, physically consistent 12-month drought progression hydrograph.
 * Generates a gradual decreasing SPEI towards the peak of the drought, then back to 0,
 * with non-zero physical deficits (P - PET) and realistic plant root-zone soil buffering.
 * Exactly matches droughtutils.ts:buildDroughtProgressionHydrograph().
 */
inline DroughtEvent BuildDroughtProgressionHydrograph(
    int rank,
    float target_peak,
    const DroughtDistributionParams& params,
    const std::vector<float>& monthly_sm,
    const std::vector<MonthlyClimatologyPoint>& annual_clim,
    double seed = 42.0,
    int custom_peak_month = -1,
    float custom_duration = -1.0f)
{
    (void)seed;
    float alpha = params.alpha;
    float beta = params.beta;
    float gamma = params.gamma;

    // Center the drought in the 12-month window:
    int default_peak_month = 5 + (((rank * 2 + 1) % 3) - 1);
    int peak_month_idx = (custom_peak_month >= 0) ? custom_peak_month : default_peak_month;

    // Duration of drought (months with SPEI <= -1.0):
    // Higher severity ranks last 5-6 months; lower ranks last ~2 months
    double rank_ratio = (rank - 1) / 19.0;
    float default_duration = static_cast<float>(5.2 - std::pow(rank_ratio, 0.6) * 3.2);
    float duration_months = (custom_duration > 0.0f) ? custom_duration : default_duration;
    float sigma = std::max(0.75f, duration_months / 2.355f);

    std::vector<MonthlyWaterBalancePoint> months(12);
    int duration = 0;
    float severity = 0.0f;
    float min_spei = 0.0f;
    float total_deficit = 0.0f;
    float sum_p = 0.0f;
    float sum_pet = 0.0f;
    int agri_duration = 0;
    float agri_severity = 0.0f;
    float min_agri_spei = 0.0f;
    float total_agri_deficit = 0.0f;
    float total_buffer_mitigated = 0.0f;

    for (int m = 0; m < 12; ++m)
    {
        float clim_pet = (m < (int)annual_clim.size()) ? annual_clim[m].pet : 60.0f;
        float clim_p   = (m < (int)annual_clim.size()) ? annual_clim[m].precipitation : 60.0f;
        float clim_d   = clim_p - clim_pet;

        // Smooth, monotonic Gaussian curve centered on peakMonth:
        float diff = static_cast<float>(m - peak_month_idx);
        float gauss_factor = std::exp(-(diff * diff) / (2.0f * sigma * sigma));
        float spei = target_peak * gauss_factor;
        spei = std::clamp(spei, -2.60f, 0.20f);

        // Deficit anomaly relative to climatology from Log-Logistic distribution:
        // Median deficit (SPEI = 0) is gamma + alpha
        float median_d = gamma + alpha;
        float d_anom = SPEIToDeficit(spei, alpha, beta, gamma) - median_d;
        float d = clim_d + d_anom;

        // Physical precipitation and PET response during drought:
        float p = std::max(0.0f, clim_p + d_anom * 0.5f);
        float pet = std::max(clim_pet, p - d);
        // Guarantee P - PET === d identically:
        p = pet + d;
        if (p < 0.0f)
        {
            p = 0.0f;
            pet = -d;
        }

        // Root-zone soil moisture buffer:
        float clim_sm = (m < (int)monthly_sm.size()) ? monthly_sm[m] : 0.25f;
        float depletion = (d < 0.0f) ? std::min(0.12f, std::abs(d) / 400.0f) : 0.0f;
        float sm = std::max(0.04f, clim_sm - depletion);
        float raw_soil_buffer = ComputeAvailableSoilWaterBuffer(sm, 0.60f);

        // Soil buffering mitigates ~20-25% of meteorological stress
        float buffer_mitigation_fraction = std::clamp(raw_soil_buffer / 160.0f, 0.15f, 0.25f);
        float agri_spei = spei;
        float buf_used = 0.0f;
        float agri_d = d;
        if (spei < 0.0f)
        {
            agri_spei = spei * (1.0f - buffer_mitigation_fraction);
            buf_used = std::abs(d_anom) * buffer_mitigation_fraction;
            agri_d = d + buf_used;
        }

        if (spei <= -1.0f)
        {
            duration++;
            severity += std::abs(spei + 1.0f);
            if (d < 0.0f) total_deficit += std::abs(d);
        }
        if (spei < min_spei) min_spei = spei;

        if (agri_spei <= -1.0f)
        {
            agri_duration++;
            agri_severity += std::abs(agri_spei + 1.0f);
            if (agri_d < 0.0f) total_agri_deficit += std::abs(agri_d);
        }
        if (agri_spei < min_agri_spei) min_agri_spei = agri_spei;
        if (spei <= -1.0f || agri_spei <= -1.0f)
        {
            total_buffer_mitigated += buf_used;
        }

        sum_p += p;
        sum_pet += pet;

        months[m].month_index = m + 1;
        months[m].month_label = "M" + std::to_string(m + 1);
        months[m].precipitation = p;
        months[m].pet = pet;
        months[m].deficit = d;
        months[m].spei = spei;
        months[m].is_drought_month = (spei <= -1.0f);
        months[m].is_severe_drought = (spei <= -1.5f);
        months[m].is_extreme_drought = (spei <= -2.0f);
        months[m].soil_moisture = sm;
        months[m].soil_buffer_mm = raw_soil_buffer;
        months[m].agri_deficit = agri_d;
        months[m].agri_spei = agri_spei;
        months[m].is_agri_drought_month = (agri_spei <= -1.0f);
        months[m].is_valid = true;
    }

    if (total_deficit == 0.0f && months[peak_month_idx].deficit < 0.0f)
    {
        total_deficit = std::abs(months[peak_month_idx].deficit);
    }
    if (total_agri_deficit == 0.0f && agri_duration > 0 && months[peak_month_idx].agri_deficit < 0.0f)
    {
        total_agri_deficit = std::abs(months[peak_month_idx].agri_deficit);
    }
    if (total_buffer_mitigated == 0.0f && total_deficit > 0.0f)
    {
        total_buffer_mitigated = std::max(0.0f, total_deficit - total_agri_deficit);
    }

    DroughtEvent ev;
    ev.rank = rank;
    ev.duration = std::max(1, duration);
    ev.severity = severity;
    ev.peak_intensity = min_spei;
    ev.total_deficit_mm = total_deficit;
    ev.mean_precip_mm = sum_p / 12.0f;
    ev.mean_pet_mm = sum_pet / 12.0f;
    ev.months = months;
    ev.agri_duration = std::max(0, agri_duration);
    ev.agri_severity = agri_severity;
    ev.agri_peak_intensity = min_agri_spei;
    ev.total_agri_deficit_mm = total_agri_deficit;
    ev.soil_buffer_mitigated_mm = total_buffer_mitigated;
    return ev;
}

/**
 * Mixes real historical event data with physical distribution scaling to produce realistic,
 * properly ranging events in terms of timing (seasonality), magnitude (return period probability),
 * and duration. Exactly matches fetchTop20DroughtEvents() blending logic in droughtutils.ts.
 */
inline std::vector<DroughtEvent> BlendRealAndSyntheticTopEvents(
    const DroughtDistributionParams& params,
    const std::vector<HistoricalEventFeature>& real_features = {},
    int count = 20)
{
    double lat = params.lat;
    double lng = params.lng;
    double seed = std::abs(lat * 1000.0 + lng * 100.0);
    auto baseline_sm = GenerateSyntheticSoilMoisture(lat, lng);
    auto annual_clim = GetAnnualClimatology(lat, lng);
    int summer_peak = (lat >= 0.0) ? 6 : 0; // July (idx 6) in NH, Jan (idx 0) in SH

    std::vector<DroughtEvent> final_events;
    std::vector<int> used_timing_peaks;

    for (int rank = 1; rank <= count; ++rank)
    {
        float target_peak = GetTargetProbabilityPeakSPEI(rank, seed);
        double rank_progress = (rank - 1) / 19.0;
        float syn_duration = static_cast<float>(5.2 - std::pow(rank_progress, 0.6) * 3.2);

        float blended_peak = target_peak;
        int blended_month = (summer_peak - 2 + ((rank * 2 + 1) % 5) + 12) % 12;
        int blended_duration = static_cast<int>(std::round(syn_duration));

        bool has_real = (rank - 1 < (int)real_features.size() && real_features[rank - 1].has_data);
        if (has_real)
        {
            const auto& rf = real_features[rank - 1];
            if (rf.peak_spei <= -2.60f)
            {
                // Temper raw COG raster clamping (-3.00 or -3.72)
                blended_peak = 0.45f * rf.peak_spei * 0.78f + 0.55f * target_peak;
            }
            else
            {
                blended_peak = 0.60f * rf.peak_spei + 0.40f * target_peak;
            }
            blended_duration = std::clamp(static_cast<int>(std::round(0.6f * rf.duration + 0.4f * syn_duration)), 2, 6);
            blended_month = rf.peak_month;
        }

        // Clamp magnitude to realistic physical limits:
        blended_peak = std::clamp(blended_peak, -2.55f, -1.02f);

        // Timing variety: If this peak month already appeared more than twice, stagger
        int same_month_count = 0;
        for (int m : used_timing_peaks) {
            if (m == blended_month) same_month_count++;
        }
        if (same_month_count >= 2)
        {
            int shift = (rank % 2 == 0 ? 1 : -1) * (same_month_count / 2);
            blended_month = (blended_month + shift + 12) % 12;
        }
        used_timing_peaks.push_back(blended_month);

        auto ev = BuildDroughtProgressionHydrograph(
            rank,
            blended_peak,
            params,
            baseline_sm,
            annual_clim,
            seed,
            blended_month,
            static_cast<float>(blended_duration)
        );
        final_events.push_back(ev);
    }

    // Sort strictly by peak intensity from Rank 1 (worst) to Rank 20 (mildest)
    std::sort(final_events.begin(), final_events.end(), [](const DroughtEvent& a, const DroughtEvent& b) {
        return a.peak_intensity < b.peak_intensity;
    });
    for (size_t idx = 0; idx < final_events.size(); ++idx)
    {
        final_events[idx].rank = static_cast<int>(idx + 1);
    }

    return final_events;
}

/**
 * Generates synthetic top drought events for a given location and distribution.
 * Implements Yevjevich's Theory of Runs on monthly precipitation & PET series.
 */
inline std::vector<DroughtEvent> GenerateSyntheticTopEvents(
    const DroughtDistributionParams& params,
    int count = 20)
{
    return BlendRealAndSyntheticTopEvents(params, {}, count);
}

// -----------------------------------------------------------------------------
// Theory of Runs Drought Evaluation
// -----------------------------------------------------------------------------

/**
 * Executes Yevjevich's Theory of Runs over a multi-month drought timeseries
 * for both Meteorological Drought (pure atmospheric SPEI) and Agricultural Drought
 * (root-zone soil moisture buffered Agri-SPEI).
 * Standard drought threshold: SPEI <= -1.0 (moderate drought threshold).
 * Severe: SPEI <= -1.5, Extreme: SPEI <= -2.0.
 */
inline CellRunHazardMetrics AnalyzeDroughtRuns(
    const std::vector<MonthlyWaterBalance>& months,
    float drought_threshold = -1.0f)
{
    CellRunHazardMetrics metrics;
    metrics.peak_intensity = 0.0f;
    metrics.agri_peak_intensity = 0.0f;

    int valid_month_count = 0;
    float sum_spei = 0.0f;
    float sum_agri_spei = 0.0f;
    float sum_sm = 0.0f;
    bool in_drought_episode = false;
    bool in_agri_drought_episode = false;

    for (const auto& m : months)
    {
        if (!m.is_valid || pcr::isMV(m.spei) || !std::isfinite(m.spei))
        {
            continue;
        }

        valid_month_count++;
        sum_spei += m.spei;

        float a_spei = (!pcr::isMV(m.agri_spei) && std::isfinite(m.agri_spei)) ? m.agri_spei : m.spei;
        sum_agri_spei += a_spei;
        sum_sm += m.soil_moisture;

        // 1. Meteorological Drought (Atmospheric SPEI)
        if (m.spei < metrics.peak_intensity)
        {
            metrics.peak_intensity = m.spei;
        }

        float cur_def = 0.0f;
        if (m.deficit < 0.0f)
        {
            cur_def = -m.deficit;
        }
        else if (m.deficit > 0.0f && m.spei <= 0.0f)
        {
            cur_def = m.deficit;
        }
        else if (m.precipitation < m.pet)
        {
            cur_def = m.pet - m.precipitation;
        }

        if (cur_def > 0.0f && m.spei <= 0.0f)
        {
            metrics.total_deficit_mm += cur_def;
        }

        if (m.spei <= drought_threshold)
        {
            metrics.drought_months++;
            metrics.total_severity += std::abs(m.spei - drought_threshold);

            if (!in_drought_episode)
            {
                in_drought_episode = true;
                metrics.event_count_drought++;
            }
        }
        else
        {
            in_drought_episode = false;
        }

        // 2. Agricultural / Soil Moisture Drought (Root-Zone Water Availability)
        if (a_spei < metrics.agri_peak_intensity)
        {
            metrics.agri_peak_intensity = a_spei;
        }

        float cur_agri_def = 0.0f;
        if (m.agri_deficit < 0.0f)
        {
            cur_agri_def = -m.agri_deficit;
        }

        if (cur_agri_def > 0.0f && a_spei <= 0.0f)
        {
            metrics.agri_total_deficit_mm += cur_agri_def;
        }

        if (a_spei <= drought_threshold)
        {
            metrics.agri_drought_months++;
            metrics.agri_total_severity += std::abs(a_spei - drought_threshold);

            if (!in_agri_drought_episode)
            {
                in_agri_drought_episode = true;
                metrics.agri_event_count++;
            }
        }
        else
        {
            in_agri_drought_episode = false;
        }
    }

    if (valid_month_count > 0)
    {
        metrics.mean_spei = sum_spei / valid_month_count;
        metrics.agri_mean_spei = sum_agri_spei / valid_month_count;
        metrics.mean_soil_moisture = sum_sm / valid_month_count;
    }

    return metrics;
}

/**
 * Applies physical climate and adaptation adjustments to a single monthly step,
 * utilizing continuous standardized linear scaling (deltaSPEI = deltaD / sigma_D)
 * exactly matching simulateDroughtScenario() in droughtutils.ts.
 */
inline MonthlyWaterBalance ApplyScenarioAndAdaptation(
    const MonthlyWaterBalance& baseline,
    const DroughtDistributionParams& dist,
    const DroughtScenarioModifiers& mod,
    float* inout_soil_storage = nullptr)
{
    MonthlyWaterBalance adapted;
    adapted.is_valid = baseline.is_valid;
    if (!baseline.is_valid || pcr::isMV(baseline.precipitation) || pcr::isMV(baseline.pet))
    {
        pcr::setMV(adapted.precipitation);
        pcr::setMV(adapted.pet);
        pcr::setMV(adapted.deficit);
        pcr::setMV(adapted.spei);
        pcr::setMV(adapted.soil_moisture);
        pcr::setMV(adapted.soil_buffer_mm);
        pcr::setMV(adapted.agri_deficit);
        pcr::setMV(adapted.agri_spei);
        adapted.is_valid = false;
        return adapted;
    }

    // 1. Precipitation with Climate Scaling and Agricultural Irrigation: P* = max(0, P * rain_mult + irrig)
    float p_mod = std::max(0.0f, baseline.precipitation * mod.rain_multiplier + mod.irrigation_supply_mm);

    // 2. PET with Thermodynamic Warming, Wind Scaling, and Evaporation Reduction: PET* = max(5.0, PET * net_pet_factor)
    float pet_clim_mult = mod.pet_multiplier;
    if (pet_clim_mult <= 0.0f)
    {
        pet_clim_mult = ComputePETClimateMultiplier(mod.delta_temp_c, mod.wind_multiplier);
    }
    float et_fraction = std::clamp(mod.et_reduction_fraction, 0.0f, 0.50f);
    float net_pet_factor = pet_clim_mult * std::max(0.20f, 1.0f - et_fraction);
    float pet_mod = std::max(5.0f, baseline.pet * net_pet_factor);

    // 3. Adjusted Climatic Water Balance Deficit
    float delta_d = (p_mod - baseline.precipitation) - (pet_mod - baseline.pet);
    float d_mod = baseline.deficit + delta_d;

    // 4. Standard deviation of climatic water balance for Log-Logistic distribution: sigma_D = 1.81 * alpha / beta
    // Continuous scaling prevents non-linear -3.00 CDF floor collapse
    float safe_beta = std::clamp(dist.beta, 0.8f, 4.0f);
    float safe_alpha = std::max(15.0f, dist.alpha);
    float sigma_d = std::max(20.0f, (1.81f * safe_alpha) / safe_beta);

    float delta_spei = delta_d / sigma_d;
    float spei_mod = std::clamp(baseline.spei + delta_spei, -3.0f, 3.0f);

    // 5. Agricultural Drought Deficit with Soil Moisture Storage and Retention Buffer
    float b_sm = baseline.soil_moisture;
    float b_buffer = (baseline.soil_buffer_mm > 0.0f) ? baseline.soil_buffer_mm : ComputeAvailableSoilWaterBuffer(b_sm, mod.root_zone_depth_m);

    float moisture_gain = ((p_mod - baseline.precipitation) / 450.0f) + (mod.soil_retention_buffer_mm / 500.0f);
    float s_sm = std::clamp(b_sm + moisture_gain, 0.04f, 0.55f);
    float s_buffer = ComputeAvailableSoilWaterBuffer(s_sm, mod.root_zone_depth_m);

    float additional_buffer = (s_buffer - b_buffer) + mod.soil_retention_buffer_mm;
    float delta_agri_d = delta_d + additional_buffer;
    float d_agri_mod = baseline.agri_deficit + delta_agri_d;

    float delta_agri_spei = delta_agri_d / sigma_d;
    float agri_spei_mod = std::clamp(baseline.agri_spei + delta_agri_spei, -3.0f, 3.0f);

    float current_storage = s_buffer + mod.soil_retention_buffer_mm;
    if (inout_soil_storage != nullptr)
    {
        *inout_soil_storage = current_storage;
    }

    adapted.month_index = baseline.month_index;
    adapted.month_label = baseline.month_label;
    adapted.precipitation = p_mod;
    adapted.pet = pet_mod;
    adapted.deficit = d_mod;
    adapted.spei = spei_mod;
    adapted.is_drought_month = (spei_mod <= -1.0f);
    adapted.is_severe_drought = (spei_mod <= -1.5f);
    adapted.is_extreme_drought = (spei_mod <= -2.0f);
    adapted.soil_moisture = s_sm;
    adapted.soil_buffer_mm = current_storage;
    adapted.agri_deficit = d_agri_mod;
    adapted.agri_spei = agri_spei_mod;
    adapted.is_agri_drought_month = (agri_spei_mod <= -1.0f);
    return adapted;
}

// -----------------------------------------------------------------------------
// Climate Scenarios & Adaptation Stress-Test Engine
// -----------------------------------------------------------------------------

/**
 * Executes the Physical Adaptation & Climate Change Simulation Engine:
 * Adjusts monthly water balance D* = (P* + Irrigation) - (PET* * (1 - ET_reduction)) + Retention
 * and evaluates SPEI* and SPEI_agri* using continuous standardized perturbation scaling.
 * Exactly matches droughtutils.ts:simulateDroughtScenario().
 */
inline std::vector<EvaluatedEventScenario> SimulateDroughtScenario(
    const std::vector<DroughtEvent>& events,
    const DroughtDistributionParams& params,
    const DroughtScenarioModifiers& modifiers)
{
    float alpha = params.alpha;
    float beta = params.beta;
    float gamma = params.gamma;

    float pet_climate_mult = modifiers.pet_multiplier;
    if (pet_climate_mult <= 0.0f)
    {
        pet_climate_mult = ComputePETClimateMultiplier(modifiers.delta_temp_c, modifiers.wind_multiplier);
    }
    float et_fraction = std::clamp(modifiers.et_reduction_fraction, 0.0f, 0.50f);
    float net_pet_factor = pet_climate_mult * std::max(0.20f, 1.0f - et_fraction);

    // Standard deviation of climatic water balance for Log-Logistic distribution: sigma_D = 1.81 * alpha / beta
    float safe_beta = std::clamp(beta, 0.8f, 4.0f);
    float safe_alpha = std::max(15.0f, alpha);
    float sigma_d = std::max(20.0f, (1.81f * safe_alpha) / safe_beta);

    std::vector<EvaluatedEventScenario> evaluated;
    evaluated.reserve(events.size());

    for (const auto& ev : events)
    {
        float baseline_severity = 0.0f;
        float scenario_severity = 0.0f;
        float base_peak_spei = 0.0f;
        float scen_peak_spei = 0.0f;
        int base_duration = 0;
        int scen_duration = 0;
        float base_deficit_mm = 0.0f;
        float scen_deficit_mm = 0.0f;

        float base_agri_severity = 0.0f;
        float scen_agri_severity = 0.0f;
        float base_agri_peak = 0.0f;
        float scen_agri_peak = 0.0f;
        int base_agri_duration = 0;
        int scen_agri_duration = 0;
        float base_agri_deficit_mm = 0.0f;
        float scen_agri_deficit_mm = 0.0f;
        float soil_buffer_mitigated_mm = 0.0f;
        float scen_soil_buffer_mitigated_mm = 0.0f;

        EvaluatedEventScenario ev_scen;
        ev_scen.rank = ev.rank;
        ev_scen.months.resize(ev.months.size());

        for (size_t m_idx = 0; m_idx < ev.months.size(); ++m_idx)
        {
            const auto& m = ev.months[m_idx];
            float b_p = m.precipitation;
            float b_pet = m.pet;
            float b_d = m.deficit;
            float b_spei = m.spei;
            if (b_spei <= -3.2f || !std::isfinite(b_spei))
            {
                b_spei = DeficitToSPEI(b_d, alpha, beta, gamma);
            }
            b_spei = std::clamp(b_spei, -3.0f, 3.0f);

            float b_sm = m.soil_moisture;
            float b_buffer = (m.soil_buffer_mm > 0.0f) ? m.soil_buffer_mm : ComputeAvailableSoilWaterBuffer(b_sm, modifiers.root_zone_depth_m);
            float b_agri_deficit = (m.agri_deficit != 0.0f) ? m.agri_deficit : b_d;
            float b_agri_spei = m.agri_spei;
            if (b_agri_spei <= -3.2f || !std::isfinite(b_agri_spei))
            {
                b_agri_spei = DeficitToSPEI(b_agri_deficit, alpha, beta, gamma);
            }
            b_agri_spei = std::clamp(b_agri_spei, -3.0f, 3.0f);

            // Adjusted Precipitation: P* = P * rainMultiplier + irrigation
            float s_p = std::max(0.0f, b_p * modifiers.rain_multiplier + modifiers.irrigation_supply_mm);
            // Adjusted Potential Evapotranspiration
            float s_pet = std::max(5.0f, b_pet * net_pet_factor);

            // Meteorological water balance deficit change relative to baseline:
            // Evaluates to 0.0 at baseline identically
            float delta_d = (s_p - b_p) - (s_pet - b_pet);
            float s_d = b_d + delta_d;

            // Standardized meteorological SPEI perturbation:
            float delta_spei = delta_d / sigma_d;
            float s_spei = std::clamp(b_spei + delta_spei, -3.0f, 3.0f);

            // Scenario Soil Moisture & Agricultural Buffer:
            float moisture_gain = ((s_p - b_p) / 450.0f) + (modifiers.soil_retention_buffer_mm / 500.0f);
            float s_sm = std::clamp(b_sm + moisture_gain, 0.04f, 0.55f);
            float s_buffer = ComputeAvailableSoilWaterBuffer(s_sm, modifiers.root_zone_depth_m);

            // Agricultural water balance deficit change relative to baseline:
            float additional_buffer = (s_buffer - b_buffer) + modifiers.soil_retention_buffer_mm;
            float delta_agri_d = delta_d + additional_buffer;
            float s_agri_deficit = b_agri_deficit + delta_agri_d;

            // Standardized agricultural SPEI perturbation:
            float delta_agri_spei = delta_agri_d / sigma_d;
            float s_agri_spei = std::clamp(b_agri_spei + delta_agri_spei, -3.0f, 3.0f);

            bool is_drought_base = (b_spei <= -1.0f);
            bool is_drought_scen = (s_spei <= -1.0f);
            bool is_agri_drought_base = (b_agri_spei <= -1.0f);
            bool is_agri_drought_scen = (s_agri_spei <= -1.0f);

            if (is_drought_base)
            {
                base_duration++;
                baseline_severity += std::abs(b_spei + 1.0f);
                if (b_d < 0.0f) base_deficit_mm += std::abs(b_d);
            }
            if (is_drought_scen)
            {
                scen_duration++;
                scenario_severity += std::abs(s_spei + 1.0f);
                if (s_d < 0.0f) scen_deficit_mm += std::abs(s_d);
            }

            if (is_agri_drought_base)
            {
                base_agri_duration++;
                base_agri_severity += std::abs(b_agri_spei + 1.0f);
                if (b_agri_deficit < 0.0f) base_agri_deficit_mm += std::abs(b_agri_deficit);
            }
            if (is_agri_drought_scen)
            {
                scen_agri_duration++;
                scen_agri_severity += std::abs(s_agri_spei + 1.0f);
                if (s_agri_deficit < 0.0f) scen_agri_deficit_mm += std::abs(s_agri_deficit);
            }

            if (is_drought_base || is_agri_drought_base)
            {
                soil_buffer_mitigated_mm += std::max(0.0f, std::abs(b_d) - std::abs(b_agri_deficit));
            }
            if (is_drought_scen || is_agri_drought_scen)
            {
                scen_soil_buffer_mitigated_mm += std::max(0.0f, std::abs(s_d) - std::abs(s_agri_deficit));
            }

            if (b_spei < base_peak_spei) base_peak_spei = b_spei;
            if (s_spei < scen_peak_spei) scen_peak_spei = s_spei;
            if (b_agri_spei < base_agri_peak) base_agri_peak = b_agri_spei;
            if (s_agri_spei < scen_agri_peak) scen_agri_peak = s_agri_spei;

            auto& em = ev_scen.months[m_idx];
            em.month_index = m.month_index;
            em.month_label = m.month_label;
            em.baseline_precip = b_p;
            em.scenario_precip = s_p;
            em.baseline_pet = b_pet;
            em.scenario_pet = s_pet;
            em.baseline_deficit = b_d;
            em.scenario_deficit = s_d;
            em.baseline_spei = b_spei;
            em.scenario_spei = s_spei;
            em.baseline_soil_moisture = b_sm;
            em.scenario_soil_moisture = s_sm;
            em.baseline_soil_buffer_mm = b_buffer;
            em.scenario_soil_buffer_mm = s_buffer;
            em.baseline_agri_deficit = b_agri_deficit;
            em.scenario_agri_deficit = s_agri_deficit;
            em.baseline_agri_spei = b_agri_spei;
            em.scenario_agri_spei = s_agri_spei;
            em.is_drought_base = is_drought_base;
            em.is_drought_scen = is_drought_scen;
            em.is_agri_drought_base = is_agri_drought_base;
            em.is_agri_drought_scen = is_agri_drought_scen;
        }

        if (base_deficit_mm == 0.0f && ev.total_deficit_mm > 0.0f) base_deficit_mm = ev.total_deficit_mm;
        if (base_agri_deficit_mm == 0.0f && ev.total_agri_deficit_mm > 0.0f) base_agri_deficit_mm = ev.total_agri_deficit_mm;
        if (soil_buffer_mitigated_mm == 0.0f && ev.soil_buffer_mitigated_mm > 0.0f) soil_buffer_mitigated_mm = ev.soil_buffer_mitigated_mm;
        if (scen_deficit_mm == 0.0f && base_deficit_mm > 0.0f) scen_deficit_mm = base_deficit_mm;
        if (scen_agri_deficit_mm == 0.0f && base_agri_deficit_mm > 0.0f) scen_agri_deficit_mm = base_agri_deficit_mm;

        float delta_sev = scenario_severity - baseline_severity;
        float sev_change_pct = (baseline_severity > 0.0f) ? (delta_sev / baseline_severity) * 100.0f : 0.0f;
        float net_ben = std::max(0.0f, scen_deficit_mm - scen_agri_deficit_mm);

        ev_scen.baseline_severity = baseline_severity;
        ev_scen.scenario_severity = scenario_severity;
        ev_scen.severity_change_percent = sev_change_pct;
        ev_scen.baseline_peak_spei = base_peak_spei;
        ev_scen.scenario_peak_spei = scen_peak_spei;
        ev_scen.delta_peak_spei = scen_peak_spei - base_peak_spei;
        ev_scen.baseline_duration = base_duration;
        ev_scen.scenario_duration = scen_duration;
        ev_scen.baseline_deficit_mm = base_deficit_mm;
        ev_scen.scenario_deficit_mm = scen_deficit_mm;
        ev_scen.delta_deficit_mm = scen_deficit_mm - base_deficit_mm;
        ev_scen.net_benefit_mm = net_ben;
        ev_scen.baseline_agri_severity = base_agri_severity;
        ev_scen.scenario_agri_severity = scen_agri_severity;
        ev_scen.baseline_agri_peak_spei = base_agri_peak;
        ev_scen.scenario_agri_peak_spei = scen_agri_peak;
        ev_scen.baseline_agri_duration = base_agri_duration;
        ev_scen.scenario_agri_duration = scen_agri_duration;
        ev_scen.agri_duration = base_agri_duration;
        ev_scen.baseline_agri_deficit_mm = base_agri_deficit_mm;
        ev_scen.scenario_agri_deficit_mm = scen_agri_deficit_mm;
        ev_scen.soil_buffer_mitigated_mm = soil_buffer_mitigated_mm;
        ev_scen.scenario_soil_buffer_mitigated_mm = scen_soil_buffer_mitigated_mm;

        evaluated.push_back(ev_scen);
    }

    return evaluated;
}

// -----------------------------------------------------------------------------
// Site Drought Risk Score Calculation
// -----------------------------------------------------------------------------

/**
 * Calculates a consolidated 0-100 Drought Risk Score and classification
 * for site planners and non-expert stakeholders.
 * Physically calibrated based on UNEP Aridity Index (P / PET), seasonal climatic deficit,
 * and plant root-zone soil moisture buffering (FAO-56).
 * Exactly matches droughtutils.ts:calculateSiteDroughtRiskScore().
 */
inline SiteDroughtRiskScore CalculateSiteDroughtRiskScore(
    const DroughtDistributionParams& params,
    const std::vector<DroughtEvent>& top_events = {},
    const std::vector<MonthlyClimatologyPoint>& custom_annual_clim = {})
{
    auto annual_clim = (!custom_annual_clim.empty()) ? custom_annual_clim : GetAnnualClimatology(params.lat, params.lng);

    int deficit_months_count = 0;
    float total_summer_deficit = 0.0f;
    float total_annual_p = 0.0f;
    float total_annual_pet = 0.0f;

    for (const auto& m : annual_clim)
    {
        total_annual_p += m.precipitation;
        total_annual_pet += m.pet;
        if (m.is_deficit)
        {
            deficit_months_count++;
            total_summer_deficit += std::abs(m.deficit);
        }
    }

    float aridity_index = total_annual_p / std::max(1.0f, total_annual_pet);

    const DroughtEvent* worst_event = (!top_events.empty()) ? &top_events[0] : nullptr;
    int worst_duration = worst_event ? worst_event->duration : (deficit_months_count > 0 ? 5 : 2);
    float worst_peak = worst_event ? worst_event->peak_intensity : -1.8f;
    int worst_agri_duration = worst_event ? worst_event->agri_duration : (deficit_months_count > 0 ? 2 : 0);
    float worst_agri_peak = worst_event ? worst_event->agri_peak_intensity : (aridity_index < 0.5f ? -1.8f : -0.9f);

    float score_f = 0.0f;

    // Component 1: Climatological Aridity & Seasonal Deficit (0 to 45 points)
    if (aridity_index < 0.20f)
    {
        score_f += 38.0f + std::min(7.0f, (total_summer_deficit / 500.0f) * 7.0f);
    }
    else if (aridity_index < 0.55f)
    {
        score_f += 26.0f + std::min(14.0f, (total_summer_deficit / 350.0f) * 14.0f);
    }
    else if (aridity_index < 0.85f)
    {
        score_f += 16.0f + std::min(10.0f, (total_summer_deficit / 250.0f) * 10.0f);
    }
    else if (aridity_index < 1.25f)
    {
        score_f += 10.0f + std::min(8.0f, (total_summer_deficit / 120.0f) * 8.0f);
    }
    else
    {
        score_f += std::min(6.0f, (total_summer_deficit / 150.0f) * 6.0f);
    }

    // Component 2: Realized Drought Duration (0 to 30 points)
    int effective_duration = std::max(worst_agri_duration, static_cast<int>(std::round(worst_duration * 0.7f)));
    if (effective_duration <= 1)
    {
        score_f += 3.0f;
    }
    else if (effective_duration == 2)
    {
        score_f += 5.0f;
    }
    else if (effective_duration == 3)
    {
        score_f += 9.0f;
    }
    else if (effective_duration <= 5)
    {
        score_f += 12.0f + ((effective_duration - 3) / 2.0f) * 4.0f;
    }
    else
    {
        score_f += 17.0f + std::min(11.0f, (effective_duration - 5) * 2.2f);
    }

    // Component 3: Drought Event Peak Severity (0 to 25 points)
    float effective_peak_mag = std::max(std::abs(worst_agri_peak), std::abs(worst_peak) * 0.85f);
    if (effective_peak_mag <= 1.0f)
    {
        score_f += 3.0f;
    }
    else if (effective_peak_mag <= 1.4f)
    {
        score_f += 4.0f + ((effective_peak_mag - 1.0f) / 0.4f) * 3.0f;
    }
    else if (effective_peak_mag <= 2.0f)
    {
        score_f += 7.0f + ((effective_peak_mag - 1.4f) / 0.6f) * 4.0f;
    }
    else
    {
        score_f += 11.0f + std::min(7.0f, (effective_peak_mag - 2.0f) * 8.0f);
    }

    int score = static_cast<int>(std::round(std::clamp(score_f, 5.0f, 95.0f)));

    SiteDroughtRiskScore res;
    res.score = score;
    res.annual_deficit_months = deficit_months_count;
    res.typical_summer_deficit_mm = total_summer_deficit;
    res.drought_frequency_years = std::max(3, static_cast<int>(std::round(35.0f - (score / 100.0f) * 32.0f)));
    res.worst_duration_months = worst_duration;
    res.worst_agri_duration_months = worst_agri_duration;
    res.worst_peak_spei = worst_peak;
    res.worst_agri_peak_spei = worst_agri_peak;

    if (score >= 75)
    {
        res.level = "Extreme";
        res.color = "#dc2626";
        res.badge_class = "bg-rose-100 text-rose-800 border-rose-300";
        res.headline = "Extreme Drought Vulnerability";
        res.summary = "Severe structural water deficits with prolonged multi-month droughts. Emergency water storage and strict irrigation conservation are essential.";
    }
    else if (score >= 50)
    {
        res.level = "High";
        res.color = "#ea580c";
        res.badge_class = "bg-orange-100 text-orange-800 border-orange-300";
        res.headline = "High Drought Hazard";
        res.summary = "Evaporative demand significantly outstrips rainfall during extended periods. Vegetation and local reservoirs face stress during dry spells.";
    }
    else if (score >= 20)
    {
        res.level = "Moderate";
        res.color = "#d97706";
        res.badge_class = "bg-amber-100 text-amber-800 border-amber-300";
        res.headline = "Moderate Seasonal Drought Risk";
        res.summary = "Regular seasonal water deficits can cause soil moisture stress during dry years. Manageable with standard irrigation and soil water conservation.";
    }
    else
    {
        res.level = "Low";
        res.color = "#16a34a";
        res.badge_class = "bg-emerald-100 text-emerald-800 border-emerald-300";
        res.headline = "Low Overall Drought Vulnerability";
        res.summary = "Rainfall consistently keeps pace with evaporative demand. Strong soil moisture buffering prevents severe agricultural stress.";
    }

    res.recommended_mitigation = (score >= 50) ? "Drip irrigation + organic mulching + 30-day rainwater storage reserve" :
                                 ((score >= 20) ? "Soil organic matter enrichment + deficit irrigation scheduling" :
                                                  "Baseline soil cover and monitoring (natural water buffer sufficient)");

    return res;
}

/**
 * Backward-compatible overload for callers providing individual parameters.
 */
inline SiteDroughtRiskScore CalculateSiteDroughtRiskScore(
    double lat, double lng,
    const DroughtDistributionParams& params,
    int worst_duration = 7,
    float worst_peak = -2.3f)
{
    DroughtDistributionParams p = params;
    p.lat = lat;
    p.lng = lng;
    DroughtEvent dummy;
    dummy.duration = worst_duration;
    dummy.peak_intensity = worst_peak;
    dummy.agri_duration = std::max(0, worst_duration - 4);
    dummy.agri_peak_intensity = worst_peak + 0.8f;
    return CalculateSiteDroughtRiskScore(p, {dummy});
}

} // namespace FastDrought

// -----------------------------------------------------------------------------
// FastDrought Core Engine Class
// -----------------------------------------------------------------------------

class FastDroughtModel
{
public:
    FastDroughtModel() = default;

    /**
     * Generates location-tailored distribution parameters (Alpha, Beta, Gamma)
     * based on latitude/longitude climatology when COG is unavailable.
     */
    static FastDrought::DroughtDistributionParams PlaceholderDistribution(double lat, double lon)
    {
        return FastDrought::GenerateSyntheticDroughtParams(lat, lon);
    }

    /**
     * Generates a realistic discrete multi-month historical drought event progression hydrograph
     * for a given rank, matching extreme value order statistics and FAO-56 soil buffering.
     */
    static std::vector<FastDrought::MonthlyWaterBalance> GenerateSyntheticDroughtEvent(
        int rank,
        const FastDrought::DroughtDistributionParams& dist,
        int duration_months = 0,
        const std::vector<float>& monthly_sm = {},
        double lat = 45.0, double lng = 0.0)
    {
        double seed = std::abs(lat * 1000.0 + lng * 100.0);
        float target_peak = FastDrought::GetTargetProbabilityPeakSPEI(rank, seed);
        auto baseline_sm = (!monthly_sm.empty()) ? monthly_sm : FastDrought::GenerateSyntheticSoilMoisture(lat, lng);
        auto annual_clim = FastDrought::GetAnnualClimatology(lat, lng);

        auto ev = FastDrought::BuildDroughtProgressionHydrograph(
            rank,
            target_peak,
            dist,
            baseline_sm,
            annual_clim,
            seed,
            -1,
            static_cast<float>(duration_months)
        );
        return ev.months;
    }

    /**
     * Generates the Top 20 historical multi-month drought events blending real historical
     * event features (timing, duration, intensity) with synthetic order statistics and FAO-56 soil buffering.
     */
    static std::vector<FastDrought::DroughtEvent> GenerateTop20Events(
        const FastDrought::DroughtDistributionParams& dist,
        const std::vector<FastDrought::HistoricalEventFeature>& real_features = {},
        int count = 20)
    {
        return FastDrought::BlendRealAndSyntheticTopEvents(dist, real_features, count);
    }

    /**
     * Executes the physical adaptation and climate change scenario simulation engine
     * across the given set of drought events.
     */
    static std::vector<FastDrought::EvaluatedEventScenario> SimulateScenario(
        const std::vector<FastDrought::DroughtEvent>& events,
        const FastDrought::DroughtDistributionParams& dist,
        const FastDrought::DroughtScenarioModifiers& modifiers)
    {
        return FastDrought::SimulateDroughtScenario(events, dist, modifiers);
    }
};

#endif // FASTDROUGHT_H
