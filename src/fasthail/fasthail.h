//
// Created for FastHail
//

#ifndef FASTHAIL_H
#define FASTHAIL_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../../common/defines.h"
#include "../../common/fferror.h"
#include "../../common/threading/parallel.h"
#include "../../geospatial/map/map.h"

namespace FastHailModel
{
    // -------------------------------------------------------------------------
    // Dataset Constants (Cloudflare R2 / cogdata worker: era5l_hail/)
    // All COGs in era5l_hail/ are in EPSG:3857 (Float32, NODATA = -9999.0)
    // -------------------------------------------------------------------------
    inline const std::string kHailCogBaseUrl = "https://cogdata.fastflood.workers.dev/?path=";
    inline const std::string kHailGevParamsCog = "era5l_hail/hail_size_gev_parameters.tif";
    inline const std::string kHailClimateRelChangeSizeCog = "era5l_hail/hail_climate_relchange_size_48b.tif";
    inline const std::string kHailClimateRelChangeOccCog = "era5l_hail/hail_climate_relchange_occurrence_48b.tif";
    inline const std::string kHailAnnualExceedanceProbCog = "era5l_hail/hail_annual_exceedance_probability.tif";
    inline const std::string kHailMeanAnnualEventRateCog = "era5l_hail/hail_mean_annual_event_rate.tif";
    inline const std::string kHailSizeReturnPeriodYearsCog = "era5l_hail/hail_size_return_period_years.tif";
    inline const std::string kHailModelDiagnosticsCog = "era5l_hail/hail_model_diagnostics.tif";
    inline const std::string kHailAiBiasDiagnosticsCog = "era5l_hail/hail_ai_bias_corrected_diagnostics.tif";

    // Physical constants for hailstone kinematics
    inline constexpr double kHailIceDensityKgM3 = 900.0;     // Standard hard hailstone density [kg/m^3] (0.9 g/cm^3)
    inline constexpr double kSeaLevelAirDensityKgM3 = 1.225; // Standard sea-level air density [kg/m^3]
    inline constexpr double kAirScaleHeightM = 8500.0;       // Atmospheric density scale height [m]
    inline constexpr double kGravityMS2 = 9.80665;           // Gravitational acceleration [m/s^2]
    inline constexpr double kHailDragCoeff = 0.55;           // Mean drag coefficient Cd for natural hailstones
    inline constexpr double kWindCouplingFactor = 0.85;      // Horizontal storm wind coupling efficiency eta_w
    inline constexpr double kMaxHailDiameterCm = 25.0;       // Physical cap on hailstone diameter [cm]

    // -------------------------------------------------------------------------
    // GEV Return-Level Evaluation (matching era5l_hail/README.md convention)
    // Band 1: gev_shape_xi (xi)
    // Band 2: gev_location_cm (mu_AI, cm)
    // Band 3: gev_scale_cm (sigma_AI, cm)
    //
    // Let y = -ln(1 - 1/T).
    // If |xi| < 1e-8: d_T = mu_AI - sigma_AI * ln(y)
    // Else:           d_T = mu_AI + (sigma_AI / xi) * (y^(-xi) - 1)
    // -------------------------------------------------------------------------
    inline double evaluateGEVHailDiameterCm(double xi, double mu_cm, double sigma_cm, double returnPeriodYears)
    {
        if (!std::isfinite(xi) || !std::isfinite(mu_cm) || !std::isfinite(sigma_cm) || sigma_cm <= 0.0)
        {
            return 0.0;
        }

        double T = std::max(1.0001, returnPeriodYears);
        double p = 1.0 - (1.0 / T);
        double y = -std::log(p);
        if (y <= 0.0 || !std::isfinite(y))
        {
            return 0.0;
        }

        double d_cm = 0.0;
        if (std::abs(xi) < 1e-8)
        {
            d_cm = mu_cm - sigma_cm * std::log(y);
        }
        else
        {
            // Clamp xi to prevent numerical overflow in extreme tails
            double xi_clamped = std::clamp(xi, -0.75, 0.75);
            d_cm = mu_cm + (sigma_cm / xi_clamped) * (std::pow(y, -xi_clamped) - 1.0);
        }

        if (!std::isfinite(d_cm) || d_cm < 0.0)
        {
            return 0.0;
        }
        return std::min(kMaxHailDiameterCm, d_cm);
    }

    // Fallback latitude-based GEV hailstone diameter [cm] if remote COGs are unreachable
    inline double placeholderGEVHailDiameterCm(double latitudeDeg, double returnPeriodYears)
    {
        double absLat = std::abs(latitudeDeg);
        double T = std::clamp(returnPeriodYears, 1.5, 2000.0);

        // Peak convective hail belt in mid-latitudes (25-50 deg)
        double beltFactor = std::exp(-0.5 * std::pow((absLat - 37.0) / 16.0, 2.0));
        double mu_cm = 1.2 + 1.8 * beltFactor;
        double sigma_cm = 0.45 + 0.55 * beltFactor;
        double xi = -0.05;

        return evaluateGEVHailDiameterCm(xi, mu_cm, sigma_cm, T);
    }

    // Fallback 2-year thunderstorm wind speed [m/s] at 10m height if wind COGs are unreachable
    inline double placeholderStormWindSpeedMS(double latitudeDeg, double windReturnPeriodYears = 2.0)
    {
        double absLat = std::abs(latitudeDeg);
        double T = std::clamp(windReturnPeriodYears, 1.1, 1000.0);
        double base2yr = 13.5 + 3.0 * std::exp(-0.5 * std::pow((absLat - 38.0) / 18.0, 2.0));
        double scale = 1.8 * std::log(std::max(1.0, T / 2.0));
        return std::clamp(base2yr + scale, 5.0, 45.0);
    }

    // -------------------------------------------------------------------------
    // 48-Band Climate Change Indexing & Interpolation (matching era5l_hail/README.md)
    // Outer loop  - Horizon  (hIdx = 0..3): 2040, 2060, 2080, 2100 (12 bands per horizon)
    // Middle loop - Scenario (sIdx = 0..3): ssp126, ssp245, ssp370, ssp585 (3 bands per scenario)
    // Inner loop  - Quantile (qIdx = 0..2): q15, q50, q85
    // Band (1-based) = hIdx * 12 + sIdx * 3 + qIdx + 1
    // -------------------------------------------------------------------------
    inline int getClimateBand1Based(int hIdx, int sIdx, int qIdx)
    {
        int h = std::clamp(hIdx, 0, 3);
        int s = std::clamp(sIdx, 0, 3);
        int q = std::clamp(qIdx, 0, 2);
        return h * 12 + s * 3 + q + 1;
    }

    inline double interpolateHailClimateRelChange(
        double targetYear,
        double rel2040,
        double rel2060,
        double rel2080,
        double rel2100)
    {
        double yr = std::clamp(targetYear, 2020.0, 2100.0);
        if (yr <= 2020.0) return 0.0;
        if (yr <= 2040.0)
        {
            double f = (yr - 2020.0) / (2040.0 - 2020.0);
            return f * rel2040;
        }
        else if (yr <= 2060.0)
        {
            double f = (yr - 2040.0) / (2060.0 - 2040.0);
            return rel2040 + f * (rel2060 - rel2040);
        }
        else if (yr <= 2080.0)
        {
            double f = (yr - 2060.0) / (2080.0 - 2060.0);
            return rel2060 + f * (rel2080 - rel2060);
        }
        else
        {
            double f = (yr - 2080.0) / (2100.0 - 2080.0);
            return rel2080 + f * (rel2100 - rel2080);
        }
    }

    // Nearest threshold band (1..8) for d in {2, 3, 4, 5, 6, 7, 8, 10} cm
    inline int getThresholdBand1Based(double diameterThresholdCm, double* matchedThresholdCm = nullptr)
    {
        static const double thresholds[8] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
        int bestIdx = 0;
        double minDiff = std::abs(diameterThresholdCm - thresholds[0]);
        for (int i = 1; i < 8; ++i)
        {
            double diff = std::abs(diameterThresholdCm - thresholds[i]);
            if (diff < minDiff)
            {
                minDiff = diff;
                bestIdx = i;
            }
        }
        if (matchedThresholdCm)
        {
            *matchedThresholdCm = thresholds[bestIdx];
        }
        return bestIdx + 1;
    }

    // -------------------------------------------------------------------------
    // Hail Kinematics, Kinetic Energy, Impact Energy Density & Vulnerability
    // -------------------------------------------------------------------------

    // Spherical hailstone mass [kg] from diameter D [cm]
    inline double hailstoneMassKg(double diameterCm)
    {
        if (diameterCm <= 0.0 || !std::isfinite(diameterCm)) return 0.0;
        double d_m = diameterCm * 0.01;
        return (M_PI / 6.0) * kHailIceDensityKgM3 * (d_m * d_m * d_m);
    }

    // Vertical terminal fall velocity vt [m/s] as a function of diameter D [cm] and elevation z [m]
    // vt = sqrt((4 * rho_h * g * d_m) / (3 * Cd * rho_a(z)))
    inline double terminalFallVelocityMS(double diameterCm, double elevationM = 0.0)
    {
        if (diameterCm <= 0.0 || !std::isfinite(diameterCm)) return 0.0;
        double d_m = diameterCm * 0.01;
        double z = std::clamp(std::isfinite(elevationM) ? elevationM : 0.0, -500.0, 8000.0);
        double rho_air = kSeaLevelAirDensityKgM3 * std::exp(-z / kAirScaleHeightM);
        double vt = std::sqrt((4.0 * kHailIceDensityKgM3 * kGravityMS2 * d_m) /
                              (3.0 * kHailDragCoeff * std::max(0.2, rho_air)));
        return vt;
    }

    // Resultant 3D impact velocity [m/s] combining vertical fall speed vt and horizontal wind v_wind
    inline double resultantImpactVelocityMS(double vtMS, double windSpeedMS)
    {
        if (vtMS <= 0.0) return 0.0;
        double vw = std::max(0.0, std::isfinite(windSpeedMS) ? windSpeedMS : 0.0);
        double vh = kWindCouplingFactor * vw;
        return std::sqrt(vtMS * vtMS + vh * vh);
    }

    // Single-stone kinetic energy at impact [Joules]: Ek = 0.5 * m * v_impact^2
    inline double singleStoneKineticEnergyJ(double diameterCm, double vImpactMS)
    {
        if (diameterCm <= 0.0 || vImpactMS <= 0.0) return 0.0;
        double m = hailstoneMassKg(diameterCm);
        return 0.5 * m * vImpactMS * vImpactMS;
    }

    // Cumulative storm hail impact kinetic energy density [J/m^2] (Waldvogel 1978 / Hohl et al. 2002 hailpad proxy)
    // Uses an asymptotic saturation formulation respecting the exponential drop in giant hailstone
    // number concentration (Cheng & English 1983), bounding catastrophic supercells to ~5,500-6,000 J/m^2.
    inline double stormImpactEnergyDensityJM2(double diameterCm, double vtMS, double vImpactMS)
    {
        if (diameterCm <= 0.0 || vtMS <= 0.0 || vImpactMS <= 0.0) return 0.0;
        double windAmplification = (vImpactMS * vImpactMS) / std::max(1e-6, vtMS * vtMS);

        // Physical saturation model:
        // Matches Hohl et al. (2002) and Swiss/US hailpad observations for small/medium hail,
        // saturating asymptotically towards ~5,500 J/m^2 for giant hail (>8-10 cm)
        const double E_max_jm2 = 5500.0; // Maximum observed supercell storm kinetic energy density [J/m^2]
        const double d0_cm = 4.2;        // Characteristic diameter [cm] reaching ~63% of maximum
        double d_cm = std::max(0.0, diameterCm);
        double e_still_jm2 = E_max_jm2 * (1.0 - std::exp(-std::pow(d_cm / d0_cm, 2.2)));

        return e_still_jm2 * std::clamp(windAmplification, 1.0, 3.5);
    }

    // Effective hailstone diameter [cm] after accounting for wind amplification and protection reduction
    // Since Ek_still ~ D^4 (m ~ D^3 and vt^2 ~ D), D_eff = D * (Ek_effective / Ek_still)^(1/4)
    inline double effectiveHailDiameterCm(double rawDiameterCm, double ekEffectiveJ, double ekStillJ)
    {
        if (rawDiameterCm <= 0.0 || ekEffectiveJ <= 0.0 || ekStillJ <= 0.0) return 0.0;
        double ratio = std::max(0.0, ekEffectiveJ / ekStillJ);
        return std::min(kMaxHailDiameterCm, rawDiameterCm * std::pow(ratio, 0.25));
    }

    // TORRO Hailstorm Intensity Scale (H0 to H10) based on effective impact kinetic energy [J] / diameter [cm]
    inline float computeTorroScale(double effectiveDiameterCm, double ekEffectiveJ)
    {
        if (effectiveDiameterCm < 0.5 || ekEffectiveJ <= 0.05) return 0.0f;
        // TORRO scale relationship with kinetic energy / effective size
        // H0: ~0.5-1.0 cm, H1: 1.0-1.5 cm, H2: 1.5-2.0 cm, H3: 2.0-3.0 cm (~4-8 J),
        // H4: 3.0-4.0 cm (~10-20 J), H5: 4.0-5.0 cm (~20-40 J), H6: 5.0-6.5 cm (~40-80 J),
        // H7: 6.5-8.0 cm, H8: 8.0-10.0 cm, H9: 10.0-12.5 cm, H10: >12.5 cm
        double h = 0.0;
        if (effectiveDiameterCm <= 2.0)
        {
            h = (effectiveDiameterCm - 0.5) / 0.5;
        }
        else if (effectiveDiameterCm <= 5.0)
        {
            h = 3.0 + (effectiveDiameterCm - 2.0) * (2.0 / 3.0);
        }
        else
        {
            h = 5.0 + (effectiveDiameterCm - 5.0) * (5.0 / 7.5);
        }
        return static_cast<float>(std::clamp(h, 0.0, 10.0));
    }

    // Standardized roof/cladding hail damage ratio proxy [0.0 - 1.0] (logistic/Weibull fragility curve, Hohl et al. 2002 / FEMA)
    inline float computeHailDamageRatio(double ekEffectiveJ, double impactEnergyEffectiveJM2)
    {
        if (ekEffectiveJ <= 1.0 && impactEnergyEffectiveJM2 <= 20.0) return 0.0f;
        // Fragility onset around Ek ~ 5 J (D ~ 2.5 cm) / E_impact ~ 100 J/m^2,
        // 50% MDR around Ek ~ 35 J (D ~ 5.2 cm) / E_impact ~ 1500 J/m^2
        double x_ke = std::max(0.0, (ekEffectiveJ - 2.0) / 32.0);
        double d_ke = 1.0 - std::exp(-std::pow(x_ke, 1.65));

        double x_ie = std::max(0.0, (impactEnergyEffectiveJM2 - 30.0) / 1400.0);
        double d_ie = 1.0 - std::exp(-std::pow(x_ie, 1.45));

        double dmg = 0.6 * d_ke + 0.4 * d_ie;
        return static_cast<float>(std::clamp(dmg, 0.0, 1.0));
    }

    struct HailModelOutputs
    {
        std::shared_ptr<cTMapf> HAIL_SIZE_EFF = nullptr;    // Effective hailstone diameter [cm] (wind + protection adjusted)
        std::shared_ptr<cTMapf> HAIL_SIZE_RAW = nullptr;    // Meteorological GEV hailstone diameter [cm]
        std::shared_ptr<cTMapf> KINETIC_ENERGY = nullptr;   // Protected single-stone kinetic impact energy [J]
        std::shared_ptr<cTMapf> KINETIC_ENERGY_RAW = nullptr; // Unprotected single-stone kinetic impact energy [J]
        std::shared_ptr<cTMapf> IMPACT_ENERGY = nullptr;    // Protected cumulative storm impact energy density [J/m^2]
        std::shared_ptr<cTMapf> IMPACT_ENERGY_RAW = nullptr;// Unprotected cumulative storm impact energy density [J/m^2]
        std::shared_ptr<cTMapf> V_IMPACT = nullptr;         // Total 3D hailstone impact velocity [m/s]
        std::shared_ptr<cTMapf> V_TERMINAL = nullptr;       // Vertical terminal fall velocity [m/s]
        std::shared_ptr<cTMapf> WIND_LOCAL = nullptr;       // Local storm wind speed at impact height [m/s]
        std::shared_ptr<cTMapf> PROTECTION = nullptr;       // Combined adaptation & canopy protection reduction factor [0.0 - 1.0]
        std::shared_ptr<cTMapf> TORRO_SCALE = nullptr;      // TORRO Hailstorm Intensity Scale [0 - 10]
        std::shared_ptr<cTMapf> DAMAGE_RATIO = nullptr;     // Standardized roof/structure damage proxy [0.0 - 1.0]
    };

    inline int RunFastHail(
        const std::shared_ptr<cTMapf>& dem,
        const std::shared_ptr<cTMapf>& bh,
        const std::shared_ptr<cTMapf>& ch,
        const std::shared_ptr<cTMapf>& hail_size_cm,
        const std::shared_ptr<cTMapf>& wind_speed_ms,
        const std::shared_ptr<cTMapf>& adaptation_prot,
        const std::shared_ptr<cTMapf>& adaptation_thresh_j,
        double default_wind_ms,
        float canopy_interception_factor,
        HailModelOutputs& outputs)
    {
        if (!dem || !hail_size_cm)
        {
            return -1;
        }

        const int rows = dem->nrRows();
        const int cols = dem->nrCols();

        outputs.HAIL_SIZE_EFF = dem->GetCopySP0();
        outputs.HAIL_SIZE_RAW = dem->GetCopySP0();
        outputs.KINETIC_ENERGY = dem->GetCopySP0();
        outputs.KINETIC_ENERGY_RAW = dem->GetCopySP0();
        outputs.IMPACT_ENERGY = dem->GetCopySP0();
        outputs.IMPACT_ENERGY_RAW = dem->GetCopySP0();
        outputs.V_IMPACT = dem->GetCopySP0();
        outputs.V_TERMINAL = dem->GetCopySP0();
        outputs.WIND_LOCAL = dem->GetCopySP0();
        outputs.PROTECTION = dem->GetCopySP0();
        outputs.TORRO_SCALE = dem->GetCopySP0();
        outputs.DAMAGE_RATIO = dem->GetCopySP0();

        #pragma omp parallel for schedule(static)
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                float z_val = dem->data[r][c];
                float d_raw = (r < hail_size_cm->nrRows() && c < hail_size_cm->nrCols()) ? hail_size_cm->data[r][c] : 0.0f;

                if (pcr::isMV(z_val) || pcr::isMV(d_raw) || !std::isfinite(d_raw) || d_raw < 0.0f)
                {
                    d_raw = 0.0f;
                }
                if (!std::isfinite(z_val) || pcr::isMV(z_val))
                {
                    z_val = 0.0f;
                }

                float h_b = (bh && r < bh->nrRows() && c < bh->nrCols() && !pcr::isMV(bh->data[r][c]) && std::isfinite(bh->data[r][c]))
                            ? std::max(0.0f, bh->data[r][c]) : 0.0f;
                float h_c = (ch && r < ch->nrRows() && c < ch->nrCols() && !pcr::isMV(ch->data[r][c]) && std::isfinite(ch->data[r][c]))
                            ? std::max(0.0f, ch->data[r][c]) : 0.0f;

                // Base 10m storm wind speed [m/s]
                double v10 = default_wind_ms;
                if (wind_speed_ms && r < wind_speed_ms->nrRows() && c < wind_speed_ms->nrCols())
                {
                    float wv = wind_speed_ms->data[r][c];
                    if (!pcr::isMV(wv) && std::isfinite(wv) && wv >= 0.0f)
                    {
                        v10 = static_cast<double>(wv);
                    }
                }

                // Adjust wind speed for building roof height and canopy shelter
                double v_local = v10;
                if (h_b > 3.0f)
                {
                    // Logarithmic wind profile adjustment at roof height
                    const double z0 = 0.25;
                    double z_eff = std::clamp(static_cast<double>(h_b), 3.0, 80.0);
                    double log_factor = std::log((z_eff + z0) / z0) / std::log((10.0 + z0) / z0);
                    v_local *= std::clamp(log_factor, 0.85, 1.45);
                }

                // Natural tree canopy wind shelter & hail interception when canopy overtops surface/roof
                double canopy_prot = 0.0;
                if (h_c > h_b + 1.5f)
                {
                    double overtop = std::min(15.0, static_cast<double>(h_c - h_b));
                    double shelter_frac = 1.0 - std::exp(-0.15 * overtop);
                    v_local *= (1.0 - 0.5 * shelter_frac);
                    canopy_prot = std::clamp(static_cast<double>(canopy_interception_factor) * shelter_frac, 0.0, 0.75);
                }

                // User adaptation protection factor [0.0 - 1.0] and impact energy threshold [J]
                double user_prot = 0.0;
                if (adaptation_prot && r < adaptation_prot->nrRows() && c < adaptation_prot->nrCols())
                {
                    float ap = adaptation_prot->data[r][c];
                    if (!pcr::isMV(ap) && std::isfinite(ap) && ap > 0.0f)
                    {
                        user_prot = std::clamp(static_cast<double>(ap), 0.0, 1.0);
                    }
                }

                double thresh_j = 0.0;
                if (adaptation_thresh_j && r < adaptation_thresh_j->nrRows() && c < adaptation_thresh_j->nrCols())
                {
                    float tj = adaptation_thresh_j->data[r][c];
                    if (!pcr::isMV(tj) && std::isfinite(tj) && tj > 0.0f)
                    {
                        thresh_j = static_cast<double>(tj);
                    }
                }

                // Combine natural canopy protection and engineered roof/adaptation protection
                double combined_prot = 1.0 - (1.0 - user_prot) * (1.0 - canopy_prot);
                combined_prot = std::clamp(combined_prot, 0.0, 0.98);

                // Compute hailstone kinematics
                double vt = terminalFallVelocityMS(d_raw, z_val);
                double v_imp = resultantImpactVelocityMS(vt, v_local);

                double ek_still = singleStoneKineticEnergyJ(d_raw, vt);
                double ek_raw = singleStoneKineticEnergyJ(d_raw, v_imp);
                double ie_raw = stormImpactEnergyDensityJM2(d_raw, vt, v_imp);

                // Apply protection reduction and threshold absorption (e.g. UL 2218 Class 1-4 roof rating)
                double ek_eff = ek_raw * (1.0 - combined_prot);
                double ie_eff = ie_raw * (1.0 - combined_prot);

                if (thresh_j > 0.0 && ek_raw > 0.0)
                {
                    if (ek_raw <= thresh_j)
                    {
                        // Below the rated impact resistance threshold, 90% of residual energy is absorbed elastically
                        ek_eff *= 0.10;
                        ie_eff *= 0.10;
                    }
                    else
                    {
                        // Above threshold, subtract absorbed threshold energy before applying protection factor
                        double excess_ratio = std::max(0.10, (ek_raw - 0.75 * thresh_j) / ek_raw);
                        ek_eff = std::min(ek_eff, ek_raw * excess_ratio * (1.0 - 0.5 * user_prot));
                        ie_eff = std::min(ie_eff, ie_raw * excess_ratio * (1.0 - 0.5 * user_prot));
                    }
                    combined_prot = std::clamp(1.0 - (ek_eff / std::max(1e-9, ek_raw)), 0.0, 0.99);
                }

                double d_eff = effectiveHailDiameterCm(d_raw, ek_eff, ek_still);
                float torro = computeTorroScale(d_eff, ek_eff);
                float dmg_ratio = computeHailDamageRatio(ek_eff, ie_eff);

                outputs.HAIL_SIZE_RAW->data[r][c] = static_cast<float>(d_raw);
                outputs.HAIL_SIZE_EFF->data[r][c] = static_cast<float>(d_eff);
                outputs.KINETIC_ENERGY_RAW->data[r][c] = static_cast<float>(ek_raw);
                outputs.KINETIC_ENERGY->data[r][c] = static_cast<float>(ek_eff);
                outputs.IMPACT_ENERGY_RAW->data[r][c] = static_cast<float>(ie_raw);
                outputs.IMPACT_ENERGY->data[r][c] = static_cast<float>(ie_eff);
                outputs.V_TERMINAL->data[r][c] = static_cast<float>(vt);
                outputs.V_IMPACT->data[r][c] = static_cast<float>(v_imp);
                outputs.WIND_LOCAL->data[r][c] = static_cast<float>(v_local);
                outputs.PROTECTION->data[r][c] = static_cast<float>(combined_prot);
                outputs.TORRO_SCALE->data[r][c] = torro;
                outputs.DAMAGE_RATIO->data[r][c] = dmg_ratio;
            }
        }

        return 0;
    }

} // namespace FastHailModel

#endif // FASTHAIL_H

