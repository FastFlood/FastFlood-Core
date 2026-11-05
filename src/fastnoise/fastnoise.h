//
// Created for FastNoise: Traffic and Environmental Noise Pollution Model
// Based on Reken- en meetvoorschrift geluid 2012, Bijlage III (SRM 1 and SRM 2)
// and Tygron Traffic Noise Overlay standard methodology.
//

#ifndef FASTNOISE_H
#define FASTNOISE_H

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cstdint>

#include "../../common/defines.h"
#include "../../common/fferror.h"
#include "../geospatial/map/map.h"
#include "../../math/fastmath.h"

#if FF_HAS_OMP
#include "omp.h"
#endif

#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace FastNoiseModel
{
    // Fast precomputed cosine taper lookup table for angular beam wedge spreading
    // Replaces expensive transcendental calls with a cache-friendly L1 table
    struct CosWedgeLUT
    {
        float table[1025];
        CosWedgeLUT()
        {
            for (int i = 0; i <= 1024; ++i)
            {
                double xi = i / 1024.0;
                table[i] = static_cast<float>(0.5 * (1.0 + std::cos(M_PI * 0.5 * xi)));
            }
        }
        inline float get(double xi) const
        {
            if (xi <= 0.0) return 1.0f;
            if (xi >= 1.0) return 0.5f;
            double pos = xi * 1024.0;
            int idx = static_cast<int>(pos);
            float frac = static_cast<float>(pos - idx);
            return table[idx] * (1.0f - frac) + table[idx + 1] * frac;
        }
    };
    inline const CosWedgeLUT g_cosWedgeLUT;

    // Fast atomic float maximum update helper (thread-safe without locking overhead)
    inline void atomicMaxFloat(float* addr, float val)
    {
        float cur = *addr;
        if (val <= cur) return;

        #if defined(__GNUC__) || defined(__clang__)
        uint32_t val_bits;
        std::memcpy(&val_bits, &val, sizeof(float));
        uint32_t cur_bits;
        std::memcpy(&cur_bits, &cur, sizeof(float));
        while (val > cur)
        {
            if (__atomic_compare_exchange_n(reinterpret_cast<uint32_t*>(addr), &cur_bits, val_bits, true, __ATOMIC_RELAXED, __ATOMIC_RELAXED))
                break;
            std::memcpy(&cur, &cur_bits, sizeof(float));
        }
        #elif defined(_MSC_VER)
        uint32_t val_bits;
        std::memcpy(&val_bits, &val, sizeof(float));
        uint32_t cur_bits;
        std::memcpy(&cur_bits, &cur, sizeof(float));
        while (val > cur)
        {
            uint32_t prev = _InterlockedCompareExchange(reinterpret_cast<volatile long*>(addr), static_cast<long>(val_bits), static_cast<long>(cur_bits));
            if (prev == cur_bits) break;
            cur_bits = prev;
            std::memcpy(&cur, &cur_bits, sizeof(float));
        }
        #else
        #pragma omp critical
        {
            if (val > *addr) *addr = val;
        }
        #endif
    }

    // Output rasters structure
    struct NoiseModelOutputs
    {
        std::shared_ptr<cTMapf> LDEN = nullptr;        // Day-Evening-Night level in dB(A)
        std::shared_ptr<cTMapf> LDAY = nullptr;        // Day level (07:00-19:00) in dB(A)
        std::shared_ptr<cTMapf> LEVENING = nullptr;    // Evening level (19:00-23:00) in dB(A)
        std::shared_ptr<cTMapf> LNIGHT = nullptr;      // Night level (23:00-07:00) in dB(A)
        std::shared_ptr<cTMapf> SOURCES = nullptr;     // Rasterized source emission in dB(A)
        std::shared_ptr<cTMapf> KRUISPUNTEN = nullptr; // Intersection acceleration correction in dB(A)
        std::shared_ptr<cTMapf> AIRPORTS = nullptr;    // Rasterized airport noise emission in dB(A)
        std::shared_ptr<cTMapf> INDUSTRY = nullptr;    // Rasterized industrial/port noise emission in dB(A)
        std::shared_ptr<cTMapf> BARRIERS = nullptr;    // Noise screen/barrier heights in m
        std::shared_ptr<cTMapf> REFLECTIONS = nullptr; // Facade reflection boost in dB(A)
    };

    // -------------------------------------------------------------------------
    // RMG 2012 Bijlage III Road Noise Emission Calculation
    // -------------------------------------------------------------------------
    inline double calculateRoadEmissionCategory(
        double alpha, double beta, double v0,
        double speedKmH, double hourlyIntensity)
    {
        if (hourlyIntensity <= 0.0) return -999.0;
        double v = std::clamp(speedKmH, 30.0, 160.0);
        double e = alpha + beta * std::log10(v / v0) + 10.0 * std::log10(hourlyIntensity);
        return e;
    }

    // Full road emission according to RMG 2012 (light, medium, heavy vehicles)
    inline double calculateRMG2012RoadEmission(
        double speedKmH,
        double Q_lv, double Q_mv, double Q_zv,
        double c_wegdek = 0.0)
    {
        // Reference constants for DAB (Dicht Asfaltbeton):
        // lv: alpha=63.7, beta=28.0, v0=80 km/h
        // mv: alpha=70.3, beta=20.0, v0=70 km/h
        // zv: alpha=73.5, beta=18.0, v0=70 km/h
        double e_lv = calculateRoadEmissionCategory(63.7, 28.0, 80.0, speedKmH, Q_lv);
        double e_mv = calculateRoadEmissionCategory(70.3, 20.0, 70.0, speedKmH, Q_mv);
        double e_zv = calculateRoadEmissionCategory(73.5, 18.0, 70.0, speedKmH, Q_zv);

        double sum_lin = 0.0;
        if (e_lv > -900.0) sum_lin += std::pow(10.0, e_lv / 10.0);
        if (e_mv > -900.0) sum_lin += std::pow(10.0, e_mv / 10.0);
        if (e_zv > -900.0) sum_lin += std::pow(10.0, e_zv / 10.0);

        if (sum_lin <= 0.0) return -999.0;
        return 10.0 * std::log10(sum_lin) + c_wegdek;
    }

    // -------------------------------------------------------------------------
    // Kruispunt (Intersection) Detection and Optrekcorrectie Calculation
    // RMG 2012 Bijlage III Section 1.6.1 & Section 2.5.1
    // -------------------------------------------------------------------------
    inline std::shared_ptr<cTMapf> detectKruispunten(
        const std::shared_ptr<cTMapf>& roadsMap,
        float dx, float dy,
        float heavyVehiclePct = 10.0f)
    {
        if (!roadsMap) return nullptr;
        int rows = roadsMap->nrRows();
        int cols = roadsMap->nrCols();
        auto krpMap = std::shared_ptr<cTMapf>(roadsMap->GetCopy0());

        // 1. Identify intersection nodes on the road network
        // Using Rutovitz crossing number and connectivity analysis:
        // A road intersection is where 3 or more distinct branches meet.
        std::vector<std::pair<int, int>> rawJunctions;

        // 8-neighborhood offsets in clockwise cyclic order: N, NE, E, SE, S, SW, W, NW
        const int dRow[8] = {-1, -1,  0,  1,  1,  1,  0, -1};
        const int dCol[8] = { 0,  1,  1,  1,  0, -1, -1, -1};

        for (int r = 1; r < rows - 1; ++r)
        {
            for (int c = 1; c < cols - 1; ++c)
            {
                float val = roadsMap->data[r][c];
                if (pcr::isMV(val) || val <= 0.0f) continue;

                int roadNeighbors = 0;
                int nArr[8];
                for (int i = 0; i < 8; ++i)
                {
                    float nVal = roadsMap->data[r + dRow[i]][c + dCol[i]];
                    nArr[i] = (!pcr::isMV(nVal) && nVal > 0.0f) ? 1 : 0;
                    if (nArr[i]) roadNeighbors++;
                }

                if (roadNeighbors >= 3)
                {
                    // Calculate crossing transitions 0 -> 1 in cycle
                    int transitions = 0;
                    for (int i = 0; i < 8; ++i)
                    {
                        int next_i = (i + 1) % 8;
                        if (nArr[i] == 0 && nArr[next_i] == 1)
                        {
                            transitions++;
                        }
                    }

                    // True branch node has >= 3 distinct branch transitions
                    // or 4+ neighbors with 2+ transitions
                    if (transitions >= 3 || (roadNeighbors >= 4 && transitions >= 2))
                    {
                        rawJunctions.emplace_back(r, c);
                    }
                }
            }
        }

        // 2. Cluster / deduplicate junction nodes within cluster distance
        // (to condense thick lines into singular intersection points)
        std::vector<std::pair<int, int>> junctionNodes;
        int clusterDistPix = std::max(2, static_cast<int>(std::round(15.0 / std::min(dx, dy))));
        int clusterDist2 = clusterDistPix * clusterDistPix;

        std::vector<bool> merged(rawJunctions.size(), false);
        for (size_t i = 0; i < rawJunctions.size(); ++i)
        {
            if (merged[i]) continue;
            long long sumR = rawJunctions[i].first;
            long long sumC = rawJunctions[i].second;
            int count = 1;

            for (size_t j = i + 1; j < rawJunctions.size(); ++j)
            {
                if (merged[j]) continue;
                int dr = rawJunctions[i].first - rawJunctions[j].first;
                int dc = rawJunctions[i].second - rawJunctions[j].second;
                if (dr * dr + dc * dc <= clusterDist2)
                {
                    sumR += rawJunctions[j].first;
                    sumC += rawJunctions[j].second;
                    count++;
                    merged[j] = true;
                }
            }
            junctionNodes.emplace_back(static_cast<int>(sumR / count), static_cast<int>(sumC / count));
        }

        FF::coutI << "[FastNoise] Kruispunt detection identified " << junctionNodes.size()
                  << " distinct intersections on the road network." << std::endl;

        if (junctionNodes.empty())
        {
            return krpMap;
        }

        // 3. RMG 2012 Bijlage III optrekcorrectie:
        // C_kruispunt(a) = deltaL_max * (1.0 - a / 150.0) for a <= 150 m
        // deltaL_max increases with percentage of medium/heavy vehicles (p)
        const double deltaL_max = std::clamp(1.5 + 2.5 * (heavyVehiclePct / 100.0), 1.5, 4.0);
        const double maxDist = 150.0; // meters per RMG 2012 Bijlage III
        int radiusR = static_cast<int>(std::ceil(maxDist / dy));
        int radiusC = static_cast<int>(std::ceil(maxDist / dx));

        // Parallel stamp of kruispunt boost directly onto nearby road cells
        #pragma omp parallel for schedule(dynamic, 8)
        for (size_t j = 0; j < junctionNodes.size(); ++j)
        {
            int jR = junctionNodes[j].first;
            int jC = junctionNodes[j].second;

            int minR = std::max(0, jR - radiusR);
            int maxR = std::min(rows - 1, jR + radiusR);
            int minC = std::max(0, jC - radiusC);
            int maxC = std::min(cols - 1, jC + radiusC);

            for (int r = minR; r <= maxR; ++r)
            {
                double dY = (r - jR) * dy;
                for (int c = minC; c <= maxC; ++c)
                {
                    float rVal = roadsMap->data[r][c];
                    if (pcr::isMV(rVal) || rVal <= 0.0f) continue;

                    double dX = (c - jC) * dx;
                    double dist = std::sqrt(dX * dX + dY * dY);
                    if (dist <= maxDist)
                    {
                        float boost = static_cast<float>(deltaL_max * (1.0 - dist / maxDist));
                        #pragma omp atomic update
                        krpMap->data[r][c] = (krpMap->data[r][c] > boost) ? krpMap->data[r][c] : boost;
                    }
                }
            }
        }

        return krpMap;
    }

    // -------------------------------------------------------------------------
    // Facade Reflection Term C_reflectie (RMG 2012 Bijlage III Section 1.7)
    // -------------------------------------------------------------------------
    inline std::shared_ptr<cTMapf> calculateFacadeReflections(
        const std::shared_ptr<cTMapf>& roadsMap,
        const std::shared_ptr<cTMapf>& bhMap,
        float dx, float dy)
    {
        if (!roadsMap) return nullptr;
        int rows = roadsMap->nrRows();
        int cols = roadsMap->nrCols();
        auto reflMap = std::shared_ptr<cTMapf>(roadsMap->GetCopy0());
        if (!bhMap) return reflMap;

        // Street canyon detection: check if buildings flank roads within ~25m
        int searchRadiusCells = std::max(1, static_cast<int>(std::round(25.0 / std::max(dx, dy))));

        #pragma omp parallel for schedule(dynamic, 16)
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                float rVal = roadsMap->data[r][c];
                if (pcr::isMV(rVal) || rVal <= 0.0f)
                {
                    reflMap->data[r][c] = 0.0f;
                    continue;
                }

                int bldgCountEast = 0, bldgCountWest = 0;
                int bldgCountNorth = 0, bldgCountSouth = 0;

                for (int step = 1; step <= searchRadiusCells; ++step)
                {
                    if (c + step < cols && !pcr::isMV(bhMap->data[r][c + step]) && bhMap->data[r][c + step] > 2.5f) bldgCountEast++;
                    if (c - step >= 0 && !pcr::isMV(bhMap->data[r][c - step]) && bhMap->data[r][c - step] > 2.5f) bldgCountWest++;
                    if (r + step < rows && !pcr::isMV(bhMap->data[r + step][c]) && bhMap->data[r + step][c] > 2.5f) bldgCountSouth++;
                    if (r - step >= 0 && !pcr::isMV(bhMap->data[r - step][c]) && bhMap->data[r - step][c] > 2.5f) bldgCountNorth++;
                }

                // If street is flanked by buildings on opposite sides
                bool flankedX = (bldgCountEast > 0 && bldgCountWest > 0);
                bool flankedY = (bldgCountNorth > 0 && bldgCountSouth > 0);

                if (flankedX || flankedY)
                {
                    // C_reflectie = 1.5 * f_obj (up to 1.5 - 2.0 dB)
                    reflMap->data[r][c] = 1.5f;
                }
                else if (bldgCountEast > 0 || bldgCountWest > 0 || bldgCountNorth > 0 || bldgCountSouth > 0)
                {
                    reflMap->data[r][c] = 0.8f;
                }
                else
                {
                    reflMap->data[r][c] = 0.0f;
                }
            }
        }

        return reflMap;
    }

    // -------------------------------------------------------------------------
    // Main Simulation Execution: RunFastNoise
    // Implements RMG 2012 Bijlage III and Tygron Traffic Noise Overlay
    // -------------------------------------------------------------------------
    inline int RunFastNoise(
        const std::shared_ptr<cTMapf>& DEM,
        const std::shared_ptr<cTMapf>& BH,
        const std::shared_ptr<cTMapf>& CH,
        const std::shared_ptr<cTMapf>& ROADS,
        const std::shared_ptr<cTMapf>& TRANSIT,
        const std::shared_ptr<cTMapf>& COMMERCIAL,
        const std::shared_ptr<cTMapf>& BARRIERS,
        const std::shared_ptr<cTMapf>& GROUND_ABS,
        double defaultTrafficFlow,
        double defaultTrafficSpeed,
        double tempC,
        double rhPct,
        double windSpeed,
        double windDirDeg,
        double receiverHeight,
        double backgroundNoiseDBA,
        double maxDistMeters = 2500.0,
        const std::string& roadSurfaceType = "dab",
        bool enableKruispuntDetection = true,
        NoiseModelOutputs& outputs = *(new NoiseModelOutputs()),
        const std::shared_ptr<cTMapf>& AIRPORTS = nullptr,
        const std::shared_ptr<cTMapf>& INDUSTRY = nullptr,
        const std::shared_ptr<cTMapf>& ROAD_SPEED = nullptr,
        double heavyVehiclePct = 10.0,
        int numBounces = 2)
    {
        if (!DEM)
        {
            FF_ERROR("Cannot execute FastNoise without DEM");
            return -1;
        }

        if (maxDistMeters <= 0.0) maxDistMeters = 2500.0;

        int rows = DEM->nrRows();
        int cols = DEM->nrCols();
        double dx = std::abs(DEM->cellSizeX());
        double dy = std::abs(DEM->cellSizeY());
        if (dx <= 0.0) dx = std::abs(DEM->cellSize());
        if (dy <= 0.0) dy = std::abs(DEM->cellSize());
        if (dx <= 0.0) dx = 2.5;
        if (dy <= 0.0) dy = 2.5;

        FF::coutI << "=================================================================" << std::endl;
        FF::coutI << "[FastNoise] Starting noise simulation per RMG 2012 Bijlage III / Tygron ("
                  << rows << "x" << cols << " cells, res=" << dx << "x" << dy << "m)..." << std::endl;
        FF::coutI << "            Default Traffic: " << defaultTrafficFlow << " veh/day @ "
                  << defaultTrafficSpeed << " km/h (" << roadSurfaceType << ")" << std::endl;
        FF::coutI << "            Atmosphere: " << tempC << " C, " << rhPct << "% RH, Wind: "
                  << windSpeed << " m/s (" << windDirDeg << " deg)" << std::endl;
        FF::coutI << "            Assessment Height: " << receiverHeight << " m above ground" << std::endl;
        FF::coutI << "            Max Calculation Radius: " << maxDistMeters << " m" << std::endl;
        FF::coutI << "            Facade Reflection Bounces: " << numBounces << " (Specular + 22.5 deg Angular Wedges)" << std::endl;
        FF::coutI << "=================================================================" << std::endl;

        // 1. Road surface correction C_wegdek
        double c_wegdek = 0.0;
        if (roadSurfaceType == "zoab" || roadSurfaceType == "silent") c_wegdek = -2.5;
        else if (roadSurfaceType == "klinkers" || roadSurfaceType == "bricks") c_wegdek = 3.0;

        // 2. Prepare Digital Surface Model (DSM = DEM + BH + Barrier)
        auto dsmMap = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                float z = DEM->data[r][c];
                if (pcr::isMV(z))
                {
                    pcr::setMV(dsmMap->data[r][c]);
                    continue;
                }
                float bh = (BH && !pcr::isMV(BH->data[r][c])) ? std::max(0.0f, BH->data[r][c]) : 0.0f;
                float bar = (BARRIERS && !pcr::isMV(BARRIERS->data[r][c])) ? std::max(0.0f, BARRIERS->data[r][c]) : 0.0f;
                dsmMap->data[r][c] = z + bh + bar;
            }
        }

        // 3. Kruispunt detection & correction
        std::shared_ptr<cTMapf> kruispuntMap = nullptr;
        if (enableKruispuntDetection && ROADS)
        {
            kruispuntMap = detectKruispunten(ROADS, static_cast<float>(dx), static_cast<float>(dy), static_cast<float>(heavyVehiclePct));
        }
        else
        {
            kruispuntMap = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        }
        outputs.KRUISPUNTEN = kruispuntMap;

        // 4. Facade reflection
        auto reflMap = calculateFacadeReflections(ROADS, BH, static_cast<float>(dx), static_cast<float>(dy));
        outputs.REFLECTIONS = reflMap;

        // 5. Structure to hold active sound emission sources on the grid
        auto sourcesRasterDay = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        auto sourcesRasterEv  = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        auto sourcesRasterNi  = std::shared_ptr<cTMapf>(DEM->GetCopy0());

        // Derive time-of-day intensities from 24h total:
        // Day (12h): 6.6% per hour
        // Evening (4h): 3.3% per hour
        // Night (8h): 1.0% per hour
        double q_day = defaultTrafficFlow * 0.066;
        double q_ev = defaultTrafficFlow * 0.033;
        double q_ni = defaultTrafficFlow * 0.010;

        int totalSourceCells = 0;

        // Scan all source maps (ROADS, TRANSIT, COMMERCIAL, INDUSTRY, AIRPORTS)
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c])) continue;

                double e_d = -999.0;
                double e_e = -999.0;
                double e_n = -999.0;

                // Roads
                if (ROADS && !pcr::isMV(ROADS->data[r][c]) && ROADS->data[r][c] > 0.0f)
                {
                    float factor = ROADS->data[r][c];
                    // If ROADS stores custom intensity or multiplier
                    double curTrafficFlow = (factor > 10.0f) ? static_cast<double>(factor) : (defaultTrafficFlow * factor);
                    double curQDay = curTrafficFlow * 0.066;
                    double curQEv  = curTrafficFlow * 0.033;
                    double curQNi  = curTrafficFlow * 0.010;

                    // Road speed from ROAD_SPEED or default reference speed
                    double curSpeed = (ROAD_SPEED && !pcr::isMV(ROAD_SPEED->data[r][c]) && ROAD_SPEED->data[r][c] > 0.0f)
                                      ? static_cast<double>(ROAD_SPEED->data[r][c])
                                      : defaultTrafficSpeed;

                    // Dynamic vehicle composition per Dutch RMG 2012 / CROW:
                    // Motorways and high-speed trunks carry a higher fraction of heavy freight (zv = 10%)
                    // Residential and low-speed urban streets carry mostly light passenger vehicles (zv = 1%)
                    double p_lv = 0.88, p_mv = 0.08, p_zv = 0.04;
                    if (curSpeed >= 90.0)
                    {
                        p_lv = 0.80; p_mv = 0.10; p_zv = 0.10; // Motorway / Autoweg freight corridor
                    }
                    else if (curSpeed <= 40.0)
                    {
                        p_lv = 0.96; p_mv = 0.03; p_zv = 0.01; // Woonstraat / 30 km/h zone
                    }

                    double ed = calculateRMG2012RoadEmission(curSpeed, curQDay * p_lv, curQDay * p_mv, curQDay * p_zv, c_wegdek);
                    double ee = calculateRMG2012RoadEmission(curSpeed, curQEv * std::clamp(p_lv + 0.04, 0.0, 1.0), curQEv * std::max(0.01, p_mv - 0.02), curQEv * std::max(0.01, p_zv - 0.02), c_wegdek);
                    double en = calculateRMG2012RoadEmission(curSpeed, curQNi * std::max(0.50, p_lv - 0.06), curQNi * (p_mv + 0.03), curQNi * (p_zv + 0.03), c_wegdek);

                    // Add Kruispunt boost & facade reflection
                    double c_krp = (kruispuntMap && !pcr::isMV(kruispuntMap->data[r][c])) ? kruispuntMap->data[r][c] : 0.0;
                    double c_ref = (reflMap && !pcr::isMV(reflMap->data[r][c])) ? reflMap->data[r][c] : 0.0;

                    ed += (c_krp + c_ref);
                    ee += (c_krp + c_ref);
                    en += (c_krp + c_ref);

                    e_d = ed;
                    e_e = ee;
                    e_n = en;
                }

                // Transit / Rail
                if (TRANSIT && !pcr::isMV(TRANSIT->data[r][c]) && TRANSIT->data[r][c] > 0.0f)
                {
                    float rf = TRANSIT->data[r][c];
                    double e_rail_d = 72.0;
                    double e_rail_e = 68.0;
                    double e_rail_n = 64.0;
                    if (rf > 30.0f)
                    {
                        // Direct custom dB level in raster
                        e_rail_d = rf;
                        e_rail_e = rf - 4.0;
                        e_rail_n = rf - 8.0;
                    }
                    else if (rf >= 1.5f)
                    {
                        // High-speed / heavy rail corridor
                        e_rail_d = 78.0;
                        e_rail_e = 74.0;
                        e_rail_n = 70.0;
                    }
                    else if (rf <= 0.6f)
                    {
                        // Tram / light rail
                        e_rail_d = 64.0;
                        e_rail_e = 60.0;
                        e_rail_n = 54.0;
                    }

                    e_d = (e_d > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_d / 10.0) + std::pow(10.0, e_rail_d / 10.0))) : e_rail_d;
                    e_e = (e_e > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_e / 10.0) + std::pow(10.0, e_rail_e / 10.0))) : e_rail_e;
                    e_n = (e_n > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_n / 10.0) + std::pow(10.0, e_rail_n / 10.0))) : e_rail_n;
                }

                // Commercial / Retail
                if (COMMERCIAL && !pcr::isMV(COMMERCIAL->data[r][c]) && COMMERCIAL->data[r][c] > 0.0f)
                {
                    double e_comm_d = 62.0;
                    double e_comm_e = 58.0;
                    double e_comm_n = 52.0;

                    e_d = (e_d > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_d / 10.0) + std::pow(10.0, e_comm_d / 10.0))) : e_comm_d;
                    e_e = (e_e > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_e / 10.0) + std::pow(10.0, e_comm_e / 10.0))) : e_comm_e;
                    e_n = (e_n > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_n / 10.0) + std::pow(10.0, e_comm_n / 10.0))) : e_comm_n;
                }

                // Industrial / Port Logistics
                if (INDUSTRY && !pcr::isMV(INDUSTRY->data[r][c]) && INDUSTRY->data[r][c] > 0.0f)
                {
                    float indFactor = INDUSTRY->data[r][c];
                    double e_ind_d = (indFactor > 1.5f) ? 76.0 : 72.0;
                    double e_ind_e = (indFactor > 1.5f) ? 72.0 : 68.0;
                    double e_ind_n = (indFactor > 1.5f) ? 66.0 : 62.0;

                    e_d = (e_d > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_d / 10.0) + std::pow(10.0, e_ind_d / 10.0))) : e_ind_d;
                    e_e = (e_e > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_e / 10.0) + std::pow(10.0, e_ind_e / 10.0))) : e_ind_e;
                    e_n = (e_n > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_n / 10.0) + std::pow(10.0, e_ind_n / 10.0))) : e_ind_n;
                }

                // Airports / Aviation Ground & Movements
                if (AIRPORTS && !pcr::isMV(AIRPORTS->data[r][c]) && AIRPORTS->data[r][c] > 0.0f)
                {
                    float aVal = AIRPORTS->data[r][c];
                    double e_air_d = 70.0;
                    double e_air_e = 66.0;
                    double e_air_n = 60.0;

                    if (aVal >= 90.0f) // Runway
                    {
                        e_air_d = 94.0;
                        e_air_e = 90.0;
                        e_air_n = 84.0;
                    }
                    else if (aVal >= 75.0f) // Taxiway / Helipad
                    {
                        e_air_d = 80.0;
                        e_air_e = 76.0;
                        e_air_n = 70.0;
                    }
                    else if (aVal >= 60.0f) // Apron / Terminal Gate
                    {
                        e_air_d = 74.0;
                        e_air_e = 70.0;
                        e_air_n = 64.0;
                    }
                    else // General airfield / aerodrome boundary
                    {
                        e_air_d = 68.0;
                        e_air_e = 64.0;
                        e_air_n = 58.0;
                    }

                    e_d = (e_d > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_d / 10.0) + std::pow(10.0, e_air_d / 10.0))) : e_air_d;
                    e_e = (e_e > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_e / 10.0) + std::pow(10.0, e_air_e / 10.0))) : e_air_e;
                    e_n = (e_n > -900.0) ? (10.0 * std::log10(std::pow(10.0, e_n / 10.0) + std::pow(10.0, e_air_n / 10.0))) : e_air_n;
                }

                if (e_d > -900.0)
                {
                    sourcesRasterDay->data[r][c] = static_cast<float>(e_d);
                    sourcesRasterEv->data[r][c]  = static_cast<float>(e_e);
                    sourcesRasterNi->data[r][c]  = static_cast<float>(e_n);
                    totalSourceCells++;
                }
            }
        }

        outputs.SOURCES = sourcesRasterDay;
        if (AIRPORTS) outputs.AIRPORTS = AIRPORTS;
        if (INDUSTRY) outputs.INDUSTRY = INDUSTRY;
        if (BARRIERS) outputs.BARRIERS = BARRIERS;

        FF::coutI << "[FastNoise] Total noise source emission cells in domain: " << totalSourceCells << std::endl;

        // Allocate output sound level rasters
        auto mapLDay = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        auto mapLEvening = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        auto mapLNight = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        auto mapLDen = std::shared_ptr<cTMapf>(DEM->GetCopy0());

        if (totalSourceCells == 0)
        {
            FF::coutI << "[FastNoise] Warning: No noise sources detected. Outputting background level." << std::endl;
            for (int r = 0; r < rows; ++r)
            {
                for (int c = 0; c < cols; ++c)
                {
                    if (!pcr::isMV(DEM->data[r][c]))
                    {
                        mapLDay->data[r][c] = static_cast<float>(backgroundNoiseDBA);
                        mapLEvening->data[r][c] = static_cast<float>(backgroundNoiseDBA);
                        mapLNight->data[r][c] = static_cast<float>(backgroundNoiseDBA);
                        mapLDen->data[r][c] = static_cast<float>(backgroundNoiseDBA);
                    }
                }
            }
            outputs.LDAY = mapLDay;
            outputs.LEVENING = mapLEvening;
            outputs.LNIGHT = mapLNight;
            outputs.LDEN = mapLDen;
            return 0;
        }

        // Air absorption coefficient (ISO 9613-1 / RMG 2012)
        // Standard ~0.005 dB/m for 500-1000 Hz at 20 C, 70% RH
        double alphaAir = 0.005 * (1.0 + 0.02 * (20.0 - tempC) + 0.01 * (70.0 - rhPct) / 10.0);
        alphaAir = std::clamp(alphaAir, 0.002, 0.015);

        // Allocate linear acoustic energy accumulation rasters
        auto accumDay = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        auto accumEv  = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        auto accumNi  = std::shared_ptr<cTMapf>(DEM->GetCopy0());

        // Initialize source cells with their own emission energy
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                float sD = sourcesRasterDay->data[r][c];
                if (sD > 0.0f)
                {
                    float sE = sourcesRasterEv->data[r][c];
                    float sN = sourcesRasterNi->data[r][c];
                    accumDay->data[r][c] = static_cast<float>(std::pow(10.0, sD / 10.0));
                    accumEv->data[r][c]  = static_cast<float>(std::pow(10.0, sE / 10.0));
                    accumNi->data[r][c]  = static_cast<float>(std::pow(10.0, sN / 10.0));
                }
            }
        }

        // 6. Fast Multi-Directional Grid Sweeping with Wedges & Facade Reflections
        struct SweepDir
        {
            std::string name;
            int dr, dc;
            int pr, pc; // Orthogonal step vector for lateral angular wedge spreading
            double wedgeFactor; // tan(half_angle) for lateral angular wedge spreading
        };

        const std::vector<SweepDir> sweepDirs = {
            // Cardinal directions (0, 90, 180, 270 deg) - angular gap to neighbor on either side is 26.565 deg.
            // Using wedgeFactor = 0.30 (tan ~16.7 deg) ensures seamless fill across the 26.565 deg gap to adjacent knight directions.
            {"East (0 deg)",                    0,  1,   1,  0, 0.30},
            {"West (180 deg)",                  0, -1,   1,  0, 0.30},
            {"South (90 deg)",                  1,  0,   0,  1, 0.30},
            {"North (270 deg)",                -1,  0,   0,  1, 0.30},

            // Diagonal directions (45, 135, 225, 315 deg) - angular gap to neighbor is 18.435 deg.
            // Using wedgeFactor = 0.24 (tan ~13.5 deg) ensures generous overlap with knight directions.
            {"South-East (45 deg)",             1,  1,  -1,  1, 0.24},
            {"South-West (135 deg)",            1, -1,   1,  1, 0.24},
            {"North-East (315 deg)",           -1,  1,   1,  1, 0.24},
            {"North-West (225 deg)",           -1, -1,  -1,  1, 0.24},

            // Knight directions (26.6, 63.4, 116.6, 153.4, 206.6, 243.4, 296.6, 333.4 deg)
            // Neighbor on one side is cardinal (gap 26.565 deg), on other side is diagonal (gap 18.435 deg).
            // Using wedgeFactor = 0.28 (tan ~15.6 deg) reaches seamlessly into cardinal and diagonal wedges.
            {"East-South-East (26.6 deg)",      1,  2,  -2,  1, 0.28},
            {"East-North-East (333.4 deg)",    -1,  2,   2,  1, 0.28},
            {"West-South-West (153.4 deg)",     1, -2,   2,  1, 0.28},
            {"West-North-West (206.6 deg)",    -1, -2,  -2,  1, 0.28},
            {"South-South-East (63.4 deg)",     2,  1,  -1,  2, 0.28},
            {"South-South-West (116.6 deg)",    2, -1,   1,  2, 0.28},
            {"North-North-East (296.6 deg)",   -2,  1,   1,  2, 0.28},
            {"North-North-West (243.4 deg)",   -2, -1,  -1,  2, 0.28}
        };

        std::vector<double> sweepDs(sweepDirs.size());
        std::vector<double> sweepDp(sweepDirs.size());
        for (size_t i = 0; i < sweepDirs.size(); ++i)
        {
            sweepDs[i] = std::sqrt((sweepDirs[i].dr * dy) * (sweepDirs[i].dr * dy) + (sweepDirs[i].dc * dx) * (sweepDirs[i].dc * dx));
            sweepDp[i] = std::sqrt((sweepDirs[i].pr * dy) * (sweepDirs[i].pr * dy) + (sweepDirs[i].pc * dx) * (sweepDirs[i].pc * dx));
        }

        // Precompute disjoint parallel start cells along the boundary for all 16 directions
        std::vector<std::vector<std::pair<int, int>>> allStartCells(sweepDirs.size());
        for (size_t dIdx = 0; dIdx < sweepDirs.size(); ++dIdx)
        {
            int dr = sweepDirs[dIdx].dr;
            int dc = sweepDirs[dIdx].dc;
            for (int r = 0; r < rows; ++r)
            {
                for (int c = 0; c < cols; ++c)
                {
                    int prevR = r - dr;
                    int prevC = c - dc;
                    if (prevR < 0 || prevR >= rows || prevC < 0 || prevC >= cols)
                    {
                        allStartCells[dIdx].emplace_back(r, c);
                    }
                }
            }
        }

        struct SourceAccessor
        {
            const cTMapf* mapD = nullptr;
            const cTMapf* mapE = nullptr;
            const cTMapf* mapN = nullptr;
            const float* flatD = nullptr;
            const float* flatE = nullptr;
            const float* flatN = nullptr;
            int cols = 0;

            inline float getD(int r, int c) const {
                if (mapD) return mapD->data[r][c];
                if (flatD) return flatD[r * cols + c];
                return -999.0f;
            }
            inline float getE(int r, int c) const {
                if (mapE) return mapE->data[r][c];
                if (flatE) return flatE[r * cols + c];
                return -999.0f;
            }
            inline float getN(int r, int c) const {
                if (mapN) return mapN->data[r][c];
                if (flatN) return flatN[r * cols + c];
                return -999.0f;
            }
        };

        auto isBldg = [&](int r, int c) -> float {
            if (r < 0 || r >= rows || c < 0 || c >= cols) return 0.0f;
            if (pcr::isMV(DEM->data[r][c])) return 0.0f;
            return (BH && !pcr::isMV(BH->data[r][c]) && BH->data[r][c] > 2.0f) ? 1.0f : 0.0f;
        };

        int maxThreads = 1;
        #if FF_HAS_OMP
        maxThreads = omp_get_max_threads();
        #endif

        struct ReflectedHit
        {
            int r, c;
            float linD, linE, linN;
        };

        std::vector<std::vector<ReflectedHit>> currentReflHits(16);
        auto accumRefl = std::shared_ptr<cTMapf>(DEM->GetCopy0());

        FF::coutI << "[FastNoise] Evaluating acoustic propagation with 16-directional grid sweeping, angular beam wedges, and facade reflections..." << std::endl;
        auto tStart = std::chrono::high_resolution_clock::now();

        auto runSweepDirection = [&](
            int dIdx,
            const SourceAccessor& srcAcc,
            bool recordReflections,
            bool isReflectionPass,
            std::vector<std::vector<ReflectedHit>>& outReflHits,
            std::vector<float>& targetDirDay,
            std::vector<float>& targetDirEv,
            std::vector<float>& targetDirNi)
        {
            const auto& sd = sweepDirs[dIdx];
            int dr = sd.dr;
            int dc = sd.dc;
            int pr = sd.pr;
            int pc = sd.pc;
            double ds = sweepDs[dIdx];
            double dp = sweepDp[dIdx];
            const auto& startCells = allStartCells[dIdx];

            double vx = (dc * dx) / ds;
            double vy = (-dr * dy) / ds;

            auto depositCell = [&](int ir, int ic, double dLat, float linValD, float linValE, float linValN, double w_max) -> bool {
                if (ir < 0 || ir >= rows || ic < 0 || ic >= cols) return false;
                if (pcr::isMV(DEM->data[ir][ic])) return false;
                float k_bh = (BH && !pcr::isMV(BH->data[ir][ic])) ? BH->data[ir][ic] : 0.0f;
                float k_bar = (BARRIERS && !pcr::isMV(BARRIERS->data[ir][ic])) ? BARRIERS->data[ir][ic] : 0.0f;
                if (k_bh > 2.0f || k_bar > 0.5f) return false;

                int kIdx = ir * cols + ic;
                double xi = (w_max > 1e-6) ? (dLat / w_max) : 0.0;
                float W = g_cosWedgeLUT.get(xi);

                float valD = linValD * W;
                if (valD > targetDirDay[kIdx])
                {
                    targetDirDay[kIdx] = valD;
                    targetDirEv[kIdx]  = linValE * W;
                    targetDirNi[kIdx]  = linValN * W;
                }
                return true;
            };

            for (int li = 0; li < static_cast<int>(startCells.size()); ++li)
            {
                int r = startCells[li].first;
                int c = startCells[li].second;

                bool active = false;
                double curDist = 0.0;
                double srcD = 0.0, srcE = 0.0, srcN = 0.0;
                double maxObsH = 0.0, canopyDist = 0.0;

                while (r >= 0 && r < rows && c >= 0 && c < cols)
                {
                    float zDem = DEM->data[r][c];
                    if (pcr::isMV(zDem))
                    {
                        active = false;
                        r += dr;
                        c += dc;
                        continue;
                    }

                    float sD = srcAcc.getD(r, c);
                    float sE = srcAcc.getE(r, c);
                    float sN = srcAcc.getN(r, c);

                    if (sD > 0.0f)
                    {
                        active = true;
                        curDist = 0.0;
                        maxObsH = 0.0;
                        canopyDist = 0.0;
                        srcD = sD;
                        srcE = sE;
                        srcN = sN;

                        r += dr;
                        c += dc;
                        continue;
                    }
                    else if (active)
                    {
                        float bh = (BH && !pcr::isMV(BH->data[r][c])) ? BH->data[r][c] : 0.0f;

                        if (bh > 2.0f)
                        {
                            // Impact with exterior building facade
                            if (recordReflections)
                            {
                                double d_impact = curDist + ds;
                                if (d_impact <= maxDistMeters)
                                {
                                    double d_a = 10.0 * std::log10(2.0 * M_PI * std::max(1.0, d_impact));
                                    double d_l = alphaAir * d_impact;
                                    double d_b = 0.5 * std::min(4.5, 1.5 + 2.0 * std::log10(1.0 + d_impact / 50.0));
                                    double d_screen = 0.0;
                                    if (maxObsH > 0.3)
                                    {
                                        double dMid = d_impact * 0.5;
                                        double pathDiff = 2.0 * std::sqrt(dMid * dMid + maxObsH * maxObsH) - d_impact;
                                        double N = (2.0 * pathDiff) / 0.68;
                                        d_screen = std::clamp(10.0 * std::log10(3.0 + 20.0 * std::max(0.0, N)), 0.0, 25.0);
                                    }
                                    double d_veg = std::clamp(0.05 * canopyDist, 0.0, 10.0);
                                    double d_m = std::max(0.0, 1.5 * (1.0 - 2.0 * (0.75 + receiverHeight) / d_impact));
                                    double d_tot = d_a + d_l + d_b + d_screen + d_veg - d_m;

                                    double l_inc_d = srcD - d_tot;
                                    double l_inc_e = srcE - d_tot;
                                    double l_inc_n = srcN - d_tot;

                                    if (l_inc_d > 20.0)
                                    {
                                        float bw = isBldg(r - 1, c - 1) + 2.0f * isBldg(r, c - 1) + isBldg(r + 1, c - 1);
                                        float be = isBldg(r - 1, c + 1) + 2.0f * isBldg(r, c + 1) + isBldg(r + 1, c + 1);
                                        float bn = isBldg(r - 1, c - 1) + 2.0f * isBldg(r - 1, c) + isBldg(r - 1, c + 1);
                                        float bs = isBldg(r + 1, c - 1) + 2.0f * isBldg(r + 1, c) + isBldg(r + 1, c + 1);

                                        double gx = (be - bw) / dx;
                                        double gy = (bn - bs) / dy;
                                        double nx = -gx;
                                        double ny = -gy;
                                        double nlen = std::sqrt(nx * nx + ny * ny);
                                        if (nlen < 1e-4) { nx = -vx; ny = -vy; }
                                        else { nx /= nlen; ny /= nlen; }

                                        double v_dot_n = vx * nx + vy * ny;
                                        if (v_dot_n > 0.0) { nx = -nx; ny = -ny; v_dot_n = -v_dot_n; }

                                        double rx = vx - 2.0 * v_dot_n * nx;
                                        double ry = vy - 2.0 * v_dot_n * ny;
                                        double rlen = std::sqrt(rx * rx + ry * ry);
                                        if (rlen > 1e-6) { rx /= rlen; ry /= rlen; }

                                        int bestDir = -1;
                                        double maxDot = -2.0;
                                        for (int k = 0; k < 16; ++k)
                                        {
                                            double k_vx = (sweepDirs[k].dc * dx) / sweepDs[k];
                                            double k_vy = (-sweepDirs[k].dr * dy) / sweepDs[k];
                                            double dot = rx * k_vx + ry * k_vy;
                                            if (dot > maxDot) { maxDot = dot; bestDir = k; }
                                        }

                                        if (bestDir >= 0)
                                        {
                                            double reflD = std::pow(10.0, (l_inc_d - 0.97) / 10.0);
                                            double reflE = std::pow(10.0, (l_inc_e - 0.97) / 10.0);
                                            double reflN = std::pow(10.0, (l_inc_n - 0.97) / 10.0);
                                            int prevR = r - dr;
                                            int prevC = c - dc;
                                            if (prevR >= 0 && prevR < rows && prevC >= 0 && prevC < cols)
                                            {
                                                outReflHits[bestDir].push_back({
                                                    prevR, prevC,
                                                    static_cast<float>(reflD),
                                                    static_cast<float>(reflE),
                                                    static_cast<float>(reflN)
                                                });
                                            }
                                        }
                                    }
                                }
                            }

                            active = false;
                            r += dr;
                            c += dc;
                            continue;
                        }

                        curDist += ds;
                        if (curDist > maxDistMeters)
                        {
                            active = false;
                            r += dr;
                            c += dc;
                            continue;
                        }

                        float bar = (BARRIERS && !pcr::isMV(BARRIERS->data[r][c])) ? BARRIERS->data[r][c] : 0.0f;
                        float obs = bh + bar;
                        if (obs > maxObsH) maxObsH = obs;

                        if (CH && !pcr::isMV(CH->data[r][c]) && CH->data[r][c] > 2.0f)
                        {
                            canopyDist += ds;
                        }

                        double gB = 0.5;
                        if (GROUND_ABS && !pcr::isMV(GROUND_ABS->data[r][c]))
                        {
                            gB = std::clamp(static_cast<double>(GROUND_ABS->data[r][c]), 0.0, 1.0);
                        }

                        double d_a = 10.0 * std::log10(2.0 * M_PI * std::max(1.0, curDist));
                        double d_l = alphaAir * curDist;
                        double d_b = gB * std::min(4.5, 1.5 + 2.0 * std::log10(1.0 + curDist / 50.0));
                        double d_screen = 0.0;
                        if (maxObsH > 0.3)
                        {
                            double dMid = curDist * 0.5;
                            double pathDiff = 2.0 * std::sqrt(dMid * dMid + maxObsH * maxObsH) - curDist;
                            double N = (2.0 * pathDiff) / 0.68;
                            d_screen = std::clamp(10.0 * std::log10(3.0 + 20.0 * std::max(0.0, N)), 0.0, 25.0);
                        }
                        double d_veg = std::clamp(0.05 * canopyDist, 0.0, 10.0);
                        double d_m = std::max(0.0, 1.5 * (1.0 - 2.0 * (0.75 + receiverHeight) / curDist));

                        double d_tot = d_a + d_l + d_b + d_screen + d_veg - d_m;
                        double l_d = srcD - d_tot;
                        double l_e = srcE - d_tot;
                        double l_n = srcN - d_tot;

                        // Once attenuated below audibility floor (0 dB), sound can no longer recover along this ray
                        if (l_d <= 0.0)
                        {
                            active = false;
                            r += dr;
                            c += dc;
                            continue;
                        }

                        float linD = static_cast<float>(std::pow(10.0, l_d / 10.0));
                        float linE = static_cast<float>(std::pow(10.0, l_e / 10.0));
                        float linN = static_cast<float>(std::pow(10.0, l_n / 10.0));

                        // Centerline cell (k = 0, W = 1.0)
                        depositCell(r, c, 0.0, linD, linE, linN, 1.0);

                        // Angular Beam Wedge Fan Spreading
                        double w_max = sd.wedgeFactor * curDist;
                        int maxSteps = static_cast<int>(std::floor(w_max / dp));
                        int maxAllowedSteps = std::max(64, static_cast<int>(std::ceil(maxDistMeters / dp)));
                        maxSteps = std::min(maxSteps, maxAllowedSteps);

                        // Dense supercover centerline connection along ray direction (prevents diagonal skips)
                        if (dr != 0 || dc != 0)
                        {
                            int minR = std::min(r, r + dr);
                            int maxR = std::max(r, r + dr);
                            int minC = std::min(c, c + dc);
                            int maxC = std::max(c, c + dc);
                            for (int ir = minR; ir <= maxR; ++ir)
                            {
                                for (int ic = minC; ic <= maxC; ++ic)
                                {
                                    depositCell(ir, ic, 0.0, linD, linE, linN, 1.0);
                                }
                            }
                        }

                        // Lateral wedge fan spreading
                        if (pr == 0 || pc == 0)
                        {
                            // Cardinal directions (orthogonal 1D steps, no diagonal bridges needed)
                            for (int k = 1; k <= maxSteps; ++k)
                            {
                                int ir = r + k * pr;
                                int ic = c + k * pc;
                                if (!depositCell(ir, ic, k * dp, linD, linE, linN, w_max)) break;
                            }
                            for (int k = 1; k <= maxSteps; ++k)
                            {
                                int ir = r - k * pr;
                                int ic = c - k * pc;
                                if (!depositCell(ir, ic, k * dp, linD, linE, linN, w_max)) break;
                            }
                        }
                        else
                        {
                            // Diagonal and Knight directions: 2D supercover steps to fill all intermediate pixels
                            for (int k = 1; k <= maxSteps; ++k)
                            {
                                int r_prev = r + (k - 1) * pr;
                                int c_prev = c + (k - 1) * pc;
                                int r_curr = r + k * pr;
                                int c_curr = c + k * pc;

                                int minR = std::min(r_prev, r_curr);
                                int maxR = std::max(r_prev, r_curr);
                                int minC = std::min(c_prev, c_curr);
                                int maxC = std::max(c_prev, c_curr);

                                bool hitObstacle = false;
                                for (int ir = minR; ir <= maxR; ++ir)
                                {
                                    for (int ic = minC; ic <= maxC; ++ic)
                                    {
                                        double dLat = std::abs((ir - r) * pr * dy * dy + (ic - c) * pc * dx * dx) / dp;
                                        if (!depositCell(ir, ic, dLat, linD, linE, linN, w_max))
                                        {
                                            hitObstacle = true;
                                        }
                                    }
                                }
                                if (hitObstacle) break;
                            }

                            for (int k = 1; k <= maxSteps; ++k)
                            {
                                int r_prev = r - (k - 1) * pr;
                                int c_prev = c - (k - 1) * pc;
                                int r_curr = r - k * pr;
                                int c_curr = c - k * pc;

                                int minR = std::min(r_prev, r_curr);
                                int maxR = std::max(r_prev, r_curr);
                                int minC = std::min(c_prev, c_curr);
                                int maxC = std::max(c_prev, c_curr);

                                bool hitObstacle = false;
                                for (int ir = minR; ir <= maxR; ++ir)
                                {
                                    for (int ic = minC; ic <= maxC; ++ic)
                                    {
                                        double dLat = std::abs((ir - r) * pr * dy * dy + (ic - c) * pc * dx * dx) / dp;
                                        if (!depositCell(ir, ic, dLat, linD, linE, linN, w_max))
                                        {
                                            hitObstacle = true;
                                        }
                                    }
                                }
                                if (hitObstacle) break;
                            }
                        }
                    }

                    r += dr;
                    c += dc;
                }
            }
        };

        // Multi-pass simulation (Primary emission + facade reflection bounces)
        for (int pass = 0; pass < numBounces; ++pass)
        {
            bool isPrimary = (pass == 0);
            bool recordReflections = (pass + 1 < numBounces);

            if (isPrimary)
            {
                FF::coutI << "[FastNoise] Evaluating Pass 0: Primary acoustic emissions (16 directional wedges in parallel)..." << std::endl;
                std::vector<std::vector<float>> allDirDay(16, std::vector<float>(rows * cols, 0.0f));
                std::vector<std::vector<float>> allDirEv(16, std::vector<float>(rows * cols, 0.0f));
                std::vector<std::vector<float>> allDirNi(16, std::vector<float>(rows * cols, 0.0f));
                std::vector<std::vector<std::vector<ReflectedHit>>> allReflHits(16, std::vector<std::vector<ReflectedHit>>(16));

                GetTPool().parallel_for(0, 16, [&](int dIdx) {
                    SourceAccessor srcAcc;
                    srcAcc.mapD = sourcesRasterDay.get();
                    srcAcc.mapE = sourcesRasterEv.get();
                    srcAcc.mapN = sourcesRasterNi.get();
                    srcAcc.cols = cols;

                    runSweepDirection(dIdx, srcAcc, recordReflections, false,
                                      allReflHits[dIdx],
                                      allDirDay[dIdx], allDirEv[dIdx], allDirNi[dIdx]);
                });
                GetTPool().waitext();

                // Accumulate all 16 directional passes into global maps in parallel across rows
                GetTPool().parallel_for(0, rows, [&](int r) {
                    int rowOffset = r * cols;
                    for (int c = 0; c < cols; ++c) {
                        int idx = rowOffset + c;
                        float sumD = 0.0f, sumE = 0.0f, sumN = 0.0f;
                        for (int d = 0; d < 16; ++d) {
                            sumD += allDirDay[d][idx];
                            sumE += allDirEv[d][idx];
                            sumN += allDirNi[d][idx];
                        }
                        if (sumD > 0.0f) {
                            accumDay->data[r][c] += sumD;
                            accumEv->data[r][c]  += sumE;
                            accumNi->data[r][c]  += sumN;
                        }
                    }
                });
                GetTPool().waitext();

                if (recordReflections) {
                    for (int d = 0; d < 16; ++d) {
                        for (int k = 0; k < 16; ++k) {
                            if (!allReflHits[d][k].empty()) {
                                currentReflHits[k].insert(currentReflHits[k].end(), allReflHits[d][k].begin(), allReflHits[d][k].end());
                            }
                        }
                    }
                }

                auto tNow = std::chrono::high_resolution_clock::now();
                double elapsedSoFar = std::chrono::duration<double>(tNow - tStart).count();
                FF::coutI << "[FastNoise] Pass 0 completed in " << std::fixed << std::setprecision(2) << elapsedSoFar << "s." << std::endl;
            }
            else
            {
                int totalHits = 0;
                for (int k = 0; k < 16; ++k) totalHits += static_cast<int>(currentReflHits[k].size());
                if (totalHits == 0)
                {
                    FF::coutI << "[FastNoise] Pass " << pass << ": No facade reflections detected; skipping bounce." << std::endl;
                    break;
                }

                FF::coutI << "[FastNoise] Evaluating Pass " << pass << ": Reflection Bounce " << pass 
                          << " (" << totalHits << " facade sources off buildings in parallel)..." << std::endl;

                std::vector<int> activeDirs;
                for (int d = 0; d < 16; ++d) {
                    if (!currentReflHits[d].empty()) activeDirs.push_back(d);
                }

                if (activeDirs.empty()) break;

                std::vector<std::vector<ReflectedHit>> nextReflHits(16);
                std::vector<std::vector<float>> reflDirDay(16, std::vector<float>(rows * cols, 0.0f));
                std::vector<std::vector<float>> reflDirEv(16, std::vector<float>(rows * cols, 0.0f));
                std::vector<std::vector<float>> reflDirNi(16, std::vector<float>(rows * cols, 0.0f));
                std::vector<std::vector<std::vector<ReflectedHit>>> reflDirHits(16, std::vector<std::vector<ReflectedHit>>(16));

                GetTPool().parallel_for(0, static_cast<int>(activeDirs.size()), [&](int aIdx) {
                    int dIdx = activeDirs[aIdx];

                    std::vector<float> secSrcDay(rows * cols, -999.0f);
                    std::vector<float> secSrcEv(rows * cols, -999.0f);
                    std::vector<float> secSrcNi(rows * cols, -999.0f);

                    for (const auto& h : currentReflHits[dIdx])
                    {
                        int idx = h.r * cols + h.c;
                        float prevD = secSrcDay[idx];
                        if (prevD > -900.0f)
                        {
                            secSrcDay[idx] = 10.0f * std::log10(std::pow(10.0f, prevD / 10.0f) + h.linD);
                            secSrcEv[idx]  = 10.0f * std::log10(std::pow(10.0f, secSrcEv[idx] / 10.0f) + h.linE);
                            secSrcNi[idx]  = 10.0f * std::log10(std::pow(10.0f, secSrcNi[idx] / 10.0f) + h.linN);
                        }
                        else
                        {
                            secSrcDay[idx] = 10.0f * std::log10(std::max(1e-10f, h.linD));
                            secSrcEv[idx]  = 10.0f * std::log10(std::max(1e-10f, h.linE));
                            secSrcNi[idx]  = 10.0f * std::log10(std::max(1e-10f, h.linN));
                        }
                    }

                    SourceAccessor srcAcc;
                    srcAcc.flatD = secSrcDay.data();
                    srcAcc.flatE = secSrcEv.data();
                    srcAcc.flatN = secSrcNi.data();
                    srcAcc.cols = cols;

                    runSweepDirection(dIdx, srcAcc, recordReflections, true,
                                      reflDirHits[dIdx],
                                      reflDirDay[dIdx], reflDirEv[dIdx], reflDirNi[dIdx]);
                });
                GetTPool().waitext();

                // Accumulate reflections into global maps and accumRefl in parallel
                GetTPool().parallel_for(0, rows, [&](int r) {
                    int rowOffset = r * cols;
                    for (int c = 0; c < cols; ++c) {
                        int idx = rowOffset + c;
                        float sumD = 0.0f, sumE = 0.0f, sumN = 0.0f;
                        for (int d = 0; d < 16; ++d) {
                            sumD += reflDirDay[d][idx];
                            sumE += reflDirEv[d][idx];
                            sumN += reflDirNi[d][idx];
                        }
                        if (sumD > 0.0f) {
                            accumDay->data[r][c] += sumD;
                            accumEv->data[r][c]  += sumE;
                            accumNi->data[r][c]  += sumN;
                            if (accumRefl) accumRefl->data[r][c] += sumD;
                        }
                    }
                });
                GetTPool().waitext();

                if (recordReflections) {
                    for (int d = 0; d < 16; ++d) {
                        for (int k = 0; k < 16; ++k) {
                            if (!reflDirHits[d][k].empty()) {
                                nextReflHits[k].insert(nextReflHits[k].end(), reflDirHits[d][k].begin(), reflDirHits[d][k].end());
                            }
                        }
                    }
                }

                currentReflHits = std::move(nextReflHits);
            }
        }

        // 7. Post-sweep gap healing:
        // Ensure no isolated terrain pixels were skipped due to discrete grid rasterization or diagonal steps.
        {
            auto origDay = std::shared_ptr<cTMapf>(accumDay->GetCopy0());
            auto origEv  = std::shared_ptr<cTMapf>(accumEv->GetCopy0());
            auto origNi  = std::shared_ptr<cTMapf>(accumNi->GetCopy0());

            GetTPool().parallel_for(1, rows - 1, [&](int r) {
                for (int c = 1; c < cols - 1; ++c)
                {
                    if (pcr::isMV(DEM->data[r][c])) continue;
                    float bh = (BH && !pcr::isMV(BH->data[r][c])) ? BH->data[r][c] : 0.0f;
                    float bar = (BARRIERS && !pcr::isMV(BARRIERS->data[r][c])) ? BARRIERS->data[r][c] : 0.0f;
                    if (bh > 2.0f || bar > 0.5f) continue;

                    float curVal = origDay->data[r][c];

                    float nSumD = 0.0f, nSumE = 0.0f, nSumN = 0.0f;
                    int validN = 0;
                    float maxN = 0.0f;

                    for (int dr = -1; dr <= 1; ++dr)
                    {
                        for (int dc = -1; dc <= 1; ++dc)
                        {
                            if (dr == 0 && dc == 0) continue;
                            int nr = r + dr;
                            int nc = c + dc;
                            if (pcr::isMV(DEM->data[nr][nc])) continue;
                            float n_bh = (BH && !pcr::isMV(BH->data[nr][nc])) ? BH->data[nr][nc] : 0.0f;
                            float n_bar = (BARRIERS && !pcr::isMV(BARRIERS->data[nr][nc])) ? BARRIERS->data[nr][nc] : 0.0f;
                            if (n_bh > 2.0f || n_bar > 0.5f) continue;

                            float nValD = origDay->data[nr][nc];
                            if (nValD > 0.0f)
                            {
                                nSumD += nValD;
                                nSumE += origEv->data[nr][nc];
                                nSumN += origNi->data[nr][nc];
                                if (nValD > maxN) maxN = nValD;
                                validN++;
                            }
                        }
                    }

                    if (validN >= 3 && maxN > 0.0f)
                    {
                        float avgD = nSumD / validN;
                        if (curVal < 0.15f * avgD)
                        {
                            accumDay->data[r][c] = avgD;
                            accumEv->data[r][c]  = nSumE / validN;
                            accumNi->data[r][c]  = nSumN / validN;
                        }
                    }
                }
            });
            GetTPool().waitext();
        }

        // 8. Post-process: Convert accumulated linear energy + background noise to decibels
        GetTPool().parallel_for(0, rows, [&](int r) {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(mapLDay->data[r][c]);
                    pcr::setMV(mapLEvening->data[r][c]);
                    pcr::setMV(mapLNight->data[r][c]);
                    pcr::setMV(mapLDen->data[r][c]);
                    continue;
                }

                // Add background ambient noise floor
                double bgLin = std::pow(10.0, backgroundNoiseDBA / 10.0);
                double dLin = accumDay->data[r][c] + bgLin;
                double eLin = accumEv->data[r][c]  + bgLin;
                double nLin = accumNi->data[r][c]  + bgLin;

                double totDay = 10.0 * std::log10(std::max(1e-10, dLin));
                double totEv  = 10.0 * std::log10(std::max(1e-10, eLin));
                double totNi  = 10.0 * std::log10(std::max(1e-10, nLin));

                // Official Lden rating level (with +5 dB evening, +10 dB night penalties per RMG 2012 / EU Directive):
                double linLden = (12.0 / 24.0) * dLin
                               + (4.0 / 24.0)  * (eLin * std::pow(10.0, 5.0 / 10.0))
                               + (8.0 / 24.0)  * (nLin * std::pow(10.0, 10.0 / 10.0));
                double totLden = 10.0 * std::log10(std::max(1e-10, linLden));

                mapLDay->data[r][c] = static_cast<float>(totDay);
                mapLEvening->data[r][c] = static_cast<float>(totEv);
                mapLNight->data[r][c] = static_cast<float>(totNi);
                mapLDen->data[r][c] = static_cast<float>(totLden);
            }
        });
        GetTPool().waitext();

        auto tEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = tEnd - tStart;
        FF::coutI << "[FastNoise] Propagation completed in " << std::fixed << std::setprecision(2) << elapsed.count() << " seconds." << std::endl;

        outputs.LDAY = mapLDay;
        outputs.LEVENING = mapLEvening;
        outputs.LNIGHT = mapLNight;
        outputs.LDEN = mapLDen;

        // Populate outputs.REFLECTIONS with accumulated facade reflection sound field
        if (!outputs.REFLECTIONS)
        {
            outputs.REFLECTIONS = std::shared_ptr<cTMapf>(DEM->GetCopy0());
        }
        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                if (pcr::isMV(DEM->data[r][c]))
                {
                    pcr::setMV(outputs.REFLECTIONS->data[r][c]);
                    continue;
                }
                float rLin = accumRefl->data[r][c];
                outputs.REFLECTIONS->data[r][c] = (rLin > 1e-10f) ? static_cast<float>(10.0 * std::log10(rLin)) : 0.0f;
            }
        }

        // Log domain statistics
        float minLden = 1e30f, maxLden = -1e30f;
        double sumLden = 0.0;
        int validCells = 0;

        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                float val = mapLDen->data[r][c];
                if (!pcr::isMV(val) && std::isfinite(val))
                {
                    minLden = std::min(minLden, val);
                    maxLden = std::max(maxLden, val);
                    sumLden += val;
                    validCells++;
                }
            }
        }

        if (validCells > 0)
        {
            double meanLden = sumLden / validCells;
            FF::coutI << "=================================================================" << std::endl;
            FF::coutI << "[FastNoise] Domain Lden Statistics (RMG 2012 Bijlage III):" << std::endl;
            FF::coutI << "            Mean: " << std::fixed << std::setprecision(1) << meanLden << " dB(A)" << std::endl;
            FF::coutI << "            Min:  " << minLden << " dB(A)" << std::endl;
            FF::coutI << "            Max:  " << maxLden << " dB(A)" << std::endl;
            FF::coutI << "            Evaluated cells: " << validCells << " (" << rows << "x" << cols << ")" << std::endl;
            FF::coutI << "=================================================================" << std::endl;
        }

        return 0;
    }
}

#endif // FASTNOISE_H
