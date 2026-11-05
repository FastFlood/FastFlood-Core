#ifndef FASTMODELDATA_H
#define FASTMODELDATA_H

#include <memory>
#include <vector>
#include <string>
#include "map/map.h"


struct FastModelLandUseData {
    std::shared_ptr<cTMapf> lu;
    std::shared_ptr<cTMapf> manning;
    bool isLUDef = false;
};

struct FastModelSoilData {
    std::shared_ptr<cTMapf> infiltration;
    std::shared_ptr<cTMapf> thetas;
    std::shared_ptr<cTMapf> soilclay;
    std::shared_ptr<cTMapf> soilsilt;
    std::shared_ptr<cTMapf> soilsand;
    std::shared_ptr<cTMapf> soilorganic;
    std::shared_ptr<cTMapf> soilgravel;
};

struct FastModelSoilMoistureStats {
    std::shared_ptr<cTMapf> sm_mean = nullptr;       // Annual mean volumetric soil moisture [m3/m3]
    std::shared_ptr<cTMapf> sm_min = nullptr;        // Annual minimum (driest month) soil moisture [m3/m3]
    std::shared_ptr<cTMapf> sm_max = nullptr;        // Annual maximum (wettest month) soil moisture [m3/m3]
    std::shared_ptr<cTMapf> sm_var = nullptr;        // Seasonal range (max - min) [m3/m3]
    std::shared_ptr<cTMapf> sm_quantile = nullptr;   // Target quantile soil moisture (e.g. q_sm) [m3/m3]
    std::vector<std::shared_ptr<cTMapf>> sm_monthly; // 12 monthly soil moisture rasters (Jan-Dec) [m3/m3]
};

// Functions to fetch data
std::shared_ptr<FastModelLandUseData> FetchFastModelLandUseData(std::shared_ptr<cTMapf> ref_map, const std::string& token, const std::string& customurl);
std::shared_ptr<FastModelSoilData> FetchFastModelSoilData(std::shared_ptr<cTMapf> ref_map, const std::string& token, const std::string& customurl);
std::shared_ptr<cTMapf> FetchFastModelSoilMoistureData(std::shared_ptr<cTMapf> ref_map, std::shared_ptr<cTMapf> dem, double q_sm, const std::string& token, const std::string& customurl);
std::shared_ptr<FastModelSoilMoistureStats> FetchFastModelSoilMoistureMonthlyStats(std::shared_ptr<cTMapf> ref_map, std::shared_ptr<cTMapf> dem, double q_sm, const std::string& token, const std::string& customurl);
std::shared_ptr<cTMapf> FetchFastModelSoilDepthData(std::shared_ptr<cTMapf> ref_map, const std::string& token, const std::string& customurl);

#endif // FASTMODELDATA_H
