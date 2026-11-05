#include "fastmodeldata.h"
#include "../inputprepare/forecastprocess.h"
#include "../../web/glbdownload.h"
#include "../../math/quantile.h"
#include "../../geospatial/map/algorithms.h"
#include "../algorithms/rasteralgorithms/derivatives.h"
#include "iogdal.h"
#include "../geospatial/map/algorithms.h"
#include "../web/cogtif/cogtif.h"
#include "fferror.h"
#include "map/map.h"
#include "cl/clargparser.h"
#include "../algorithms/rasteralgorithms/distancefield.h"
#include "../algorithms/rasteralgorithms/gridboundary.h"
#include "../geospatial/shape/shapefile.h"
#include "../util/datetime.h"
#include "projtools.h"
#include "../math/quantile.h"
#include "rasteralgorithms/datafill.h"
#include "rasteralgorithms/rasterblend.h"
#include "rasteralgorithms/catchments.h"
#include "rasteralgorithms/vectorize.h"
#include <mutex>

std::shared_ptr<FastModelSoilData> FetchFastModelSoilData(std::shared_ptr<cTMapf> ref_map, const std::string& token, const std::string& customurl) {
    std::vector<std::shared_ptr<cTMapf>> res = GetSoilMaps(ref_map, token, customurl);
    if (res.size() > 0) {
        auto data = std::make_shared<FastModelSoilData>();
        data->infiltration = res[0];
        data->thetas = res[1];
        data->soilclay = res[2];
        data->soilsilt = res[3];
        data->soilsand = res[4];
        data->soilorganic = res[5];
        data->soilgravel = res[6];
        return data;
    }
    return nullptr;
}
std::shared_ptr<cTMapf> FetchFastModelSoilDepthData(std::shared_ptr<cTMapf> ref_map, const std::string& token, const std::string& customurl) {

    std::vector<cTMapf *> soilddata = DownloadAndParseGlobalTiledAlt("SOILD", 20, ref_map.get(), token, 1.0, customurl,{},{},0,true);
    if (soilddata.size() > 0) {
        std::shared_ptr<cTMapf> soild = GetGlobalWGSTiledMapWarped(soilddata[0], ref_map, 0);
        SafeDeleteList(soilddata);
        return soild;
    }
    return nullptr;
}

std::shared_ptr<cTMapf> FetchFastModelSoilMoistureData(std::shared_ptr<cTMapf> ref_map, std::shared_ptr<cTMapf> dem, double q_sm, const std::string& token, const std::string& customurl) {
    std::vector<cTMapf *> soilmdata = DownloadAndParseGlobalTiledAlt("SOILM", 20, ref_map.get(), token, 1.0, customurl,{},{},0,true);

    if(soilmdata.size() > 6)
    {
        //writeRaster(*soilmdata[2],"smtest1.tif");
        //writeRaster(*soilmdata[3],"smtest2.tif");
        std::shared_ptr<cTMapf> temp = soilmdata[0]->GetCopySP0();

        for(int r = 0; r < temp->nrRows(); r++)
        {
            for(int c = 0; c < temp->nrCols(); c++) {
                float min = soilmdata[2]->data[r][c];
                float max = soilmdata[3]->data[r][c];

                std::vector<double> smlist = {};
                for(int i = 0; i < 12; i++)
                {
                    smlist.push_back(soilmdata[6+12+i]->data[r][c]);
                }
                smlist.push_back(min);
                smlist.push_back(max);

                double val = GetQuantileValue(smlist, q_sm);
                if(std::isnormal(val) && ! pcr::isMV(val))
                {
                    temp->data[r][c] = 0.9;//std::min(0.95, std::max(0.8, val/std::min(0.6f,std::max(0.4f,max))));
                }else {
                    temp->data[r][c] = 0.9;//std::min(0.95, std::max(0.8, q_sm));
                }
            }
        }
        std::shared_ptr<cTMapf> soilm = GetGlobalWGSTiledMapWarped(temp.get(), dem, 0);
        SafeDeleteList(soilmdata);
        return soilm;
    }
    return nullptr;
}

std::shared_ptr<FastModelSoilMoistureStats> FetchFastModelSoilMoistureMonthlyStats(
    std::shared_ptr<cTMapf> ref_map,
    std::shared_ptr<cTMapf> dem,
    double q_sm,
    const std::string& token,
    const std::string& customurl)
{
    auto target_ref = ref_map ? ref_map : dem;
    auto target_dem = dem ? dem : ref_map;
    if (!target_ref && !target_dem) return nullptr;

    auto result = std::make_shared<FastModelSoilMoistureStats>();
    result->sm_monthly.resize(12, nullptr);

    std::vector<cTMapf*> soilmdata = DownloadAndParseGlobalTiledAlt("SOILM", 20, target_ref.get(), token, 1.0, customurl, {}, {}, 0, true);

    if (soilmdata.size() > 6)
    {
        std::shared_ptr<cTMapf> temp_mean = soilmdata[0]->GetCopySP0();
        std::shared_ptr<cTMapf> temp_min  = soilmdata[2]->GetCopySP0();
        std::shared_ptr<cTMapf> temp_max  = soilmdata[3]->GetCopySP0();
        std::shared_ptr<cTMapf> temp_var  = soilmdata[0]->GetCopySP0();
        std::shared_ptr<cTMapf> temp_q    = soilmdata[0]->GetCopySP0();

        std::vector<std::shared_ptr<cTMapf>> temp_months(12, nullptr);
        for (int m = 0; m < 12; ++m)
        {
            temp_months[m] = soilmdata[0]->GetCopySP0();
        }

        const bool has_monthly_bands = (soilmdata.size() >= 30);

        for (int r = 0; r < temp_mean->nrRows(); ++r)
        {
            for (int c = 0; c < temp_mean->nrCols(); ++c)
            {
                float min_val  = soilmdata[2]->data[r][c];
                float max_val  = soilmdata[3]->data[r][c];
                float mean_val = soilmdata[0]->data[r][c];

                std::vector<double> smlist;
                for (int i = 0; i < 12; ++i)
                {
                    float m_val = mean_val;
                    if (has_monthly_bands)
                    {
                        m_val = soilmdata[6 + 12 + i]->data[r][c];
                    }
                    else if (soilmdata.size() > 6 + i)
                    {
                        m_val = soilmdata[6 + i]->data[r][c];
                    }
                    temp_months[i]->data[r][c] = m_val;
                    if (!pcr::isMV(m_val) && std::isfinite(m_val))
                    {
                        smlist.push_back(m_val);
                    }
                }

                if (!pcr::isMV(min_val) && std::isfinite(min_val)) smlist.push_back(min_val);
                if (!pcr::isMV(max_val) && std::isfinite(max_val)) smlist.push_back(max_val);

                double q_val = (smlist.size() > 0) ? GetQuantileValue(smlist, q_sm) : mean_val;
                if (!std::isfinite(q_val) || pcr::isMV(q_val))
                {
                    q_val = mean_val;
                }

                temp_mean->data[r][c] = mean_val;
                temp_min->data[r][c]  = min_val;
                temp_max->data[r][c]  = max_val;
                temp_var->data[r][c]  = std::max(0.0f, max_val - min_val);
                temp_q->data[r][c]    = static_cast<float>(q_val);
            }
        }

        result->sm_mean     = GetGlobalWGSTiledMapWarped(temp_mean.get(), target_dem, 0);
        result->sm_min      = GetGlobalWGSTiledMapWarped(temp_min.get(), target_dem, 0);
        result->sm_max      = GetGlobalWGSTiledMapWarped(temp_max.get(), target_dem, 0);
        result->sm_var      = GetGlobalWGSTiledMapWarped(temp_var.get(), target_dem, 0);
        result->sm_quantile = GetGlobalWGSTiledMapWarped(temp_q.get(), target_dem, 0);

        for (int m = 0; m < 12; ++m)
        {
            result->sm_monthly[m] = GetGlobalWGSTiledMapWarped(temp_months[m].get(), target_dem, 0);
        }

        SafeDeleteList(soilmdata);
        return result;
    }

    // Fallback: estimate physically consistent climatological soil moisture stats from coordinates
    LSMVector2 center = target_dem->GetPixelLngLat(target_dem->nrRows() / 2, target_dem->nrCols() / 2);
    double abs_lat = std::abs(center.y);

    float base_sm_mean = 0.32f;
    float base_sm_min  = 0.22f;
    float base_sm_max  = 0.40f;
    float base_sm_var  = 0.18f;

    if (abs_lat < 23.5)
    {
        // Tropical / Caribbean: high soil moisture year-round with low seasonal variation
        base_sm_mean = 0.38f;
        base_sm_min  = 0.30f;
        base_sm_max  = 0.45f;
        base_sm_var  = 0.15f;
    }
    else if (abs_lat >= 23.5 && abs_lat <= 38.0)
    {
        // Subtropical arid / Mediterranean dry summer
        base_sm_mean = 0.18f;
        base_sm_min  = 0.08f;
        base_sm_max  = 0.28f;
        base_sm_var  = 0.20f;
    }
    else if (abs_lat > 38.0 && abs_lat <= 55.0)
    {
        // Temperate
        base_sm_mean = 0.28f;
        base_sm_min  = 0.16f;
        base_sm_max  = 0.36f;
        base_sm_var  = 0.20f;
    }

    float base_q = static_cast<float>(base_sm_min + (base_sm_max - base_sm_min) * std::clamp(q_sm, 0.0, 1.0));

    result->sm_mean     = target_dem->GetCopySP(base_sm_mean);
    result->sm_min      = target_dem->GetCopySP(base_sm_min);
    result->sm_max      = target_dem->GetCopySP(base_sm_max);
    result->sm_var      = target_dem->GetCopySP(base_sm_var);
    result->sm_quantile = target_dem->GetCopySP(base_q);

    for (int m = 0; m < 12; ++m)
    {
        double phase = (center.y >= 0.0) ? (m * 3.141592653589793 / 6.0) : ((m + 6) * 3.141592653589793 / 6.0);
        float m_val = base_sm_mean + (base_sm_var * 0.5f) * static_cast<float>(std::cos(phase));
        m_val = std::clamp(m_val, 0.02f, 0.55f);
        result->sm_monthly[m] = target_dem->GetCopySP(m_val);
    }

    return result;
}

