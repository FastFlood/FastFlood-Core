#pragma once

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include "../../common/defines.h"
#include "../../common/fferror.h"
#include "../../geospatial/map/map.h"
#include "../../geospatial/shape/shapefile.h"
#include "../../geospatial/shape/shapelayer.h"
#include "../../geospatial/shape/feature.h"
#include "../../geospatial/shape/polygon.h"
#include "../../geospatial/shape/ring.h"
#include "../../io/iogdalshape.h"
#include "../../io/iogdal.h"
#include "../../algorithms/rasteralgorithms/rasterize.h"
#include "../../web/gendownload.h"
#include "fastrisk_vulnerability.h"
#include "fastrisk_raster.h"

namespace FastRisk {

struct PopulationExposureResult {
    std::vector<float> depths;
    std::vector<double> count_impacted;
    std::vector<double> count_total;
    ShapeFile* grid_shape = nullptr;
};

// ============================================================================
// 1. DATA DOWNLOADING / RETRIEVAL
// ============================================================================

inline ShapeFile* DownloadGlobalBuildings(cTMapf* ref, const std::string& token) {
    if (!ref) return nullptr;
    FF::coutI << "Downloading global building footprints from FastFlood PMTiles..." << std::endl;
    ShapeFile* data_build = DownloadAndParseGlobalPMTiles(
        ref,
        "https://pmtiles.fastflood.workers.dev/goog_msft_osm.pmtiles/{z}/{x}/{y}.mvt",
        {15},
        token, 1.0, 1.0
    );
    return data_build;
}

inline ShapeFile* DownloadGlobalOSM(cTMapf* ref, const std::string& token) {
    if (!ref) return nullptr;
    FF::coutI << "Downloading global OSM infrastructure from FastFlood PMTiles..." << std::endl;
    ShapeFile* data_osm = DownloadAndParseGlobalPMTiles(
        ref,
        "https://pmtiles.fastflood.workers.dev/osm_20240318.pmtiles/{z}/{x}/{y}.mvt",
        {15},
        token, 1.0, 1.0
    );
    return data_osm;
}

inline ShapeFile* DownloadGlobalRoads(cTMapf* ref, const std::string& token, ShapeFile** out_full_osm = nullptr) {
    ShapeFile* data_osm = DownloadGlobalOSM(ref, token);
    if (!data_osm) return nullptr;

    ShapeFile* roads = GetVectorLayerByName(data_osm, "roads");
    if (out_full_osm) {
        *out_full_osm = data_osm;
    } else {
        delete data_osm;
    }
    return roads;
}

inline ShapeFile* DownloadGlobalAdmin(cTMapf* ref, const std::string& token) {
    if (!ref) return nullptr;
    FF::coutI << "Downloading administrative boundaries from FastFlood PMTiles..." << std::endl;
    ShapeFile* data_adm = DownloadAndParseGlobalPMTiles(
        ref,
        "https://pmtiles.fastflood.workers.dev/adm.pmtiles/{z}/{x}/{y}.mvt",
        {10},
        token, 1.0, 1.0
    );
    if (data_adm && data_adm->GetLayerCount() > 0) {
        ShapeLayer* l_adm = data_adm->GetLayer(0);
        for (int i = l_adm->GetAttributeCount() - 1; i > 0; --i) {
            if (!l_adm->GetAttributeName(i).starts_with("NAME")) {
                l_adm->RemoveAttribute(i);
            }
        }
    }
    return data_adm;
}

inline ShapeFile* DownloadGlobalAgriPlots(cTMapf* ref, const std::string& token) {
    if (!ref) return nullptr;
    FF::coutI << "Extracting agricultural plots from OSM and ESA WorldCover..." << std::endl;

    ShapeFile* agri_sf = new ShapeFile();
    ShapeLayer* agri_layer = new ShapeLayer();
    agri_layer->SetName("agriplots");
    agri_layer->AddAttribute("type", ATTRIBUTE_TYPE_STRING);
    agri_layer->AddAttribute("source", ATTRIBUTE_TYPE_STRING);

    // 1. Check OSM landuse layer if available
    ShapeFile* data_osm = DownloadGlobalOSM(ref, token);
    if (data_osm) {
        ShapeFile* landuse = GetVectorLayerByName(data_osm, "landuse");
        if (landuse && landuse->GetLayerCount() > 0) {
            ShapeLayer* l = landuse->GetLayer(0);
            int class_idx = l->GetAttributeIndex("class");
            if (class_idx < 0) class_idx = l->GetAttributeIndex("type");
            if (class_idx < 0) class_idx = l->GetAttributeIndex("landuse");

            for (int i = 0; i < l->GetFeatureCount(); ++i) {
                Feature* f = l->GetFeature(i);
                std::string cls = (class_idx >= 0) ? l->GetAttributeValue(class_idx, i) : "";
                std::transform(cls.begin(), cls.end(), cls.begin(), ::tolower);

                bool is_agri = (cls.find("farm") != std::string::npos ||
                                cls.find("agri") != std::string::npos ||
                                cls.find("crop") != std::string::npos ||
                                cls.find("orchard") != std::string::npos ||
                                cls.find("vineyard") != std::string::npos ||
                                cls.find("meadow") != std::string::npos);

                if (is_agri || class_idx < 0) {
                    Feature* f_new = f->Copy();
                    agri_layer->AddFeature(f_new);
                }
            }
            delete landuse;
        }
        delete data_osm;
    }

    // 2. Supplement with ESA WorldCover (class 40 = Cropland) if vector is sparse
    if (agri_layer->GetFeatureCount() == 0) {
        std::vector<cTMapf*> bands_lu = GetAndWarpRTMap(ref, "lu", {12, 11, 10, 9, 8, 7, 6}, token, "nearest", false, 1e-30, "", true);
        if (!bands_lu.empty() && bands_lu[0]) {
            cTMapf* lu = bands_lu[0];
            int step = std::max(1, static_cast<int>(std::round(50.0 / std::abs(lu->cellSizeX()))));
            for (int r = 0; r < lu->Rows(); r += step) {
                for (int c = 0; c < lu->Cols(); c += step) {
                    int val = static_cast<int>(std::round(lu->data[r][c]));
                    if (val == 40) { // 40 = Cropland
                        double x0 = lu->west() + c * lu->cellSizeX();
                        double y0 = lu->north() + r * lu->cellSizeY();
                        double x1 = x0 + std::min(step, lu->Cols() - c) * lu->cellSizeX();
                        double y1 = y0 + std::min(step, lu->Rows() - r) * lu->cellSizeY();

                        Feature* f = new Feature();
                        LSMPolygon* poly = new LSMPolygon();
                        Ring* ring = new Ring();
                        ring->AddVertex(x0, y0);
                        ring->AddVertex(x1, y0);
                        ring->AddVertex(x1, y1);
                        ring->AddVertex(x0, y1);
                        ring->AddVertex(x0, y0);
                        poly->SetExtRing(ring);
                        f->AddShape(poly);
                        agri_layer->AddFeature(f);
                    }
                }
            }
            delete bands_lu[0];
        }
    }

    agri_sf->AddLayer(agri_layer);
    agri_sf->SetProjection(ref->projectionstring());
    FF::coutI << "Total agricultural plots identified: " << agri_layer->GetFeatureCount() << std::endl;
    return agri_sf;
}

// ============================================================================
// 2. POPULATION EXPOSURE MODELING
// ============================================================================

inline PopulationExposureResult ComputePopulationExposure(
    cTMapf* I1,
    ShapeFile* adm_shape,
    double impact_threshold,
    const std::string& token,
    const std::string& popdataout_file = ""
) {
    PopulationExposureResult res;
    if (!I1) return res;

    // Generate a 10x10 gridded shape over the intensity extent
    ShapeFile* sgrid = new ShapeFile();
    ShapeLayer* sgridl = new ShapeLayer();
    sgridl->SetName("grid");
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            Feature* f = new Feature();
            LSMPolygon* p = new LSMPolygon();
            Ring* ring = new Ring();
            ring->AddVertex(I1->west() + (static_cast<float>(i)) * I1->cellSizeX() * (static_cast<float>(I1->Cols())) / 10.0f,
                            I1->north() + (static_cast<float>(j)) * I1->cellSizeY() * (static_cast<float>(I1->Rows())) / 10.0f);
            ring->AddVertex(I1->west() + (static_cast<float>(i + 1)) * I1->cellSizeX() * (static_cast<float>(I1->Cols())) / 10.0f,
                            I1->north() + (static_cast<float>(j)) * I1->cellSizeY() * (static_cast<float>(I1->Rows())) / 10.0f);
            ring->AddVertex(I1->west() + (static_cast<float>(i + 1)) * I1->cellSizeX() * (static_cast<float>(I1->Cols())) / 10.0f,
                            I1->north() + (static_cast<float>(j + 1)) * I1->cellSizeY() * (static_cast<float>(I1->Rows())) / 10.0f);
            ring->AddVertex(I1->west() + (static_cast<float>(i)) * I1->cellSizeX() * (static_cast<float>(I1->Cols())) / 10.0f,
                            I1->north() + (static_cast<float>(j + 1)) * I1->cellSizeY() * (static_cast<float>(I1->Rows())) / 10.0f);
            ring->AddVertex(I1->west() + (static_cast<float>(i)) * I1->cellSizeX() * (static_cast<float>(I1->Cols())) / 10.0f,
                            I1->north() + (static_cast<float>(j)) * I1->cellSizeY() * (static_cast<float>(I1->Rows())) / 10.0f);
            p->SetExtRing(ring);
            f->AddShape(p);
            sgridl->AddFeature(f);
        }
    }
    sgrid->AddLayer(sgridl);
    sgrid->SetProjection(I1->projectionstring());
    res.grid_shape = sgrid;

    // Fetch WorldPop 2025 population data
    std::vector<cTMapf*> bands_pop = GetAndWarpRTMap(I1, "wpop_2025", {7, 8, 9, 10, 11}, token, "bilinear", true, 1e-30);

    float dx = static_cast<float>(I1->CellSizeXApproximateMeters());
    std::shared_ptr<cTMapf> impop = I1->GetCopySP0();

    res.depths = {
        static_cast<float>(1.0 * impact_threshold),
        static_cast<float>(2.0 * impact_threshold),
        static_cast<float>(5.0 * impact_threshold),
        static_cast<float>(10.0 * impact_threshold),
        static_cast<float>(20.0 * impact_threshold),
        static_cast<float>(50.0 * impact_threshold)
    };
    res.count_impacted.assign(res.depths.size(), 0.0);
    res.count_total.assign(res.depths.size(), 0.0);

    if (!bands_pop.empty() && bands_pop[0]) {
        for (int r = 0; r < bands_pop[0]->Rows(); ++r) {
            for (int c = 0; c < bands_pop[0]->Cols(); ++c) {
                float val = bands_pop[0]->data[r][c];
                if (pcr::isMV(val) || !std::isnormal(val) || val < 0.0f) {
                    bands_pop[0]->data[r][c] = 0.0f;
                }
            }
        }

        if (!popdataout_file.empty()) {
            writeRaster(*bands_pop[0], popdataout_file);
        }

        for (size_t i = 0; i < res.depths.size(); ++i) {
            float depthh = res.depths[i];
            float tot_pop_i = 0.0f;
            float tot_pop_i_imp = 0.0f;

            for (int r = 0; r < bands_pop[0]->Rows(); ++r) {
                for (int c = 0; c < bands_pop[0]->Cols(); ++c) {
                    float pop_cell = bands_pop[0]->data[r][c] * dx * dx / (100.0f * 100.0f);
                    tot_pop_i += pop_cell;
                    float imp_val = (I1->data[r][c] > depthh) ? pop_cell : 0.0f;
                    impop->data[r][c] = imp_val;
                    tot_pop_i_imp += imp_val;
                }
            }

            std::string attr_name = "pop_i1_" + std::to_string(i);
            if (adm_shape) {
                RasterizeGetShapePropsTot(adm_shape, attr_name, impop.get());
            }
            if (sgrid) {
                RasterizeGetShapePropsTot(sgrid, attr_name, impop.get());
            }

            res.count_impacted[i] = tot_pop_i_imp;
            res.count_total[i] = tot_pop_i;
        }

        delete bands_pop[0];
    }

    std::cout << "Thresholded impact numbers population: " << std::endl;
    for (size_t i = 0; i < res.depths.size(); ++i) {
        std::cout << "Threshold " << res.depths[i] << " count "
                  << static_cast<long long>(res.count_impacted[i])
                  << " out of total " << static_cast<long long>(res.count_total[i]) << std::endl;
    }

    return res;
}

// ============================================================================
// 3. SHAPE PROP SAMPLING & VECTOR EXPORT
// ============================================================================

inline void SampleIntensityOnVectors(ShapeFile* sf, const std::string& attrib_name, cTMapf* intensity_map) {
    if (sf && intensity_map) {
        RasterizeGetShapeProps(sf, attrib_name, intensity_map);
    }
}

inline std::vector<int> ComputeThresholdCounts(
    ShapeFile* sf,
    const std::string& attrib_name,
    const std::vector<float>& thresholds
) {
    std::vector<int> counts(thresholds.size(), 0);
    if (!sf) return counts;

    sf->ForEachShape([&](ShapeLayer* l, Feature* f, int j, LSMShape* s) {
        int id = l->GetAttributeIndex(attrib_name);
        if (id >= 0) {
            std::string val_str = l->GetAttributeValue(id, j);
            if (!val_str.empty()) {
                try {
                    double d = std::stod(val_str);
                    for (size_t i = 0; i < thresholds.size(); ++i) {
                        if (d > thresholds[i]) {
                            counts[i]++;
                        }
                    }
                } catch (...) {}
            }
        }
    });

    return counts;
}

inline void SaveVectorOutputs(
    ShapeFile* sf,
    const std::string& base_filename,
    const std::string& layer_name,
    bool export_pmtiles,
    int minzoom = 10,
    int maxzoom = 16
) {
    if (!sf || base_filename.empty()) return;

    // Export GPKG
    std::string gpkg_path = base_filename + ".gpkg";
    FF::coutI << "Saving layer '" << layer_name << "' to " << gpkg_path << std::endl;
    SaveVector(sf, gpkg_path, "GPKG", true);

    // Export PMTiles if requested
    if (export_pmtiles) {
        std::string pmtiles_path = base_filename + ".pmtiles";
        FF::coutI << "Saving layer '" << layer_name << "' to PMTiles: " << pmtiles_path << std::endl;
        SaveVectorPMTiles(sf, pmtiles_path, layer_name, minzoom, maxzoom, 50000);
    }
}

} // namespace FastRisk
