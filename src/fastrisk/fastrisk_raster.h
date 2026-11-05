#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <mutex>

#include "../../common/defines.h"
#include "../../common/fferror.h"
#include "../../geospatial/map/map.h"
#include "../../geospatial/shape/shapefile.h"
#include "../../geospatial/shape/boundingbox.h"
#include "../../geospatial/geoprojection.h"
#include "../../web/cogtif/cogtif.h"
#include "../../io/iogdal.h"

namespace FastRisk {

struct RasterStats {
    double min_val = 0.0;
    double max_val = 0.0;
    double mean_val = 0.0;
    int count = 0;
};

// Generic, multi-source raster reader supporting local GeoTIFFs, virtual rasters,
// and remote Cloud-Optimized GeoTIFFs (COG) with intelligent size-based caching
// and spatial patch caching for high-performance exposure and risk modeling.
class FastRiskRasterReader {
public:
    struct CachedPatch {
        double minx = 0.0;
        double maxx = 0.0;
        double miny = 0.0;
        double maxy = 0.0;
        int band = 1;
        std::shared_ptr<cTMapf> raster;

        bool Contains(double x, double y) const {
            return (x >= minx && x <= maxx && y >= miny && y <= maxy);
        }

        bool Encloses(double b_minx, double b_miny, double b_maxx, double b_maxy) const {
            return (b_minx >= minx && b_maxx <= maxx && b_miny >= miny && b_maxy <= maxy);
        }
    };

    FastRiskRasterReader() = default;
    ~FastRiskRasterReader() = default;

    bool Initialize(const std::string& input_path, const std::string& token = "") {
        m_path = input_path;
        m_token = token;
        m_raster = nullptr;
        m_cog = nullptr;
        m_cached_patches.clear();

        // Recognize COG URLs: either explicitly prefixed with "cog:" or standard http/https COG URLs
        if (input_path.rfind("cog:", 0) == 0 || input_path.rfind("http://", 0) == 0 || input_path.rfind("https://", 0) == 0) {
            m_is_cog = true;
            m_cog_url = (input_path.rfind("cog:", 0) == 0) ? input_path.substr(4) : input_path;
            FF::coutI << "Initializing Cloud-Optimized GeoTIFF (COG) reader: " << m_cog_url << std::endl;

            m_cog = std::make_unique<COG::COGTiff>(m_cog_url, m_token);
            if (!m_cog->Initialize()) {
                ff_error = "Failed to initialize COGTiff from url: " + m_cog_url;
                FF::coutE << ff_error << std::endl;
                return false;
            }

            m_width = m_cog->GetWidth();
            m_height = m_cog->GetHeight();
            m_bands = m_cog->GetBands();
            m_projection = m_cog->GetProjection();
            m_cog->GetGeoTransform(m_geoTransform);

            FF::coutI << "COG metadata verified: " << m_width << "x" << m_height
                      << ", bands=" << m_bands << ", dx=" << m_geoTransform[1]
                      << ", dy=" << m_geoTransform[5] << std::endl;
            return true;
        } else {
            m_is_cog = false;
            FF::coutI << "Loading local raster dataset: " << input_path << std::endl;
            cTMapf* loaded = load_map(input_path);
            if (loaded) {
                m_raster = std::shared_ptr<cTMapf>(loaded);
            } else {
                m_raster = nullptr;
            }
            if (!m_raster) {
                ff_error = "Error while reading raster: " + input_path;
                FF::coutE << ff_error << std::endl;
                return false;
            }

            m_width = m_raster->Cols();
            m_height = m_raster->Rows();
            m_bands = 1;
            m_projection = m_raster->projectionstring();
            m_geoTransform[0] = m_raster->west();
            m_geoTransform[1] = m_raster->cellSizeX();
            m_geoTransform[2] = 0.0;
            m_geoTransform[3] = m_raster->north();
            m_geoTransform[4] = 0.0;
            m_geoTransform[5] = m_raster->cellSizeY();

            FF::coutI << "Local raster loaded: " << m_width << "x" << m_height
                      << ", dx=" << m_geoTransform[1] << ", dy=" << m_geoTransform[5] << std::endl;
            return true;
        }
    }

    // Optimization strategy according to user requirements:
    // - If local file: always in memory.
    // - If COG:
    //   - If total elements at risk > 100 OR grid size <= 2500x2500: fetch into memory.
    //   - If < 100 elements at risk and grid size > 2500x2500: on-demand window sampling with spatial cache!
    bool PrepareForExposure(int total_elements_at_risk, const BoundingBox* bbox = nullptr) {
        if (!m_is_cog) {
            return (m_raster != nullptr);
        }

        bool is_small_grid = (m_width <= 2500 && m_height <= 2500);
        bool is_large_portfolio = (total_elements_at_risk > 100);

        if (is_small_grid) {
            FF::coutI << "COG grid size (" << m_width << "x" << m_height
                      << ") <= 2500x2500: pre-fetching full raster into memory..." << std::endl;
            double ulx = m_geoTransform[0];
            double uly = m_geoTransform[3];
            double brx = ulx + m_width * m_geoTransform[1];
            double bry = uly + m_height * m_geoTransform[5];
            m_raster = m_cog->GetRegion(ulx, uly, brx, bry, m_width, m_height, 0, 1);
            return (m_raster != nullptr);
        } else if (is_large_portfolio && bbox != nullptr) {
            double dx = std::abs(m_geoTransform[1]);
            double dy = std::abs(m_geoTransform[5]);
            int out_w = std::min(m_width, static_cast<int>(std::ceil(std::abs(bbox->GetMaxX() - bbox->GetMinX()) / dx)));
            int out_h = std::min(m_height, static_cast<int>(std::ceil(std::abs(bbox->GetMaxY() - bbox->GetMinY()) / dy)));

            if (out_w <= 3500 && out_h <= 3500) {
                FF::coutI << "Large portfolio (" << total_elements_at_risk
                          << " elements) within extent (" << out_w << "x" << out_h
                          << "): pre-fetching bounding box region..." << std::endl;
                m_raster = m_cog->GetRegion(bbox->GetMinX(), bbox->GetMaxY(), bbox->GetMaxX(), bbox->GetMinY(), out_w, out_h, -1, 1);
                return (m_raster != nullptr);
            }
        }

        FF::coutI << "Selective element count (" << total_elements_at_risk
                  << ") on large grid (" << m_width << "x" << m_height
                  << "): activating fast on-demand spatial tile cache." << std::endl;
        return true;
    }

    std::shared_ptr<cTMapf> GetRaster(int band = 1) {
        if (!m_raster && m_is_cog && m_cog) {
            double ulx = m_geoTransform[0];
            double uly = m_geoTransform[3];
            double brx = ulx + m_width * m_geoTransform[1];
            double bry = uly + m_height * m_geoTransform[5];
            int req_w = std::min(m_width, 2500);
            int req_h = std::min(m_height, 2500);
            m_raster = m_cog->GetRegion(ulx, uly, brx, bry, req_w, req_h, -1, band);
        }
        return m_raster;
    }

    bool HasFullRaster() const {
        return (m_raster != nullptr);
    }

    int GetWidth() const { return m_width; }
    int GetHeight() const { return m_height; }
    int GetBands() const { return m_bands; }
    double GetDx() const { return m_geoTransform[1]; }
    double GetDy() const { return m_geoTransform[5]; }
    double GetOriginX() const { return m_geoTransform[0]; }
    double GetOriginY() const { return m_geoTransform[3]; }
    std::string GetProjection() const { return m_projection; }
    void GetGeoTransform(double* gt) const { std::memcpy(gt, m_geoTransform, 6 * sizeof(double)); }

    // Coordinate conversions
    void WorldToPixel(double wx, double wy, int& col, int& row) const {
        col = static_cast<int>(std::floor((wx - m_geoTransform[0]) / m_geoTransform[1]));
        row = static_cast<int>(std::floor((wy - m_geoTransform[3]) / m_geoTransform[5]));
    }

    void PixelToWorld(int col, int row, double& wx, double& wy) const {
        wx = m_geoTransform[0] + (col + 0.5) * m_geoTransform[1];
        wy = m_geoTransform[3] + (row + 0.5) * m_geoTransform[5];
    }

    // Point sampling
    double SamplePoint(double x, double y, int band = 1) {
        if (m_raster) {
            float val = m_raster->Sample(x, y);
            return (pcr::isMV(val) || std::isnan(val) || val <= -9999.0f) ? 0.0 : static_cast<double>(val);
        }

        // Check spatial patch cache first
        std::shared_ptr<cTMapf> patch = FindOrCreatePatch(x, y, x, y, 100.0, band);
        if (patch) {
            float val = patch->Sample(x, y);
            return (pcr::isMV(val) || std::isnan(val) || val <= -9999.0f) ? 0.0 : static_cast<double>(val);
        }

        if (m_cog) {
            auto pix = m_cog->GetPixel(x, y, band);
            if (pix && pix->Rows() > 0 && pix->Cols() > 0) {
                float val = pix->data[0][0];
                return (pcr::isMV(val) || std::isnan(val) || val <= -9999.0f) ? 0.0 : static_cast<double>(val);
            }
        }
        return 0.0;
    }

    // Window sampling with full summary statistics
    RasterStats SampleWindow(double x, double y, double radius_m, int band = 1) {
        RasterStats stats;

        cTMapf* target_map = m_raster.get();
        std::shared_ptr<cTMapf> patch = nullptr;

        if (!target_map) {
            patch = FindOrCreatePatch(x - radius_m, y - radius_m, x + radius_m, y + radius_m, radius_m * 2.0, band);
            target_map = patch.get();
        }

        if (target_map) {
            double dx = std::abs(target_map->cellSizeX());
            double dy = std::abs(target_map->cellSizeY());
            int rx = std::max(1, static_cast<int>(std::ceil(radius_m / dx)));
            int ry = std::max(1, static_cast<int>(std::ceil(radius_m / dy)));

            float x_pix = (x - target_map->west()) / target_map->cellSizeX();
            float y_pix = (y - target_map->north()) / target_map->cellSizeY();
            int c_center = static_cast<int>(std::floor(x_pix));
            int r_center = static_cast<int>(std::floor(y_pix));

            double sum = 0.0;
            bool first = true;

            for (int r = std::max(0, r_center - ry); r <= std::min(target_map->Rows() - 1, r_center + ry); ++r) {
                for (int c = std::max(0, c_center - rx); c <= std::min(target_map->Cols() - 1, c_center + rx); ++c) {
                    float v = target_map->data[r][c];
                    if (!pcr::isMV(v) && !std::isnan(v) && v > -9999.0f) {
                        double vd = static_cast<double>(v);
                        if (first) {
                            stats.min_val = vd;
                            stats.max_val = vd;
                            first = false;
                        } else {
                            if (vd < stats.min_val) stats.min_val = vd;
                            if (vd > stats.max_val) stats.max_val = vd;
                        }
                        sum += vd;
                        stats.count++;
                    }
                }
            }
            if (stats.count > 0) {
                stats.mean_val = sum / stats.count;
            }
        }

        return stats;
    }

    double SampleWindowMax(double x, double y, double radius_m, int band = 1) {
        return SampleWindow(x, y, radius_m, band).max_val;
    }

    double SampleWindowMean(double x, double y, double radius_m, int band = 1) {
        return SampleWindow(x, y, radius_m, band).mean_val;
    }

    // Geometry / Shape sampling
    double SampleShape(LSMShape* shape, GeoProjection shape_proj, const std::string& stat_type = "max", int band = 1) {
        if (!shape) return 0.0;

        BoundingBox bbox = shape->GetBoundingBox();
        GeoProjection raster_proj = GeoProjection::FromString(m_projection);
        if (shape_proj.IsValid() && raster_proj.IsValid() && shape_proj.GetEPSG() != raster_proj.GetEPSG()) {
            GeoCoordTransformer* tr = GeoCoordTransformer::GetCoordTransformer(shape_proj, raster_proj);
            if (tr) {
                LSMVector2 pmin = tr->Transform(LSMVector2(bbox.GetMinX(), bbox.GetMinY()));
                LSMVector2 pmax = tr->Transform(LSMVector2(bbox.GetMaxX(), bbox.GetMaxY()));
                bbox = BoundingBox(std::min(pmin.x, pmax.x), std::max(pmin.x, pmax.x),
                                   std::min(pmin.y, pmax.y), std::max(pmin.y, pmax.y));
                delete tr;
            }
        }

        cTMapf* target_map = m_raster.get();
        std::shared_ptr<cTMapf> patch = nullptr;

        if (!target_map) {
            patch = FindOrCreatePatch(bbox.GetMinX(), bbox.GetMinY(), bbox.GetMaxX(), bbox.GetMaxY(), 100.0, band);
            target_map = patch.get();
        }

        if (target_map) {
            int c1 = static_cast<int>(std::floor((bbox.GetMinX() - target_map->west()) / target_map->cellSizeX()));
            int c2 = static_cast<int>(std::ceil((bbox.GetMaxX() - target_map->west()) / target_map->cellSizeX()));
            int r1 = static_cast<int>(std::floor((bbox.GetMaxY() - target_map->north()) / target_map->cellSizeY()));
            int r2 = static_cast<int>(std::ceil((bbox.GetMinY() - target_map->north()) / target_map->cellSizeY()));

            int r_start = std::max(0, std::min(r1, r2));
            int r_end = std::min(target_map->Rows() - 1, std::max(r1, r2));
            int c_start = std::max(0, std::min(c1, c2));
            int c_end = std::min(target_map->Cols() - 1, std::max(c1, c2));

            double max_val = 0.0;
            double sum_val = 0.0;
            int count = 0;

            for (int r = r_start; r <= r_end; ++r) {
                for (int c = c_start; c <= c_end; ++c) {
                    float v = target_map->data[r][c];
                    if (!pcr::isMV(v) && !std::isnan(v) && v > -9999.0f) {
                        double vd = static_cast<double>(v);
                        if (vd > max_val) max_val = vd;
                        sum_val += vd;
                        count++;
                    }
                }
            }

            if (stat_type == "mean" && count > 0) {
                return sum_val / count;
            }
            return max_val;
        }

        return 0.0;
    }

    double SampleShapeMax(LSMShape* shape, GeoProjection shape_proj, int band = 1) {
        return SampleShape(shape, shape_proj, "max", band);
    }

    double SampleShapeMean(LSMShape* shape, GeoProjection shape_proj, int band = 1) {
        return SampleShape(shape, shape_proj, "mean", band);
    }

    // Extract an arbitrary region
    std::shared_ptr<cTMapf> GetRegion(double ulx, double uly, double brx, double bry, int out_w, int out_h, int band = 1) {
        if (m_raster) {
            MaskedRaster<float> raster_data(out_h, out_w, uly, ulx, (brx - ulx) / out_w, (bry - uly) / out_h);
            for (int r = 0; r < out_h; ++r) {
                for (int c = 0; c < out_w; ++c) {
                    double wx = ulx + (c + 0.5) * ((brx - ulx) / out_w);
                    double wy = uly + (r + 0.5) * ((bry - uly) / out_h);
                    raster_data[r][c] = m_raster->Sample(wx, wy);
                }
            }
            return std::make_shared<cTMapf>(std::move(raster_data), m_projection, "RasterRegion");
        }
        if (m_cog) {
            return m_cog->GetRegion(ulx, uly, brx, bry, out_w, out_h, -1, band);
        }
        return nullptr;
    }

private:
    std::shared_ptr<cTMapf> FindOrCreatePatch(double minx, double miny, double maxx, double maxy, double margin_m, int band) {
        std::lock_guard<std::mutex> lock(m_cache_mutex);

        // Check existing cached patches
        for (const auto& patch : m_cached_patches) {
            if (patch.band == band && patch.Encloses(minx, miny, maxx, maxy)) {
                return patch.raster;
            }
        }

        // Not found, fetch new patch with buffer margin
        if (!m_cog) return nullptr;

        double dx = std::abs(m_geoTransform[1]);
        double dy = std::abs(m_geoTransform[5]);
        double buffer = std::max(margin_m, dx * 30.0); // at least 30 pixels buffer

        double fetch_minx = minx - buffer;
        double fetch_maxx = maxx + buffer;
        double fetch_miny = miny - buffer;
        double fetch_maxy = maxy + buffer;

        int out_w = std::clamp(static_cast<int>(std::ceil((fetch_maxx - fetch_minx) / dx)), 16, 1024);
        int out_h = std::clamp(static_cast<int>(std::ceil((fetch_maxy - fetch_miny) / dy)), 16, 1024);

        auto new_patch = m_cog->GetRegion(fetch_minx, fetch_maxy, fetch_maxx, fetch_miny, out_w, out_h, -1, band);
        if (new_patch) {
            // Keep cache size bounded (max 16 regional patches)
            if (m_cached_patches.size() >= 16) {
                m_cached_patches.erase(m_cached_patches.begin());
            }
            CachedPatch cp;
            cp.minx = fetch_minx;
            cp.maxx = fetch_maxx;
            cp.miny = fetch_miny;
            cp.maxy = fetch_maxy;
            cp.band = band;
            cp.raster = new_patch;
            m_cached_patches.push_back(cp);
            return new_patch;
        }

        return nullptr;
    }

    std::string m_path;
    std::string m_token;
    bool m_is_cog = false;
    std::string m_cog_url;

    std::unique_ptr<COG::COGTiff> m_cog;
    std::shared_ptr<cTMapf> m_raster;
    std::vector<CachedPatch> m_cached_patches;
    std::mutex m_cache_mutex;

    int m_width = 0;
    int m_height = 0;
    int m_bands = 1;
    double m_geoTransform[6] = {0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    std::string m_projection;
};

} // namespace FastRisk
