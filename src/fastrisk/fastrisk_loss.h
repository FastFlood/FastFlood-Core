#pragma once

#include <string>
#include <vector>
#include <memory>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "../../common/defines.h"
#include "../../geospatial/map/map.h"
#include "../../geospatial/shape/shapefile.h"
#include "../../geospatial/shape/shapelayer.h"
#include "../../geospatial/shape/feature.h"
#include "fastrisk_vulnerability.h"
#include "fastrisk_raster.h"

namespace FastRisk {

struct SinglePropertyLossResult {
    double structure_area_m2 = 0.0;
    double grounds_area_m2 = 0.0;
    double structure_rcv = 0.0;
    double grounds_rcv = 0.0;
    double total_property_rcv = 0.0;

    double depth_at_grade_m = 0.0;
    double first_floor_height_m = 0.3;
    double depth_rel_ffe_m = 0.0;
    double grounds_depth_m = 0.0;

    double structural_damage_pct = 0.0;
    double grounds_damage_pct = 0.0;
    double total_damage_pct = 0.0;

    double structural_loss = 0.0;
    double grounds_loss = 0.0;
    double total_loss = 0.0;

    std::string access_status = "Normal";
    std::string access_notes = "Unrestricted vehicle and pedestrian access";
    std::string mon_unit = "EUR";
};

struct PortfolioRiskSummary {
    std::string hazard_type = "flood";
    std::string vuln_curve_name = "JRC";
    std::string mon_unit = "EUR";

    int total_buildings = 0;
    int impacted_buildings = 0;
    double buildings_rcv = 0.0;
    double buildings_loss = 0.0;

    int total_roads = 0;
    int impacted_roads = 0;
    double roads_rcv = 0.0;
    double roads_loss = 0.0;

    int total_agriplots = 0;
    int impacted_agriplots = 0;
    double agriplots_rcv = 0.0;
    double agriplots_loss = 0.0;

    double total_loss = 0.0;
    double total_rcv = 0.0;

    void PrintSummary() const {
        std::cout << "\n========================================================\n";
        std::cout << "               FASTRISK ASSESSMENT SUMMARY\n";
        std::cout << "========================================================\n";
        std::cout << "Hazard Peril:           " << hazard_type << "\n";
        std::cout << "Vulnerability Model:    " << vuln_curve_name << "\n";
        std::cout << "Currency:               " << mon_unit << "\n\n";

        if (total_buildings > 0) {
            double dmg_pct = (buildings_rcv > 0) ? (buildings_loss / buildings_rcv) * 100.0 : 0.0;
            std::cout << "BUILDINGS:\n";
            std::cout << "  Total Count:          " << total_buildings << "\n";
            std::cout << "  Impacted Count:       " << impacted_buildings << "\n";
            std::cout << "  Replacement Value:    " << mon_unit << " " << std::fixed << std::setprecision(2) << buildings_rcv << "\n";
            std::cout << "  Estimated Loss:       " << mon_unit << " " << buildings_loss << " (" << dmg_pct << "%)\n\n";
        }

        if (total_roads > 0) {
            double dmg_pct = (roads_rcv > 0) ? (roads_loss / roads_rcv) * 100.0 : 0.0;
            std::cout << "ROADS & TRANSPORTATION:\n";
            std::cout << "  Total Count:          " << total_roads << "\n";
            std::cout << "  Impacted Count:       " << impacted_roads << "\n";
            std::cout << "  Replacement Value:    " << mon_unit << " " << roads_rcv << "\n";
            std::cout << "  Estimated Loss:       " << mon_unit << " " << roads_loss << " (" << dmg_pct << "%)\n\n";
        }

        if (total_agriplots > 0) {
            double dmg_pct = (agriplots_rcv > 0) ? (agriplots_loss / agriplots_rcv) * 100.0 : 0.0;
            std::cout << "AGRICULTURE & CROPLAND:\n";
            std::cout << "  Total Count:          " << total_agriplots << "\n";
            std::cout << "  Impacted Count:       " << impacted_agriplots << "\n";
            std::cout << "  Crop/Land Value:      " << mon_unit << " " << agriplots_rcv << "\n";
            std::cout << "  Estimated Loss:       " << mon_unit << " " << agriplots_loss << " (" << dmg_pct << "%)\n\n";
        }

        double total_dmg_pct = (total_rcv > 0) ? (total_loss / total_rcv) * 100.0 : 0.0;
        std::cout << "--------------------------------------------------------\n";
        std::cout << "TOTAL EXPOSURE VALUE:   " << mon_unit << " " << total_rcv << "\n";
        std::cout << "TOTAL ESTIMATED LOSS:   " << mon_unit << " " << total_loss << " (" << total_dmg_pct << "%)\n";
        std::cout << "========================================================\n\n";
    }

    std::string ToJSON() const {
        std::ostringstream ss;
        ss << "{\n";
        ss << "  \"hazard_type\": \"" << hazard_type << "\",\n";
        ss << "  \"vuln_curve\": \"" << vuln_curve_name << "\",\n";
        ss << "  \"mon_unit\": \"" << mon_unit << "\",\n";
        ss << "  \"buildings\": {\n";
        ss << "    \"total_count\": " << total_buildings << ",\n";
        ss << "    \"impacted_count\": " << impacted_buildings << ",\n";
        ss << "    \"total_rcv\": " << buildings_rcv << ",\n";
        ss << "    \"total_loss\": " << buildings_loss << "\n";
        ss << "  },\n";
        ss << "  \"roads\": {\n";
        ss << "    \"total_count\": " << total_roads << ",\n";
        ss << "    \"impacted_count\": " << impacted_roads << ",\n";
        ss << "    \"total_rcv\": " << roads_rcv << ",\n";
        ss << "    \"total_loss\": " << roads_loss << "\n";
        ss << "  },\n";
        ss << "  \"agriculture\": {\n";
        ss << "    \"total_count\": " << total_agriplots << ",\n";
        ss << "    \"impacted_count\": " << impacted_agriplots << ",\n";
        ss << "    \"total_rcv\": " << agriplots_rcv << ",\n";
        ss << "    \"total_loss\": " << agriplots_loss << "\n";
        ss << "  },\n";
        ss << "  \"total_rcv\": " << total_rcv << ",\n";
        ss << "  \"total_loss\": " << total_loss << "\n";
        ss << "}";
        return ss.str();
    }
};

// ============================================================================
// ACCESSIBILITY & GEOMETRY HELPERS
// ============================================================================

inline std::pair<std::string, std::string> EvaluateAccessibility(double intensity, const std::string& peril = "flood") {
    std::string p = peril;
    std::transform(p.begin(), p.end(), p.begin(), ::tolower);

    if (p.find("wind") != std::string::npos) {
        // Wind speed at 2m (m/s)
        if (intensity < 17.0) { // < Beaufort 8
            return {"Normal", "Normal vehicle travel and pedestrian conditions"};
        } else if (intensity < 24.5) { // Beaufort 8-9 (Severe gale)
            return {"Caution", "High crosswinds; caution for high-sided vehicles"};
        } else if (intensity < 32.7) { // Beaufort 10-11 (Storm / Violent storm)
            return {"Impaired", "Severe gale/storm winds; debris hazard, access restricted"};
        } else { // >= 32.7 (Hurricane force)
            return {"Severed", "Hurricane force winds; impassable roadway, extreme hazard"};
        }
    } else if (p.find("heat") != std::string::npos || p.find("temp") != std::string::npos) {
        // Max ambient temperature (°C)
        if (intensity < 38.0) {
            return {"Normal", "Normal driving and road surface conditions"};
        } else if (intensity < 43.0) {
            return {"Caution", "Elevated temperatures; vehicle cooling stress"};
        } else if (intensity < 48.0) {
            return {"Impaired", "Extreme heat; road asphalt softening, heavy vehicle restrictions"};
        } else {
            return {"Severed", "Severe heat emergency; asphalt bleeding / pavement rutting / blowout hazard"};
        }
    } else if (p.find("slide") != std::string::npos || p.find("debris") != std::string::npos || p.find("mud") != std::string::npos) {
        // Landslide / debris flow deposit depth (m)
        if (intensity < 0.05) {
            return {"Normal", "Unrestricted roadway access"};
        } else if (intensity < 0.15) {
            return {"Caution", "Minor mud/debris accumulation; reduce speed"};
        } else if (intensity < 0.40) {
            return {"Impaired", "Moderate slide debris; low-clearance vehicles blocked"};
        } else {
            return {"Severed", "Heavy slide deposit (>0.4m); road blocked / structurally severed"};
        }
    } else {
        // Flood depth (m)
        if (intensity < 0.08) { // < ~3 inches / 0.25 ft
            return {"Normal", "Unrestricted vehicle and pedestrian access"};
        } else if (intensity < 0.15) { // 3 - 6 inches / 0.25 - 0.50 ft
            return {"Caution", "Minor standing water (<6 in); passenger vehicle caution"};
        } else if (intensity < 0.30) { // 6 - 12 inches / 0.50 - 1.00 ft
            return {"Impaired", "Moderate flood (6-12 in); low-clearance vehicles stalled, access restricted"};
        } else { // > 12 inches
            return {"Severed", "Severe flood (>12 in); impassable roadway, vehicle washout hazard"};
        }
    }
}

inline double CalculateFeatureAreaM2(Feature* f, GeoProjection proj, double center_lat = 0.0) {
    if (!f) return 0.0;
    double raw_area = f->Area();
    if (raw_area <= 0.0) return 0.0;

    int epsg = proj.GetEPSG();
    if (epsg == 3857) {
        // Web Mercator area scale correction: cos(lat)^2
        double lat_rad = center_lat * (M_PI / 180.0);
        double cos_val = std::cos(lat_rad);
        return raw_area * (cos_val * cos_val);
    } else if (epsg == 4326 || raw_area < 1e-4) {
        // Latitude/longitude in degrees: ~111,139 m per degree
        double lat_rad = center_lat * (M_PI / 180.0);
        double m_per_deg_lat = 111139.0;
        double m_per_deg_lon = 111139.0 * std::cos(lat_rad);
        return raw_area * m_per_deg_lat * m_per_deg_lon;
    }
    // Projected CRS already in meters (e.g. UTM)
    return raw_area;
}

inline double CalculateFeatureLengthM(Feature* f, GeoProjection proj, double center_lat = 0.0) {
    if (!f) return 0.0;
    double raw_len = f->Length();
    if (raw_len <= 0.0) return 0.0;

    int epsg = proj.GetEPSG();
    if (epsg == 3857) {
        double lat_rad = center_lat * (M_PI / 180.0);
        return raw_len * std::cos(lat_rad);
    } else if (epsg == 4326 || raw_len < 1.0) {
        return raw_len * 111139.0;
    }
    return raw_len;
}

// ============================================================================
// 1. SINGLE PROPERTY LOSS ESTIMATION
// ============================================================================

inline SinglePropertyLossResult ComputeSinglePropertyLoss(
    ShapeFile* structure_sf,
    ShapeFile* grounds_sf,
    FastRiskRasterReader& hazard_reader,
    const VulnerabilityCurve& bldg_vcurve,
    const VulnerabilityCurve& grounds_vcurve,
    double bldg_val_per_m2,
    double grounds_val_per_m2,
    double ffh_m,
    const std::string& mon_unit,
    const std::string& hazard_type = "flood"
) {
    SinglePropertyLossResult res;
    res.first_floor_height_m = ffh_m;
    res.mon_unit = mon_unit;

    std::string htype = hazard_type;
    std::transform(htype.begin(), htype.end(), htype.begin(), ::tolower);
    bool is_flood = (htype.find("heat") == std::string::npos &&
                     htype.find("temp") == std::string::npos &&
                     htype.find("wind") == std::string::npos);

    // 1. Structure Assessment
    if (structure_sf && structure_sf->GetLayerCount() > 0) {
        ShapeLayer* l_bldg = structure_sf->GetLayer(0);
        if (l_bldg->GetFeatureCount() > 0) {
            Feature* f_bldg = l_bldg->GetFeature(0);
            BoundingBox bbox = f_bldg->GetBoundingBox();
            double center_lat = (bbox.GetMinY() + bbox.GetMaxY()) / 2.0;
            if (std::abs(center_lat) > 90.0) {
                // If projected coordinates (like 3857), approximate latitude
                center_lat = 42.0;
            }

            res.structure_area_m2 = CalculateFeatureAreaM2(f_bldg, structure_sf->GetProjection(), center_lat);
            if (res.structure_area_m2 <= 0.0) {
                res.structure_area_m2 = 180.0; // fallback ~1900 sqft
            }
            res.structure_rcv = res.structure_area_m2 * bldg_val_per_m2;

            // Sample hazard intensity on structure footprint
            if (f_bldg->GetShapeCount() > 0) {
                res.depth_at_grade_m = hazard_reader.SampleShapeMax(f_bldg->GetShape(0), structure_sf->GetProjection());
            } else {
                res.depth_at_grade_m = hazard_reader.SamplePoint((bbox.GetMinX() + bbox.GetMaxX()) / 2.0, (bbox.GetMinY() + bbox.GetMaxY()) / 2.0);
            }

            if (is_flood) {
                if (res.depth_at_grade_m >= 0.01) {
                    res.depth_rel_ffe_m = res.depth_at_grade_m - res.first_floor_height_m;

                    // Evaluate building damage
                    double living_dmg_pct = bldg_vcurve.EvaluateDamagePct(res.depth_rel_ffe_m);
                    VulnerabilityCurve fnd_curve = GetHAZUSFoundationCurve();
                    double foundation_dmg_pct = fnd_curve.EvaluateDamagePct(res.depth_at_grade_m);

                    res.structural_damage_pct = std::max(foundation_dmg_pct, living_dmg_pct);
                    res.structural_loss = (res.structural_damage_pct / 100.0) * res.structure_rcv;
                }
            } else {
                res.depth_rel_ffe_m = res.depth_at_grade_m;
                res.structural_damage_pct = bldg_vcurve.EvaluateDamagePct(res.depth_at_grade_m);
                res.structural_loss = (res.structural_damage_pct / 100.0) * res.structure_rcv;
            }

            // Set attributes on structure layer
            l_bldg->AddAttribute("depth_m", ATTRIBUTE_TYPE_FLOAT64);
            l_bldg->AddAttribute("ffh_m", ATTRIBUTE_TYPE_FLOAT64);
            l_bldg->AddAttribute("depth_ffe_m", ATTRIBUTE_TYPE_FLOAT64);
            l_bldg->AddAttribute("dmg_pct", ATTRIBUTE_TYPE_FLOAT64);
            l_bldg->AddAttribute("loss", ATTRIBUTE_TYPE_FLOAT64);
            l_bldg->AddAttribute("rcv", ATTRIBUTE_TYPE_FLOAT64);
            l_bldg->AddAttribute("area_m2", ATTRIBUTE_TYPE_FLOAT64);
            l_bldg->AddAttribute("hazard", ATTRIBUTE_TYPE_STRING);
            l_bldg->AddAttribute("mon_unit", ATTRIBUTE_TYPE_STRING);
            l_bldg->AddAttribute("vuln_curve", ATTRIBUTE_TYPE_STRING);

            l_bldg->SetFeatureAttribute(0, "depth_m", res.depth_at_grade_m);
            l_bldg->SetFeatureAttribute(0, "ffh_m", res.first_floor_height_m);
            l_bldg->SetFeatureAttribute(0, "depth_ffe_m", res.depth_rel_ffe_m);
            l_bldg->SetFeatureAttribute(0, "dmg_pct", res.structural_damage_pct);
            l_bldg->SetFeatureAttribute(0, "loss", res.structural_loss);
            l_bldg->SetFeatureAttribute(0, "rcv", res.structure_rcv);
            l_bldg->SetFeatureAttribute(0, "area_m2", res.structure_area_m2);
            l_bldg->SetFeatureAttribute(0, "hazard", hazard_type);
            l_bldg->SetFeatureAttribute(0, "mon_unit", res.mon_unit);
            l_bldg->SetFeatureAttribute(0, "vuln_curve", bldg_vcurve.name);
        }
    }

    // 2. Grounds / Site Assessment
    if (grounds_sf && grounds_sf->GetLayerCount() > 0) {
        ShapeLayer* l_grounds = grounds_sf->GetLayer(0);
        if (l_grounds->GetFeatureCount() > 0) {
            Feature* f_grounds = l_grounds->GetFeature(0);
            BoundingBox bbox = f_grounds->GetBoundingBox();
            double center_lat = (bbox.GetMinY() + bbox.GetMaxY()) / 2.0;
            if (std::abs(center_lat) > 90.0) center_lat = 42.0;

            res.grounds_area_m2 = CalculateFeatureAreaM2(f_grounds, grounds_sf->GetProjection(), center_lat);
            if (res.grounds_area_m2 <= 0.0) {
                res.grounds_area_m2 = res.structure_area_m2 * 3.0; // fallback parcel
            }
            res.grounds_rcv = res.grounds_area_m2 * grounds_val_per_m2;

            if (f_grounds->GetShapeCount() > 0) {
                res.grounds_depth_m = hazard_reader.SampleShapeMax(f_grounds->GetShape(0), grounds_sf->GetProjection());
            } else {
                res.grounds_depth_m = hazard_reader.SamplePoint((bbox.GetMinX() + bbox.GetMaxX()) / 2.0, (bbox.GetMinY() + bbox.GetMaxY()) / 2.0);
            }

            auto [stat, notes] = EvaluateAccessibility(res.grounds_depth_m, hazard_type);
            res.access_status = stat;
            res.access_notes = notes;

            res.grounds_damage_pct = grounds_vcurve.EvaluateDamagePct(res.grounds_depth_m);
            // Continuous grounds damage is defined as percentage of grounds value (or building RCV fallback)
            double basis_val = (res.grounds_rcv > 0) ? res.grounds_rcv : res.structure_rcv;
            res.grounds_loss = (res.grounds_damage_pct / 100.0) * basis_val;

            // Set attributes on grounds layer
            l_grounds->AddAttribute("site_depth_m", ATTRIBUTE_TYPE_FLOAT64);
            l_grounds->AddAttribute("dmg_pct", ATTRIBUTE_TYPE_FLOAT64);
            l_grounds->AddAttribute("loss", ATTRIBUTE_TYPE_FLOAT64);
            l_grounds->AddAttribute("grounds_rcv", ATTRIBUTE_TYPE_FLOAT64);
            l_grounds->AddAttribute("area_m2", ATTRIBUTE_TYPE_FLOAT64);
            l_grounds->AddAttribute("hazard", ATTRIBUTE_TYPE_STRING);
            l_grounds->AddAttribute("access_stat", ATTRIBUTE_TYPE_STRING);
            l_grounds->AddAttribute("access_note", ATTRIBUTE_TYPE_STRING);
            l_grounds->AddAttribute("mon_unit", ATTRIBUTE_TYPE_STRING);

            l_grounds->SetFeatureAttribute(0, "site_depth_m", res.grounds_depth_m);
            l_grounds->SetFeatureAttribute(0, "dmg_pct", res.grounds_damage_pct);
            l_grounds->SetFeatureAttribute(0, "loss", res.grounds_loss);
            l_grounds->SetFeatureAttribute(0, "grounds_rcv", res.grounds_rcv);
            l_grounds->SetFeatureAttribute(0, "area_m2", res.grounds_area_m2);
            l_grounds->SetFeatureAttribute(0, "hazard", hazard_type);
            l_grounds->SetFeatureAttribute(0, "access_stat", res.access_status);
            l_grounds->SetFeatureAttribute(0, "access_note", res.access_notes);
            l_grounds->SetFeatureAttribute(0, "mon_unit", res.mon_unit);
        }
    }

    res.total_property_rcv = res.structure_rcv + res.grounds_rcv;
    res.total_loss = res.structural_loss + res.grounds_loss;
    res.total_damage_pct = (res.total_property_rcv > 0) ? (res.total_loss / res.total_property_rcv) * 100.0 : 0.0;

    return res;
}

// ============================================================================
// 2. PORTFOLIO / TOTAL AREA LOSS ESTIMATION
// ============================================================================

inline void ComputeBuildingsLoss(
    ShapeFile* buildings_sf,
    const VulnerabilityCurve& vcurve,
    double bldg_val_per_m2,
    double ffh_m,
    const std::string& mon_unit,
    PortfolioRiskSummary& summary,
    const std::string& hazard_type = "flood"
) {
    if (!buildings_sf) return;

    std::string htype = hazard_type;
    std::transform(htype.begin(), htype.end(), htype.begin(), ::tolower);
    bool is_flood = (htype.find("heat") == std::string::npos &&
                     htype.find("temp") == std::string::npos &&
                     htype.find("wind") == std::string::npos);

    for (int l = 0; l < buildings_sf->GetLayerCount(); ++l) {
        ShapeLayer* layer = buildings_sf->GetLayer(l);

        layer->AddAttribute("depth_m", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("depth_ffe_m", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("dmg_ratio", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("dmg_pct", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("loss", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("rcv", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("hazard", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("mon_unit", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("vuln_curve", ATTRIBUTE_TYPE_STRING);

        int i1_idx = layer->GetAttributeIndex("i1");

        for (int i = 0; i < layer->GetFeatureCount(); ++i) {
            Feature* f = layer->GetFeature(i);
            summary.total_buildings++;

            double depth_m = 0.0;
            if (i1_idx >= 0) {
                std::string s_val = layer->GetAttributeValue(i1_idx, i);
                if (!s_val.empty()) {
                    try { depth_m = std::stod(s_val); } catch (...) {}
                }
            }
            if (depth_m < 0.0 || std::isnan(depth_m)) depth_m = 0.0;

            BoundingBox bbox = f->GetBoundingBox();
            double center_lat = (bbox.GetMinY() + bbox.GetMaxY()) / 2.0;
            if (std::abs(center_lat) > 90.0) center_lat = 42.0;
            double area_m2 = CalculateFeatureAreaM2(f, buildings_sf->GetProjection(), center_lat);
            if (area_m2 <= 0.0) area_m2 = 150.0; // fallback area

            double rcv = area_m2 * bldg_val_per_m2;
            summary.buildings_rcv += rcv;

            double depth_ffe_m = is_flood ? (depth_m - ffh_m) : depth_m;
            double dmg_ratio = 0.0;
            double loss = 0.0;

            if (is_flood) {
                if (depth_m >= 0.01) {
                    summary.impacted_buildings++;
                    dmg_ratio = vcurve.EvaluateDamageRatio(std::max(0.0, depth_ffe_m));
                    loss = dmg_ratio * rcv;
                    summary.buildings_loss += loss;
                }
            } else {
                dmg_ratio = vcurve.EvaluateDamageRatio(depth_m);
                if (dmg_ratio > 0.0) {
                    summary.impacted_buildings++;
                    loss = dmg_ratio * rcv;
                    summary.buildings_loss += loss;
                }
            }

            layer->SetFeatureAttribute(i, "depth_m", depth_m);
            layer->SetFeatureAttribute(i, "depth_ffe_m", depth_ffe_m);
            layer->SetFeatureAttribute(i, "dmg_ratio", dmg_ratio);
            layer->SetFeatureAttribute(i, "dmg_pct", dmg_ratio * 100.0);
            layer->SetFeatureAttribute(i, "loss", loss);
            layer->SetFeatureAttribute(i, "rcv", rcv);
            layer->SetFeatureAttribute(i, "hazard", hazard_type);
            layer->SetFeatureAttribute(i, "mon_unit", mon_unit);
            layer->SetFeatureAttribute(i, "vuln_curve", vcurve.name);
        }
    }
}

inline void ComputeRoadsLoss(
    ShapeFile* roads_sf,
    const VulnerabilityCurve& vcurve,
    double road_val_per_m2,
    const std::string& mon_unit,
    PortfolioRiskSummary& summary,
    const std::string& hazard_type = "flood"
) {
    if (!roads_sf) return;

    for (int l = 0; l < roads_sf->GetLayerCount(); ++l) {
        ShapeLayer* layer = roads_sf->GetLayer(l);

        layer->AddAttribute("depth_m", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("dmg_ratio", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("dmg_pct", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("loss", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("asset_val", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("hazard", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("access_stat", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("access_note", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("mon_unit", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("vuln_curve", ATTRIBUTE_TYPE_STRING);

        int i1_idx = layer->GetAttributeIndex("i1");

        for (int i = 0; i < layer->GetFeatureCount(); ++i) {
            Feature* f = layer->GetFeature(i);
            summary.total_roads++;

            double depth_m = 0.0;
            if (i1_idx >= 0) {
                std::string s_val = layer->GetAttributeValue(i1_idx, i);
                if (!s_val.empty()) {
                    try { depth_m = std::stod(s_val); } catch (...) {}
                }
            }
            if (depth_m < 0.0 || std::isnan(depth_m)) depth_m = 0.0;

            BoundingBox bbox = f->GetBoundingBox();
            double center_lat = (bbox.GetMinY() + bbox.GetMaxY()) / 2.0;
            if (std::abs(center_lat) > 90.0) center_lat = 42.0;

            // Width of polylines is assumed 25 meters as per user specification!
            double area_m2 = 0.0;
            if (f->GetShapeCount() > 0 && (f->GetShape(0)->GetType() == LISEM_SHAPE_LINE || f->GetShape(0)->GetType() == LISEM_SHAPE_MULTILINE)) {
                double len_m = CalculateFeatureLengthM(f, roads_sf->GetProjection(), center_lat);
                area_m2 = len_m * 25.0; // 25 meters standard polyline road buffer width
            } else {
                area_m2 = CalculateFeatureAreaM2(f, roads_sf->GetProjection(), center_lat);
            }
            if (area_m2 <= 0.0) area_m2 = 500.0;

            double asset_val = area_m2 * road_val_per_m2;
            summary.roads_rcv += asset_val;

            auto [stat, note] = EvaluateAccessibility(depth_m, hazard_type);
            double dmg_ratio = vcurve.EvaluateDamageRatio(depth_m);
            double loss = 0.0;

            if (dmg_ratio > 0.0) {
                summary.impacted_roads++;
                loss = dmg_ratio * asset_val;
                summary.roads_loss += loss;
            }

            layer->SetFeatureAttribute(i, "depth_m", depth_m);
            layer->SetFeatureAttribute(i, "dmg_ratio", dmg_ratio);
            layer->SetFeatureAttribute(i, "dmg_pct", dmg_ratio * 100.0);
            layer->SetFeatureAttribute(i, "loss", loss);
            layer->SetFeatureAttribute(i, "asset_val", asset_val);
            layer->SetFeatureAttribute(i, "hazard", hazard_type);
            layer->SetFeatureAttribute(i, "access_stat", stat);
            layer->SetFeatureAttribute(i, "access_note", note);
            layer->SetFeatureAttribute(i, "mon_unit", mon_unit);
            layer->SetFeatureAttribute(i, "vuln_curve", vcurve.name);
        }
    }
}

inline void ComputeAgriLoss(
    ShapeFile* agri_sf,
    const VulnerabilityCurve& vcurve,
    double agri_val_per_m2,
    const std::string& mon_unit,
    PortfolioRiskSummary& summary,
    const std::string& hazard_type = "flood"
) {
    if (!agri_sf) return;

    for (int l = 0; l < agri_sf->GetLayerCount(); ++l) {
        ShapeLayer* layer = agri_sf->GetLayer(l);

        layer->AddAttribute("depth_m", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("dmg_ratio", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("dmg_pct", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("loss", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("crop_val", ATTRIBUTE_TYPE_FLOAT64);
        layer->AddAttribute("hazard", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("mon_unit", ATTRIBUTE_TYPE_STRING);
        layer->AddAttribute("vuln_curve", ATTRIBUTE_TYPE_STRING);

        int i1_idx = layer->GetAttributeIndex("i1");

        for (int i = 0; i < layer->GetFeatureCount(); ++i) {
            Feature* f = layer->GetFeature(i);
            summary.total_agriplots++;

            double depth_m = 0.0;
            if (i1_idx >= 0) {
                std::string s_val = layer->GetAttributeValue(i1_idx, i);
                if (!s_val.empty()) {
                    try { depth_m = std::stod(s_val); } catch (...) {}
                }
            }
            if (depth_m < 0.0 || std::isnan(depth_m)) depth_m = 0.0;

            BoundingBox bbox = f->GetBoundingBox();
            double center_lat = (bbox.GetMinY() + bbox.GetMaxY()) / 2.0;
            if (std::abs(center_lat) > 90.0) center_lat = 42.0;
            double area_m2 = CalculateFeatureAreaM2(f, agri_sf->GetProjection(), center_lat);
            if (area_m2 <= 0.0) area_m2 = 1000.0;

            double crop_val = area_m2 * agri_val_per_m2;
            summary.agriplots_rcv += crop_val;

            double dmg_ratio = vcurve.EvaluateDamageRatio(depth_m);
            double loss = 0.0;

            if (dmg_ratio > 0.0) {
                summary.impacted_agriplots++;
                loss = dmg_ratio * crop_val;
                summary.agriplots_loss += loss;
            }

            layer->SetFeatureAttribute(i, "depth_m", depth_m);
            layer->SetFeatureAttribute(i, "dmg_ratio", dmg_ratio);
            layer->SetFeatureAttribute(i, "dmg_pct", dmg_ratio * 100.0);
            layer->SetFeatureAttribute(i, "loss", loss);
            layer->SetFeatureAttribute(i, "crop_val", crop_val);
            layer->SetFeatureAttribute(i, "hazard", hazard_type);
            layer->SetFeatureAttribute(i, "mon_unit", mon_unit);
            layer->SetFeatureAttribute(i, "vuln_curve", vcurve.name);
        }
    }
}

} // namespace FastRisk
