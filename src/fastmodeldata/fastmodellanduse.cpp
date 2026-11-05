#include "fastmodeldata.h"
#include "../inputprepare/forecastprocess.h"
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

std::shared_ptr<FastModelLandUseData> FetchFastModelLandUseData(std::shared_ptr<cTMapf> ref_map, const std::string& token, const std::string& customurl) {
    std::vector<std::shared_ptr<cTMapf>> res = GetLUManningMaps(ref_map, token, customurl);
    if (res.size() > 0) {
        auto data = std::make_shared<FastModelLandUseData>();
        data->lu = res[0];
        data->manning = res[1];
        data->isLUDef = true;
        return data;
    }
    return nullptr;
}
