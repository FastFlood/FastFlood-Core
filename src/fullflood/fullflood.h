
#ifndef FULLFLOOD_3CASCHELPER_H
#define FULLFLOOD_H
/*#include <stdio.h>
#include "map.h"
#include "parallel.h"
#include "emscripten.h"

#include "slope.h"

#include "defines.h"
#include "hillshade.h"
#include "filter.h"
#include "ldd.h"
#include "derivatives.h"
#include "mapstats.h"

#include "spread.h"

#include "maxwave.h"
#include "spreadflow.h"
#include "rasterize.h"
#include "soilstats.h"
#include "demproc.h"
#include "demchannel.h"
#include "oceanfill.h"
#include "infiltration.h"
#include "upscale.h"
#include "flowwidthq.h"*/

//#include "../threading/tp_simple.h"
//ThreadPool  p;


#include "../util/lisemmath.h"
#include "fullflood_state.h"
#include "fullflood_kinematic.h"
#include "fullflood_riemann.h"
#include "fullflood_3caschelper.h"
#include "fullflood_mmsolve.h"
#include "fullflood_courant.h"
#include "fullflood_gridadapt.h"
#include "fullflood_cache.h"
#include "fullflood_gridresample.h"

#include "fullflood_svsolve.h"

#include "fullflood_step.h"
#include "fullflood_init.h"
#include "fullflood_clean.h"



inline FF_FULLFLOODT * InitFullFloodModel(cTMapf * dem, cTMapf * N, cTMapf * ksat, cTMapf * rain, double duration, bool do_sdt, bool do_cache)
{
    int in_fluxiters = 100.0;
    double delta = 0.01;

    float data_chanw = 0;
    float data_chanh = 0;

    float in_dx = dem->CellSizeXApproximateMeters();

    cTMapf *stormdraindata = nullptr;
    FF_FULLFLOODT *FULLFLOOD = new FF_FULLFLOODT();
    FULLFLOOD->SetIterations(in_fluxiters);

    FULLFLOOD->SetCoreOptions(do_sdt >= 1,
                                do_cache >= 1,
                                data_chanw != 0 && data_chanh != 0);
    FULLFLOOD->SetCoreMapsAndAllocate(
              dem
            , N
            ,ksat
            , rain
               ,duration
               ,nullptr
               , nullptr
               , nullptr,
               delta,
                in_dx,0,0,0,0);
    FULLFLOOD->SetStormDrains(stormdraindata);
    FULLFLOOD->SetWHInitial(0.0);
    /*if(n_culverts > 0)
    {
        FULLFLOOD.SetCulverts(n_culverts,culvert_q, culvert_fraction,
                                 culverts_source_row, culverts_source_col,
                                     culverts_target_row, culverts_target_col);
    }
    if(do_q_neg)
    {
        FULLFLOOD.SetReservoirs(shapes, points,slengths,stypes,respasscap,cap);
    }
    if(do_ocean > 0.5)
    {

        FULLFLOOD.SetOceanBoundary(oclevel);
    }
    FULLFLOOD.SetCustomReturnPointers(ret,retvx,retvy,retq,ret_whmax,ret_vmax,ret_at);*/



    return FULLFLOOD;
}

inline void StepFullFlood(FF_FULLFLOODT * FULLFLOOD, double dt, double courant, double courant_forced, double scale_req, bool use_openmp, bool do_flow)
{

    do_algo_fullflood_step_onstate(*FULLFLOOD,dt,courant,-1.0,courant_forced,-1,scale_req,use_openmp,do_flow);

}

/*void main()
{
     FullFloodState @state = CreateFullFloodModel(dem_lidar.tif,
         dem_lidar.tif * 0.0 + 0.1,
                  dem_lidar.tif * 0.0 + 1.0,
                           dem_lidar.tif * 0.0 + 10.0,3);

                       StepFullFlood(state,1.0);
                   }*/

#endif

