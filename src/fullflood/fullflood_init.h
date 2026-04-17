#ifndef FULLFLOOD_INIT_H
#define FULLFLOOD_INIT_H

#include "../algorithms/rasteralgorithms/rasterize.h"
#include "../algorithms/rasteralgorithms/demchannel.h"
#include "../algorithms/rasteralgorithms/dinfinite.h"
#include "../algorithms/rasteralgorithms/monotonicreconstruct.h"
#include "../algorithms/rasteralgorithms/spreadflow.h"
#include "../algorithms/rasteralgorithms/slope.h"
#include "../algorithms/rasteralgorithms/spread.h"
#include "../algorithms/rasteralgorithms/derivatives.h"
#include "../algorithms/rasteralgorithms/demproc.h"
#include "../common/defines.h"

extern "C"
{

inline void EMSCRIPTEN_KEEPALIVE do_algo_fullflood_init(
    float*ret,float*retvx,float*retvy,float*retq,
    float*data_dem,float*data_man, float*data_rain,
    float * data_chanw, float *data_chanh,
    int w, int h,
    float in_dx,
    float duration,
    int do_depr, int correct_hq, int pressure_iters, float quality_p, int in_fluxiters, float delta,
    float diffusivity_coefficient, float man_chan,
    float * discharge_out, float * discharge_out2,
    int do_ocean, float oclevel,
    int do_q_neg,int shapes, float * points, float * slengths, float * stypes, float * respasscap, float * cap,
    int do_qdetail_out, float * ret_b, float * ret_l, float * ret_avg_vel, float * ret_ss_comp,
    float conc_vel_mult,
    int do_qbase, float qbase_max, float * ret_qbase,
    int do_chancor, float min_drainchancor, float * dem_channelized,
    int n_culverts, float * culvert_q, float * culvert_fraction,float * culverts_source_row, float * culverts_source_col,
    float * culverts_target_row, float * culverts_target_col,
    int do_new_crosssections, int do_new_manning, float * stormdraindata,
    int do_output_vel, float * vel_x, float * vel_y, float chan_a, float chan_b, float chan_c, float chan_d,

    float * ret_whmax, float * ret_vmax, float * ret_at,
        float wh_initial = 0.0, int do_cache_act = 1, int do_sdt_act = 1,
        int ff_input_cores = 8,
        float * gradx1 = 0, float * gradx2 = 0, float * grady1=0, float * grady2=0
)
{

#ifdef FASTFLOOD_WASM
    //omp_set_dynamic(0);
    GetTPool().set_active_threads(ff_input_cores);
    FF::cout << "set threads " << ff_input_cores <<std::endl;
    //omp_set_num_threads(std::stoi(ff_input_cores));
#endif
    FF::cout << " init fullflood, do sdt: " << do_sdt_act << " do_caching: " << do_cache_act << " do channel: " << (data_chanw != 0 && data_chanh != 0) << " dims: " << w << " " <<h<< std::endl;

    FF_FULLFLOOD.SetIterations(in_fluxiters);

    FF_FULLFLOOD.SetCoreOptions(do_sdt_act >= 1,
                                do_cache_act >= 1,
                                data_chanw != 0 && data_chanh != 0);
   FF_FULLFLOOD.SetCoreMapsAndAllocate(
              new cTMapf(data_dem,w,h,0.0,0.0, in_dx,- in_dx)
            , new cTMapf(data_man,w,h,0.0,0.0, in_dx,- in_dx)
            , new cTMapf(0,w,h,0.0,0.0, in_dx,- in_dx)
            , new cTMapf(data_rain,w,h,0.0,0.0, in_dx,- in_dx)
               ,duration
               , data_chanw? new cTMapf(data_chanw,w,h,0.0,0.0, in_dx,- in_dx): nullptr
                , data_chanh?new cTMapf(data_chanh,w,h,0.0,0.0, in_dx,- in_dx) : nullptr,
               dem_channelized? new cTMapf(dem_channelized,w,h,0.0,0.0, in_dx,- in_dx) : nullptr,
               delta,in_dx,
               new cTMapf(gradx1,w,h,0.0,0.0, in_dx,- in_dx),
               new cTMapf(gradx2,w,h,0.0,0.0, in_dx,- in_dx),
               new cTMapf(grady1,w,h,0.0,0.0, in_dx,- in_dx),
               new cTMapf(grady2,w,h,0.0,0.0, in_dx,- in_dx)
               );
    FF_FULLFLOOD.SetStormDrains(stormdraindata != nullptr? new cTMapf(stormdraindata,w,h,0.0,0.0, in_dx,- in_dx) : nullptr);
    FF_FULLFLOOD.SetWHInitial(wh_initial);
    if(n_culverts > 0)
    {
        FF_FULLFLOOD.SetCulverts(n_culverts,culvert_q, culvert_fraction,
                                 culverts_source_row, culverts_source_col,
                                     culverts_target_row, culverts_target_col);
    }
    if(do_q_neg)
    {
        FF_FULLFLOOD.SetReservoirs(shapes, points,slengths,stypes,respasscap,cap);
    }
    if(do_ocean > 0.5)
    {

        FF_FULLFLOOD.SetOceanBoundary(oclevel);
    }
    FF_FULLFLOOD.SetCustomReturnPointers(ret,retvx,retvy,retq,ret_whmax,ret_vmax,ret_at);

}
}

#endif
