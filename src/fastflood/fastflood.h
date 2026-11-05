#ifndef FASTFLOOD_H
#define FASTFLOOD_H

#include <stdio.h>
#include "../geospatial/map/map.h"
#include "../common/threading/parallel.h"

#include "../util/lsm_vector3.h"
#include "../util/lsm_vector2.h"
#include "../algorithms/rasteralgorithms/slope.h"
#include "../io/csv.h"
#include "../common/defines.h"
#include "../algorithms/rasteralgorithms/hillshade.h"
#include "../algorithms/filter/filter.h"
#include "../algorithms/rasteralgorithms/ldd.h"
#include "../algorithms/rasteralgorithms/derivatives.h"
#include "../util/mapstats.h"
#include "../algorithms/rasteralgorithms/mapstats.h"

#include "../algorithms/rasteralgorithms/spread.h"

#include "maxwave.h"
#include "../algorithms/rasteralgorithms/spreadflow.h"
#include "../util/rasterize.h"
#include "../algorithms/soil/soilstats.h"
#include "../algorithms/rasteralgorithms/demproc.h"
#include "../algorithms/rasteralgorithms/demchannel.h"
#include "oceanfill.h"
#include "../algorithms/soil/infiltration.h"
#include "upscale.h"
#include "flowwidthq.h"

#include "velcor.h"

#include "../algorithms/rasteralgorithms/distancefield.h"



extern "C"
{
	



/*void EMSCRIPTEN_KEEPALIVE do_algo_fastflood(
 *
    float*ret,                      // OUT REQ-> pointer to float array to contain wh results in meter
    float*data_dem,                 // IN REQ-> pointer to elevation input map (meter)
    float*data_man,                 // IN REQ-> pointer to manning surface roughness coefficient map (metric)
    float*data_rain,                // IN REQ-> pointer to rainfall input map (mm/h)

    float * data_chanw,             // IN OPT-> pointer to channel width map (meter)
    float *data_chanh,              // IN OPT-> pointer to channel depth map (meter)

    int w, int h,                   // IN REQ-> width and height of the 2d input rasters

    float in_dx,                    // IN REQ-> cellsize in meter

    float duration,                 // IN REQ-> event duration in hours

    int do_depr,                    // IN REQ-> do not do hydro corr if >= 1, default = 0
    int correct_hq,                 // IN REQ-> do correct discahrge if >= 1, default = 0

    int pressure_iters,             // IN REQ-> do not do hydro corr if >= 1, default = 40
    float quality_p,                // IN REQ->
    int in_fluxiters,               // IN REQ-> number of flux sweep iters, default = 40

    float delta,                    // IN REQ-> delta for terrain correction for flow network, default = 0.001
    float diffusivity_coefficient,  // IN REQ-> diffusivity coefficient, default = 0.2
    float man_chan,                 // IN REQ-> mannings value for channels, default = 0.02

    float * discharge_out,          // OUT OPT-> Output discharge float array (m3/s)
    float * discharge_out2,         // OUT OPT-> Output discharge float array (m3/s)

    int do_ocean,                   // IN REQ-> do ocean boundary condition if >= 1, default = 0
    float level,                    // IN REQ-> ocean boundary condition height in meters, default = 0

    int do_q_neg,
    int shapes,
    float * points,
    float * slengths,
    float * stypes,
    float * respasscap,
    float * cap,

    int do_qdetail_out,
    float * ret_b,
    float * ret_l,
    float * ret_avg_vel,
    float * ret_ss_comp,

    float conc_vel_mult,

    int do_qbase,
    float qbase_max,
    float * ret_qbase,

    int do_chancor,
    float min_drainchancor,

    float * dem_channelized,

    int n_culverts,
    float * culvert_q,
    float * culvert_fraction,
    float * culverts_source_row,
    float * culverts_source_col,
    float * culverts_target_row,
    float * culverts_target_col,

    int do_new_crosssections,
    int do_new_manning,

    float * stormdraindata,

    int do_output_vel,
    float * vel_x,
    float * vel_y )
{

}*/


inline void ff_do_algo_fastflood(
    float*ret,
    float*data_dem,float*data_man, float*data_rain, float * data_chanw, float *data_chanh,
    int w, int h, 
    float in_dx, 
    float duration, 
    int do_depr, int correct_hq, int pressure_iters, float quality_p, int in_fluxiters, float delta, 
    float diffusivity_coefficient, float man_chan, 
    float * discharge_out, float * discharge_out2, 
    int do_ocean, float level, 
    int do_q_neg,
    int shapes,
    float * points,
    float * slengths,
    float * stypes,
    float * respasscap,
    float * cap,
    int do_qdetail_out, float * ret_b, float * ret_l, float * ret_avg_vel, float * ret_ss_comp,
    float conc_vel_mult,
    int do_qbase, float qbase_max, float * ret_qbase,
    int do_chancor, float min_drainchancor, float * dem_channelized,
    int n_culverts,
    float * culvert_q,
    float * culvert_fraction,
    float * culverts_source_row,
    float * culverts_source_col,
    float * culverts_target_row,
    float * culverts_target_col,
    int do_new_crosssections, int do_new_manning, float * stormdraindata,
    int do_output_vel, float * vel_x, float * vel_y, bool export_chanwh,
    int n_sources,
    float * val_sources,
    float * row_sources,
    float * col_sources,
    float timeslice,
    float * dem_mon_org,
    float * da_source_pix,
    int noqfregion_rmin,
    int noqfregion_rmax,
    int noqfregion_cmin,
    int noqfregion_cmax,
    float * demnoshapes,
    int do_wave,
    float wave_h,
    float wave_freq,
    float * culvert_qtot = nullptr,
    float * culvert_qmaxmodel = nullptr,
    float * reservoir_qtot = nullptr,
    float * reservoir_final_vols = nullptr,
    double north = 0.0,
    double west = 0.0,
    double cellsize_x = 0.0,
    double cellsize_y = 0.0
        )
{
	
    FF::cout << "start fastflood function " << data_dem << " cols: " << w << " rows: " << h << " dx: " << in_dx << " duration: " << duration << std::endl;
    if(export_chanwh)
    {
        FF::cout << "output wh with channel" <<std::endl;
    }else
    {
        FF::cout << "output wh without channel" <<std::endl;
    }

    double c_north = north;
    double c_west = west;
    double c_cellsize_x = cellsize_x == 0.0 ? (double)in_dx : cellsize_x;
    double c_cellsize_y = cellsize_y == 0.0 ? (double)-in_dx : cellsize_y;

    if ((std::fabs(c_cellsize_x) < 0.5 && std::fabs(c_north) <= 90.0) || in_dx < 0.5f)
    {
        double avg_lat = c_north + 0.5 * (double)h * (cellsize_y != 0.0 ? cellsize_y : -std::fabs(c_cellsize_x));
        double cs = (cellsize_x != 0.0) ? std::fabs(cellsize_x) : (double)in_dx;
        double est_dx = cs * 111321.0 * std::cos(avg_lat * 3.141592653589793 / 180.0);
        if (in_dx < 0.5f)
        {
            in_dx = (float)est_dx;
        }
        c_cellsize_x = (double)in_dx;
        c_cellsize_y = (double)-in_dx;
    }
    else if (in_dx >= 0.5f)
    {
        c_cellsize_x = (double)in_dx;
        c_cellsize_y = (double)-in_dx;
    }

    cTMapf * RAININ = new cTMapf(da_source_pix,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
    cTMapf * DEM_ORG_MON = new cTMapf(dem_mon_org,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
    cTMapf * DEM_ORG = new cTMapf(data_dem,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);


    //writeRaster(*RAININ,"rainintest.tif");
    cTMapf * DEM_NOSHAPES = nullptr;
    if(demnoshapes)
    {
        FF::cout << "dem noshapes present" << std::endl;
        DEM_NOSHAPES = new cTMapf(demnoshapes,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
    }

    cTMapf * N = new cTMapf(data_man,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
    cTMapf * Rain_ORG = new cTMapf(data_rain,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
   
	bool corr = correct_hq > 0;
    bool do_depressions = do_depr > 0;
	bool do_output_discharge = discharge_out != 0;
	
    cTMapf * DEM_ORGCHAN = nullptr;
    if(dem_channelized != 0)
    {
        DEM_ORGCHAN = new cTMapf(dem_channelized,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
    }
    
    cTMapf * SDC = nullptr;
    cTMapf * SDCTOT = nullptr;
    bool do_stormdrain = false;
    if(stormdraindata != 0)
    {
        do_stormdrain = true;
        FF::cout <<" do stormdrain " <<std::endl;
        SDC = new cTMapf(stormdraindata,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
        SDCTOT = DEM_ORG->GetCopy0();
    }
	int piters = pressure_iters;
	int itersp = (int)(quality_p *15);
	
	bool do_chan = data_chanw != 0;

	//printf("do flood algorithm : %i %i %i %i %f do chan? %i culverts: %i",  do_depressions, corr, piters, iterations,duration, do_chan?1:0,n_culverts);
	
    int n_pix = 0;
    cTMapf * QBASE = nullptr;
    cTMapf * QBASES = nullptr;
    bool qbase_set = false;

    cTMapf * CULVERTID = nullptr;
    

    std::vector<float> CULVERT_QTOT;
    std::vector<float> CULVERT_Q;
    std::vector<float> CULVERT_FRAC;
    std::vector<float> culvert_ffactors;
    std::vector<int> CULVERT_RS;
    std::vector<int> CULVERT_CS;
    std::vector<int> CULVERT_RT;
    std::vector<int> CULVERT_CT;

    std::cout << "Number of culverts: " << n_culverts << std::endl;
    if(n_culverts > 0 || do_stormdrain)
    {

        printf("do culverts ");


        for(int i = 0; i <n_culverts; i++)
        {
            FF::cout << "culvert " << i << " " << culvert_q[i] << std::endl;
            int r_source = (int)(culverts_source_row[i]);
            int c_source = (int)(culverts_source_col[i]);
            int r_target = (int)(culverts_target_row[i]);
            int c_target = (int)(culverts_target_col[i]);

            CULVERT_RS.push_back(r_source);
            CULVERT_CS.push_back(c_source);
            CULVERT_RT.push_back(r_target);
            CULVERT_CT.push_back(c_target);
            
            CULVERT_Q.push_back(culvert_q[i]);
            //FF::cout << " culvert: " << i << " " <<  culvert_q[i] << " " <<(culvert_fraction[i]/100.0f) << " locs: " << r_source << " " << c_source << "   "  << r_target << " "<<c_target <<  std::endl;
            CULVERT_FRAC.push_back(std::max(0.0f,std::min(1.0f,(culvert_fraction[i]/100.0f))));
            CULVERT_QTOT.push_back(0.0);
        }
        CULVERTID = DEM_ORG->GetCopy0();
         for (int r = 0; r < DEM_ORG->Rows(); r++) {
            for (int c = 0; c < DEM_ORG->Cols(); c++) {

                CULVERTID->data[r][c] = -1.0;

                for(int i = 0; i <n_culverts; i++)
                {
                    int r_source = (int)(culverts_source_row[i]);
                    int c_source = (int)(culverts_source_col[i]);
                    int r_target = (int)(culverts_target_row[i]);
                    int c_target = (int)(culverts_target_col[i]);

                    if((r == r_source && c == c_source) || (r == r_target && c == c_target))
                    {
                        //FF::cout << "added culvert to culvert ID map" << std::endl;
                        CULVERTID->data[r][c] = i;
                    }
                }
            }
         }



    }
    
	
	//keep track of previous total (compensated) discharge
    //and of previous water height
    FF::cout << "Has specific dem channel? " << (DEM_ORGCHAN != nullptr) << std::endl;
    cTMapf * DEM = DEM_ORGCHAN? DEM_ORGCHAN->GetCopy():DEM_ORG->GetCopy();
    cTMapf * WH_TOT = DEM->GetCopy0();
    cTMapf * SLOPE = DEM->GetCopy0();
    cTMapf * WH_PREV = DEM->GetCopy0();
    cTMapf * QSS_PREV = DEM->GetCopy0();
    cTMapf * Rain = Rain_ORG->GetCopy();


    cTMapf *WHO = nullptr;
    cTMapf *WHODEM = nullptr;
    float rain_temp_total = true;
    if(do_ocean > 0 || do_wave > 0 )
    {
        printf("ocean level %g , wave level %g", level, wave_h);
        if (do_wave) {

            WHO = do_fill_ocean_wave_decay(DEM_ORG,N,wave_h,wave_freq,level,in_fluxiters, 1.35);
        }else {
            WHO = do_fill_ocean(DEM_ORG,level,in_fluxiters);
        }

        //check if we should do rainfall at all?
        rain_temp_total = false;
        for (int r = 0; r < DEM->Rows(); r++) {
            for (int c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(Rain_ORG->data[r][c])) {

                    if(Rain_ORG->data[r][c]>0.0)
                    {
                        rain_temp_total = true;
                        break;
                    }

                }
            }
            if(rain_temp_total)
            {
                break;
            }
        }
        if(rain_temp_total)
        {
            WHODEM = DEM->GetCopy();
        }


    }
    

    cTMapf * ReservoirID = nullptr;
    std::vector<int> reservoirids;
    std::vector<float> reservoirvols;
    std::vector<float> reservoircvols;
    std::vector<float> reservoir_inflow;
    std::vector<int> reservoirlink;
    std::vector<float> reservoirpassvols;
    if(do_q_neg || n_culverts > 0 || do_stormdrain)
    {
        FF::cout << "do q neg " <<std::endl;
        ReservoirID = DEM->GetCopy();

        for (int r = 0; r < DEM->Rows(); r++) {
            for (int c = 0; c < DEM->Cols(); c++) {
                pcr::setMV(ReservoirID->data[r][c]);
            }
        }

        //we rasterize the negative discharge shapes and store an ID to refer to the storage and throughflow total
        if(shapes > 0)
        {
            rasterize_reservoirs(ReservoirID, shapes,points, slengths,stypes);
            FF::cout << "shapes n " << shapes << std::endl;
            //create a list with total capacity
            for(int i = 0; i < shapes; i++)
            {
                FF::cout << "shape " << i << " " << cap[i] << " " << respasscap[i] << std::endl;
                reservoirids.push_back(i);
                reservoirvols.push_back(cap[i]);
                reservoircvols.push_back(0.0);
                reservoir_inflow.push_back(0.0);
                float qh = respasscap[i];
                reservoirpassvols.push_back(qh);
            }
        }
    }

    cTMapf * CHANW = nullptr;
    cTMapf * CHANH = nullptr;
    cTMapf * ChanStorage = nullptr;
    if(do_chan)
	{
        CHANW = new cTMapf(data_chanw,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
        CHANH = new cTMapf(data_chanh,w,h,c_north,c_west, c_cellsize_x,c_cellsize_y);
		ChanStorage = DEM->GetCopy0();
	}
	if(piters <= 1)
    {
        piters = 1;
    }
    if(do_ocean && !rain_temp_total)
    {
        piters = 0;
    }

    
    cTMapf * q_detail_w =nullptr; //DEM->GetCopy0();
    if(do_qdetail_out > 0 && piters > 1)
    {
        q_detail_w =DEM->GetCopy0();

        for(int r = 0; r < DEM->Rows(); r++) {
				for(int c = 0; c < DEM->Cols(); c++)
				{
					int index = r* DEM->Cols() + c;
					if(!pcr::isMV(DEM->data[r][c]))
					{
                        q_detail_w->data[r][c] = 0.0;
                        
                    }
                }
        }


    }

    if(do_qdetail_out > 0)
    {
         for(int r = 0; r < DEM->Rows(); r++) {
				for(int c = 0; c < DEM->Cols(); c++)
				{
					int index = r* DEM->Cols() + c;
					if(!pcr::isMV(DEM->data[r][c]))
					{   
                        ret_b[index] = 0.0;
                        ret_l[index] = 0.0;
                        ret_avg_vel[index] = 0.0;
                        ret_ss_comp[index] = 0.0;
                    }
                }
         }
    }


#ifdef FASTFLOOD_WASM
	emscripten_run_script("increment_progress_fastflood()");
#endif

    float * local_ss_comp = nullptr;
    float * local_ss_weight = nullptr;
    bool use_local_ss = (do_qdetail_out <= 0) && (n_culverts > 0) && (culvert_qtot != nullptr);
    if (use_local_ss)
    {
        int total_pixels = DEM->Rows() * DEM->Cols();
        local_ss_comp = new float[total_pixels]();
        if (piters > 1)
        {
            local_ss_weight = new float[total_pixels]();
        }
    }
	
	for(int ip = 0; ip < piters; ip++)
    {
        if(itersp < 0)
        {
            itersp = 40;
        }
        //first we get the right gradient maps

        FF::cout << "piters: " << ip <<std::endl;



        cTMapf * DEMD = nullptr;

        //then do the elevation model fix
        cTMapf * DEM2 = nullptr;

        cTMapf * FlowWL = nullptr;
        cTMapf * FlowWR = nullptr;
        if(do_output_vel)
        {
            FlowWL = DEM->GetCopy0();
            FlowWR = DEM->GetCopy0();
        }


        cTMapf* dmo =DEM_ORGCHAN? DEM_ORGCHAN:DEM_ORG;
		if(do_chan)
		{
			//#pragma omp parallel for collapse(2)
			for (int r = 0; r < DEM->Rows(); r++) {
				for (int c = 0; c < DEM->Cols(); c++) {
					if (!pcr::isMV(DEM->data[r][c])) {

                        DEM->data[r][c] = dmo->data[r][c] + WH_TOT->data[r][c] - (CHANH->data[r][c]*CHANW->data[r][c]/in_dx) + (ChanStorage->data[r][c] /(in_dx*in_dx));

					}
				}
			}
			
		}else{
			//#pragma omp parallel for collapse(2)
			for (int r = 0; r < DEM->Rows(); r++) {
				for (int c = 0; c < DEM->Cols(); c++) {
					if (!pcr::isMV(DEM->data[r][c])) {

                        DEM->data[r][c] = dmo->data[r][c] + WH_TOT->data[r][c];

					}
				}
			}
		}
        

        float dx = in_dx;
        FF::cout << "in_fluxiters " << in_fluxiters <<std::endl;
        int iters= in_fluxiters * std::max(1.0f,(float)std::max(DEM->nrCols(),DEM->nrRows())/1000.0f);
        cTMapf * GradX1 = SlopeX1(DEM); GradX1->setcellSizeX(in_dx); GradX1->setcellSizeY(-in_dx);
        cTMapf * GradY1 = SlopeY1(DEM); GradY1->setcellSizeX(in_dx); GradY1->setcellSizeY(-in_dx);
        cTMapf * GradX2 = SlopeX2(DEM); GradX2->setcellSizeX(in_dx); GradX2->setcellSizeY(-in_dx);
        cTMapf * GradY2 = SlopeY2(DEM); GradY2->setcellSizeX(in_dx); GradY2->setcellSizeY(-in_dx);
        cTMapf * one = DEM->GetCopy1();
	    FF::cout << "real fluxiters " << iters << " " << dx << std::endl;
        //get a map with seed points for corrected dem
        //MaskedRaster<float> raster_data3(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
        cTMapf *Seeds = new cTMapf(nullptr, DEM->data.nr_cols(), DEM->data.nr_rows(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data3),DEM->projection(),"");

        int r = 0,c = 0;
        int seeds = 0;
        float mindem = MapMinimumRed(DEM);

        //#pragma omp parallel for collapse(2)
        for (r = 0; r < DEM->Rows(); r++) {
            for (c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(DEM->data[r][c])) {

                    if(ip == 0)
                    {
                        Rain->data[r][c] = (1.0/((float)(piters)))*(dx*dx*Rain_ORG->data[r][c])/3600000.0;
                    }
                    //if dem < smallest dem, or if it is an edge cell, make it a seed for the new dem
                    if(OUTORMV(DEM,r,c) || OUTORMV(DEM,r,c+1) || OUTORMV(DEM,r,c-1) || OUTORMV(DEM,r+1,c) || OUTORMV(DEM,r-1,c) || (DEM->data[r][c] < mindem + 1e-6))
                    {
                        seeds ++;
                        Seeds->data[r][c] =1.0;
                    }else
                    {
                        Seeds->data[r][c] =0.0;
                    }
                }else
                {
                    Seeds->data[r][c] =0.0;
                }
            }
        }
#ifdef FASTFLOOD_WASM

        FF::cout << "Start initial processing threads "<< std::endl;
		emscripten_run_script("increment_progress_fastflood()");
#endif
        cTMapf * DEM2CHAN = nullptr;

        cTMapf * DEMDM = nullptr;
        cTMapf * DEMDID = nullptr;
        //get depressions

        std::vector<int> depressionids;
        std::vector<float> depressionvols;
        std::vector<float> depressionfeatureinteract;
        std::vector<float> depressioncvols;

        std::vector<std::thread> workers1;
		//alter N based on channel presence
		
		if(do_chan)
		{
			
			workers1.push_back(std::thread([N,DEM,CHANW,CHANH,man_chan,in_dx,iters,delta]()
			{
				
				for (int r = 0; r < DEM->Rows(); r++) {
						 for (int c = 0; c < DEM->Cols(); c++) {
							 if (!pcr::isMV(DEM->data[r][c])) {
								 
								 float width = CHANW->data[r][c];
								 width = std::min(in_dx,width);
                                 N->data[r][c] = std::max(0.00001f,(width * man_chan + std::max(0.00001f,(in_dx - width))*N->data[r][c])/in_dx);
							 }
						 }
				}
				

			}));
		}
		
		
        workers1.push_back(std::thread([&DEM2,&DEM2CHAN, ip,DEM_ORGCHAN,DEM,Seeds,Rain,N,GradX1,GradX2,GradY1,GradY2,iters,delta]()
        {
			
            //allways need this one
            {
                FF::cout << "start fastflood DEM correction " << delta << " " << iters << std::endl;

                DEM2 = SpreadMonotonicReconstruct(DEM,/*0.00001*/delta,iters);//SpreadDirectionalAbsMaxMD(Seeds,DEM,GradX1,GradX2,GradY1,GradY2);

            }
            

        }));
        workers1.push_back(std::thread([&DEM2,&DEM2CHAN, ip,DEM_ORGCHAN,DEM,Seeds,Rain,N,GradX1,GradX2,GradY1,GradY2,iters,delta]()
        {
			if(DEM_ORGCHAN != nullptr && ip == 0)
            {
                FF::cout << "using channelized dem "<< std::endl;
                DEM2CHAN = SpreadMonotonicReconstruct(DEM_ORGCHAN,/*0.00001*/delta,iters);//SpreadDirectionalAbsMaxMD(Seeds,DEM,GradX1,GradX2,GradY1,GradY2);

            }

        }));
        if(do_depressions)
        {
			printf("do clump etc.. \n");
            workers1.push_back(std::thread([dx,DEM_ORG,&depressioncvols,DEM_NOSHAPES,&depressionfeatureinteract,&depressionvols,&depressionids,&DEMDID,DEM,&DEMD,&DEMDM,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
            {
				
				printf("a thread2 started \n");
                //DEMD will contain the depression depth
                //MaskedRaster<float> raster_data5();
                DEMDM = new cTMapf(nullptr,  DEM->data.nr_cols(),DEM->data.nr_rows(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data5),DEM->projection(),"");

                DEMD = SpreadDepressionFind(DEM_ORG);
				
                int depressioncells = 0;
                 for (int r = 0; r < DEM->Rows(); r++) {
                     for (int c = 0; c < DEM->Cols(); c++) {
                         if (!pcr::isMV(DEM->data[r][c])) {
                             float valn = DEMD->data[r][c]-DEM_ORG->data[r][c];

                             if(valn > 0.05)
                             {
                                 DEMDM->data[r][c] = DEMD->data[r][c];
                                 depressioncells++;
                                 DEMD->data[r][c] = valn;
                             }else
                             {
                                 pcr::setMV(DEMDM->data[r][c]);
                                 DEMD->data[r][c] = 0.0;
                             }
                         }else
                         {
                             pcr::setMV(DEMDM->data[r][c]);
                         }
                     }
                 }
				printf("a thread2 clump \n");
                 DEMDID = Clump(DEMDM,0.1,true);

				 printf("a thread2 clumpd \n");

                 //get depression volumes
                 float cellarea = dx * dx;

                 depressioncells = 0;
                 for(int r = 0; r < DEMDID->data.nr_rows();r++)
                 {
                     for(int c = 0; c < DEMDID->data.nr_cols();c++)
                     {
                         if(!pcr::isMV(DEMDID->data[r][c]))
                         {
                             depressioncells ++;
                             bool found = false;
                             int class_current = (int)DEMDID->data[r][c];
                             //FF::cout << "add class " << class_current << " " << depressionids.size() <<  std::endl;

                             if(class_current > (((int)depressionids.size())-5))
                             {
                                 //FF::cout << "add class " << class_current << " " << depressionids.size() <<  std::endl;
                                 int to_add = (5 + class_current - depressionids.size());
                                 //FF::cout << "to_add " << to_add << std::endl;
                                 int oldsize = depressionids.size();
                                 for(int idr = 0; idr < to_add; idr++)
                                 {
                                     //FF::cout << "add" << std::endl;
                                     depressionids.push_back(oldsize + idr);
                                     depressionvols.push_back(0.0);
                                     depressionfeatureinteract.push_back(0.0);
                                 }

                             }
                             {
                                //FF::cout << "addr " << depressionvols.size() << " " << class_current << std::endl;
                                 depressionvols.at(class_current) += cellarea * DEMD->data[r][c];

                                 if(DEM_NOSHAPES)
                                 {
                                    int r3 = r; int c3 = c;
                                    if(r3 >= 0 && c3 >=0 && r3 < DEM_ORG->nrRows() && c3 < DEM_ORG->nrCols())
                                    {
                                        depressionfeatureinteract.at(class_current) += std::fabs(DEM_NOSHAPES->data[r3][c3] - DEM_ORG->data[r3][c3]);
                                    }r3 = r-1; c3 = c-1;
                                    if(r3 >= 0 && c3 >=0 && r3 < DEM_ORG->nrRows() && c3 < DEM_ORG->nrCols())
                                    {
                                        depressionfeatureinteract.at(class_current) += std::fabs(DEM_NOSHAPES->data[r3][c3] - DEM_ORG->data[r3][c3]);
                                    }r3 = r+1; c3 = c-1;
                                    if(r3 >= 0 && c3 >=0 && r3 < DEM_ORG->nrRows() && c3 < DEM_ORG->nrCols())
                                    {
                                        depressionfeatureinteract.at(class_current) += std::fabs(DEM_NOSHAPES->data[r3][c3] - DEM_ORG->data[r3][c3]);
                                    }r3 = r-1; c3 = c+1;
                                    if(r3 >= 0 && c3 >=0 && r3 < DEM_ORG->nrRows() && c3 < DEM_ORG->nrCols())
                                    {
                                        depressionfeatureinteract.at(class_current) += std::fabs(DEM_NOSHAPES->data[r3][c3] - DEM_ORG->data[r3][c3]);
                                    }r3 = r+1; c3 = c+1;
                                    if(r3 >= 0 && c3 >=0 && r3 < DEM_ORG->nrRows() && c3 < DEM_ORG->nrCols())
                                    {
                                        depressionfeatureinteract.at(class_current) += std::fabs(DEM_NOSHAPES->data[r3][c3] - DEM_ORG->data[r3][c3]);
                                    }
                                 }
                             }
                         }

                     }
                 }

                 for(int i =0; i < depressionvols.size(); i++)
                 {
                     if(depressionvols.at(i) > 50.0*50.0*2.5) //cutoff point for what are real depressions
                     {
                         if(DEM_NOSHAPES)
                         {
                             if(depressionfeatureinteract.at(i) < 0.01)
                             {
                                 depressionvols.at(i) = 0.0;
                             }
                         }else
                         {
                             depressionvols.at(i) = 0.0;
                         }
                     }
                     depressioncvols.push_back(0.0);
                 }
                 for(int r = 0; r < DEMDID->data.nr_rows();r++)
                 {
                     for(int c = 0; c < DEMDID->data.nr_cols();c++)
                     {
                         if(!pcr::isMV(DEMDID->data[r][c]))
                         {
                             bool found = false;
                             int class_current = (int)DEMDID->data[r][c];

                             if(depressionvols.at(class_current) < 1e-12f)
                             {
                                 if(DEM_NOSHAPES)
                                 {
                                     if(depressionfeatureinteract.at(class_current) < 0.01)
                                     {
                                        DEMD->data[r][c] = 0.0;
                                        pcr::setMV(DEMDID->data[r][c]);
                                     }
                                 }else
                                 {
                                    DEMD->data[r][c] = 0.0;
                                    pcr::setMV(DEMDID->data[r][c]);
                                 }
                             }

                         }
                     }
                 }


                 FF::cout << "found depressions: " << depressioncvols.size()<< std::endl;
                 if(depressioncvols.size() == 0)
                 {
                    depressionids.push_back(0);
                    depressionvols.push_back(0.0);
                    depressioncvols.push_back(0.0);

                 }

            }));
        }


        if((do_q_neg || do_stormdrain) && !do_depressions)
        {
            FF::cout << "do depressions empty for neg q" <<std::endl;
            workers1.push_back(std::thread([dx,&depressioncvols,&depressionvols,&depressionids,&DEMDID,DEM,&DEMD,&DEMDM,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
            {
                DEMDM = new cTMapf(nullptr,  DEM->data.nr_cols(),DEM->data.nr_rows(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data5),DEM->projection(),"");
                DEMD = DEM->GetCopy();
                DEMDID =  DEM->GetCopy0();

                 int depressioncells = 0;
                 for(int r = 0; r < DEMDID->data.nr_rows();r++)
                 {
                     for(int c = 0; c < DEMDID->data.nr_cols();c++)
                     {
                        pcr::setMV(DEMDID->data[r][c]);

                     }
                 }

                depressionids.push_back(0);
                depressionvols.push_back(0.0);
                            
                depressioncvols.push_back(0.0);


            }));
        }

#ifdef FASTFLOOD_WASM
        FF::cout << "wait for threads to complete "<< std::endl;
		emscripten_run_script("increment_progress_fastflood()");
#endif
       // Looping every thread via for_each
       // Then, joins one by one, and this works like barrier
       std::for_each(workers1.begin(), workers1.end(), [](std::thread &t)
       {
           t.join();
       });





        //get new normalized flux network
        //#pragma omp parallel for collapse(2)
        if(dem_channelized != 0 && ip == 0)
        {
            //writeRaster(*DEM2CHAN,"dem_ff_fill.tif");
            std::cout << "channelized dem for grads" << std::endl;

            for (r = 0; r < DEM->Rows(); r++) {
                for (c = 0; c < DEM->Cols(); c++) {
                    if (!pcr::isMV(DEM->data[r][c])) {


                        float gx1 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_L));
                        float gx2 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_R));
                        float gy1 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_L));
                        float gy2 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_R));

                        float gtot = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));

                        if(do_depressions)
                        {
                            if(pcr::isMV(DEMDID->data[r][c]))
                            {
                                float gx12 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_L));
                                float gx22 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_R));
                                float gy12 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_L));
                                float gy22 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_R));
                                SLOPE->data[r][c] = std::max(1e-12f,std::fabs(gx12) + std::fabs(gx22) + std::fabs(gy12) + std::fabs(gy22));
                            }else
                            {
                                SLOPE->data[r][c] = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));
                            }
                        }else
                        {
                            SLOPE->data[r][c] = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));
                        }

                        GradX1->data[r][c] = gx1/gtot;
                        GradX2->data[r][c] = gx2/gtot;
                        GradY1->data[r][c] = gy1/gtot;
                        GradY2->data[r][c] = gy2/gtot;

                    } else
                    {
                        pcr::setMV(GradX1->data[r][c]);
                        pcr::setMV(GradX2->data[r][c]);
                        pcr::setMV(GradY1->data[r][c]);
                        pcr::setMV(GradY2->data[r][c]);
                    }
                }
            }
        }else
        {

            //writeRaster(*DEM2CHAN,"dem_ff_fill.tif");
            std::cout << "no channelized dem for grads "<< std::endl;

            for (int r = 0; r < DEM->Rows(); r++) {
                for (int c = 0; c < DEM->Cols(); c++) {
                    if (!pcr::isMV(DEM->data[r][c])) {

                        float gx1 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM2,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_L));
                        float gx2 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM2,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_R));
                        float gy1 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM2,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_L));
                        float gy2 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM2,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_R));

                        float gtot = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));


                        if(do_depressions)
                        {
                            if(pcr::isMV(DEMDID->data[r][c]))
                            {
                                float gx12 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_L));
                                float gx22 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_R));
                                float gy12 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_L));
                                float gy22 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_R));
                                SLOPE->data[r][c] = std::max(1e-12f,std::fabs(gx12) + std::fabs(gx22) + std::fabs(gy12) + std::fabs(gy22));
                            }else
                            {
                                SLOPE->data[r][c] = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));
                            }
                        }else
                        {
                            SLOPE->data[r][c] = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));
                        }
                        GradX1->data[r][c] = gx1/gtot;
                        GradX2->data[r][c] = gx2/gtot;
                        GradY1->data[r][c] = gy1/gtot;
                        GradY2->data[r][c] = gy2/gtot;

                    } else
                    {
                        pcr::setMV(GradX1->data[r][c]);
                        pcr::setMV(GradX2->data[r][c]);
                        pcr::setMV(GradY1->data[r][c]);
                        pcr::setMV(GradY2->data[r][c]);
                    }
                }
            }
        }


        if(do_chancor >=1 )
        {
             ApplyConditionalLDD4Correction(DEM_ORG,GradX1,GradX2,GradY1,GradY2,iters,dx,min_drainchancor,RAININ);
        }

         //here we can solve the splitting issue, we need to do the following:
         DoRiverSplitSearchAndPatch(DEM_ORG,DEM_ORG,DEM_ORG_MON,GradX1,GradX2,GradY1,GradY2,n_sources,val_sources,row_sources,col_sources,in_dx,iters);
         //should we do this beforehand? to improve output for the channel stuff?
         //and include the splitting in there as well?


        //do flow accumulation
        cTMapf * q;

        //MaskedRaster<float> raster_data4();
        cTMapf *RainN = new cTMapf(nullptr, DEM->data.nr_cols(), DEM->data.nr_rows(),DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data4),DEM->projection(),"");
        //MaskedRaster<float> raster_data5(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
        cTMapf *Vel = new cTMapf(nullptr, DEM->data.nr_cols(),  DEM->data.nr_rows(),DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data5),DEM->projection(),"");

        cTMapf * accN = nullptr;
        cTMapf * accqreal = nullptr;
        cTMapf * accaccqreal = nullptr;
        cTMapf * acc1 = nullptr;
        cTMapf * accacc1 = nullptr;
        cTMapf * accvel = nullptr;
        cTMapf * accsources = nullptr;
        cTMapf * accs = nullptr;

        std::cout << "Number of boundary conditions: " << n_sources << std::endl;
        if(n_sources > 0)
        {

            accsources = DEM->GetCopy0();
            for(int i = 0; i < n_sources; i++)
            {
                FF::cout << "in da:  "  << i << "  "  << val_sources[i] << "  " << row_sources[i] << "  " << col_sources[i]<<std::endl;
                int r = (int)(row_sources[i]);
                int c = (int)(col_sources[i]);
                float val = std::max(0.0f,val_sources[i]) /((in_dx * in_dx)/1000000.0); //convert km2 to pixels

                if(r >-1 && c> -1 && r < DEM->Rows() && c < DEM->Cols())
                {
                    FF::cout << "add " << val << std::endl;
                    accsources->data[r][c] += val;
                }
            }

            accs = SpreadFlowMD2(accsources,GradX1,GradX2,GradY1,GradY2,iters);
        }


        // vector container stores threads
        std::vector<std::thread> workers;
		
#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        if(do_qbase > 0  && qbase_set == false)
        {
            acc1 = SpreadFlowMD2(one,GradX1,GradX2,GradY1,GradY2,iters);

            QBASES = DEM->GetCopy0();

            float pix_acc_max = 0.0;

            if (accs) {
                for(int r = 0; r < DEM->Rows(); r++) {
                    for(int c = 0; c < DEM->Cols(); c++)
                    {
                        //int index = r* DEM->Cols() + c;
                        if(!pcr::isMV(DEM->data[r][c]))
                        {
                            n_pix ++;
                            pix_acc_max = std::max(accs->data[r][c] + acc1->data[r][c],pix_acc_max);
                        }
                    }
                }
            }else {
                for(int r = 0; r < DEM->Rows(); r++) {
                    for(int c = 0; c < DEM->Cols(); c++)
                    {
                        //int index = r* DEM->Cols() + c;
                        if(!pcr::isMV(DEM->data[r][c]))
                        {
                            n_pix ++;
                            pix_acc_max = std::max(acc1->data[r][c],pix_acc_max);
                        }
                    }
                }
            }


            std::cout << "baseflow active: "<< qbase_max << " " << pix_acc_max << std::endl;
            //baseflow source per pixel
            float val_qbase_dist = qbase_max/std::max(1.0f,pix_acc_max);
            if (accsources) {
                for(int r = 0; r < DEM->Rows(); r++) {
                    for(int c = 0; c < DEM->Cols(); c++)
                    {
                        //int index = r* DEM->Cols() + c;
                        if(!pcr::isMV(DEM->data[r][c]))
                        {
                            QBASES->data[r][c] = val_qbase_dist + accsources->data[r][c] * val_qbase_dist;

                        }
                    }
                }
            }else {
                for(int r = 0; r < DEM->Rows(); r++) {
                    for(int c = 0; c < DEM->Cols(); c++)
                    {
                        //int index = r* DEM->Cols() + c;
                        if(!pcr::isMV(DEM->data[r][c]))
                        {
                            QBASES->data[r][c] = val_qbase_dist;

                        }
                    }
                }
            }


            QBASE = SpreadFlowMD2(QBASES,GradX1,GradX2,GradY1,GradY2,iters);

            //writeRaster(*QBASE,"qb.tif");
            qbase_set = true;
            
        }

        


        if(!do_depressions && !do_q_neg && !do_stormdrain)
        {
            
            /*q = SpreadFlowMD2(Rain,GradX1,GradX2,GradY1,GradY2,iters);
            for (int r = 0; r < DEM->Rows(); r++) {
                    for (int c = 0; c < DEM->Cols(); c++) {
                        if (!pcr::isMV(DEM->data[r][c])) {

                            RainN->data[r][c] = Rain->data[r][c] * N->data[r][c];
                        }
                    }
                }

                accN = SpreadFlowMD2(RainN,GradX1,GradX2,GradY1,GradY2,iters);

                acc1 = SpreadFlowMD2(one,GradX1,GradX2,GradY1,GradY2,iters);
				
                accacc1 = SpreadFlowMD2(acc1,GradX1,GradX2,GradY1,GradY2,iters);
					for (int r = 0; r < DEM->Rows(); r++) {
                    for (int c = 0; c < DEM->Cols(); c++) {
                        if (!pcr::isMV(DEM->data[r][c])) {

                            //9.68*(accuflux.tif/1e6)**0.32
                            float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));

                            Vel->data[r][c] = conc_vel_mult*acc1->data[r][c] * std::pow(q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                        }
                    }
                }
				
                accvel = SpreadFlowMD2(Vel,GradX1,GradX2,GradY1,GradY2,iters);*/

            workers.push_back(std::thread([&q,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
            {
				
                FF::cout << "do q calc " << iters << std::endl;
                q = SpreadFlowMD2(Rain,GradX1,GradX2,GradY1,GradY2,iters);

            }));

            workers.push_back(std::thread([&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
            {
				
				for (int r = 0; r < DEM->Rows(); r++) {
                    for (int c = 0; c < DEM->Cols(); c++) {
                        if (!pcr::isMV(DEM->data[r][c])) {

                            RainN->data[r][c] = Rain->data[r][c] * N->data[r][c];
                        }
                    }
                }

                accN = SpreadFlowMD2(RainN,GradX1,GradX2,GradY1,GradY2,iters);


            }));
            

            workers.push_back(std::thread([one,DEM,&acc1,&accacc1,N,GradX1,GradX2,GradY1,GradY2,iters]()
            {
				
				
                if(acc1 == nullptr)
                {
                    acc1 = SpreadFlowMD2(one,GradX1,GradX2,GradY1,GradY2,iters);
                }
				
                accacc1 = SpreadFlowMD2(acc1,GradX1,GradX2,GradY1,GradY2,iters);
				
				

            }));

            std::for_each(workers.begin(), workers.end(), [](std::thread &t)
                {
                    t.join();
                });
                workers.clear();

            workers.push_back(std::thread([&Vel,&acc1,&QBASE,qbase_set,QSS_PREV,conc_vel_mult, &q,SLOPE,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
            {

				 if(qbase_set)
                    {
                        for (int r = 0; r < DEM->Rows(); r++) {
                            for (int c = 0; c < DEM->Cols(); c++) {
                                if (!pcr::isMV(DEM->data[r][c])) {

                                    //9.68*(accuflux.tif/1e6)**0.32
                                    float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                    Vel->data[r][c] = conc_vel_mult*q->data[r][c]  *std::pow(QBASE->data[r][c] + QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                }
                            }
                        }
                    }else
                    {
                        for (int r = 0; r < DEM->Rows(); r++) {
                            for (int c = 0; c < DEM->Cols(); c++) {
                                if (!pcr::isMV(DEM->data[r][c])) {

                                    //9.68*(accuflux.tif/1e6)**0.32
                                    float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                    Vel->data[r][c] = conc_vel_mult*q->data[r][c] *std::pow(QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                }
                            }
                        }

                    }

                accvel = SpreadFlowMD2(Vel,GradX1,GradX2,GradY1,GradY2,iters);


            }));

            workers.push_back(std::thread([&DEMDID,conc_vel_mult,qbase_set,&QBASE,QSS_PREV,&SLOPE,&acc1, &q,&depressioncvols,&depressionvols,&Vel,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters, &accaccqreal]()
            {

                accaccqreal = SpreadFlowMD2(q,GradX1,GradX2,GradY1,GradY2,iters);

            }));

        }else
        {
			
			if(do_depressions && (!do_q_neg) && n_culverts == 0 && ! do_stormdrain)
            {
                workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,&q,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
                {
                    FF::cout << "do q calc depr " << iters << std::endl;

                    q = SpreadFlowMD2(Rain,GradX1,GradX2,GradY1,GradY2,iters);
                    //writeRaster(*q,"qtestorg.tif");
                }));
               
                    
                    workers.push_back(std::thread([accsources,duration,&ReservoirID,&accqreal,&DEMDID,&depressioncvols,&depressionvols,&q,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters,dx]()
                    {
                        ReservoirID = nullptr;
                        accqreal =SpreadFlowMD2Depression(Rain,GradX1,GradX2,GradY1,GradY2,iters,DEMDID,depressionvols,depressioncvols,true,duration * 60.0 * 60.0,accsources,dx);


                    }));

                

                workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
                {
                    
                    for (int r = 0; r < DEM->Rows(); r++) {
                        for (int c = 0; c < DEM->Cols(); c++) {
                            if (!pcr::isMV(DEM->data[r][c])) {

                                RainN->data[r][c] = Rain->data[r][c] * N->data[r][c];
                            }
                        }
                    }
                    
                    accN = SpreadFlowMD2(RainN,GradX1,GradX2,GradY1,GradY2,iters);


                }));

                workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,one,DEM,&acc1,&accacc1,N,GradX1,GradX2,GradY1,GradY2,iters]()
                {
                    
                    if(acc1 == nullptr)
                    {   
                        acc1 = SpreadFlowMD2(one,GradX1,GradX2,GradY1,GradY2,iters);
                    }
                    accacc1 = SpreadFlowMD2(acc1,GradX1,GradX2,GradY1,GradY2,iters);
                    

                }));
                
                std::for_each(workers.begin(), workers.end(), [](std::thread &t)
                {
                    t.join();
                });
                workers.clear();

                workers.push_back(std::thread([&DEMDID,conc_vel_mult,qbase_set,&QBASE,QSS_PREV,&acc1,&q,SLOPE,&depressioncvols,&depressionvols,&Vel,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
                {
                    if(qbase_set)
                    {
                        for (int r = 0; r < DEM->Rows(); r++) {
                            for (int c = 0; c < DEM->Cols(); c++) {
                                if (!pcr::isMV(DEM->data[r][c])) {

                                    //9.68*(accuflux.tif/1e6)**0.32
                                    float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                    Vel->data[r][c] = conc_vel_mult*q->data[r][c]  *std::pow(QBASE->data[r][c] + QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                }
                            }
                        }
                    }else
                    {
                        for (int r = 0; r < DEM->Rows(); r++) {
                            for (int c = 0; c < DEM->Cols(); c++) {
                                if (!pcr::isMV(DEM->data[r][c])) {

                                    //9.68*(accuflux.tif/1e6)**0.32
                                    float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                    Vel->data[r][c] = conc_vel_mult*q->data[r][c]   *std::pow(QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                }
                            }
                        }

                    }
                    accvel = SpreadFlowMD2(Vel,GradX1,GradX2,GradY1,GradY2,iters);
                    


                }));

                workers.push_back(std::thread([&DEMDID,conc_vel_mult,qbase_set,&QBASE,QSS_PREV,&SLOPE,&acc1, &q,&depressioncvols,&depressionvols,&Vel,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters,&accaccqreal]()
                {

                    accaccqreal = SpreadFlowMD2(q,GradX1,GradX2,GradY1,GradY2,iters);

                }));

            }else{

                
                 if(n_culverts > 0 || do_stormdrain)
                {
                     std::vector<float> culvert_ffactors(CULVERT_Q.size(), 1.0f);
                     if(acc1 == nullptr)
                     {
                         acc1 = SpreadFlowMD2Culvert(one,GradX1,GradX2,GradY1,GradY2,iters, CULVERTID,CULVERT_FRAC,CULVERT_RS,CULVERT_CS, CULVERT_RT, CULVERT_CT);
                     }
                     for(size_t i = 0; i < CULVERT_Q.size(); i++)
                     {
                         int r = CULVERT_RS[i];
                         int c = CULVERT_CS[i];
                         float da = 0.0f;
                         if(r >= 0 && r < DEM->Rows() && c >= 0 && c < DEM->Cols())
                         {
                             float val = acc1->data[r][c];
                             if(!pcr::isMV(val))
                             {
                                 da = val * dx * dx / 1000000.0f;
                             }
                         }
                         culvert_ffactors[i] = 1.0f + std::pow(da / 100.0f, 1.28f);
                     }

                      workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,&q,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters, &CULVERTID,&CULVERT_FRAC,&CULVERT_RS,&CULVERT_CS, &CULVERT_RT, &CULVERT_CT]()
                    {
                          FF::cout << "do q calc depr culv" << iters << std::endl;

                        q = SpreadFlowMD2Culvert(Rain,GradX1,GradX2,GradY1,GradY2,iters, CULVERTID,CULVERT_FRAC,CULVERT_RS,CULVERT_CS, CULVERT_RT, CULVERT_CT);
                        //writeRaster(*q,"dischargetest.tif");
                        //writeRaster(*CULVERTID,"culvertstest.tif");

                    }));


                    workers.push_back(std::thread([duration,&ReservoirID,&accqreal,&DEMDID, &SDC, &SDCTOT, &depressioncvols,&depressionvols,&q,&accN, &CULVERTID,&reservoircvols,&reservoir_inflow,&reservoirpassvols, &CULVERT_Q, &culvert_ffactors, &CULVERT_FRAC,&CULVERT_QTOT,&reservoirvols,&CULVERT_RS,&CULVERT_CS, &CULVERT_RT, &CULVERT_CT, DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
                    {
                        FF::cout << "do culvert flow acc" << std::endl;
                        accqreal = SpreadFlowMD2DepressionReservoirCulvert(Rain,GradX1,GradX2,GradY1,GradY2,iters,DEMDID,depressionvols,depressioncvols,true,ReservoirID, reservoircvols,reservoir_inflow,reservoirpassvols,reservoirvols,CULVERTID, CULVERT_Q, culvert_ffactors, CULVERT_FRAC,CULVERT_QTOT,CULVERT_RS,CULVERT_CS, CULVERT_RT, CULVERT_CT,duration * 60.0 * 60.0, SDC, SDCTOT);


                    }));

                    workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters, &CULVERTID,&CULVERT_FRAC,&CULVERT_RS,&CULVERT_CS, &CULVERT_RT, &CULVERT_CT]()
                    {
                        
                        for (int r = 0; r < DEM->Rows(); r++) {
                            for (int c = 0; c < DEM->Cols(); c++) {
                                if (!pcr::isMV(DEM->data[r][c])) {

                                    RainN->data[r][c] = Rain->data[r][c] * N->data[r][c];
                                }
                            }
                        }
                        
                        accN = SpreadFlowMD2Culvert(RainN,GradX1,GradX2,GradY1,GradY2,iters, CULVERTID,CULVERT_FRAC,CULVERT_RS,CULVERT_CS, CULVERT_RT, CULVERT_CT);


                    }));


                    workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,one,DEM,&acc1,&accacc1,N,GradX1,GradX2,GradY1,GradY2,iters, &CULVERTID,&CULVERT_FRAC,&CULVERT_RS,&CULVERT_CS, &CULVERT_RT, &CULVERT_CT]()
                    {
                        
                        if(acc1 == nullptr)
                        {   
                            acc1 = SpreadFlowMD2Culvert(one,GradX1,GradX2,GradY1,GradY2,iters, CULVERTID,CULVERT_FRAC,CULVERT_RS,CULVERT_CS, CULVERT_RT, CULVERT_CT);
                        }
                        accacc1 = SpreadFlowMD2Culvert(acc1,GradX1,GradX2,GradY1,GradY2,iters, CULVERTID,CULVERT_FRAC,CULVERT_RS,CULVERT_CS, CULVERT_RT, CULVERT_CT);
                    
                    }));
                    
                    std::for_each(workers.begin(), workers.end(), [](std::thread &t)
                    {
                        t.join();
                    });
                    workers.clear();
                    workers.push_back(std::thread([&DEMDID,conc_vel_mult,qbase_set,&QBASE,QSS_PREV,&SLOPE,&acc1, &q,&depressioncvols,&depressionvols,&Vel,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters, &CULVERTID,&CULVERT_FRAC,&CULVERT_RS,&CULVERT_CS, &CULVERT_RT, &CULVERT_CT]()
                    {
                        if(qbase_set)
                        {
                            for (int r = 0; r < DEM->Rows(); r++) {
                                for (int c = 0; c < DEM->Cols(); c++) {
                                    if (!pcr::isMV(DEM->data[r][c])) {

                                        //9.68*(accuflux.tif/1e6)**0.32
                                        float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                        Vel->data[r][c] = conc_vel_mult*q->data[r][c]  *std::pow(QBASE->data[r][c] + QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                    }
                                }
                            }
                        }else
                        {
                            for (int r = 0; r < DEM->Rows(); r++) {
                                for (int c = 0; c < DEM->Cols(); c++) {
                                    if (!pcr::isMV(DEM->data[r][c])) {

                                        //9.68*(accuflux.tif/1e6)**0.32
                                        float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                        Vel->data[r][c] = conc_vel_mult*q->data[r][c]  *std::pow(QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                    }
                                }
                            }

                        }
                        
                        accvel = SpreadFlowMD2Culvert(Vel,GradX1,GradX2,GradY1,GradY2,iters, CULVERTID,CULVERT_FRAC,CULVERT_RS,CULVERT_CS, CULVERT_RT, CULVERT_CT);
                        

                    }));

                    workers.push_back(std::thread([&DEMDID,conc_vel_mult,qbase_set,&QBASE,QSS_PREV,&SLOPE,&acc1, &q,&depressioncvols,&depressionvols,&Vel,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters,&accaccqreal]()
                    {

                        accaccqreal = SpreadFlowMD2(q,GradX1,GradX2,GradY1,GradY2,iters);

                    }));


                }else
                {
                    workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,&q,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
                    {
                        FF::cout << "do q calc depr2" << iters << std::endl;

                        q = SpreadFlowMD2(Rain,GradX1,GradX2,GradY1,GradY2,iters);


                    }));


                    workers.push_back(std::thread([duration,&accqreal,&DEMDID,&depressioncvols,&depressionvols,&q,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters,&ReservoirID,& reservoircvols,&reservoir_inflow,&reservoirpassvols,&reservoirvols]()
                    {
                        FF::cout << "do non-culvert flow acc with reservoir and depression" << std::endl;
                        accqreal = SpreadFlowMD2DepressionReservoir(Rain,GradX1,GradX2,GradY1,GradY2,iters,DEMDID,depressionvols,depressioncvols,true,ReservoirID, reservoircvols,reservoir_inflow,reservoirpassvols,reservoirvols,duration * 60.0 * 60.0);


                    }));

                       workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,&accN,DEM,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
                    {
                        
                        for (int r = 0; r < DEM->Rows(); r++) {
                            for (int c = 0; c < DEM->Cols(); c++) {
                                if (!pcr::isMV(DEM->data[r][c])) {

                                    RainN->data[r][c] = Rain->data[r][c] * N->data[r][c];
                                }
                            }
                        }
                        
                        accN = SpreadFlowMD2(RainN,GradX1,GradX2,GradY1,GradY2,iters);


                    }));


                    workers.push_back(std::thread([&DEMDID,&depressioncvols,&depressionvols,one,DEM,&acc1,&accacc1,N,GradX1,GradX2,GradY1,GradY2,iters]()
                    {
                        
                        if(acc1 == nullptr)
                        {   
                            acc1 = SpreadFlowMD2(one,GradX1,GradX2,GradY1,GradY2,iters);
                        }
                        accacc1 = SpreadFlowMD2(acc1,GradX1,GradX2,GradY1,GradY2,iters);
                    
                    }));
                    
                    std::for_each(workers.begin(), workers.end(), [](std::thread &t)
                    {
                        t.join();
                    });
                    workers.clear();
                    workers.push_back(std::thread([&DEMDID,conc_vel_mult,qbase_set,&QBASE,QSS_PREV,&SLOPE,&acc1, &q,&depressioncvols,&depressionvols,&Vel,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters]()
                    {
                        if(qbase_set)
                        {
                            for (int r = 0; r < DEM->Rows(); r++) {
                                for (int c = 0; c < DEM->Cols(); c++) {
                                    if (!pcr::isMV(DEM->data[r][c])) {

                                        //9.68*(accuflux.tif/1e6)**0.32
                                        float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                        Vel->data[r][c] = conc_vel_mult*q->data[r][c]  * std::pow(QBASE->data[r][c] + QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                    }
                                }
                            }
                        }else
                        {
                            for (int r = 0; r < DEM->Rows(); r++) {
                                for (int c = 0; c < DEM->Cols(); c++) {
                                    if (!pcr::isMV(DEM->data[r][c])) {

                                        //9.68*(accuflux.tif/1e6)**0.32
                                        float slope = SLOPE->data[r][c];//0.5*std::max(1e-12f,std::fabs(GradX1->data[r][c]) + std::fabs(GradX2->data[r][c]) + std::fabs(GradY1->data[r][c]) + std::fabs(GradY2->data[r][c]));
                                        Vel->data[r][c] = conc_vel_mult*q->data[r][c]  * std::pow(QSS_PREV->data[r][c] + q->data[r][c],0.25) * sqrt(slope)/std::max(0.001f,N->data[r][c]);
                                    }
                                }
                            }

                        }
                        
                        accvel = SpreadFlowMD2(Vel,GradX1,GradX2,GradY1,GradY2,iters);
                        

                    }));

                    workers.push_back(std::thread([&DEMDID,conc_vel_mult,qbase_set,&QBASE,QSS_PREV,&SLOPE,&acc1, &q,&depressioncvols,&depressionvols,&Vel,DEM,&accvel,RainN,Rain,N,GradX1,GradX2,GradY1,GradY2,iters,&accaccqreal]()
                    {

                        accaccqreal = SpreadFlowMD2(q,GradX1,GradX2,GradY1,GradY2,iters);

                    }));
                }

             




            }
            
        }


        // Looping every thread via for_each
        // Then, joins one by one, and this works like barrier
        std::for_each(workers.begin(), workers.end(), [](std::thread &t)
        {
            t.join();
        });


        FF::cout << "workers done" <<std::endl;        
		
		
		//printf("depressed workers done \n");
#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        
        cTMapf * Qin = nullptr;



        if(do_new_crosssections > 0 )
        {
            Qin = DEM->GetCopy0();
        }

	    /*writeRaster(*DEM2,"dem2test.tif");
	    writeRaster(*GradX1,"gradx1.tif");
	    writeRaster(*GradX2,"gradx2.tif");
	    writeRaster(*GradY1,"grady1.tif");
	    writeRaster(*GradY2,"grady2.tif");
	    writeRaster(*Rain,"rain.tif");
	    writeRaster(*q,"qtest.tif");

	    writeRaster(*accqreal,"qrtest.tif");*/


        //get average distance, N and velocity
        //#pragma omp parallel for collapse(2)
        for (r = 0; r < DEM->Rows(); r++) {
            for (c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(DEM->data[r][c])) {

                    RainN->data[r][c] = accN->data[r][c]/std::max(0.00001f,q->data[r][c]);
                    Vel->data[r][c] = accvel->data[r][c]/std::max(0.00001f,accaccqreal->data[r][c]);
                    accacc1->data[r][c] = dx * accacc1->data[r][c]/std::max(1.0f,acc1->data[r][c]);

                    //get steady-state compensation
                    //Map h =  pow( q* N /max(0.001,sqrt(max(0.001,Slope(dem)))),q * 0.0 + (5.0/6.0));
                    float i_dx = dx;
                    float i_L = accacc1->data[r][c] * 3.0/2.0;//initial guess for iteration based on average distance of circular segment to center of circle
                    float i_n = acc1->data[r][c];
                    float i_b = 1.0; //start iteration with linear
                    float i_Lav = accacc1->data[r][c];//average travel distance to outlet

                    bool do_break = false;
                    //this should typically only take a small amount of iterations (<5)
                    for(int i = 0; i < 10; i++)
                    {

                        float Lnew = ((2.0 + i_b) *i_Lav)/(1.0 + i_b); // N
                        float bnew = (-std::log(i_dx/i_L) + std::log(1.0/i_n))/std::log(i_dx/i_L); // 

                        if(((Lnew - i_L) < 0.001) && ((bnew - i_b) < 0.001))
                        {
                            do_break = true;
                        }
                        i_L = Lnew;
                        i_b = bnew;
                        if(do_break)
                        {
                            break;
                        }
                    }

                    if(!std::isnormal(i_b) || !std::isnormal(i_L))
                    {
                        i_b = 1.0;
                        i_L = accacc1->data[r][c] * 3.0/2.0;
                    }

                    i_b = std::max(0.00001f,i_b);
                    i_L = std::max(dx,i_L);

                    //our velocity was calculated using power law depth estimate, so should we rescale based on the actual q estimate?

                    float ss_comp_x = 2.5 * duration * 60.0 * 60.0 * Vel->data[r][c]/i_L; //relative event duration
                    float ss_comp;
                    if(ss_comp_x < 1.0)
                    {
						
						float term_diffusivity = std::min(1.0,std::exp(-diffusivity_coefficient*std::log((1.0/std::min(1.0f,std::max(0.0001f,ss_comp_x))))));
					
                        ss_comp = std::min(1.0,(1.0 - std::pow(1.0 - ss_comp_x,1.0 + i_b))) * term_diffusivity; //1 - (1 - x)^(1 + b)
                    }else
                    {
                        ss_comp = 1.0;
                    }
                    if(n_sources > 0)
                    {
                        ss_comp = (accs->data[r][c] * 1.0 + acc1->data[r][c] * ss_comp)/std::max(1e-12f,accs->data[r][c]  + acc1->data[r][c]);
                    }
 

                    if(do_qdetail_out > 0)
                    {
                        float weighth = 1.0;
                        if(piters > 1)
                        {

                            if(do_depressions || do_q_neg)
                            {
                                float w_cell = accqreal->data[r][c];
                                if (CULVERTID != nullptr && CULVERTID->data[r][c] >= 0.0f)
                                {
                                    w_cell = std::max(w_cell, q->data[r][c]);
                                }
                                q_detail_w->data[r][c] += w_cell;
                                weighth = w_cell;
                            }else{
                                q_detail_w->data[r][c] += q->data[r][c];
                                weighth = q->data[r][c];
                            }

                        }

                        int index = r* DEM->Cols() + c;
                        ret_b[index] += weighth * i_b;
                        ret_l[index] += weighth * i_L;
                        ret_avg_vel[index] += weighth * Vel->data[r][c];
                        ret_ss_comp[index] += weighth * ss_comp;
                        
                    }

                    if (use_local_ss)
                    {
                        float weighth = 1.0;
                        if (piters > 1)
                        {
                            if (do_depressions || do_q_neg)
                            {
                                float w_cell = accqreal->data[r][c];
                                if (CULVERTID != nullptr && CULVERTID->data[r][c] >= 0.0f)
                                {
                                    w_cell = std::max(w_cell, q->data[r][c]);
                                }
                                weighth = w_cell;
                            }
                            else
                            {
                                weighth = q->data[r][c];
                            }
                            int index = r * DEM->Cols() + c;
                            local_ss_weight[index] += weighth;
                        }
                        int index = r * DEM->Cols() + c;
                        local_ss_comp[index] += weighth * ss_comp;
                    }

                    float SlopeX = UF2D_Derivative(DEM,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_LR);
                    float SlopeY = UF2D_Derivative(DEM,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_LR);
                    float Slope = std::fabs(SlopeX) + std::fabs(SlopeY);

                    //better estimate for dx in flow width?
                    //I think we can!
                    float q_add_base = 0.0;
                    if(qbase_set && ip == 0)
                    {
                        q_add_base = QBASE->data[r][c];
                    }

                    if(do_new_crosssections > 0 )
                    {
                        if(!do_depressions  && !do_q_neg && !do_stormdrain)
                        {
                            float qhere = ss_comp*q->data[r][c];
                            if(timeslice > 0.0)
                            {
                                qhere = GetDischargeCurveAtTime(qhere,ss_comp,i_L,Vel->data[r][c],i_b,0.0,duration,timeslice);
                            }

                            Qin->data[r][c] = q_add_base + qhere;
                            QSS_PREV->data[r][c] += ss_comp * q->data[r][c];
                        }else
                        {
                            float qhere = ss_comp*accqreal->data[r][c];
                            if(timeslice > 0.0)
                            {
                                qhere = GetDischargeCurveAtTime(qhere,ss_comp,i_L,Vel->data[r][c],i_b,0.0,duration,timeslice);
                            }

                            Qin->data[r][c] = q_add_base + qhere;
                            QSS_PREV->data[r][c] += ss_comp * accqreal->data[r][c];
                        }
                    }else
                    {
                        //fill with h, but compensate for any previous flow (due to velocity increase, q increases non-linearly with h)
                        if(!do_depressions  && !do_q_neg && !do_stormdrain)
                        {
                            float qhere = ss_comp*q->data[r][c];
                            if(timeslice > 0.0)
                            {
                                qhere = GetDischargeCurveAtTime(qhere,ss_comp,i_L,Vel->data[r][c],i_b,0.0,duration,timeslice);
                            }

                            float q_old = QSS_PREV->data[r][c];
                            float h_old = std::pow((1.0/(dx)) * ss_comp * q_old * N->data[r][c]/std::max(0.001f,Slope),3.0/5.0);
                            one->data[r][c] = std::max(0.0,-h_old + std::pow((1.0/(dx)) * (qhere + ss_comp * (q_add_base + q_old)) * N->data[r][c]/std::max(0.001f,Slope),3.0/5.0));

                            WH_PREV->data[r][c] += one->data[r][c];
                            QSS_PREV->data[r][c] += ss_comp * q->data[r][c];
                        }else
                        {
                            float qhere = ss_comp*accqreal->data[r][c];
                            if(timeslice > 0.0)
                            {
                                qhere = GetDischargeCurveAtTime(qhere,ss_comp,i_L,Vel->data[r][c],i_b,0.0,duration,timeslice);
                            }
                            float q_old = QSS_PREV->data[r][c];
                            float h_old = std::pow((1.0/(dx)) * ss_comp * q_old * N->data[r][c]/std::max(0.001f,Slope),3.0/5.0);
                            one->data[r][c] = std::max(0.0,-h_old + std::pow((1.0/(dx)) *(qhere + ss_comp * (q_add_base + q_old)) * N->data[r][c]/std::max(0.001f,Slope),3.0/5.0));
                            WH_PREV->data[r][c] += one->data[r][c];
                            QSS_PREV->data[r][c] += ss_comp *accqreal->data[r][c];
                        }


                    }
                    
                  
                }
            }
        }

        //estimate flow width based on logarithmic cascaded cross-sectional estimates from elevation dataset
        //min width is still dx to prevent local elevation noise from ruining the quality of the estimate
        //DEM already includes multi-step with WH_prev
        //increasing Q? linearly interpolate between steps to find actual width?
        //volume transform? traditional application of method to find width of re-distribution?
        //Trapezoid steps? Exact solutions for trapezoids to find the perfect H/W value related to the height
        //This would give us an actual volume to re-distribute, instead of using manning
        //this would avoid non-linear correction as well
        //I need to visualize these results for the paper though

        cTMapf * FlowH_Pre = nullptr;

        int iters_cascade = 6; // 3^6 = 729 pixels width potentially, that should be anough for most use-cases

        if(do_new_crosssections > 0 )
        {
            //do a backwater effect propagation first

            //writeRaster(*SLOPE,"slope_chan.tif");

            std::vector<cTMapf *> resfwpd = GetFlowCrossSectionalVolumeFromRiverDischargePreDistribute2(DEM_ORG,Qin,N, do_chan?CHANW:0,do_chan? CHANH:0, SLOPE, iters_cascade, quality_p, FlowWL,FlowWR);

            FF::cout << "copy results " << std::endl;
            FlowH_Pre = resfwpd[0];

            //writeRaster(*DEM_ORG,"dem_org.tif");
            //writeRaster(*FlowH_Pre,"wh_pre.tif");
            //pre-redistribution using the exponentially cascading structure
            //FlowH_Pre = PreDistribute(DEM,FlowA,FlowHCS, do_chan?CHANW:0,do_chan? CHANH:0, iters_cascade);

            for (r = 0; r < DEM->Rows(); r++) {
                for (c = 0; c < DEM->Cols(); c++) {
                    if (!pcr::isMV(DEM->data[r][c])) {

                        //our old re-distribution scheme starts with all water in one gricell, so we divide by dx
                        //remove any volume we already distributed with the pre-distribution scheme
                        one->data[r][c] = FlowH_Pre->data[r][c];//(FlowA->data[r][c]/dx);
                        WH_PREV->data[r][c] += one->data[r][c];
                        
                    }
                }
            }

            FF::cout << "delete results " << std::endl;
            SAFE_DELETE(Qin);
            
            for(int i = 0; i < resfwpd.size(); i++)
            {
                SAFE_DELETE(resfwpd[i]);
            }

            FF::cout << "done" << std::endl;
        }



        // based on these volumes, do the re-distribution


         for (r = 0; r < DEM->Rows(); r++) {
            for (c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(DEM->data[r][c])) {


                    if(!std::isnormal(one->data[r][c]))
                    {
                        //one->data[r][c] = 0.0;
                        //WH_PREV->data[r][c] = 0.0;
                    }
                    if(do_depressions )
                    { 

                        //is this cell in a depression?
                        if(!pcr::isMV(DEMDID->data[r][c]))
                        {
                            int id = DEMDID->data[r][c];
                            //FF::cout << "access depression at " << id << std::endl;
                            if(!(id >= depressioncvols.size() || id < 0))
                            {
                                if(depressionvols.at(id) > 1e-12f)
                                {
                                    //to what extent is this depression filled?
                                    float frac = depressioncvols.at(id)/ std::max(0.01f,depressionvols.at(id));
                                    one->data[r][c] += DEMD->data[r][c] * frac;
                                    WH_PREV->data[r][c] += DEMD->data[r][c] * frac;

                                }
                            }
                            
                        }
                    }
					
					if(do_chan)
					{
						//potential storage in channel?
						float store_max = CHANW->data[r][c] * CHANH->data[r][c] * in_dx;
						float store_left = std::max(0.0f,store_max -ChanStorage->data[r][c]);
						
						float store_lefth = store_left/(in_dx*in_dx);
						float moveh = std::min(store_lefth,one->data[r][c]);
						//move water from our current one map (containing inverted water height) into the channel
						if (!export_chanwh)
						{
						    ChanStorage->data[r][c] += moveh * in_dx*in_dx;
						    one->data[r][c] -= moveh;
						}

					}

                }
            }
        }


            //writeRaster(*one,"wh_pre_minchan.tif");

#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        cTMapf * hfinal;
    //itersp = 0;
    //itersp = 0;

        if(!do_ocean)
        {

            if(corr)
            {
                if(!do_depressions)
                {
                    //if we dont do depressions, we better use the corrected dem as it improves throughflow in low resolution rough terrain
                   hfinal = DiffusiveQMaxWave(DEM_ORG,one,itersp,0.15);
                }else
                {
                    hfinal = DiffusiveQMaxWave(DEM_ORG,one ,itersp,0.15);
                }
            }else
            {
                if (n_culverts > 0)
                {
                    hfinal = DiffusiveMaxWaveCulverts(DEM_ORG, one, itersp, 0.15,
                                                      CULVERTID, CULVERT_Q, CULVERT_FRAC, CULVERT_QTOT,
                                                      CULVERT_RS, CULVERT_CS, CULVERT_RT, CULVERT_CT,
                                                      culvert_ffactors);
                }
                else if(!do_depressions)
                {
                    //if we dont do depressions, we better use the corrected dem as it improves throughflow in low resolution rough terrain
                   hfinal = DiffusiveMaxWave(DEM_ORG,one ,itersp,0.15);
                }else
                {
                    hfinal = DiffusiveMaxWave(DEM_ORG,one,itersp,0.15);
                }
            }
        }else
        {
            if(!do_depressions)
            {
                   for (int r = 0; r < DEM->Rows(); r++) {
                    for (int c = 0; c < DEM->Cols(); c++) {
                        if (!pcr::isMV(DEM->data[r][c])) {

                            WHODEM->data[r][c] = DEM_ORG->data[r][c] + WHO->data[r][c];

                        }
                    }
                }
            }else
            {
                for (int r = 0; r < DEM->Rows(); r++) {
                    for (int c = 0; c < DEM->Cols(); c++) {
                        if (!pcr::isMV(DEM->data[r][c])) {

                            WHODEM->data[r][c] = DEM_ORG->data[r][c] + WHO->data[r][c];

                        }
                    }
                }
            }

            //water redistribution for pressure advection
            
            if(corr)
            {
                if(!do_depressions)
                {
                    //if we dont do depressions, we better use the corrected dem as it improves throughflow in low resolution rough terrain
                   hfinal = DiffusiveQMaxWave(WHODEM,one,itersp,0.15);
                }else
                {
                    hfinal = DiffusiveQMaxWave(WHODEM,one ,itersp,0.15);
                }
            }else
            {
                if (n_culverts > 0) {
                    hfinal = DiffusiveMaxWaveCulverts(WHODEM, one, itersp, 0.15,
                                                      CULVERTID, CULVERT_Q, CULVERT_FRAC, CULVERT_QTOT,
                                                      CULVERT_RS, CULVERT_CS, CULVERT_RT, CULVERT_CT,
                                                      culvert_ffactors);
                }else {

                    hfinal = DiffusiveMaxWave(WHODEM,one ,itersp,0.15);
                }
            }
        }


        
#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        //return flow height
        //get average distance, N and velocity
        //#pragma omp parallel for collapse(2)
        for(r = 0; r < DEM->Rows(); r++) {
            for(c = 0; c < DEM->Cols(); c++)
            {
                if(!pcr::isMV(DEM->data[r][c]))
                {

                    WH_TOT->data[r][c] += hfinal->data[r][c];
                }
            }
        }

        //writeRaster(*WH_TOT,"whpost.tif");
            
        if(do_output_vel)
        {   


            FF::cout << "do vel output " << std::endl;
            cTMapf * whdem;
            if(!do_depressions)
            {
                whdem  = DEM2->GetCopy();
            }else
            {
                whdem  = DEM->GetCopy();
            }
                    
            for(int r = 0; r < DEM->Rows(); r++) {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                    
                    if(!pcr::isMV(whdem->data[r][c]))
                    {
                        whdem->data[r][c] += WH_TOT->data[r][c];
                    }
                }
            }

            //smooth the slopes for velocity directions
            cTMapf * whdem2 =whdem;
            whdem = WindowAverageWR(whdem, WH_TOT,in_dx/50.0);

            //get river specific gradient
            //track encountered average discharge, vs current local discharge
            //then use this to interpolate to the floodplain if required
            
            cTMapf * sxg = DEM->GetCopy();
            cTMapf * syg = DEM->GetCopy();
            for(int r = 0; r < DEM->Rows(); r++) {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                    
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        sxg->data[r][c] = GradX1->data[r][c] + GradX2->data[r][c];
                        syg->data[r][c] = GradY1->data[r][c] + GradY2->data[r][c];

                        

                    }else
                    {
                        sxg->data[r][c] = 0.0;
                        syg->data[r][c] = 0.0;
                    }
                }
            }

            //cTMapf * gradx_avw = WindowDistributeWR(sxg, FlowW, QSS_PREV,in_dx*2.0);
            //cTMapf * grady_avw = WindowDistributeWR(syg, FlowW, QSS_PREV,in_dx*2.0);
 
            cTMapf * gradx_avw = WindowDistributeWRDir(sxg,sxg,syg,FlowWL,FlowWR, QSS_PREV,in_dx*0.5);
            cTMapf * grady_avw = WindowDistributeWRDir(syg,sxg,syg,FlowWL,FlowWR, QSS_PREV,in_dx*0.5);
 
            for(int r = 0; r < DEM->Rows(); r++) {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                      int index = r* DEM->Cols() + c;


                     if(!pcr::isMV(DEM->data[r][c]))
                    {
                        float  dem =whdem->data[r][c];
                        float dem_xl = dem;
                        float dem_xr = dem;
                        float dem_yl = dem;
                        float dem_yr = dem;

                        //slope
                        if(c > 0)
                        {
                            dem_xl = pcr::isMV(DEM->data[r][c-1])?dem:whdem->data[r][c-1];
                        }else {
                            dem_xl = dem;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            dem_xr = pcr::isMV(DEM->data[r][c+1])?dem:whdem->data[r][c+1];
                        }else {
                            dem_xr = dem;
                        }
                        if(r > 0)
                        {
                            dem_yl = pcr::isMV(DEM->data[r-1][c])?dem:whdem->data[r-1][c];
                        }else {
                            dem_yl = dem;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            dem_yr = pcr::isMV(DEM->data[r+1][c])?dem:whdem->data[r+1][c];
                        }else {
                            dem_yr = dem;
                        }

                        float minx =std::min(std::min(dem,dem_xl),dem_xr);
                        float miny =std::min(std::min(dem,dem_yl),dem_yr);

                        float sx = (minx - dem)/dx;
                        float sy = (miny - dem)/dx;

                        float s = std::sqrt(sx*sx + sy*sy);

                       
                        //local diffusive wave estimate
                        
                        float wh =  WH_TOT->data[r][c];
                        float n = N->data[r][c];

                        float v =  std::sqrt(s) *std::pow(wh,2.0/3.0)/(N->data[r][c]);

                        float st = std::sqrt(sx*sx + sy*sy);
                        if(st <1e-12)
                        {
                            st = 1e-12;
                        }

                        sx =(dem_xl < dem_xr ? -1.0 : 1.0) * std::fabs(sx)/(st);
                        sy =(dem_yl < dem_yr ? -1.0 : 1.0) * std::fabs(sy)/(st);

                        float sx_q = gradx_avw->data[r][c];
                        float sy_q = grady_avw->data[r][c];
                        float st_q = std::sqrt(sx_q*sx_q + sy_q*sy_q);
                        if(st_q <1e-12)
                        {
                            sx_q = sx;
                            sy_q = sy;
                            st_q = 1e-12;
                        }else
                        {
                            
                            sx_q = (sx_q)/(st_q);
                            sy_q = (sy_q)/(st_q);
                        }
 
                        float h_us = 1.0*std::max(0.05f*in_dx,WH_TOT->data[r][c]);
                        sx = (9.81f * st * sx + h_us * sx_q)/std::max(1e-12f,9.81f*st + h_us);
                        sy = (9.81f * st * sy + h_us * sy_q)/std::max(1e-12f,9.81f * st + h_us);

                        vel_x[index] = v * sx;
                        vel_y[index] = v * sy;

                    }else
                    {
                        pcr::setMV(vel_x[index]);
                        pcr::setMV(vel_y[index]);
                    }
                }
            }

            SAFE_DELETE(whdem);
            SAFE_DELETE(whdem2);
            SAFE_DELETE(gradx_avw);
            SAFE_DELETE(grady_avw);
            //SAFE_DELETE(sq;
            SAFE_DELETE(FlowWL);
            SAFE_DELETE(FlowWR);
        }

        SAFE_DELETE(RainN);
        SAFE_DELETE(Vel);
        SAFE_DELETE(accN);
        SAFE_DELETE(acc1);
        SAFE_DELETE(accacc1);
        SAFE_DELETE(accaccqreal);
        SAFE_DELETE(accvel);
        if(n_sources > 0)
        {
            SAFE_DELETE(accs);
            SAFE_DELETE(accsources);
        }
        if(DEM2CHAN != nullptr)
        {
        SAFE_DELETE(DEM2CHAN);

        }
        SAFE_DELETE(Seeds);
        SAFE_DELETE(q);
        SAFE_DELETE(one);
        SAFE_DELETE(GradX1);
        SAFE_DELETE(GradY1);
        SAFE_DELETE(GradX2);
        SAFE_DELETE(GradY2);

        if(do_depressions || do_q_neg)
        {
            SAFE_DELETE(accqreal);
            SAFE_DELETE(DEMDM);
            SAFE_DELETE(DEMDID);
        }

        SAFE_DELETE(hfinal);


        SAFE_DELETE(DEMD);

        SAFE_DELETE(DEM2);
    }

    //now finally we fill hydraulic depressions
    //those depressions generated uniquely by the presence of water flow in the landscape.
    //we still have the original depressions in the DEM
    //recalculate based on water top surface
    //and apply water levels through filling depressions above some criteria
    //if(false)
    {

        //TODO:
        //this needs some rethinking,
        //there are edge-cases with extremely large depressions (central south africa, or for example middle east has some)
        //Here, extremely large hydrualic depressions could exist.
        //we could just get the size of the depressions and estimate filling physics based

        cTMapf* DEM3 = SpreadDepressionFind(DEM_ORG);

        //writeRaster(*DEM_ORG,"wh_demo.tif");
        //writeRaster(*DEM3,"wh_dem3.tif");

        cTMapf* DEM4 = DEM_ORG->GetCopy();
        for (int r = 0; r < DEM->Rows(); r++) {
            for (int c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(DEM->data[r][c])) {
                    DEM4->data[r][c] += WH_TOT->data[r][c];
                }
            }
        }
        //writeRaster(*DEM_ORG,"wh_demo.tif");

        //writeRaster(*DEM4,"wh_demws.tif");
        cTMapf* DEMDHYDRO = SpreadDepressionFind(DEM4);
        //writeRaster(*DEMDHYDRO,"wh_demhy.tif");

        for (int r = 0; r < DEM->Rows(); r++) {
            for (int c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(DEM->data[r][c])) {
                    float valn = DEMDHYDRO->data[r][c]-DEM4->data[r][c];

                    if(valn > 0.15 && (DEM3->data[r][c] - DEM_ORG->data[r][c]) <0.15)
                    {
                        WH_TOT->data[r][c] += DEMDHYDRO->data[r][c] - (DEM_ORG->data[r][c] + WH_TOT->data[r][c]);
                    }else if(valn > 0.15 && (DEMDHYDRO->data[r][c] - DEM3->data[r][c]) >= 0.15)
                    {
                        WH_TOT->data[r][c] += DEMDHYDRO->data[r][c] - (DEM_ORG->data[r][c] + WH_TOT->data[r][c]);
                    }
                }else
                {
                }
            }
        }


        delete DEM3;
        delete DEMDHYDRO;
    }

    //writeRaster(*WH_TOT,"wh_tot_postfill.tif");
    if(do_qdetail_out > 0 && piters > 1)
    {
        for(int r = 0; r < DEM->Rows(); r++) {
				for(int c = 0; c < DEM->Cols(); c++)
				{
					int index = r* DEM->Cols() + c;
					if(!pcr::isMV(DEM->data[r][c]))
					{
                        double w = std::max(1e-12f,q_detail_w->data[r][c]);
                        ret_b[index] = ret_b[index]/w;
                        ret_l[index] = ret_l[index]/w;
                        ret_avg_vel[index] = ret_avg_vel[index]/w;
                        ret_ss_comp[index] = ret_ss_comp[index]/w;
                    }
				}
        }
        

        delete q_detail_w;
    }

    if (use_local_ss && piters > 1)
    {
        for (int r = 0; r < DEM->Rows(); r++) {
            for (int c = 0; c < DEM->Cols(); c++)
            {
                int index = r * DEM->Cols() + c;
                if (!pcr::isMV(DEM->data[r][c]))
                {
                    double w = std::max(1e-12f, local_ss_weight[index]);
                    local_ss_comp[index] = local_ss_comp[index] / w;
                }
            }
        }
    }

    FF::cout << "do output discharge? " << do_output_discharge << std::endl;
	if(do_output_discharge)
	{
		for(int r = 0; r < DEM->Rows(); r++) {
				for(int c = 0; c < DEM->Cols(); c++)
				{
						int index = r* DEM->Cols() + c;
					if(!pcr::isMV(DEM->data[r][c]))
					{
						
						discharge_out[index] = QSS_PREV->data[r][c];
					    if (QBASE) {

					        discharge_out[index] = QSS_PREV->data[r][c] + QBASE->data[r][c];
					    }

					}else{
							
						pcr::setMV(discharge_out[index]);
					}
				}
		}
		
		
	}

    

    if (culvert_qmaxmodel != nullptr && n_culverts > 0)
    {
        for (int i = 0; i < n_culverts; i++)
        {
            float qmax = 0.0f;
            if (DEM != nullptr && culverts_source_row != nullptr && culverts_source_col != nullptr)
            {
                int r_source = (int)(culverts_source_row[i]);
                int c_source = (int)(culverts_source_col[i]);
                if (r_source >= 0 && r_source < DEM->Rows() && c_source >= 0 && c_source < DEM->Cols())
                {
                    float val = QSS_PREV->data[r_source][c_source];
                    if (QBASE)
                    {
                        val += QBASE->data[r_source][c_source];
                    }
                    if (!pcr::isMV(val))
                    {
                        qmax = val;
                    }
                }
            }
            culvert_qmaxmodel[i] = qmax;
        }
    }

    delete QSS_PREV;
    delete WH_PREV;
	

    //writeRaster(*ChanStorage,"chanstore.tif");
	if(do_chan)
	{
        if(export_chanwh)
        {
            FF::cout << "output wh with channel" <<std::endl;
            for(int r = 0; r < DEM->Rows(); r++) {
                    for(int c = 0; c < DEM->Cols(); c++)
                    {
                            int index = r* DEM->Cols() + c;
                        if(!pcr::isMV(DEM->data[r][c]))
                        {
//(ChanStorage->data[r][c] /(in_dx*in_dx)) +
                            ret[index] = std::max(0.0f, WH_TOT->data[r][c]);
                        }else{

                            pcr::setMV(ret[index]);
                        }
                    }
            }
        }else
        {
            for(int r = 0; r < DEM->Rows(); r++) {
                    for(int c = 0; c < DEM->Cols(); c++)
                    {
                            int index = r* DEM->Cols() + c;
                        if(!pcr::isMV(DEM->data[r][c]))
                        {

                            ret[index] = std::max(0.0f,WH_TOT->data[r][c]);
                        }else{

                            pcr::setMV(ret[index]);
                        }
                    }
            }
        }

	}else{
		for(int r = 0; r < DEM->Rows(); r++) {
				for(int c = 0; c < DEM->Cols(); c++)
				{
						int index = r* DEM->Cols() + c;
					if(!pcr::isMV(DEM->data[r][c]))
					{
						
                        ret[index] =  std::max(0.0f,WH_TOT->data[r][c]);
					}else{
							
						pcr::setMV(ret[index]);
					}
				}
		}
		
		
	}

    if(qbase_set == true)
    {
        delete QBASES;

        if(ret_qbase != 0)
        {
            for(int r = 0; r < DEM->Rows(); r++) {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                    int index = r* DEM->Cols() + c;
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        
                        ret_qbase[index] = QBASE->data[r][c];
                    }
                }
            }
            
        }

        delete QBASE;
    }


    if(do_q_neg || n_culverts > 0 || do_stormdrain)
    {
        delete ReservoirID;

        //check final reservoir fills and return those
        for(int i = 0; i < shapes; i++)
        {
            reservoircvols.at(i) = std::max(0.0f,reservoircvols.at(i)-reservoirpassvols.at(i));
            if (reservoir_final_vols != nullptr)
            {
                reservoir_final_vols[i] = reservoircvols.at(i);
            }
            if (reservoir_qtot != nullptr)
            {
                reservoir_qtot[i] = reservoir_inflow.at(i);
            }
        }

    }

    if(do_ocean > 0)
    {
        for(int r = 0; r < DEM->Rows(); r++) {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                    int index = r* DEM->Cols() + c;
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        
                        ret[index] += WHO->data[r][c];
                    }
                }
        }

        delete WHO;
    }
    if (culvert_qtot != nullptr && n_culverts > 0)
    {
        float * ss_comp_array = nullptr;
        if (do_qdetail_out > 0 && ret_ss_comp != nullptr)
        {
            ss_comp_array = ret_ss_comp;
        }
        else if (use_local_ss)
        {
            ss_comp_array = local_ss_comp;
        }

        for (int i = 0; i < n_culverts; i++)
        {
            float qtot = CULVERT_QTOT[i];
            float ss_val = 1.0f;
            if (ss_comp_array != nullptr && DEM != nullptr && culverts_source_row != nullptr && culverts_source_col != nullptr)
            {
                int r_source = (int)(culverts_source_row[i]);
                int c_source = (int)(culverts_source_col[i]);
                if (r_source >= 0 && r_source < DEM->Rows() && c_source >= 0 && c_source < DEM->Cols())
                {
                    int index = r_source * DEM->Cols() + c_source;
                    ss_val = ss_comp_array[index];
                    qtot *= ss_val;
                }
            }
            culvert_qtot[i] = std::min(culvert_q[i],qtot);
        }
    }

    if (local_ss_comp != nullptr)
    {
        delete[] local_ss_comp;
    }
    if (local_ss_weight != nullptr)
    {
        delete[] local_ss_weight;
    }

    if(n_culverts > 0 || do_stormdrain)
    {
        delete CULVERTID;
    }
		

        


        /*for(int r = 0; r < DEM->Rows(); r++) {
				for(int c = 0; c < DEM->Cols(); c++)
				{
						int index = r* DEM->Cols() + c;
					if(!pcr::isMV(DEM->data[r][c]))
					{
						
						ret[index] = SDCTOT->data[r][c];
					}else{
							
						pcr::setMV(ret[index]);
					}
				}
		}*/


	if(do_chan)
	{
		delete ChanStorage;
	}
    if(SDCTOT != nullptr)
    {
        delete SDCTOT;
    }
    delete DEM;
    delete SLOPE;
	
	delete WH_TOT;

    DEM_ORG = nullptr;
    DEM = nullptr;
	

    FF::cout << "End of fastflood function "<< std::endl;


}
}

#endif
