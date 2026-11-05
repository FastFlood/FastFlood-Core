#include <stdio.h>
#include "map/map.h"
#include "threading/parallel.h"
#include "emscripten.h"


#include "defines.h"
#include "hillshade.h"
#include "filter.h"
#include "ldd.h"
#include "derivatives.h"
#include "mapstats.h"

#include "spread.h"
#include "spreadflow.h"

#include "flowwidth.h"
#include "dinfinite.h"
#include "maxwave.h"
#include "rasterize.h"
#include "soilstats.h"
#include "demproc.h"
#include "demchannel.h"
#include "oceanfill.h"
#include "infiltration.h"
#include "upscale.h"
#include "soildepth.h"
#include "stability.h"
#include "runout2.h"
#include "pedotransfer.h"

int main() {
  printf("Loaded wasm module for fastrocks algorithm\n");
  return 0;
}


extern "C"
{
	

void EMSCRIPTEN_KEEPALIVE do_algo_soildepth(
    float*ret,
    float*data_dem,
    float*data_prod,
    float*data_thick,
    float*data_mov,
    float scale,
    float in_dx,
    float delta,
    int iters,
    int w, int h)
{

  cTMapf * DEM = new cTMapf(data_dem,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * ONE= DEM->GetCopy1();
  cTMapf * PROD= new cTMapf(data_prod,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * THICK= new cTMapf(data_thick,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * MOV = new cTMapf(data_mov,w,h,0.0,0.0, in_dx,- in_dx);
  std::vector<cTMapf *> dinf =  GetDInfiniteNetwork(DEM,delta,iters);

  cTMapf * ACC = AccufluxDInfinite(dinf,ONE,iters);

    //one version with and one without its own local flow accumulatio
    //one works better for small accumulated areas, and one works better for large accumulated areas
  cTMapf * SD1 = SoilDepthSolve(dinf,DEM, PROD,THICK,MOV, nullptr, 15,true);
  cTMapf * SDNACC = SoilDepthSolve(dinf,DEM, PROD,THICK,MOV, nullptr, 15,false);

    //Try to estimate the width of the valley locally for the soil depth to be spread around on
  cTMapf * FlowWidth = GetFlowWidthFromRiverVolume(DEM,SD1,5);

  cTMapf * RATIO = FlowWidth->GetCopy0();
  cTMapf * SDTEMP = FlowWidth->GetCopy0();
  cTMapf * DEMTEMP = FlowWidth->GetCopy0();
  #pragma omp parallel for collapse(2)
  for (int r = 0; r < DEM->Rows(); r++)
  {
      for(int c = 0; c < DEM->Cols(); c++)
      {
        if(!pcr::isMV(DEM->data[r][c]))
        {
          if(ACC->data[r][c] < 100.0)
          {
            FlowWidth->data[r][c] = std::min(in_dx,FlowWidth->data[r][c]);
          }

         RATIO->data[r][c] = std::max(0.0f,FlowWidth->data[r][c] - in_dx)/std::max(0.001f,FlowWidth->data[r][c]);
         SDTEMP->data[r][c] = RATIO->data[r][c] * SD1->data[r][c];
         DEMTEMP->data[r][c] = DEM->data[r][c] + (1.0 - RATIO->data[r][c]) * SD1->data[r][c];
        }
      }
  }

  cTMapf * SDTEMP2 = DiffusiveMaxWave(DEMTEMP, SDTEMP,150,0.1);

  #pragma omp parallel for collapse(2)
  for (int r = 0; r < DEM->Rows(); r++)
  {
      for(int c = 0; c < DEM->Cols(); c++)
      {
        if(!pcr::isMV(DEM->data[r][c]))
        {
         
         SDTEMP2->data[r][c] = std::max(SDNACC->data[r][c],(1.0f - RATIO->data[r][c]) * SD1->data[r][c] + SDTEMP2->data[r][c]);

            int index = r* DEM->Cols() + c;
            ret[index] = scale * SDTEMP2->data[r][c] ;
                    
        }
      }
  }
   
   delete SDTEMP2;
   delete SDTEMP;
   delete SD1;
   delete ONE;
   delete ACC;
   delete FlowWidth;
   delete RATIO;
   delete SDNACC;
   delete dinf[0];
   delete dinf[1];
   delete dinf[2];
   delete dinf[3];

}


void EMSCRIPTEN_KEEPALIVE do_algo_stability(
    float*ret, //stabililty
    float*ret2, //failure depth
    float*data_dem,
    float*data_coh,
    float*data_ifa,
    float*data_dens,
    float*data_sd,
    float*data_wc,
    float*data_seismic,
    float scale,
    int do_seismic,
    int do_failure,
    int iters,
    float in_dx,
    int w, int h,
    float coh_mult,
    float ifa_mult)
{

FF::cout << "do stability "<< iters << " " << in_dx << " " << w << " " << h <<  std::endl;
  cTMapf * DEM = new cTMapf(data_dem,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * COH = new cTMapf(data_coh,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * IFA = new cTMapf(data_ifa,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * DENS = new cTMapf(data_dens,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * SD = new cTMapf(data_sd,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * THETA = new cTMapf(data_wc,w,h,0.0,0.0, in_dx,- in_dx);

  cTMapf * PGA = nullptr;
  
  std::vector<cTMapf *> res;

  if( data_seismic != nullptr)
  {
      PGA = new cTMapf(data_seismic,w,h,0.0,0.0, in_dx,- in_dx);

      //do stability with seismic

      res = SlopeStabilityISS(DEM,SD,COH,IFA,DENS,THETA,PGA, 1.0,do_failure >= 1, coh_mult, ifa_mult);
  }else
  {


      res = SlopeStabilityIFM(DEM,SD,COH,IFA,DENS,THETA,1.0,iters,do_failure >= 1, coh_mult, ifa_mult);

  }

  if(res.size() == 1)
  {

    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++)
    {
        for(int c = 0; c < DEM->Cols(); c++)
        {
          if(!pcr::isMV(DEM->data[r][c]))
          {
             int index = r* DEM->Cols() + c;
              ret[index] = res[0]->data[r][c] ;
            
          }else
          {
            
               int index = r* DEM->Cols() + c;
             pcr::setMV(ret[index]);
          }
        }
    }

    

  }else if(res.size() == 2)
  {
  #pragma omp parallel for collapse(2)
      for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
               int index = r* DEM->Cols() + c;
                  ret[index] = res[0]->data[r][c] ;
                  ret2[index] = res[1]->data[r][c] ;
              
            }else
            {
              
               int index = r* DEM->Cols() + c;
                pcr::setMV(ret[index]);
                pcr::setMV(ret2[index]);
            }
          }
      }


    }


    


}


void EMSCRIPTEN_KEEPALIVE do_algo_runout(
    float*ret, //maximum flow height
    float*ret2, //maximum flow velocity
    float*ret3, //deposition
    float*ret4, //erosion

    float*data_dem,

    float*data_sourcef,
    float*data_sources,
    float*data_dens,
    float*data_ifa,
    float*data_size,

    float*data_flood,
    float*data_n,

    float scale,
    float delta,
    float iters,
    float in_dx,

    int do_depr,
    int do_erosion,

    float erosion_mult,
    float *data_coh,
    float *data_sd,
    float *data_ifab,
    
    float *ret_avg_flow_length,

    int w, 
    int h)
{


    bool do_depressions = do_depr >= 1?true:false;

  FF::cout << "do runout " << in_dx << " " << w << " " << h << " " << iters << std::endl;
  FF::cout << "do erosion " << do_erosion << " " << erosion_mult << std::endl;
  cTMapf * DEM = new cTMapf(data_dem,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * F = new cTMapf(data_sourcef,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * S = new cTMapf(data_sources,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * HFLOOD = nullptr;
  if(data_flood != 0)
  {
      HFLOOD = new cTMapf(data_flood,w,h,0.0,0.0, in_dx,- in_dx);
  }

  cTMapf * DEPOSITION= nullptr;
  cTMapf * DEPOSITION2= nullptr;
  cTMapf * DEPOSITIONNS = nullptr;
  if(ret3 != 0)
  {
    DEPOSITIONNS = DEM->GetCopy0();
    DEPOSITION = DEM->GetCopy0();
    DEPOSITION2 = DEM->GetCopy0();
  }
  cTMapf * ROCKSIZE = new cTMapf(data_size,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * IFA = new cTMapf(data_ifa,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * DENS = new cTMapf(data_dens,w,h,0.0,0.0, in_dx,- in_dx);
  cTMapf * N = new cTMapf(data_n,w,h,0.0,0.0, in_dx,- in_dx);

  cTMapf * COH = nullptr;
  cTMapf * IFAB  = nullptr;
  cTMapf * SD = nullptr;
  cTMapf * ERODE = nullptr;
  if(do_erosion)
  {
    ERODE = DEM->GetCopy0();
    COH = new cTMapf(data_coh,w,h,0.0,0.0, in_dx,- in_dx);
    IFAB = new cTMapf(data_ifab,w,h,0.0,0.0, in_dx,- in_dx);
    SD = new cTMapf(data_sd,w,h,0.0,0.0, in_dx,- in_dx);
  }

FF::cout << "test1 "<< std::endl;

  cTMapf * SourceONE = DEM->GetCopy1();
  cTMapf * SourceH = DEM->GetCopy1();
 
  cTMapf * DEMTEMP = DEM->GetCopy0();
  if(data_flood == 0)
  {
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++)
        {
            for(int c = 0; c < DEM->Cols(); c++)
            {
              if(!pcr::isMV(DEM->data[r][c]))
              {
                DEMTEMP->data[r][c]=DEM->data[r][c];
              }
            }
        }
  }else
  {
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++)
        {
            for(int c = 0; c < DEM->Cols(); c++)
            {
              if(!pcr::isMV(DEM->data[r][c]))
              {
                DEMTEMP->data[r][c]=DEM->data[r][c] +HFLOOD->data[r][c];
              }
            }
        }
  }

  std::vector<std::thread> workers1;


     std::vector<cTMapf *> dinf;
     
      workers1.push_back(std::thread([in_dx,&dinf,delta,DEM, iters]()
      {
        dinf =  GetDInfiniteNetwork(DEM,delta,iters);
      }));


    FF::cout << "test2 "<< std::endl;

    float total_failure_pixels = 0.0;
    float total_failure_vol = 0.0;
    float total_acc_travel_vol = 0.0;

      workers1.push_back(std::thread([&total_failure_pixels, &total_failure_vol, &total_acc_travel_vol, in_dx,&SourceH, &SourceONE, F, S,DEM, iters]()
      {
        float dx = in_dx;
        //#pragma omp parallel for collapse(2)
              for (int r = 0; r < DEM->Rows(); r++)
              {
                  for(int c = 0; c < DEM->Cols(); c++)
                  {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                      SourceONE->data[r][c] = (F->data[r][c] + S->data[r][c])> 0.001? 1.0:0.0;
                      total_failure_pixels = total_failure_pixels + (F->data[r][c] + S->data[r][c])> 0.001? 1.0:0.0;
                      SourceH->data[r][c] = (F->data[r][c] + S->data[r][c]);
                      total_failure_vol += (F->data[r][c] + S->data[r][c]);
                    }
                  }
              }
      }));;


FF::cout << "test2 "<< std::endl;


  //do channel?

    cTMapf * DEMD = nullptr;
    cTMapf * DEMDM = nullptr;
    cTMapf * DEMDID = nullptr;
    //get depressions

    std::vector<int> depressionids;
    std::vector<float> depressionvols;
    std::vector<float> depressioncvols;

		//alter N based on channel presence

  //detect reservoirs?
  if(do_depressions)
  {
//printf("do clump etc.. \n");
      workers1.push_back(std::thread([in_dx,&depressioncvols,&depressionvols,&depressionids,&DEMDID,DEM,&DEMD,&DEMDM,N,iters]()
      {

        float dx = in_dx;
  
  //printf("a thread2 started \n");
          //DEMD will contain the depression depth
          //MaskedRaster<float> raster_data5();
          DEMDM = new cTMapf(nullptr,  DEM->data.nr_cols(),DEM->data.nr_rows(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data5),DEM->projection(),"");

          DEMD = SpreadDepressionFind(DEM);
  
          int depressioncells = 0;
            for (int r = 0; r < DEM->Rows(); r++) {
                for (int c = 0; c < DEM->Cols(); c++) {
                    if (!pcr::isMV(DEM->data[r][c])) {
                        float valn = DEMD->data[r][c]-DEM->data[r][c];

                        if(valn > 0.0001)
                        {
                            DEMDM->data[r][c] = DEMD->data[r][c];
                            depressioncells++;

                        }else
                        {
                            pcr::setMV(DEMDM->data[r][c]);
                        }
                        DEMD->data[r][c] = valn;
                    }else
                    {
                        pcr::setMV(DEMDM->data[r][c]);
                    }
                }
            }
  //printf("a thread2 clump \n");
            DEMDID = Clump(DEMDM,0.1,true);

  //printf("a thread2 clumpd \n");

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
                            int oldsize = depressionids.size();
                            for(int i = 0; i < (5 + class_current - depressionids.size()); i++)
                            {
                                //FF::cout << "add" << std::endl;
                                depressionids.push_back(oldsize + i);
                                depressionvols.push_back(0.0);
                            }

                        }
                        {
                            depressionvols.at(class_current) += cellarea * DEMD->data[r][c];
                        }
                    }

                }
            }

            
            for(int i =0; i < depressionvols.size(); i++)
            {
                depressioncvols.push_back(0.0);
            }

            FF::cout << "found depressions: " << depressioncvols.size()<< " " << depressioncells <<  std::endl;
            if(depressioncvols.size() == 0)
            {
              depressionids.push_back(0);
              depressionvols.push_back(0.0);
              depressioncvols.push_back(0.0);

            }

      }));
  }

  std::for_each(workers1.begin(), workers1.end(), [](std::thread &t)
       {
           t.join();
       });

  //routing

  cTMapf * ACC1 = AccufluxDInfinite(dinf,SourceONE,iters);
cTMapf * ACCH = AccufluxDInfinite(dinf,SourceH,iters);
cTMapf * ACCACC1 = AccufluxDInfinite(dinf,ACC1,iters);
cTMapf * ACCHAVG = DEM->GetCopy0();
cTMapf * ACCDIST = DEM->GetCopy0();
#pragma omp parallel for collapse(2)
      for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              ACCHAVG->data[r][c] = ACCH->data[r][c]/std::max(ACC1->data[r][c],0.001f);
              ACCDIST->data[r][c] = in_dx * ACCACC1->data[r][c]/std::max(ACC1->data[r][c],0.001f);
              
            }
          }
      }

  //statistics
FF::cout << "test3 "<< std::endl;
   FF::cout << "test dx " << DEM->cellSizeX() << std::endl;

  cTMapf * ACCDISTVAR = DInfiniteVariance(dinf,ACCDIST,SourceONE,ACC1,iters);

  FF::cout << "test31 "<< std::endl;

  #pragma omp parallel for collapse(2)
      for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              ACCDISTVAR->data[r][c] = std::sqrt(12.0 * ACCDISTVAR->data[r][c]);
            }
          }
      }

FF::cout << "test32 "<< std::endl;

      cTMapf * LSH = DEM->GetCopy0();
#pragma omp parallel for collapse(2)
      for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              LSH->data[r][c] = ACCH->data[r][c]/ std::max(in_dx,ACCDISTVAR->data[r][c]);
            }
          }
      }
FF::cout << "test4 "<< std::endl;


  //properties

  cTMapf * FW =GetFlowWidthFromRiverVolume(DEM,LSH,5);

  //momentum and mass model solver

  std::vector<cTMapf *> results;
  
  if(!do_erosion)
  {
    results = DinfiniteLandslideMomentumDepositionSolve2(dinf,DEM,N,S,F,DENS,ACCH,FW,IFA,ROCKSIZE, HFLOOD, ACCDISTVAR,DEMDID,depressionvols,depressioncvols,iters);
  }else
  {
     results = DinfiniteLandslideMomentumDepositionErosionSolve2(dinf,DEM,N,S,F,DENS,ACCH,FW,IFA,ROCKSIZE, HFLOOD, ACCDISTVAR,DEMDID,depressionvols,depressioncvols,COH,IFAB,SD,ERODE,erosion_mult,iters);
  }
    //{materialf, materials, velavg,velacc,deposition,momentum};

  //get average travel length
  //get number of pixels of failure, and total volume of failure
  //calculate average failure volume
  //calculate the total amount of volume traveling accumulatively 
  for (int r = 0; r < DEM->Rows(); r++)
    {
        for(int c = 0; c < DEM->Cols(); c++)
        {
          if(!pcr::isMV(DEM->data[r][c]))
          {
            total_acc_travel_vol +=(results[0]->data[r][c] + results[1]->data[r][c])/(in_dx*in_dx);
          }
        }
    }

  //total accumulative travel volume/failure volume = average length of travel
  FF::cout << "average travel length:  " << total_acc_travel_vol/std::max(1e-10f,total_failure_vol) << std::endl;


  if(ret_avg_flow_length != nullptr)
  {
    ret_avg_flow_length[0] = in_dx * total_acc_travel_vol/std::max(1e-10f,total_failure_vol);
    ret_avg_flow_length[1] = in_dx * in_dx * total_failure_vol;
  }
      
  //flowwidth

  cTMapf * FHM = DEM->GetCopy0();
  #pragma omp parallel for collapse(2)
      for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              FHM->data[r][c] = ((results[0]->data[r][c] + results[1]->data[r][c])/(in_dx * std::max(in_dx,ACCDISTVAR->data[r][c])));
            }
          }
      }
FF::cout << "test5 "<< std::endl;

  cTMapf * FW2 = GetFlowWidthFromRiverVolume(DEMTEMP,FHM,5);

  std::vector<cTMapf *> results_widthdir = GetFlowWidthDirectionFromRiverVolume(DEMTEMP,FHM,5);
  //ratio
  cTMapf * RATIO = DEM->GetCopy0();
  cTMapf * HTEMP = DEM->GetCopy0();

  float scale_subcell = 1.0;
  
  scale_subcell = std::min(1.0f,(float)(std::max(0.1,std::exp(-0.075*in_dx)))); //plot exp(-0.075*x) * (1.0-exp(-0.25x)) + exp(-0.25x)


  #pragma omp parallel for collapse(2)
for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              results_widthdir[1]->data[r][c] = -results_widthdir[1]->data[r][c];
               results_widthdir[2]->data[r][c] = -results_widthdir[2]->data[r][c];
              
              float scale_subcell_here = scale_subcell;
              float frac_ratio_subcell = (std::min(0.25f * in_dx,FHM->data[r][c])/FHM->data[r][c]);
              scale_subcell_here = std::max(0.0f,std::min(1.0f,frac_ratio_subcell * scale_subcell_here + (1.0f-frac_ratio_subcell)));
              RATIO->data[r][c]=std::min(1.0f,std::max(0.0f,scale_subcell_here*results_widthdir[0]->data[r][c] - in_dx)/std::max(0.0001f,results_widthdir[0]->data[r][c]));
              DEMTEMP->data[r][c]=DEMTEMP->data[r][c] + (1.0-RATIO->data[r][c]) * FHM->data[r][c];
              HTEMP->data[r][c]= RATIO->data[r][c] * FHM->data[r][c];

              if(DEPOSITION != nullptr)
              {
                  DEPOSITIONNS->data[r][c] = 0.0;//(1.0-RATIO->data[r][c]) * results[4]->data[r][c]/(in_dx*in_dx);
                  DEPOSITION->data[r][c] =  results[4]->data[r][c]/(in_dx*in_dx);//(RATIO->data[r][c])
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
                          //to what extent is this depression filled?
                          float frac = std::min(1.0f,std::max(0.0f,depressioncvols.at(id)/ std::max(0.01f,depressionvols.at(id))));
                          HTEMP->data[r][c] += DEMD->data[r][c] * frac;
                         if(DEPOSITION != nullptr)
                          {
                              DEPOSITION->data[r][c] += DEMD->data[r][c] * frac;
                          }
                      }
                      
                  }
              }

            }
          }
      }
  //maxwave
FF::cout << "test6 "<< std::endl;

//we should also spread the deposition, as the mass travels along with the full width of the mass movement, and 
//not just along the primary channels
cTMapf * FHCORR = nullptr;
cTMapf * FHCORR2 = nullptr;
if(DEPOSITION == nullptr)
{
  FHCORR = DiffusiveMaxWavePrefDir(DEMTEMP,HTEMP, results_widthdir[1], results_widthdir[2],75,0.1);

  FF::cout << "test7 "<< std::endl;
  FHCORR2 = DiffusiveMaxWave(DEMTEMP,FHCORR,50,0.1);
}else
{
  //FHCORR = HTEMP->GetCopy();
  //FHCORR2 = HTEMP->GetCopy();
  FHCORR = DiffusiveMaxWavePrefDirDeposition(DEMTEMP,HTEMP, results_widthdir[1], results_widthdir[2],75,0.1,DEPOSITION,DEPOSITION2);

  FF::cout << "test7 with deposition "<< std::endl;
  FHCORR2 = DiffusiveMaxWaveDeposition(DEMTEMP,FHCORR,50,0.1,DEPOSITION2,DEPOSITION);


}
 

  //final solution

for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              int index = r* DEM->Cols() + c;

                        
                ret[index] = (1.0-RATIO->data[r][c]) * FHM->data[r][c] + FHCORR2->data[r][c];
                ret2[index] = results[2]->data[r][c];
            }
          }
      }

      if(DEPOSITION != nullptr)
      {
          for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              if(!pcr::isMV(DEM->data[r][c]))
              {
                int index = r* DEM->Cols() + c;

        
                ret3[index] = DEPOSITION->data[r][c] + DEPOSITIONNS->data[r][c];
              }
            }
          }
      }
    }

    if(do_erosion)
    {
       for (int r = 0; r < DEM->Rows(); r++)
      {
          for(int c = 0; c < DEM->Cols(); c++)
          {
            if(!pcr::isMV(DEM->data[r][c]))
            {
              if(!pcr::isMV(DEM->data[r][c]))
              {
                int index = r* DEM->Cols() + c;

        
                ret4[index] = ERODE->data[r][c]/(in_dx * std::max(in_dx,results_widthdir[0]->data[r][c]));
              }
            }
          }
      }
    }

FF::cout << "test8 "<< std::endl;

    if(do_depressions)
    {
      delete DEMD;
      delete DEMDM;
        delete DEMDID;
    }

    if(DEPOSITION != nullptr)
    {
      delete DEPOSITION;
      delete DEPOSITIONNS;
      delete DEPOSITION2;
    }

     if(do_erosion)
    {
      delete COH;
      delete IFAB;
      delete SD;
      delete ERODE;
    }

    if(FHCORR != nullptr)
    {
      delete FHCORR;
    }
    if(FHCORR2 != nullptr)
    {
      delete FHCORR2;
    }
      delete RATIO;
      delete DEMTEMP;
      delete HTEMP;
      delete FHM;
      delete FW2;
      delete FW;
      delete ACCDIST;
      delete ACCDISTVAR;
      delete SourceONE;
      delete SourceH;
      delete ACC1;
      delete ACCH;
      delete ACCACC1;
      delete ACCHAVG;
      delete dinf[0];
      delete dinf[1];
      delete dinf[2];
      delete dinf[3];
      delete results_widthdir[0];
      delete results_widthdir[1];
      delete results_widthdir[2];
      delete results[0];
      delete results[1];
      delete results[2];
      delete results[3];
      delete results[4];
      delete results[5];


}



}



