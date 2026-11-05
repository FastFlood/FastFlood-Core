#ifndef FF_FULLFLOODT_H
#define FF_FULLFLOODT_H

#include "../geospatial/map/map.h"
#include "../common/object.h"
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
#include "../algorithms/rasteralgorithms/lddresample.h"
#include "../fullflood/fullflood_gridresample.h"

extern "C"
{


class FF_FULLFLOODT : public Object
{

public:


  double qout = 0.0;

  //channel info
  bool do_channels = false;
  cTMapf * CHANNELLDD = nullptr;
  cTMapf * CHANNELW = nullptr;
  cTMapf * CHANNELD = nullptr;
  cTMapf * SLOPE = nullptr;
  cTMapf * DEM_ORGCHAN = nullptr;
  cTMapf * DEM2CHAN = nullptr;

  //Channel is run at R3 level, and interactions with 2D solver are as well
  cTMapf * R3CHANNELLDD = nullptr;
  cTMapf * R3CHANNELW = nullptr;
  cTMapf * R3CHANNELD = nullptr;
  cTMapf * R3SLOPE = nullptr;

  cTMapf * R3CHANNELALPHA = nullptr;
  cTMapf * R3CHANNELQN = nullptr;
  cTMapf * R3CHANNELINFIL = nullptr;

  cTMapf * R3CHANNELH = nullptr;
  cTMapf * R3CHANNELHN = nullptr;

  cTMapf * GradX1 = nullptr;
  cTMapf * GradX2 = nullptr;
  cTMapf * GradY1 = nullptr;
  cTMapf * GradY2 = nullptr;
  //culverts
  bool do_culverts = false;
  cTMapf * CULVERTID = nullptr;
  std::vector<float> CULVERT_QTOT;
  std::vector<float> CULVERT_Q;
  std::vector<float> CULVERT_FRAC;
  std::vector<int> CULVERT_RS;
  std::vector<int> CULVERT_CS;
  std::vector<int> CULVERT_RT;
  std::vector<int> CULVERT_CT;


  cTMapf * FlowTSUPERMINNEXT = nullptr;

  //ocean
  bool do_ocean = false;
  float ocean_level = 0.0;

  //baseflow
  bool do_baseflow = false;
  float baseflow_vol_in_m3sec = false;

  //reservoirs
  bool do_reservoirs = false;
  cTMapf * ReservoirID = nullptr;
  std::vector<int> reservoirids;
  std::vector<float> reservoirvols;
  std::vector<float> reservoircvols;
  std::vector<int> reservoirlink;
  std::vector<float> reservoirpassvols;


  //stormdrains
  bool do_stormdrains = false;
  cTMapf * SDC = nullptr;
  cTMapf * SDCTOT = nullptr;

  cTMapf * WHMax = nullptr;
  cTMapf * VMax = nullptr;
  cTMapf * ArrivalTime = nullptr;
  cTMapf * QCHMax = nullptr;
  std::vector<cTMapf *> WHMax_History = {};
  std::vector<cTMapf *> VxMax_History = {};
  std::vector<cTMapf *> VyMax_History = {};

  cTMapf * QDRAIN = nullptr;
  cTMapf * KSAT = nullptr;
  cTMapf * INFIL = nullptr;
  cTMapf * FlowT = nullptr;
  cTMapf * FlowTN = nullptr;
  cTMapf * FlowDT = nullptr;
  cTMapf * FlowDTL = nullptr;
  cTMapf * FlowH = nullptr;
  cTMapf * FlowUX = nullptr;
  cTMapf * FlowUY = nullptr;
  cTMapf * FlowHN = nullptr;
  cTMapf * FlowUXN = nullptr;
  cTMapf * FlowUYN = nullptr;
  cTMapf * DEM = nullptr;
  cTMapf * SGA = nullptr;
  cTMapf * MAN = nullptr;
  cTMapf * RAIN = nullptr;
  cTMapf * CRIT = nullptr;

  cTMapf * R3FlowT = nullptr;
  cTMapf * R3FlowTN = nullptr;
  cTMapf * R3FlowDT = nullptr;
  cTMapf * R3FlowDTL = nullptr;
  cTMapf * R3FlowH = nullptr;
  cTMapf * R3FlowUX = nullptr;
  cTMapf * R3FlowUY = nullptr;
  cTMapf * R3FlowHN = nullptr;
  cTMapf * R3FlowUXN = nullptr;
  cTMapf * R3FlowUYN = nullptr;
  cTMapf * R3DEM = nullptr; //average dem in subgrid
  cTMapf * R3DEMO = nullptr; //average dem in subgrid
  cTMapf * R3SGA = nullptr; //is sub grid active?
  cTMapf * R3DEMM = nullptr; //maximum dem in subgrid
  cTMapf * R3MAN = nullptr; //average manning
  cTMapf * R3RAIN = nullptr; //average manning
  cTMapf * R3KSAT = nullptr; //average manning
  cTMapf * R3CRIT = nullptr;

  cTMapf * R9FlowT = nullptr;
  cTMapf * R9FlowTN = nullptr;
  cTMapf * R9FlowDT = nullptr;
  cTMapf * R9FlowDTL = nullptr;
  cTMapf * R9FlowH = nullptr;
  cTMapf * R9FlowUX = nullptr;
  cTMapf * R9FlowUY = nullptr;
  cTMapf * R9FlowHN = nullptr;
  cTMapf * R9FlowUXN = nullptr;
  cTMapf * R9FlowUYN = nullptr;
  cTMapf * R9DEM = nullptr; //average dem in subgrid
  cTMapf * R9DEMO = nullptr; //average dem in subgrid
  cTMapf * R9SGA = nullptr; //is sub grid active?

  cTMapf * R9DEMM = nullptr; //maximum dem in subgrid
  cTMapf * R9MAN = nullptr; //average manning
  cTMapf * R9RAIN = nullptr; //average manning
  cTMapf * R9KSAT = nullptr; //average manning
  cTMapf * R9CRIT = nullptr;

  cTMapf * R27FlowT = nullptr;
  cTMapf * R27FlowTN = nullptr;
  cTMapf * R27FlowDT = nullptr;
  cTMapf * R27FlowDTL = nullptr;
  cTMapf * R27FlowH = nullptr;
  cTMapf * R27FlowUX = nullptr;
  cTMapf * R27FlowUY = nullptr;
  cTMapf * R27FlowHN = nullptr;
  cTMapf * R27FlowUXN = nullptr;
  cTMapf * R27FlowUYN = nullptr;
  cTMapf * R27DEM = nullptr; //average dem in subgrid
  cTMapf * R27DEMO = nullptr; //average dem in subgrid
  cTMapf * R27SGA = nullptr; //is sub grid active?

  cTMapf * R27DHDTN = nullptr;
  cTMapf * R27DUXDTN = nullptr;
  cTMapf * R27DUYDTN = nullptr;
  cTMapf * R9DHDTN = nullptr;
  cTMapf * R9DUXDTN = nullptr;
  cTMapf * R9DUYDTN = nullptr;
  cTMapf * R3DHDTN = nullptr;
  cTMapf * R3DUXDTN = nullptr;
  cTMapf * R3DUYDTN = nullptr;
  cTMapf * DHDTN = nullptr;
  cTMapf * DUXDTN = nullptr;
  cTMapf * DUYDTN = nullptr;


  cTMapf * R27DHDT = nullptr;
  cTMapf * R27DUXDT = nullptr;
  cTMapf * R27DUYDT = nullptr;
  cTMapf * R9DHDT = nullptr;
  cTMapf * R9DUXDT = nullptr;
  cTMapf * R9DUYDT = nullptr;
  cTMapf * R3DHDT = nullptr;
  cTMapf * R3DUXDT = nullptr;
  cTMapf * R3DUYDT = nullptr;
  cTMapf * DHDT = nullptr;
  cTMapf * DUXDT = nullptr;
  cTMapf * DUYDT = nullptr;

  cTMapf * HLLX1X = nullptr;
  cTMapf * HLLX1Y = nullptr;
  cTMapf * HLLX1Z = nullptr;
  cTMapf * HLLX2X = nullptr;
  cTMapf * HLLX2Y = nullptr;
  cTMapf * HLLX2Z = nullptr;
  cTMapf * HLLY1X = nullptr;
  cTMapf * HLLY1Y = nullptr;
  cTMapf * HLLY1Z = nullptr;
  cTMapf * HLLY2X = nullptr;
  cTMapf * HLLY2Y = nullptr;
  cTMapf * HLLY2Z = nullptr;

  cTMapf * SAVX1X = nullptr;
  cTMapf * SAVX1Y = nullptr;
  cTMapf * SAVX2X = nullptr;
  cTMapf * SAVX2Y = nullptr;
  cTMapf * SAVY1X = nullptr;
  cTMapf * SAVY1Y = nullptr;
  cTMapf * SAVY2X = nullptr;
  cTMapf * SAVY2Y = nullptr;


  cTMapf * R3HLLX1X = nullptr;
  cTMapf * R3HLLX1Y = nullptr;
  cTMapf * R3HLLX1Z = nullptr;
  cTMapf * R3HLLX2X = nullptr;
  cTMapf * R3HLLX2Y = nullptr;
  cTMapf * R3HLLX2Z = nullptr;
  cTMapf * R3HLLY1X = nullptr;
  cTMapf * R3HLLY1Y = nullptr;
  cTMapf * R3HLLY1Z = nullptr;
  cTMapf * R3HLLY2X = nullptr;
  cTMapf * R3HLLY2Y = nullptr;
  cTMapf * R3HLLY2Z = nullptr;

  cTMapf * R3SAVX1X = nullptr;
  cTMapf * R3SAVX1Y = nullptr;
  cTMapf * R3SAVX2X = nullptr;
  cTMapf * R3SAVX2Y = nullptr;
  cTMapf * R3SAVY1X = nullptr;
  cTMapf * R3SAVY1Y = nullptr;
  cTMapf * R3SAVY2X = nullptr;
  cTMapf * R3SAVY2Y = nullptr;


  cTMapf * R9HLLX1X = nullptr;
  cTMapf * R9HLLX1Y = nullptr;
  cTMapf * R9HLLX1Z = nullptr;
  cTMapf * R9HLLX2X = nullptr;
  cTMapf * R9HLLX2Y = nullptr;
  cTMapf * R9HLLX2Z = nullptr;
  cTMapf * R9HLLY1X = nullptr;
  cTMapf * R9HLLY1Y = nullptr;
  cTMapf * R9HLLY1Z = nullptr;
  cTMapf * R9HLLY2X = nullptr;
  cTMapf * R9HLLY2Y = nullptr;
  cTMapf * R9HLLY2Z = nullptr;

  cTMapf * R9SAVX1X = nullptr;
  cTMapf * R9SAVX1Y = nullptr;
  cTMapf * R9SAVX2X = nullptr;
  cTMapf * R9SAVX2Y = nullptr;
  cTMapf * R9SAVY1X = nullptr;
  cTMapf * R9SAVY1Y = nullptr;
  cTMapf * R9SAVY2X = nullptr;
  cTMapf * R9SAVY2Y = nullptr;

  cTMapf * R27HLLX1X = nullptr;
  cTMapf * R27HLLX1Y = nullptr;
  cTMapf * R27HLLX1Z = nullptr;
  cTMapf * R27HLLX2X = nullptr;
  cTMapf * R27HLLX2Y = nullptr;
  cTMapf * R27HLLX2Z = nullptr;
  cTMapf * R27HLLY1X = nullptr;
  cTMapf * R27HLLY1Y = nullptr;
  cTMapf * R27HLLY1Z = nullptr;
  cTMapf * R27HLLY2X = nullptr;
  cTMapf * R27HLLY2Y = nullptr;
  cTMapf * R27HLLY2Z = nullptr;

  cTMapf * R27SAVX1X = nullptr;
  cTMapf * R27SAVX1Y = nullptr;
  cTMapf * R27SAVX2X = nullptr;
  cTMapf * R27SAVX2Y = nullptr;
  cTMapf * R27SAVY1X = nullptr;
  cTMapf * R27SAVY1Y = nullptr;
  cTMapf * R27SAVY2X = nullptr;
  cTMapf * R27SAVY2Y = nullptr;

  cTMapf * R27SGHO = nullptr;
  cTMapf * R27SGUXO = nullptr;
  cTMapf * R27SGUYO = nullptr;

  cTMapf * R27SGCA = nullptr; //Super-grid cache active
  cTMapf * R27DEMM = nullptr; //maximum dem in subgrid
  cTMapf * R27MAN = nullptr; //average manning
  cTMapf * R27RAIN = nullptr; //average manning
  cTMapf * R27KSAT = nullptr; //average manning
  cTMapf * R27CRIT = nullptr;

  cTMapf * R27FlowDTLN = nullptr;
  cTMapf * R9FlowDTLN = nullptr;
  cTMapf * R3FlowDTLN = nullptr;
  cTMapf * FlowDTLN = nullptr;

  cTMapf * R27FlowDTC = nullptr;
  cTMapf * R9FlowDTC = nullptr;
  cTMapf * R3FlowDTC = nullptr;
  cTMapf * FlowDTC = nullptr;

  //usually at the highest resolution
  //we can take a larger timestep as well, due to the usage of FlowKinematic
  cTMapf * CHAN_LDD = nullptr;
  cTMapf * CHAN_H = nullptr;
  cTMapf * CHAN_W = nullptr;
  cTMapf * CHAN_D = nullptr;
  float chan_manning = 0.03;


  cTMapf * Cache_HOLD = nullptr;
  cTMapf * Cache_UXOLD = nullptr;
  cTMapf * Cache_UYOLD = nullptr;

  double rain_total = 0.0;

  bool do_cache = true;
  bool do_sdt = true;
  double cache_rel_change_crit = 0.25; //0.5% change

  int step = 0;
  float courant = 0.1;
  float dt_min = 0.01;
  float dt_user = 1.0;
  double time_sec = 0.0;
  float duration = 1.0;

  //return data
  float *ret = nullptr;
  float *retvx = nullptr;
  float *retvy = nullptr;
  float *retq = nullptr;
  float *ret_hmax = nullptr;
  float *ret_vmax = nullptr;
  float *ret_at = nullptr;
  int iters = 10;

  int setting_in_iters = 1;

  inline void Assign()
  {

  }

  inline void SetCoreOptions(bool in_do_sdt, bool in_do_cache, bool in_do_channel)
  {
    do_sdt = in_do_sdt;
    do_cache = in_do_cache;
    do_channels = in_do_channel;
  }

  inline void SetIterations(int x)
  {
     setting_in_iters = x;
  }
  inline void SetCoresToUse(int x)
  {

  }
  inline void SetCoreCriteria(double crit_refinement, double crit_refinement_chan, double crit_cutoff_channel)
  {

  }

  //this transfers ownership!!!!
  inline void SetCoreMapsAndAllocate(cTMapf * in_dem, cTMapf * in_N, cTMapf * in_KSAT, cTMapf * in_RAIN, float in_duration, cTMapf * chanw, cTMapf * chanh, cTMapf * dem_channelized, float delta, float in_dx, cTMapf* GradX1_i, cTMapf* GradX2_i, cTMapf*GradY1_i, cTMapf*GradY2_i)
  {
        DEM= in_dem;
        MAN = in_N;
        KSAT = in_KSAT;
        RAIN = in_RAIN;

        duration = in_duration;
        iters = setting_in_iters  * 2.0 * std::max(10.0,std::max(DEM->nrCols(),DEM->nrRows())/1000.0);

        SGA = DEM->GetCopy0();
        FlowH = DEM->GetCopy0();
        FlowUX = DEM->GetCopy0();
        FlowUY = DEM->GetCopy0();
        FlowHN = DEM->GetCopy0();
        FlowUXN = DEM->GetCopy0();
        FlowUYN = DEM->GetCopy0();
        CRIT = DEM->GetCopy1();

        step = 0;


        R3DEM = DEM->GetCopy0Cascade(3);
        R3DEMO = R3DEM;
        R3DEMM = R3DEM->GetCopy0();
        R3MAN = R3DEM->GetCopy0();
        R3KSAT = R3DEM->GetCopy0();
        R3RAIN = R3DEM->GetCopy0();
        R3SGA = R3DEM->GetCopy0();
        R3FlowH = R3DEM->GetCopy0();
        R3FlowUX = R3DEM->GetCopy0();
        R3FlowUY = R3DEM->GetCopy0();
        R3FlowHN = R3DEM->GetCopy0();
        R3FlowUXN = R3DEM->GetCopy0();
        R3FlowUYN = R3DEM->GetCopy0();
        R3CRIT = R3DEM->GetCopy1();

        R9DEM = DEM->GetCopy0Cascade(9);
        R9DEMM = R9DEM->GetCopy0();
        R9DEMO = R9DEM;
        R9SGA = R9DEM->GetCopy0();
        R9MAN = R9DEM->GetCopy0();
        R9KSAT = R9DEM->GetCopy0();
        R9RAIN = R9DEM->GetCopy0();
        R9FlowH = R9DEM->GetCopy0();
        R9FlowUX = R9DEM->GetCopy0();
        R9FlowUY = R9DEM->GetCopy0();
        R9FlowHN = R9DEM->GetCopy0();
        R9FlowUXN = R9DEM->GetCopy0();
        R9FlowUYN = R9DEM->GetCopy0();
        R9CRIT = R9DEM->GetCopy1();

        R27DEM = DEM->GetCopy0Cascade(27);
        R27DEMO = R27DEM;
        R27DEMM = R27DEM->GetCopy0();
        R27SGA = R27DEM->GetCopy0();
        R27MAN = R27DEM->GetCopy0();
        R27KSAT = R27DEM->GetCopy0();
        R27RAIN = R27DEM->GetCopy0();
        R27FlowH = R27DEM->GetCopy0();
        R27FlowUX = R27DEM->GetCopy0();
        R27FlowUY = R27DEM->GetCopy0();
        R27FlowHN = R27DEM->GetCopy0();
        R27FlowUXN = R27DEM->GetCopy0();
        R27FlowUYN = R27DEM->GetCopy0();
        R27CRIT = R27DEM->GetCopy1();


        FF::cout << "Resolutions - R27: " << R27DEM->nrRows() << " " << R27DEM->nrCols() <<std::endl;
        FF::cout << "Resolutions - R9: " << R9DEM->nrRows() << " " << R9DEM->nrCols() <<std::endl;
        FF::cout << "Resolutions - R3: " << R3DEM->nrRows() << " " << R3DEM->nrCols() <<std::endl;
        FF::cout << "Resolutions - R1: " << DEM->nrRows() << " " << DEM->nrCols() <<std::endl;

        R27DEM = SpreadMonotonicReconstruct(R27DEMO,0.00001,45);//SpreadDirectionalAbsMaxMD(Seeds,DEM,GradX1,GradX2,GradY1,GradY2);
        R9DEM = SpreadMonotonicReconstruct(R9DEMO,0.00001,45);//SpreadDirectionalAbsMaxMD(Seeds,DEM,GradX1,GradX2,GradY1,GradY2);
        R3DEM = SpreadMonotonicReconstruct(R3DEMO,0.00001,45);//SpreadDirectionalAbsMaxMD(Seeds,DEM,GradX1,GradX2,GradY1,GradY2);

        if(do_cache)
        {
            //for now, we only do the caching at R27
            //might be usefull to do at R9 instead?
            //Not sure both might be best but how
            R27SGHO = R27DEM->GetCopy0();
            R27SGUXO = R27DEM->GetCopy0();
            R27SGUYO = R27DEM->GetCopy0();

            R27SGCA = R27DEM->GetCopy0();

            Cache_HOLD = DEM->GetCopy0();
            Cache_UXOLD = DEM->GetCopy0();
            Cache_UYOLD = DEM->GetCopy0();
        }

        if(do_sdt)
        {
            FlowT = DEM->GetCopy0();
            FlowDT = DEM->GetCopy0();
            R3FlowT = R3DEM->GetCopy0();
            R3FlowDT = R3DEM->GetCopy0();
            R9FlowT = R9DEM->GetCopy0();
            R9FlowDT = R9DEM->GetCopy0();
            R27FlowT = R27DEM->GetCopy0();
            R27FlowDT = R27DEM->GetCopy0();
            FlowTN = DEM->GetCopy0();
            FlowDTL = DEM->GetCopy0();
            R3FlowTN = R3DEM->GetCopy0();
            R3FlowDTL = R3DEM->GetCopy0();
            R9FlowTN = R9DEM->GetCopy0();
            R9FlowDTL = R9DEM->GetCopy0();
            R27FlowTN = R27DEM->GetCopy0();
            R27FlowDTL = R27DEM->GetCopy0();

            R27FlowDTLN = R27DEM->GetCopy0();
            R9FlowDTLN = R9DEM->GetCopy0();
            R3FlowDTLN = R3DEM->GetCopy0();
            FlowDTLN = DEM->GetCopy0();

            R27FlowDTC = R27DEM->GetCopy0();
            R9FlowDTC = R9DEM->GetCopy0();
            R3FlowDTC = R3DEM->GetCopy0();
            FlowDTC = DEM->GetCopy0();

            FlowTSUPERMINNEXT = R27DEM->GetCopy0();


            R27DHDT = R27DEM->GetCopy0();
            R27DUXDT = R27DEM->GetCopy0();
            R27DUYDT = R27DEM->GetCopy0();
            R9DHDT = R9DEM->GetCopy0();
            R9DUXDT = R9DEM->GetCopy0();
            R9DUYDT = R9DEM->GetCopy0();
            R3DHDT = R3DEM->GetCopy0();
            R3DUXDT = R3DEM->GetCopy0();
            R3DUYDT = R3DEM->GetCopy0();
            DHDT = DEM->GetCopy0();
            DUXDT = DEM->GetCopy0();
            DUYDT = DEM->GetCopy0();

            R27DHDTN = R27DEM->GetCopy0();
            R27DUXDTN = R27DEM->GetCopy0();
            R27DUYDTN = R27DEM->GetCopy0();
            R9DHDTN = R9DEM->GetCopy0();
            R9DUXDTN = R9DEM->GetCopy0();
            R9DUYDTN = R9DEM->GetCopy0();
            R3DHDTN = R3DEM->GetCopy0();
            R3DUXDTN = R3DEM->GetCopy0();
            R3DUYDTN = R3DEM->GetCopy0();
            DHDTN = DEM->GetCopy0();
            DUXDTN = DEM->GetCopy0();
            DUYDTN = DEM->GetCopy0();


            HLLX1X = DEM->GetCopy0();
            HLLX1Y = DEM->GetCopy0();
            HLLX1Z = DEM->GetCopy0();
            HLLX2X = DEM->GetCopy0();
            HLLX2Y = DEM->GetCopy0();
            HLLX2Z = DEM->GetCopy0();
            HLLY1X = DEM->GetCopy0();
            HLLY1Y = DEM->GetCopy0();
            HLLY1Z = DEM->GetCopy0();
            HLLY2X = DEM->GetCopy0();
            HLLY2Y = DEM->GetCopy0();
            HLLY2Z = DEM->GetCopy0();

            R3HLLX1X = R3DEM->GetCopy0();
            R3HLLX1Y = R3DEM->GetCopy0();
            R3HLLX1Z = R3DEM->GetCopy0();
            R3HLLX2X = R3DEM->GetCopy0();
            R3HLLX2Y = R3DEM->GetCopy0();
            R3HLLX2Z = R3DEM->GetCopy0();
            R3HLLY1X = R3DEM->GetCopy0();
            R3HLLY1Y = R3DEM->GetCopy0();
            R3HLLY1Z = R3DEM->GetCopy0();
            R3HLLY2X = R3DEM->GetCopy0();
            R3HLLY2Y = R3DEM->GetCopy0();
            R3HLLY2Z = R3DEM->GetCopy0();

            R9HLLX1X = R9DEM->GetCopy0();
            R9HLLX1Y = R9DEM->GetCopy0();
            R9HLLX1Z = R9DEM->GetCopy0();
            R9HLLX2X = R9DEM->GetCopy0();
            R9HLLX2Y = R9DEM->GetCopy0();
            R9HLLX2Z = R9DEM->GetCopy0();
            R9HLLY1X = R9DEM->GetCopy0();
            R9HLLY1Y = R9DEM->GetCopy0();
            R9HLLY1Z = R9DEM->GetCopy0();
            R9HLLY2X = R9DEM->GetCopy0();
            R9HLLY2Y = R9DEM->GetCopy0();
            R9HLLY2Z = R9DEM->GetCopy0();

            R27HLLX1X = R27DEM->GetCopy0();
            R27HLLX1Y = R27DEM->GetCopy0();
            R27HLLX1Z = R27DEM->GetCopy0();
            R27HLLX2X = R27DEM->GetCopy0();
            R27HLLX2Y = R27DEM->GetCopy0();
            R27HLLX2Z = R27DEM->GetCopy0();
            R27HLLY1X = R27DEM->GetCopy0();
            R27HLLY1Y = R27DEM->GetCopy0();
            R27HLLY1Z = R27DEM->GetCopy0();
            R27HLLY2X = R27DEM->GetCopy0();
            R27HLLY2Y = R27DEM->GetCopy0();
            R27HLLY2Z = R27DEM->GetCopy0();


            R27SAVX1X = R27DEM->GetCopy0();
            R27SAVX1Y = R27DEM->GetCopy0();
            R27SAVX2X = R27DEM->GetCopy0();
            R27SAVX2Y = R27DEM->GetCopy0();
            R27SAVY1X = R27DEM->GetCopy0();
            R27SAVY1Y = R27DEM->GetCopy0();
            R27SAVY2X = R27DEM->GetCopy0();
            R27SAVY2Y = R27DEM->GetCopy0();

            R9SAVX1X = R9DEM->GetCopy0();
            R9SAVX1Y = R9DEM->GetCopy0();
            R9SAVX2X = R9DEM->GetCopy0();
            R9SAVX2Y = R9DEM->GetCopy0();
            R9SAVY1X = R9DEM->GetCopy0();
            R9SAVY1Y = R9DEM->GetCopy0();
            R9SAVY2X = R9DEM->GetCopy0();
            R9SAVY2Y = R9DEM->GetCopy0();

            R3SAVX1X = R3DEM->GetCopy0();
            R3SAVX1Y = R3DEM->GetCopy0();
            R3SAVX2X = R3DEM->GetCopy0();
            R3SAVX2Y = R3DEM->GetCopy0();
            R3SAVY1X = R3DEM->GetCopy0();
            R3SAVY1Y = R3DEM->GetCopy0();
            R3SAVY2X = R3DEM->GetCopy0();
            R3SAVY2Y = R3DEM->GetCopy0();

            SAVX1X = DEM->GetCopy0();
            SAVX1Y = DEM->GetCopy0();
            SAVX2X = DEM->GetCopy0();
            SAVX2Y = DEM->GetCopy0();
            SAVY1X = DEM->GetCopy0();
            SAVY1Y = DEM->GetCopy0();
            SAVY2X = DEM->GetCopy0();
            SAVY2Y = DEM->GetCopy0();


        }


        //do a full-resolution flow accumulation for upstream pixel count
        //this is a least-required indication of which areas we might want to have at higher resolution in advance


        //set the actual elevation and manning data
        //FF::cout << "do 3-cascade average " <<std::endl;
        CreateMap3CascadeAvMax(DEM,MAN,RAIN,KSAT,R3DEM,R3DEMM,R3MAN,R3RAIN,R3KSAT);
        //FF::cout << "do 9-cascade average " <<std::endl;
        CreateMap3CascadeAvMax(R3DEM,R3MAN,R3RAIN,R3KSAT,R9DEM,R9DEMM,R9MAN,R9RAIN,R9KSAT);
        //FF::cout << "do 27-cascade average " <<std::endl;
        CreateMap3CascadeAvMax(R9DEM,R9MAN,R9RAIN,R9KSAT,R27DEM,R27DEMM,R27MAN,R27RAIN,R27KSAT);



        //FF::cout << "rain av3 " << MapAverageRedRP(R3RAIN) << std::endl;
        //FF::cout << "rain av9 " << MapAverageRedRP(R9RAIN) << std::endl;
        //FF::cout << "rain av27 " << MapAverageRedRP(R27RAIN) << std::endl;

        if(dem_channelized != nullptr)
        {
            FF::cout << "Channelized dem" <<std::endl;
            DEM_ORGCHAN = dem_channelized;
        }else{
            DEM_ORGCHAN = DEM->GetCopy();
        }

        CHANNELLDD = DEM->GetCopy0();
        SLOPE = DEM->GetCopy0();
        GradX1 = DEM->GetCopy0();
        GradX2 = DEM->GetCopy0();
        GradY1 = DEM->GetCopy0();
        GradY2 = DEM->GetCopy0();



        DEM2CHAN = SpreadMonotonicReconstruct(DEM_ORGCHAN,/*0.00001*/delta,iters);//SpreadDirectionalAbsMaxMD(Seeds,DEM,GradX1,GradX2,GradY1,GradY2);
        FF::cout << "do channnels2b" <<std::endl;

            for (int r = 0; r < DEM->Rows(); r++) {
                for (int c = 0; c < DEM->Cols(); c++) {
                    if (!pcr::isMV(DEM->data[r][c])) {

                        float gx1 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_L));
                        float gx2 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_R));
                        float gy1 =std::min(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_L));
                        float gy2 =std::max(0.0,-UF2D_Derivative_scaled_filldelta(DEM2CHAN,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_R));

                        float gtot = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));

                        SLOPE->data[r][c] = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));

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

        {
            cTMapf * rtemp = DEM->GetCopy1();
            //do flow accumulation
            FF::cout << "do channnels4" <<std::endl;
            double cellareakm = (in_dx * in_dx)/1000000.0;

            LDD4FromGradients(GradX1_i, GradX2_i, GradY1_i, GradY2_i, CHANNELLDD);

         QDRAIN = SpreadFlowMD2(rtemp,GradX1_i,GradX2_i,GradY1_i,GradY2_i,iters);

         //now set the refinement critera, based on the drainage area in pixels, and the maximum local rainfall (in case of discharge boundary condition)
         //very important for final efficiency, indicates what level of detail should be done at which locations, although water heights are going
         //to impact similarly. The thresholds are set here!
         //criteria is
         for (int r = 0; r < DEM->Rows(); r++) {
             for (int c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(DEM->data[r][c])) {

                    CRIT->data[r][c] = (1.0/
                                                     (1.0 + 25 *std::log(std::max(1.0,(QDRAIN->data[r][c]-10000.0))))
                                                     )
                            *
                            (1.0/(1.0+std::max(0.0,(RAIN->data[r][c] -1000.0)/1000.0)));
                 }else
                {
                    CRIT->data[r][c] = 1.0;
                }
             }
         }

         CreateMap3CascadeMin(CRIT,R3CRIT);
         CreateMap3CascadeMin(R3CRIT,R9CRIT);
         CreateMap3CascadeMin(R9CRIT,R27CRIT);

         //reserve the most detailed level though for real high level flow
         for (int r = 0; r < DEM->Rows(); r++) {
             for (int c = 0; c < DEM->Cols(); c++) {
                if (!pcr::isMV(DEM->data[r][c])) {

                    CRIT->data[r][c] = 1.0;
                 }
             }
         }


         //channels
         //Create channel LDD
         //use channel dims as a mask?
         if( (chanw != 0 && chanh != 0))
         {
             do_channels = true;

             CHANNELW = chanw;
             CHANNELD = chanh;

             R3CHANNELLDD = R3DEM->GetCopy0();
             R3CHANNELW = R3DEM->GetCopy0();
             R3CHANNELD = R3DEM->GetCopy0();
             R3SLOPE = R3DEM->GetCopy0();
             R3CHANNELH = R3DEM->GetCopy0();
             R3CHANNELHN = R3DEM->GetCopy0();

             R3CHANNELALPHA = R3DEM->GetCopy0();
             R3CHANNELQN = R3DEM->GetCopy0();
             R3CHANNELINFIL =R3DEM->GetCopy0();


             FF::cout << "do channnels3" <<std::endl;

             ResampleLDD3Scale(R3DEM,CHANNELLDD,CHANNELW,CHANNELD,QDRAIN,R3CHANNELLDD,R3CHANNELW,R3CHANNELD);

             //writeRaster(*CHANNELW,"chw1.tif");

             //writeRaster(*R3CHANNELW,"chw3.tif");

             FF::cout << "do channnels5" <<std::endl;
                 //get R3 average slopes
                 CreateMap3CascadeAv(SLOPE, R3SLOPE);

             FF::cout << "do channnels6" <<std::endl;
                 //delete temporary maps for flow accumulation
                 delete rtemp;

             }






         }


             FF::cout << "Channels done" <<std::endl;


           WHMax = DEM->GetCopy0();
           VMax = DEM->GetCopy0();
           ArrivalTime = DEM->GetCopy0();
           QCHMax = DEM->GetCopy0();




  }

  inline void SetCulverts(    int n_culverts, float * culvert_q, float * culvert_fraction,float * culverts_source_row, float * culverts_source_col,
                              float * culverts_target_row, float * culverts_target_col)
  {


      if(n_culverts > 0 )
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
          CULVERTID = DEM->GetCopy0();
           for (int r = 0; r < DEM->Rows(); r++) {
              for (int c = 0; c < DEM->Cols(); c++) {

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


  }
  inline void SetReservoirs(int shapes, float * points, float * slengths, float * stypes, float * respasscap, float * cap)
  {


      //reservoirs
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
              float qh = respasscap[i];
              reservoirpassvols.push_back(qh);
          }
      }


  }
  inline void SetStormDrains(cTMapf * stormdraindata)
  {
      //stormdrains
      bool do_stormdrain = false;
      if(stormdraindata != nullptr)
      {
          do_stormdrain = true;
          FF::cout <<" do stormdrain " <<std::endl;
          SDC = stormdraindata;
          SDCTOT = DEM->GetCopy0();
      }

  }
  inline void SetOceanBoundary(float boundary)
  {

      //ocean
      if(do_ocean > 0.5)
      {
          do_ocean = true;
          ocean_level = boundary;
      }
  }
  inline void SetWHInitial(float wh_initial)
  {

      //start with everywhere low resolution

     if(wh_initial > 1e-12f)
         {
           FF::cout << "set initial wh " << wh_initial << std::endl;
          for (int r = 0; r < R27DEM->Rows(); r++) {
               for (int c = 0; c < R27DEM->Cols(); c++) {
                   if (!pcr::isMV(R27DEM->data[r][c])) {
                     R27FlowH->data[r][c] = wh_initial;

                   }
               }
          }
     }
  }
  inline void SetCustomReturnPointers(float *in_ret, float * in_retvx,float * in_retvy,
                                      float * in_retq, float * in_ret_whmax, float * in_ret_vmax,
                                      float * in_ret_at)
  {
      ret =in_ret;
      retvx =in_retvx;
      retvy =in_retvy;
      retq =in_retq;
      ret_hmax =in_ret_whmax;
      ret_vmax =in_ret_vmax;
      ret_at =in_ret_at;
  }
  inline void Initialize()
  {

  }
  inline void Clear()
  {

  }
  inline ~FF_FULLFLOODT()
  {
        Clear();
  }
};

inline FF_FULLFLOODT FF_FULLFLOOD;

}

inline FF_FULLFLOODT * FullFloodStateFactory()
{
    return new FF_FULLFLOODT();
}

#endif
