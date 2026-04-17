#ifndef FULLFLOOD_STEP_H
#define FULLFLOOD_STEP_H

extern "C"
{

inline void EMSCRIPTEN_KEEPALIVE do_algo_fullflood_step_onstate(FF_FULLFLOODT &FULLFLOOD,float _dt, float courant, float dt_min, float courant_forced, int steps_per_history, float scale_req, bool use_openmp =true, bool do_flow = true)
{
    FF::cout << "step fullflood"  <<std::endl;

    float dx = FULLFLOOD.DEM->cellSizeX();
    if(dt_min < 0.0)
    {
        dt_min = 1.0 * dx/10.0;
    }

  std::vector<std::vector<cTMapf*>> flowstates =
  {
      //[0]->DEM
      //[1]->DEMM
      //[2]->SGA
      //[3]->MAN
      //[4]->H
      //[5]->U
      //[6]->V
      //[7]->HN
      //[8]->UN
      //[9]->VN
      //[10]->RAIN
      //[11]->Crit

      //[12]->T
      //[13]->DT

      //[12]->SGHOx1 //super grid original height
      //[13]->SGUXOx1 //super grid original x-vel
      //[14]->SGUYOx1 //super grid original y-vel

      //[12]->SGHOx2 //super grid original height
      //[13]->SGUXOx2 //super grid original x-vel
      //[14]->SGUYOx2 //super grid original y-vel

      //[12]->SGHOy1 //super grid original height
      //[13]->SGUXOy1 //super grid original x-vel
      //[14]->SGUYOy1 //super grid original y-vel

      //[12]->SGHOy2 //super grid original height
      //[13]->SGUXOy2 //super grid original x-vel
      //[14]->SGUYOy2 //super grid original y-vel

      //[15]->SGFHLx1 //super grid height flux limiting factor
      //[16]->R27SGFUXLx1 //super grid x-vel flux limiting factor
      //[17]->R27SGFUYLx1 //super grid y-vel flux limiting factor
      //[18]->R27SGFAYx1 //super-grid active
      //[19]->R27SGFAXx1 //super-grid active
      //[20]->R27SGFMXx1 //super-grid active
      //[21]->R27SGFMYx1 //super-grid active

      //[22]->SGFHLx2 //super grid height flux limiting factor
      //[23]->R27SGFUXLx2 //super grid x-vel flux limiting factor
      //[24]->R27SGFUYLx2 //super grid y-vel flux limiting factor
      //[25]->R27SGFAYx2 //super-grid active
      //[26]->R27SGFAXx2 //super-grid active
      //[27]->R27SGFMXx2 //super-grid active
      //[28]->R27SGFMYx2 //super-grid active

      //[29]->SGFHLy1 //super grid height flux limiting factor
      //[30]->R27SGFUXLy1 //super grid x-vel flux limiting factor
      //[31]->R27SGFUYLy1 //super grid y-vel flux limiting factor
      //[32]->R27SGFAYy1 //super-grid active
      //[33]->R27SGFAXy1 //super-grid active
      //[34]->R27SGFMXy1 //super-grid active
      //[35]->R27SGFMYy1 //super-grid active

      //[36]->SGFHLy2 //super grid height flux limiting factor
      //[37]->R27SGFUXLy2 //super grid x-vel flux limiting factor
      //[38]->R27SGFUYLy2 //super grid y-vel flux limiting factor
      //[39]->R27SGFAYy2 //super-grid active
      //[40]->R27SGFAXy2 //super-grid active
      //[41]->R27SGFMXy2 //super-grid active
      //[42]->R27SGFMYy2 //super-grid active



    {

     FULLFLOOD.R27DEM,
     FULLFLOOD.R27DEMM,
     FULLFLOOD.R27SGA,
     FULLFLOOD.R27MAN,
     FULLFLOOD.R27FlowH,
     FULLFLOOD.R27FlowUX,
     FULLFLOOD.R27FlowUY,
     FULLFLOOD.R27FlowHN,
     FULLFLOOD.R27FlowUXN,
     FULLFLOOD.R27FlowUYN,
     FULLFLOOD.R27RAIN,
     FULLFLOOD.R27CRIT,

     FULLFLOOD.R27FlowT, //12
     FULLFLOOD.R27FlowDT,

     FULLFLOOD.R27DHDT, //14
     FULLFLOOD.R27DUXDT,
     FULLFLOOD.R27DUYDT,


     FULLFLOOD.R27SGCA, //17
     FULLFLOOD.R27SGHO,
     FULLFLOOD.R27SGUXO,
     FULLFLOOD.R27SGUYO,

    FULLFLOOD.R27HLLX1X, //21
    FULLFLOOD.R27HLLX1Y,
    FULLFLOOD.R27HLLX1Z,
    FULLFLOOD.R27HLLX2X,
    FULLFLOOD.R27HLLX2Y,
    FULLFLOOD.R27HLLX2Z,
    FULLFLOOD.R27HLLY1X,
    FULLFLOOD.R27HLLY1Y,
    FULLFLOOD.R27HLLY1Z,
    FULLFLOOD.R27HLLY2X,
    FULLFLOOD.R27HLLY2Y,
    FULLFLOOD.R27HLLY2Z,

    FULLFLOOD.R27FlowTN,      //33
    FULLFLOOD.R27FlowDTL,

    FULLFLOOD.R27DHDTN, //35
    FULLFLOOD.R27DUXDTN,
    FULLFLOOD.R27DUYDTN,

    FULLFLOOD.R27SAVX1X, //38
    FULLFLOOD.R27SAVX1Y,
    FULLFLOOD.R27SAVX2X,
    FULLFLOOD.R27SAVX2Y,
    FULLFLOOD.R27SAVY1X,
    FULLFLOOD.R27SAVY1Y,
    FULLFLOOD.R27SAVY2X,
    FULLFLOOD.R27SAVY2Y,
    FULLFLOOD.R27FlowDTLN, //46
    FULLFLOOD.R27FlowDTC,
    FULLFLOOD.FlowTSUPERMINNEXT // 48

      },

    {FULLFLOOD.R9DEM,FULLFLOOD.R9DEMM,FULLFLOOD.R9SGA,FULLFLOOD.R9MAN,
     FULLFLOOD.R9FlowH,FULLFLOOD.R9FlowUX,FULLFLOOD.R9FlowUY,
     FULLFLOOD.R9FlowHN,FULLFLOOD.R9FlowUXN,FULLFLOOD.R9FlowUYN,
     FULLFLOOD.R9RAIN,FULLFLOOD.R9CRIT,
     FULLFLOOD.R9FlowT, FULLFLOOD.R9FlowDT,
     FULLFLOOD.R9DHDT,
     FULLFLOOD.R9DUXDT,
     FULLFLOOD.R9DUYDT,
     nullptr,
     nullptr,
     nullptr,
     nullptr,

     FULLFLOOD.R9HLLX1X, //21
     FULLFLOOD.R9HLLX1Y,
     FULLFLOOD.R9HLLX1Z,
     FULLFLOOD.R9HLLX2X,
     FULLFLOOD.R9HLLX2Y,
     FULLFLOOD.R9HLLX2Z,
     FULLFLOOD.R9HLLY1X,
     FULLFLOOD.R9HLLY1Y,
     FULLFLOOD.R9HLLY1Z,
     FULLFLOOD.R9HLLY2X,
     FULLFLOOD.R9HLLY2Y,
     FULLFLOOD.R9HLLY2Z,

     FULLFLOOD.R9FlowTN,      //33
     FULLFLOOD.R9FlowDTL,

     FULLFLOOD.R9DHDTN, //35
     FULLFLOOD.R9DUXDTN,
     FULLFLOOD.R9DUYDTN,

     FULLFLOOD.R9SAVX1X, //38
     FULLFLOOD.R9SAVX1Y,
     FULLFLOOD.R9SAVX2X,
     FULLFLOOD.R9SAVX2Y,
     FULLFLOOD.R9SAVY1X,
     FULLFLOOD.R9SAVY1Y,
     FULLFLOOD.R9SAVY2X,
     FULLFLOOD.R9SAVY2Y,

     FULLFLOOD.R9FlowDTLN, //46
     FULLFLOOD.R9FlowDTC,
     nullptr // 48
      },

    {FULLFLOOD.R3DEM,FULLFLOOD.R3DEMM,FULLFLOOD.R3SGA,FULLFLOOD.R3MAN,
     FULLFLOOD.R3FlowH,FULLFLOOD.R3FlowUX,FULLFLOOD.R3FlowUY,
     FULLFLOOD.R3FlowHN,FULLFLOOD.R3FlowUXN,FULLFLOOD.R3FlowUYN,
     FULLFLOOD.R3RAIN,FULLFLOOD.R3CRIT,
     FULLFLOOD.R3FlowT, FULLFLOOD.R3FlowDT,
     FULLFLOOD.R3DHDT,
     FULLFLOOD.R3DUXDT,
     FULLFLOOD.R3DUYDT,
     nullptr,
     nullptr,
     nullptr,
     nullptr,

     FULLFLOOD.R3HLLX1X, //21
     FULLFLOOD.R3HLLX1Y,
     FULLFLOOD.R3HLLX1Z,
     FULLFLOOD.R3HLLX2X,
     FULLFLOOD.R3HLLX2Y,
     FULLFLOOD.R3HLLX2Z,
     FULLFLOOD.R3HLLY1X,
     FULLFLOOD.R3HLLY1Y,
     FULLFLOOD.R3HLLY1Z,
     FULLFLOOD.R3HLLY2X,
     FULLFLOOD.R3HLLY2Y,
     FULLFLOOD.R3HLLY2Z,

     FULLFLOOD.R3FlowTN,      //33
     FULLFLOOD.R3FlowDTL,

     FULLFLOOD.R3DHDTN, //35
     FULLFLOOD.R3DUXDTN,
     FULLFLOOD.R3DUYDTN,

     FULLFLOOD.R3SAVX1X, //38
     FULLFLOOD.R3SAVX1Y,
     FULLFLOOD.R3SAVX2X,
     FULLFLOOD.R3SAVX2Y,
     FULLFLOOD.R3SAVY1X,
     FULLFLOOD.R3SAVY1Y,
     FULLFLOOD.R3SAVY2X,
     FULLFLOOD.R3SAVY2Y,

     FULLFLOOD.R3FlowDTLN, //46
     FULLFLOOD.R3FlowDTC,
     nullptr // 48
      },

    {FULLFLOOD.DEM,
     FULLFLOOD.DEM,
     FULLFLOOD.SGA,
     FULLFLOOD.MAN,

     FULLFLOOD.FlowH, //4
     FULLFLOOD.FlowUX, //5
     FULLFLOOD.FlowUY, //6

     FULLFLOOD.FlowHN, //7
     FULLFLOOD.FlowUXN, //8
     FULLFLOOD.FlowUYN, //9

     FULLFLOOD.RAIN, //10
     FULLFLOOD.CRIT, //11

     FULLFLOOD.FlowT, //12
     FULLFLOOD.FlowDT, //13

     FULLFLOOD.DHDT, //14
     FULLFLOOD.DUXDT, //15
     FULLFLOOD.DUYDT, //16

     nullptr,
     nullptr,
     nullptr,
     nullptr,

     FULLFLOOD.HLLX1X, //21
     FULLFLOOD.HLLX1Y, //22
     FULLFLOOD.HLLX1Z, //23

     FULLFLOOD.HLLX2X, //24
     FULLFLOOD.HLLX2Y, //25
     FULLFLOOD.HLLX2Z, //26

     FULLFLOOD.HLLY1X, //27
     FULLFLOOD.HLLY1Y, //28
     FULLFLOOD.HLLY1Z, //29

     FULLFLOOD.HLLY2X, //30
     FULLFLOOD.HLLY2Y, //31
     FULLFLOOD.HLLY2Z, //32

     FULLFLOOD.FlowTN, //33
     FULLFLOOD.FlowDTL,//34

     FULLFLOOD.DHDTN, //35
     FULLFLOOD.DUXDTN, //36
     FULLFLOOD.DUYDTN, //37

     FULLFLOOD.SAVX1X, //38
     FULLFLOOD.SAVX1Y,
     FULLFLOOD.SAVX2X, //40
     FULLFLOOD.SAVX2Y,
     FULLFLOOD.SAVY1X, //42
     FULLFLOOD.SAVY1Y,
     FULLFLOOD.SAVY2X, //44
     FULLFLOOD.SAVY2Y,

     FULLFLOOD.FlowDTLN, //46
     FULLFLOOD.FlowDTC,
     nullptr // 48
      },
  };




  FF::cout << "Resolutions - R27: " << FULLFLOOD.R27DEM->nrRows() << " " << FULLFLOOD.R27DEM->nrCols() <<std::endl;
  FF::cout << "Resolutions - R9: " << FULLFLOOD.R9DEM->nrRows() << " " << FULLFLOOD.R9DEM->nrCols() <<std::endl;
  FF::cout << "Resolutions - R3: " << FULLFLOOD.R3DEM->nrRows() << " " << FULLFLOOD.R3DEM->nrCols() <<std::endl;
  FF::cout << "Resolutions - R1: " << FULLFLOOD.DEM->nrRows() << " " << FULLFLOOD.DEM->nrCols() <<std::endl;

  //calculate rainfall factor
  float rain_fac = 0.0;
  float dur_half = 0.5 * FULLFLOOD.duration;
  float time_hour = FULLFLOOD.time_sec/3600.0;

  //TODO: implement shape
  if (time_hour < FULLFLOOD.duration)
  {
      rain_fac = 1.0;
  }else
  {
      rain_fac = 0.0;
  }

  //FF::cout << "rain curve " << rain_fac<< std::endl;

  //FF::cout << "Flow SubGriddify " << 27.0 * FULLFLOOD.DEM->cellSizeX() << " "  << FULLFLOOD.DEM->cellSizeX() <<std::endl;
  //only touches active cells
  FlowSubgriddify(_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,courant,scale_req);

  if((FULLFLOOD.do_cache && FULLFLOOD.step % 2 == 0) && do_flow)
  {
      FlowGridCachify(_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,FULLFLOOD.Cache_HOLD,FULLFLOOD.Cache_UXOLD,FULLFLOOD.Cache_UYOLD,courant,scale_req,FULLFLOOD.cache_rel_change_crit,FULLFLOOD.do_ocean,FULLFLOOD.ocean_level,rain_fac);
      FlowGridUnCachify(_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,FULLFLOOD.Cache_HOLD,FULLFLOOD.Cache_UXOLD,FULLFLOOD.Cache_UYOLD,courant,scale_req,FULLFLOOD.cache_rel_change_crit,FULLFLOOD.do_ocean,FULLFLOOD.ocean_level,rain_fac);
  }

  float currenttime = _dt * ((float)(FULLFLOOD.step));

   FULLFLOOD.step = FULLFLOOD.step +1;
  if((FULLFLOOD.step % 10 == 0 ) && do_flow)
  {
    //only touches active cells
    FlowUnSubgriddify(_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,courant,scale_req);
  }
  //FF::cout << "Flow Core Kernel " << std::endl;
  if(FULLFLOOD.do_channels)
  {
    //Touches all 1D cells in R3
    //and all active cells with channel in them
      /*writeRaster(*FULLFLOOD.R3CHANNELLDD,"ldd3.tif");
      writeRaster(*FULLFLOOD.R3CHANNELW,"lddw3.tif");
      writeRaster(*FULLFLOOD.R3CHANNELD,"lddd3.tif");
      writeRaster(*FULLFLOOD.R3SLOPE,"ldds3.tif");
      writeRaster(*FULLFLOOD.R3CHANNELH,"lddh3.tif");*/
    flow_1D2DKinematic(_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),3.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,FULLFLOOD.R3CHANNELLDD, FULLFLOOD.R3CHANNELW, FULLFLOOD.R3CHANNELD, FULLFLOOD.R3SLOPE, FULLFLOOD.R3CHANNELH, FULLFLOOD.R3CHANNELHN, FULLFLOOD.R3CHANNELALPHA, FULLFLOOD.R3CHANNELQN,FULLFLOOD.R3CHANNELINFIL,nullptr,nullptr,nullptr,nullptr,courant, dt_min,courant_forced);

  }
  //FF::cout << "do 2d steps" <<std::endl;

#ifdef FASTFLOOD_WASM
  if(use_openmp)
  {
      flow_saintvenantdtandflow3Casc(FULLFLOOD.time_sec,_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,nullptr,nullptr,nullptr,nullptr,courant, dt_min,courant_forced,FULLFLOOD.do_ocean,FULLFLOOD.ocean_level,rain_fac,do_flow);
  }else
  {
      flow_saintvenantdtandflow3Casc_threadpool(FULLFLOOD.time_sec,_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,nullptr,nullptr,nullptr,nullptr,courant, dt_min,courant_forced,FULLFLOOD.do_ocean,FULLFLOOD.ocean_level,rain_fac,do_flow);
  }
#else
  flow_saintvenantdtandflow3Casc_threadpool(FULLFLOOD.time_sec,_dt,27.0 * FULLFLOOD.DEM->cellSizeX(),flowstates,nullptr,nullptr,nullptr,nullptr,courant, dt_min,courant_forced,FULLFLOOD.do_ocean,FULLFLOOD.ocean_level,rain_fac,do_flow);
#endif
  //FF::cout << "Convert For Visualize " <<std::endl;

  FULLFLOOD.time_sec += _dt;

  //convert the cascading structure to full raster output

  FULLFLOOD.WHMax_History = {};
  FULLFLOOD.VxMax_History = {};
  FULLFLOOD.VyMax_History = {};

  float rain = 0.0;

  float hmax = 0.0;
  float level_avg = 0.0;
  double vol_total = 0.0;
  float qout = 0.0;
  double pix_total = 0.0;
  for (int r = 0; r < FULLFLOOD.DEM->Rows(); r++) {
    for (int c = 0; c < FULLFLOOD.DEM->Cols(); c++) {

                  int index = r* FULLFLOOD.FlowH->Cols() + c;

                if (!pcr::isMV(FULLFLOOD.DEM->data[r][c])) {

                  float ratio_cache = 1.0;


                  float critreq =1.0;
                  float dtc = 0.0;

                  float dth = 0.0;
                  if(FULLFLOOD.FlowDT != nullptr)
                  {
                      dth = FULLFLOOD.FlowDT->data[r][c];
                      dtc += FULLFLOOD.FlowT->data[r][c];
                  }
                  float level = 0;
                  float channelh = 0.0;
                  float channelfrac = 0.0;
                  float channelq = 0.0;

                    float h = FULLFLOOD.FlowH->data[r][c];
                    float vx = FULLFLOOD.FlowUX->data[r][c];
                    float vy = FULLFLOOD.FlowUY->data[r][c];

                    int r_ip = r/3;
                    int c_ip = c/3;
                    if(r_ip > -1 && c_ip > -1 && r_ip < FULLFLOOD.R3FlowH->nrRows() && c_ip < FULLFLOOD.R3FlowH->nrCols())
                    {
                        critreq = FULLFLOOD.R3CRIT->data[r_ip][c_ip];
                        if(FULLFLOOD.FlowDT != nullptr)
                        {
                            dtc += FULLFLOOD.R3FlowT->data[r_ip][c_ip];
                        }
                        h += FULLFLOOD.R3FlowH->data[r_ip][c_ip];
                        vx += FULLFLOOD.R3FlowUX->data[r_ip][c_ip];
                        vy += FULLFLOOD.R3FlowUY->data[r_ip][c_ip];

                        if(FULLFLOOD.do_channels)
                        {
                            if(FULLFLOOD.R3CHANNELLDD->data[r_ip][c_ip] == 5)
                            {
                                channelq =FULLFLOOD.R3CHANNELQN->data[r_ip][c_ip];
                                qout += FULLFLOOD.R3CHANNELQN->data[r_ip][c_ip];
                            }
                            channelh = FULLFLOOD.R3CHANNELH->data[r_ip][c_ip];
                            channelfrac = std::min(1.0f,FULLFLOOD.CHANNELD->data[r][c]/std::max(1e-12f,FULLFLOOD.R3CHANNELD->data[r_ip][c_ip]));
                        }
                        if(FULLFLOOD.R3SGA->data[r_ip][c_ip] > 0.5)
                        {
                          level = 10;
                          if(FULLFLOOD.FlowDT != nullptr)
                          {
                              dth = FULLFLOOD.R3FlowDT->data[r_ip][c_ip];
                          }
                        }

                        if(FULLFLOOD.do_channels)
                        {
                            h = h + FULLFLOOD.R3CHANNELH->data[r_ip][c_ip] * channelfrac;
                        }
                    }
                    r_ip = r/9;
                    c_ip = c/9;
                    if(r_ip > -1 && c_ip > -1 && r_ip < FULLFLOOD.R9FlowH->nrRows() && c_ip < FULLFLOOD.R9FlowH->nrCols())
                    {

                        rain= FULLFLOOD.R9RAIN->data[r_ip][c_ip];
                        if(FULLFLOOD.FlowDT != nullptr)
                        {
                            dtc += FULLFLOOD.R9FlowT->data[r_ip][c_ip];
                        }
                        h += FULLFLOOD.R9FlowH->data[r_ip][c_ip];
                        vx += FULLFLOOD.R9FlowUX->data[r_ip][c_ip];
                        vy += FULLFLOOD.R9FlowUY->data[r_ip][c_ip];
                        if(FULLFLOOD.R9SGA->data[r_ip][c_ip] > 0.5 && level < 9 )
                        {
                          level = 1;
                          if(FULLFLOOD.FlowDT != nullptr)
                          {
                              dth = FULLFLOOD.R9FlowDT->data[r_ip][c_ip];
                          }
                        }
                    }

                    float is_cache = 0.0;

                    r_ip = r/27;
                    c_ip = c/27;
                    if(r_ip > -1 && c_ip > -1 && r_ip < FULLFLOOD.R27FlowH->nrRows() && c_ip < FULLFLOOD.R27FlowH->nrCols())
                    {

                        if(FULLFLOOD.do_cache)
                        {

                            if(FULLFLOOD.R27SGCA != nullptr)
                            {
                               is_cache = FULLFLOOD.R27SGCA->data[r_ip][c_ip];
                            }
                        }

                        if(is_cache)
                        {
                            ratio_cache = FULLFLOOD.R27FlowH->data[r_ip][c_ip]/std::max(1e-12f,FULLFLOOD.R27SGHO->data[r_ip][c_ip]);

                            h = h * ratio_cache;
                            vx = vx * FULLFLOOD.R27FlowUX->data[r_ip][c_ip]/(FULLFLOOD.R27SGUXO->data[r_ip][c_ip] < 0.0?std::min(-1e-12f,FULLFLOOD.R27SGUXO->data[r_ip][c_ip]):std::max(1e-12f,FULLFLOOD.R27SGUXO->data[r_ip][c_ip]));
                            vy = vy * FULLFLOOD.R27FlowUY->data[r_ip][c_ip]/(FULLFLOOD.R27SGUYO->data[r_ip][c_ip] < 0.0?std::min(-1e-12f,FULLFLOOD.R27SGUYO->data[r_ip][c_ip]):std::max(1e-12f,FULLFLOOD.R27SGUYO->data[r_ip][c_ip]));
                        }else
                        {
                            h += FULLFLOOD.R27FlowH->data[r_ip][c_ip];
                            vx += FULLFLOOD.R27FlowUX->data[r_ip][c_ip];
                            vy += FULLFLOOD.R27FlowUY->data[r_ip][c_ip];
                        }
                        if(FULLFLOOD.FlowDT != nullptr)
                        {
                            dtc += FULLFLOOD.R27FlowT->data[r_ip][c_ip];
                        }
                        if(FULLFLOOD.R27SGA->data[r_ip][c_ip] > 0.5 && level < 0.9 )
                        {
                          level = 0.5;
                          if(FULLFLOOD.FlowDT != nullptr)
                          {
                              dth = FULLFLOOD.R27FlowDT->data[r_ip][c_ip];

                          }
                        }
                    }

                    FULLFLOOD.WHMax->data[r][c] = std::max(FULLFLOOD.WHMax->data[r][c],h);
                    if(FULLFLOOD.do_channels)
                    {
                        FULLFLOOD.QCHMax->data[r][c] = std::max(  FULLFLOOD.QCHMax->data[r][c],FULLFLOOD.CHANNELD->data[r][c] * channelfrac);
                    }

                    FULLFLOOD.VMax->data[r][c]  = std::max(FULLFLOOD.VMax->data[r][c], std::sqrt(vx*vx+vy*vy));

                    if(h > 0.1)
                    {
                        if(FULLFLOOD.ArrivalTime->data[r][c] == 0.0f)
                        {
                            FULLFLOOD.ArrivalTime->data[r][c] = currenttime;
                        }
                    }

                  if(!pcr::isMV(FULLFLOOD.FlowH->data[r][c]))
                  {

                    FULLFLOOD.ret[index] =h + channelfrac *channelh;

                    if(FULLFLOOD.retvx != nullptr)
                    {
                        FULLFLOOD.retvx[index] = vx;//;//channelh;
                        FULLFLOOD.retvy[index] = vy;//;//channelh;


                        if(FULLFLOOD.do_cache )
                        {
                            if(is_cache)
                            {
                                //FF::cout << "export ratio cache: "  << ratio_cache << " " << r_ip << "  " << c_ip << " " << FULLFLOOD.R27FlowH->data[r_ip][c_ip] << " " << std::max(1e-12f,FULLFLOOD.R27SGHO->data[r_ip][c_ip]) << "  "  << h <<  "  " << FULLFLOOD.R27FlowH->data[r_ip][c_ip]/std::max(1e-12f,FULLFLOOD.R27SGHO->data[r_ip][c_ip]) << std::endl;

                                if(std::isnan(h) || !std::isfinite(h) || std::isnan(vx) || !std::isfinite(vx)|| std::isnan(vy) || !std::isfinite(vy))
                                {
                                    //FF::cout << "export ratio cache NAN found : "  << ratio_cache << " " << FULLFLOOD.R27FlowH->data[r_ip][c_ip] << " " << std::max(1e-12f,FULLFLOOD.R27SGHO->data[r_ip][c_ip]) << "  "  << h <<  std::endl;
                                }
                            }
                            FULLFLOOD.retq[index] = h * dx * std::fabs(vx) + h*dx * std::fabs(vy) + channelfrac * channelq;//std::sqrt(vx*vx+vy*vy);//FULLFLOOD.R27SGCA->data[r_ip][c_ip] > 0.5f? FULLFLOOD.R27FlowH->data[r_ip][c_ip]/std::max(1e-12f,FULLFLOOD.R27SGHO->data[r_ip][c_ip]):1.0;//h * dx * vx + h*dx * vy;//channelh;
                        }else
                        {
                            FULLFLOOD.retq[index] = h * dx * std::fabs(vx) + h*dx * std::fabs(vy) + channelfrac * channelq;//std::sqrt(vx*vx+vy*vy);//std::sqrt(vx*vx+vy*vy);//h * dx * vx + h*dx * vy;
                        }
                    }
                    vol_total += h*dx*dx;
                    pix_total += 1.0*dx*dx;
                    level_avg += level;

                    FULLFLOOD.ret_at[index] =FULLFLOOD.ArrivalTime->data[r][c];
                    FULLFLOOD.ret_hmax[index] =FULLFLOOD.WHMax->data[r][c];
                    FULLFLOOD.ret_vmax[index] =FULLFLOOD.VMax->data[r][c];

                    //are we at model edge, or at mv edge?
                    if(r == 0 || r == FULLFLOOD.DEM->nrRows() || c == 0 || c == FULLFLOOD.DEM->nrCols())
                    {
                        qout += h* std::fabs(vx) * dx + h* std::fabs(vy) * dx;
                    }else
                    {
                        if(OUTORMV(FULLFLOOD.DEM,r-1,c)||
                                OUTORMV(FULLFLOOD.DEM,r+1,c)||
                                OUTORMV(FULLFLOOD.DEM,r,c+1)||
                                OUTORMV(FULLFLOOD.DEM,r,c-1))
                        {
                            qout += h* std::fabs(vx) * dx + h* std::fabs(vy) * dx;
                        }
                    }

                    hmax = std::max(hmax, h);


                  }else{
                    pcr::setMV(FULLFLOOD.ret[index]);
                  }
                }else{
                    pcr::setMV(FULLFLOOD.ret[index]);
                  }
            }
  }

  FULLFLOOD.qout = qout;
#ifdef FASTFLOOD_WASM
    emscripten_run_script("increment_progress_fullflood()");
#endif

    FF::cout << "FastFlood Full hours/VOl/RAIN/MB/RefRa - " << FULLFLOOD.time_sec/3600.0 << " " << vol_total/std::max(1.0,pix_total) << " " << FULLFLOOD.rain_total/std::max(1.0,pix_total) << " " << vol_total/FULLFLOOD.rain_total << "  " << level_avg/std::max(1.0,pix_total) << std::endl;
}



inline void EMSCRIPTEN_KEEPALIVE do_algo_fullflood_step(float _dt, float courant, float dt_min, float courant_forced, int steps_per_history, float scale_req, bool use_openmp =true, bool do_flow = true)
{
    do_algo_fullflood_step_onstate(FF_FULLFLOOD,_dt,courant,dt_min,courant_forced,steps_per_history,scale_req,use_openmp,do_flow);

}

}
#endif

