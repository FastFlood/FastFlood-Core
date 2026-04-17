#ifndef FULLFLOOD_SVSOLVE_H
#define FULLFLOOD_SVSOLVE_H


#include "../common/defines.h"


#define OMP_COLLAPSE
//no collapse, otherwise, put collapse(2)

#define OMP_MODE dynamic,1
//dynamic means work stealing, 1 item at a time, otherwise for pre-distributed chunks, set #define OMP_MODE static

static inline void flow_saintvenantdtandflow3Casc(double tglob_in,  double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates, cTMapf * QX1 = nullptr, cTMapf * QX2 = nullptr,cTMapf * QY1 = nullptr,cTMapf * QY2 = nullptr, float courant = 0.1, float mindt = 0.0, float courant_forced = 0.25f, bool do_ocean = false, float ocean_level = 0.0, float rain_fac = 1.0, bool do_flow = true)
{
    //tglob_in;

    /*std::vector<std::vector<cTMapf*>> flowstates =
      {
        {FF_FULLFLOOD.R27DEM,FF_FULLFLOOD.R27DEMM,FF_FULLFLOOD.R27SGA,FF_FULLFLOOD.R27MAN,FF_FULLFLOOD.R27FlowH,FF_FULLFLOOD.R27FlowUX,FF_FULLFLOOD.R27FlowUY,FF_FULLFLOOD.R27FlowHN,FF_FULLFLOOD.R27FlowUXN,FF_FULLFLOOD.R27FlowUYN},
        {FF_FULLFLOOD.R9DEM,FF_FULLFLOOD.R9DEMM,FF_FULLFLOOD.R3SGA,FF_FULLFLOOD.R9MAN,FF_FULLFLOOD.R9FlowH,FF_FULLFLOOD.R9FlowUX,FF_FULLFLOOD.R9FlowUY,FF_FULLFLOOD.R9FlowHN,FF_FULLFLOOD.R9FlowUXN,FF_FULLFLOOD.R9FlowUYN},
        {FF_FULLFLOOD.R3DEM,FF_FULLFLOOD.R9DEMM,FF_FULLFLOOD.R3SGA,FF_FULLFLOOD.R3MAN,FF_FULLFLOOD.R3FlowH,FF_FULLFLOOD.R3FlowUX,FF_FULLFLOOD.R3FlowUY,FF_FULLFLOOD.R3FlowHN,FF_FULLFLOOD.R3FlowUXN,FF_FULLFLOOD.R3FlowUYN},
        {FF_FULLFLOOD.DEM,FF_FULLFLOOD.DEM,FF_FULLFLOOD.SGA,FF_FULLFLOOD.MAN,FF_FULLFLOOD.FlowH,FF_FULLFLOOD.FlowUX,FF_FULLFLOOD.FlowUY,FF_FULLFLOOD.FlowHN,FF_FULLFLOOD.FlowUXN,FF_FULLFLOOD.FlowUYN},
      }*/

    FF::cout << "saint venant solve: " << _dt << " "  << mindt << " " << courant << std::endl;
    cTMapf * DEM = flowstates[0][0];

    int dim0 = DEM->nrCols();
    int dim1 = DEM->nrRows();


    double dt = 1e6;

    double dtlast = 1e6;
    float dt_min = 1e6;
    float dt_store = _dt;
    bool stop = false;
    #pragma omp parallel shared(stop,dt_min,dt,dtlast)
    {
        double t = 0;
        double t_end = _dt;
        double tglob = 0.0;
        int iter = 0;

        if(flowstates[0][21] != nullptr)
        {
            auto fi = [&flowstates](int l, int r, int c, auto &f){

                if(r >= flowstates[l][0]->nrRows() || c >= flowstates[l][0]->nrCols())
                {
                    return;
                }
                flowstates[l][12]->Drc = 0.0;
                flowstates[l][33]->Drc = 0.0;
                flowstates[l][47]->Drc = 0.0;

                if( l < flowstates.size()-1)
                {
                    for(int n = 0; n < 3; n++)
                    {
                        for(int k = 0; k < 3; k++)
                        {
                            f(l+1,r * 3 + n,c *3 + k,f);
                        }
                    }

                }
            };

          #pragma omp for OMP_COLLAPSE schedule(OMP_MODE)
          for(int r = 0; r < flowstates[0][4]->nrRows(); r++)
          {
              for(int c = 0; c < flowstates[0][4]->nrCols(); c++)
              {

                  flowstates[0][48]->Drc = 0.0;
                  fi(0,r,c,fi);
              }

          }
        }



        while(true)
        {

            #pragma omp barrier
            #pragma omp master
            {
                //FF::cout << "time : " << t << " " << t_end  << " " << std::endl;

                dt_min = 1e6;
                dt = 1e6;
            }
            #pragma omp barrier

            if(!(t < t_end))
            {
                stop = true;
            }
            #pragma omp barrier
            if(stop)
            {
                break;
            }

            double dt_min_this = 1e6;
            double dtleft = std::max(0.0,_dt-t);
             //get dt
            #pragma omp for OMP_COLLAPSE reduction(min:dt_min,dt) schedule(OMP_MODE)
            for(int r = 0; r < DEM->nrRows(); r++)
            {
                for (int c = 0; c < DEM->nrCols(); c++)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        if(flowstates[0][21] != nullptr)
                        {
                            //if we didnt do the cell last time
                            //if((flowstates[0][47]->data[r][c] <= 1e-12f))
                            {
                                //check if this cell requires an update
                                if(tglob + 1.5 * dtlast < flowstates[0][48]->data[r][c])
                                {
                                    //FF::cout << "o"<< std::endl;
                                    continue;
                                    //or if any neighbor requires an upgrade?
                                }
                            }
                        }

                        //FF::cout << "i"<< std::endl;
                        float* ref_tnextmin = nullptr;
                        if(flowstates[0][21] != nullptr)
                        {
                            flowstates[0][48]->data[r][c] = t_end;
                            ref_tnextmin = &(flowstates[0][48]->data[r][c]);
                        }

                        //for all pixels in this supergrid pixel
                        double dt_req = GetMinDT3Casc(tglob, courant,_dt,dx, flowstates,r,c,0,mindt,_dt-t,ref_tnextmin);

                        double dth = std::min(_dt-t,std::max((double)mindt,dt_req));
                        dt_min = std::min((double)dt_min,dth);
                        dt = std::min(_dt-t,std::max((double)mindt,std::min((double)dt_min,dt)));

                        if(flowstates[0][21] != nullptr)
                        {
                            flowstates[0][13]->data[r][c] = dt_req;
                        }
                    }
                }
            }

            //wait
            #pragma omp barrier
            #pragma omp master
            {
                dt = std::min(1.5 * dtlast,dt);
              //FF::cout << " dt: " << _dt << " " << tglob << " "  << dt << " " << dt_min  << " " << dtleft << std::endl;
            }

    #pragma omp barrier
            //if sdt, set the actual timestepping flags
            if(flowstates[0][21] != nullptr)
            {
                #pragma omp for OMP_COLLAPSE schedule(OMP_MODE)
                for(int r = 0; r < DEM->nrRows(); r++)
                {
                    for (int c = 0; c < DEM->nrCols(); c++)
                    {
                        if(!pcr::isMV(DEM->data[r][c]))
                        {
                            float* ref_sp = nullptr;
                            if(flowstates[0][21] != nullptr)
                            {
                                flowstates[0][47]->data[r][c] = 0.0;
                                ref_sp = &(flowstates[0][47]->data[r][c]);
                            }


                            if(flowstates[0][21] != nullptr)
                            {
                                //if we didnt do the cell last time
                                //if(!(flowstates[0][47]->data[r][c] <= 1e-12f))
                                {
                                    //check if this cell requires an update
                                    if(tglob +dt < flowstates[0][48]->data[r][c])
                                    {
                                        if(tglob +dt < flowstates[0][48]->data[std::max(0,std::min(flowstates[0][48]->nrRows()-1,r+1))][c])
                                        {
                                            if(tglob +dt < flowstates[0][48]->data[std::max(0,std::min(flowstates[0][48]->nrRows()-1,r-1))][c])
                                            {
                                                if(tglob +dt < flowstates[0][48]->data[r][std::max(0,std::min(flowstates[0][48]->nrCols()-1,c+1))])
                                                {
                                                    if(tglob +dt < flowstates[0][48]->data[r][std::max(0,std::min(flowstates[0][48]->nrCols()-1,c-1))])
                                                    {
                                                        continue;
                                                        //or if any neighbor requires an upgrade?
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                                SetMinDT3Casc(tglob, _dt,dt,flowstates, r, c, 0, dt_min, dtleft,ref_sp);
                        }
                    }
                }
            }


            //wait
            #pragma omp barrier

            #pragma omp for OMP_COLLAPSE schedule(OMP_MODE)
            for(int r = 0; r < DEM->nrRows(); r++)
            {
                for (int c = 0; c < DEM->nrCols(); c++)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        if(flowstates[0][21] != nullptr)
                        {
                            if(flowstates[0][47]->data[r][c] < 0.5f)
                            {
                                continue;
                            }
                        }
                        //flowstates[0][7]->data[r][c] += rain_fac * flowstates[0][10]->data[r][c] * dt/(3600000.0);
                        DoMassMomentum3Casc(tglob,courant,dt,dx, flowstates,r,c,0, do_ocean, ocean_level,rain_fac);
                    }
                }
            }

            //wait
        #pragma omp barrier
            //for(int i = 0; i < flowstates.size(); i++)
            {

                auto f = [&flowstates](int l, int r, int c, auto &f){

                    if(r >= flowstates[l][0]->nrRows() || c >= flowstates[l][0]->nrCols())
                    {
                        return;
                    }
                    if(l == 0 && flowstates[l][17] != nullptr)
                    {
                        if(flowstates[l][17]->data[r][c] > 0.5f)
                        {
                            flowstates[l][4]->Drc = flowstates[l][7]->Drc; //H
                            flowstates[l][5]->Drc = flowstates[l][8]->Drc; //Ux
                            flowstates[l][6]->Drc = flowstates[l][9]->Drc; //Uy
                            if(flowstates[l][21] != nullptr)
                            {
                                flowstates[l][14]->data[r][c] = flowstates[l][35]->data[r][c];
                                flowstates[l][15]->data[r][c] = flowstates[l][36]->data[r][c];
                                flowstates[l][16]->data[r][c] = flowstates[l][37]->data[r][c];

                                flowstates[l][12]->data[r][c] = flowstates[l][33]->data[r][c];
                                flowstates[l][34]->data[r][c] = flowstates[l][46]->data[r][c];
                            }
                            return;
                        }
                    }
                    if(flowstates[l][2]->data[r][c] > 0.5 && l < flowstates.size()-1)
                    {
                        for(int n = 0; n < 3; n++)
                        {
                            for(int k = 0; k < 3; k++)
                            {
                                f(l+1,r * 3 + n,c *3 + k,f);
                            }
                        }

                    }else
                    {
                        flowstates[l][4]->Drc = flowstates[l][7]->Drc; //H
                        flowstates[l][5]->Drc = flowstates[l][8]->Drc; //Ux
                        flowstates[l][6]->Drc = flowstates[l][9]->Drc; //Uy
                        if(flowstates[l][21] != nullptr)
                        {
                            flowstates[l][14]->data[r][c] = flowstates[l][35]->data[r][c];
                            flowstates[l][15]->data[r][c] = flowstates[l][36]->data[r][c];
                            flowstates[l][16]->data[r][c] = flowstates[l][37]->data[r][c];
                            flowstates[l][12]->data[r][c] = flowstates[l][33]->data[r][c];
                            flowstates[l][34]->data[r][c] = flowstates[l][46]->data[r][c];
                        }
                    }
                };

              #pragma omp for OMP_COLLAPSE schedule(OMP_MODE)
              for(int r = 0; r < flowstates[0][4]->nrRows(); r++)
              {
                  for(int c = 0; c < flowstates[0][4]->nrCols(); c++)
                  {
                      if(!pcr::isMV(DEM->data[r][c]))
                      {
                          if(flowstates[0][21] != nullptr)
                          {
                              if(flowstates[0][47]->data[r][c] < 0.5f)
                              {
                                  continue;
                              }
                          }

                          f(0,r,c,f);
                      }
                  }
              }
            }

            //wait


            tglob = tglob +dt;
            t = t+dt;
            dtlast = dt;

            iter ++;

        }

        #pragma omp master
        {
            FF::cout << "sub iteration done, time: " << t << " steps: " << iter <<std::endl;
        }
    }
    return;
}
template<typename F>
static inline void Parallelize_Loop(int start, int end,const F &f)
{
    /*for(int i = start; i < end; i++)
    {
        f(i,i+1);
    }*/


    //parallel_for(start, end, [&] (int i) { f(i))};


    /*int n_threads = 4;// std::thread::hardware_concurrency();

    //make blocks
    std::vector<std::thread> threads;

    int n = end-start;
    int n_per = n/n_threads;

    //FF::cout << n << " " <<n_per << " " << n_threads << std::endl;

    for(int i = 0; i < n_threads; i++)
    {
        int start_h = start + i*n_per;
        int end_h = std::min(end,start + (i+1) * n_per);
        if(i == n_threads - 1)
        {
            end_h = end;
        }
        //FF::cout << " i " << start_h << " " << end_h << std::endl;
        threads.push_back(std::thread([start_h, end_h,f]()
        {
            //FF::cout << " call f " << start_h << " " << end_h << std::endl;
            f(start_h, end_h);
        }));
    }
    std::for_each(threads.begin(), threads.end(), [](std::thread &t)
    {
        //FF::cout << " join " <<std::endl;
        t.join();
    });*/


}

inline double g_mindt_thread[256];

static inline void flow_saintvenantdtandflow3Casc_threadpool(double tglob_in, double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates, cTMapf * QX1 = nullptr, cTMapf * QX2 = nullptr,cTMapf * QY1 = nullptr,cTMapf * QY2 = nullptr, float courant = 0.1, float mindt = 0.0, float courant_forced = 0.25f, bool do_ocean = false, float ocean_level = 0.0, float rain_fac = 1.0, bool do_flow = true)
{

    cTMapf * DEM = flowstates[0][0];
    std::mutex m;

    int dim0 = DEM->nrCols();
    int dim1 = DEM->nrRows();

    double dt = 1e6;

    double tglob = 0.0;
    double dtlast = 1e6;
    float dt_min = 1e6;
    float dt_store = _dt;
    bool stop = false;

    double t = 0;
    double t_end = _dt;

    int iter = 0;


    std::vector<double> g_mindt_row;
    int n_rows = DEM->nrRows();
    for(int i=0; i < DEM->nrRows();i++)
    {
        g_mindt_row.push_back(1e6);
    }

    if(flowstates[0][21] != nullptr)
    {
        auto fi = [&flowstates](int l, int r, int c, auto &f){

            if(r >= flowstates[l][0]->nrRows() || c >= flowstates[l][0]->nrCols())
            {
                return;
            }
            flowstates[l][12]->Drc = 0.0;
            flowstates[l][33]->Drc = 0.0;
            flowstates[l][47]->Drc = 0.0;



            if( l < flowstates.size()-1)
            {
                for(int n = 0; n < 3; n++)
                {
                    for(int k = 0; k < 3; k++)
                    {
                        f(l+1,r * 3 + n,c *3 + k,f);
                    }
                }

            }
        };

      GetTPool().parallel_for(0,DEM->nrRows(),[&](int j){
              int r = j;
              for(int c = 0; c < flowstates[0][4]->nrCols(); c++)
              {
              flowstates[0][48]->Drc = 0.0;
              fi(0,r,c,fi);
              }
          });


    }




    while(true)
    {
        dt_min = 1e6;
        dt = 1e6;

        if(!(t < t_end))
        {
            stop = true;
        }
        if(stop)
        {
            break;
        }

        FF::cout << "tglob: "<<tglob << " " << t << std::endl;

        double dt_min_this = 1e6;
        double dtleft = std::max(0.0,_dt-t);

        for(int i = 0; i < n_rows; i++)
        {
            g_mindt_row[i] = 1e6;
        }

         //get dt
        GetTPool().parallel_for(0,DEM->nrRows(),[&](int j){

            int r = j;
            for(int c = 0; c < flowstates[0][4]->nrCols(); c++)
            {

            double dtleft = std::max(0.0,_dt-t);
                float dt_min_here = dt_min;
                double dt_here = dt;

                    if(!pcr::isMV(DEM->data[r][c]))
                    {

                        if(flowstates[0][21] != nullptr)
                        {
                            //if we didnt do the cell last time
                            //if((flowstates[0][47]->data[r][c] <= 1e-12f))
                            {
                                //check if this cell requires an update
                                if(tglob + 1.5 * dtlast < flowstates[0][48]->data[r][c])
                                {
                                    //FF::cout << "o"<< std::endl;
                                    continue;
                                    //or if any neighbor requires an upgrade?

                                }
                            }
                        }

                        //FF::cout << "i"<< std::endl;
                        float* ref_tnextmin = nullptr;
                        if(flowstates[0][21] != nullptr)
                        {
                            flowstates[0][48]->data[r][c] = t_end;
                            ref_tnextmin = &(flowstates[0][48]->data[r][c]);
                        }

                        //for all pixels in this supergrid pixel
                        double dt_req = GetMinDT3Casc(tglob, courant,_dt,dx, flowstates,r,c,0,mindt,_dt-t,ref_tnextmin);

                        dt_min_here = std::min(dtleft,std::max((double)mindt,std::min((double)dt_min_here,dt_req)));
                        dt_here = std::min(dtleft,std::max((double)mindt,std::min((double)dt_min_here,dt_here)));

                        if(flowstates[0][21] != nullptr)
                        {
                            flowstates[0][13]->data[r][c] = dt_req;
                        }
                    }
                    //now combine back all of the dt's
                    //Guard access with mutex
                    //m.lock();
                    //dt_min = std::min(dt_min,dt_min_here);
                    //dt = std::min(dt,dt_here);
                    //m.unlock();

                    g_mindt_row[r] = std::min(g_mindt_row[r],(double)dt_min_here);
                    //do a per thread min first
                    //std::this_thread::get_id();
                    //g_mindt_thread[tid] = std::min(g_mindt_thread[tid],dt_here);

                }


            });
        GetTPool().waitext();


        for(int i = 0; i < n_rows; i++)
        {
            dt = std::min(g_mindt_row[i],dt);
        }
        dt = std::min(1.5 * dtlast,dt);
        dt_min = dt;


        FF::cout << "dt  " << dt << std::endl;

        if(flowstates[0][21] != nullptr)
        {
            GetTPool().parallel_for(0,DEM->nrRows(),[&](int j){
                //for(int j = start; j < end; j++)
                {
                    //FF::cout << "do row " << j << std::endl;
                    int r = j;
                    for(int c = 0; c < flowstates[0][4]->nrCols(); c++)
                    {

                        //FF::cout << "rc " << r << c << std::endl;
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        float* ref_sp = nullptr;
                        if(flowstates[0][21] != nullptr)
                        {
                            flowstates[0][47]->data[r][c] = 0.0;
                            ref_sp = &(flowstates[0][47]->data[r][c]);
                        }


                        if(flowstates[0][21] != nullptr)
                        {

                            //if we didnt do the cell last time
                            //if(!(flowstates[0][47]->data[r][c] <= 1e-12f))
                            {
                                if(r == 15 && c==15)
                                {
                                    //FF::cout << "check set dt " << tglob << " " << dt << " " << tglob +dt << "  sgtimemin: " << flowstates[0][48]->data[r][c] <<std::endl;
                                }

                                //check if this cell requires an update
                                if(tglob +dt < flowstates[0][48]->data[r][c])
                                {
                                    if(tglob +dt < flowstates[0][48]->data[std::max(0,std::min(flowstates[0][48]->nrRows()-1,r+1))][c])
                                    {
                                        if(tglob +dt < flowstates[0][48]->data[std::max(0,std::min(flowstates[0][48]->nrRows()-1,r-1))][c])
                                        {
                                            if(tglob +dt < flowstates[0][48]->data[r][std::max(0,std::min(flowstates[0][48]->nrCols()-1,c+1))])
                                            {
                                                if(tglob +dt < flowstates[0][48]->data[r][std::max(0,std::min(flowstates[0][48]->nrCols()-1,c-1))])
                                                {
                                                    continue;
                                                    //or if any neighbor requires an upgrade?
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                            SetMinDT3Casc(tglob, _dt,dt,flowstates, r, c, 0, dt_min, dtleft,ref_sp);
                    }
                    }

                }
            });
            GetTPool().waitext();
        }


            GetTPool().parallel_for(0,DEM->nrRows() ,[&](int j){
                //for(int j = start; j < end; j++)
                {
                    int r = j;
                    for(int c = 0; c < flowstates[0][4]->nrCols(); c++)
                    {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        if(flowstates[0][21] != nullptr)
                        {
                            if(flowstates[0][47]->data[r][c] < 0.5f)
                            {
                                continue;
                            }
                        }

                        //flowstates[0][7]->data[r][c] += rain_fac * flowstates[0][10]->data[r][c] * dt/(3600000.0);
                        DoMassMomentum3Casc(tglob, courant,dt,dx, flowstates,r,c,0, do_ocean, ocean_level, rain_fac);
                    }
                    }

                }
            });
            GetTPool().waitext();

        {
                auto f = [&flowstates](int l, int r, int c, auto &f){

                    if(r >= flowstates[l][0]->nrRows() || c >= flowstates[l][0]->nrCols())
                    {
                        return;
                    }
                    if(l == 0 && flowstates[l][17] != nullptr)
                    {
                        if(flowstates[l][17]->data[r][c] > 0.5f)
                        {
                            flowstates[l][4]->Drc = flowstates[l][7]->Drc; //H
                            flowstates[l][5]->Drc = flowstates[l][8]->Drc; //Ux
                            flowstates[l][6]->Drc = flowstates[l][9]->Drc; //Uy
                            if(flowstates[l][21] != nullptr)
                            {
                                flowstates[l][14]->data[r][c] = flowstates[l][35]->data[r][c];
                                flowstates[l][15]->data[r][c] = flowstates[l][36]->data[r][c];
                                flowstates[l][16]->data[r][c] = flowstates[l][37]->data[r][c];

                                flowstates[l][12]->data[r][c] = flowstates[l][33]->data[r][c];
                                flowstates[l][34]->data[r][c] = flowstates[l][46]->data[r][c];
                            }
                            return;
                        }
                    }
                    if(flowstates[l][2]->data[r][c] > 0.5 && l < flowstates.size()-1)
                    {
                        for(int n = 0; n < 3; n++)
                        {
                            for(int k = 0; k < 3; k++)
                            {
                                f(l+1,r * 3 + n,c *3 + k,f);
                            }
                        }

                    }else
                    {
                        flowstates[l][4]->Drc = flowstates[l][7]->Drc; //H
                        flowstates[l][5]->Drc = flowstates[l][8]->Drc; //Ux
                        flowstates[l][6]->Drc = flowstates[l][9]->Drc; //Uy
                        if(flowstates[l][21] != nullptr)
                        {
                            flowstates[l][14]->data[r][c] = flowstates[l][35]->data[r][c];
                            flowstates[l][15]->data[r][c] = flowstates[l][36]->data[r][c];
                            flowstates[l][16]->data[r][c] = flowstates[l][37]->data[r][c];
                            flowstates[l][12]->data[r][c] = flowstates[l][33]->data[r][c];
                            flowstates[l][34]->data[r][c] = flowstates[l][46]->data[r][c];
                        }
                    }
                };

                GetTPool().parallel_for(0,flowstates[0][2]->nrRows(),[&](int j){
                    //calling a recusrive lambda within our lambda, things are getting interesting
                    int r = j;
                    for(int c = 0; c < flowstates[0][4]->nrCols(); c++)
                    {
                        if(!pcr::isMV(DEM->data[r][c]))
                        {
                            if(flowstates[0][21] != nullptr)
                            {
                                if(flowstates[0][47]->data[r][c] < 0.5f)
                                {
                                    continue;
                                }
                            }

                            f(0,r,c,f);
                        }
                    }
                });
                GetTPool().waitext();
        }

        //wait
        t = t+dt;
        tglob = tglob + dt;
        dtlast = dt;
        iter ++;

        //FF::cout << " t: " << t << " " << dt << std::endl;


    }
    FF::cout << "threadpool flow calc " << iter <<std::endl;


    return;
}


inline void flow_1D2DKinematic(double _dt,float dx, float dx_chan, std::vector<std::vector<cTMapf *>> &flowstates, cTMapf * CHAN_LDD, cTMapf * CHAN_W, cTMapf * CHAN_D, cTMapf * CHAN_S, cTMapf * CHAN_H, cTMapf * CHAN_HN, cTMapf * CHAN_ALPHA, cTMapf * CHAN_QN, cTMapf * CHAN_INFIL,cTMapf * QX1 = nullptr, cTMapf * QX2 = nullptr,cTMapf * QY1 = nullptr,cTMapf * QY2 = nullptr, float courant = 0.1, float mindt = 0.0, float courant_forced = 0.25f)
{
    //CHANNEL MAPS ARE DEFINED in R3 SPACE!

 //FF::cout << "Flow Kinematic 1D2d" << std::endl;
    //do channel in/out flow
    for(int r = 0; r < CHAN_LDD->nrRows(); r++)
    {
        for(int c = 0; c < CHAN_LDD->nrCols(); c++)
        {
            if(!pcr::isMV(CHAN_LDD->data[r][c]) && CHAN_W->data[r][c] > 1e-12f && CHAN_D->data[r][c] > 1e-12f)
            {

                float h_surf_av = 0.0;
                float level = flowstates.size()-2;
                float n_pix = 0.0;
                float dx_h = dx_chan;

                    float frac = 1.0;
                //convert to water height in relevant super/subgrid

                //check r3 grid
                bool do_r1 = false;
                bool do_r2 = false;
                bool do_rs = false;
                int l_rs = 0;
                int r_rs = 0;
                int c_rs = 0;

               //FF::cout<<"channel check r1"<< std::endl;

                if(flowstates.size() >= 2)
                {
                    int l = flowstates.size()-2;
                    if(flowstates[l][2]->data[r][c] < 0.5 && flowstates.size() == 2)
                    {
                        do_r2 = true;
                        if(!pcr::isMV(flowstates[l][0]->data[r][c]))
                        {
                            h_surf_av = flowstates[l][4]->data[r][c];
                            n_pix += 1.0;

                        }
                    }else if(flowstates[l][2]->data[r][c] > 0.5)
                    {
                        do_r1 = true;
                    }else{
                        int cs = c/3;
                        int rs = r/3;
                        if(cs > -1 && rs > -1 && cs < flowstates[l-1][2]->Cols() && rs < flowstates[l-1][2]->Rows())
                        {
                            if(flowstates[l-1][2]->data[rs][cs] < 0.5)
                            {
                                do_rs = true;

                            }else{
                                do_r2 = true;
                                if(!pcr::isMV(flowstates[l][0]->data[r][c]))
                                {
                                    h_surf_av = flowstates[l][4]->data[r][c];
                                    n_pix += 1.0;

                                }
                            }
                        }else{
                            do_r2 = true;
                            if(!pcr::isMV(flowstates[l][0]->data[r][c]))
                            {
                                h_surf_av = flowstates[l][4]->data[r][c];
                                n_pix += 1.0;

                            }
                        }

                    }
                }else{
                    do_r1 = true;
                }

                //FF::cout<<"channel check r2"<< std::endl;
                if(do_r1)
                {
                    int l = flowstates.size()-1;
                    dx_h = dx_h/3.0;
                    for(int r2 = 0; r2 < 3; r2++)
                    {
                        for(int c2 = 0; c2 < 3; c2++)
                        {
                            int r1 = r *3 + r2;
                            int c1 = c *3 + c2;
                            if((c1 > -1 && r1 > -1 && c1 < flowstates[l][0]->Cols() && r1 < flowstates[l][0]->Rows()))
                            {
                                 if(!pcr::isMV(flowstates[l][0]->data[r1][c1]))
                                {
                                    h_surf_av += flowstates[l][4]->data[r1][c1];
                                    n_pix += 1.0;
                                }
                             }
                        }
                    }

                }else if(!do_r2){

                //FF::cout<<"channel check rs"<< std::endl;
                    int l = flowstates.size()-3;
                    int r_3 =r;
                    int c_3 =c;
                    while(l > -1){
                        r_3 = r_3/3;
                        c_3 = c_3/3;
                        frac = frac/9.0;
                        int r_33 = r_3/3;
                        int c_33 = c_3/3;
                        dx_h = dx_h*3.0;
                        //FF::cout << " try " << l <<std::endl;

                        if(!(c_3 > -1 && r_3 > -1 && c_3 < flowstates[l][2]->Cols() && r_3 < flowstates[l][2]->Rows()))
                        {
                            //FF::cout << " aob1 " << l << " " << c_3 << " " << flowstates[l][2]->Cols()  << " " << r_3 << " " << flowstates[l][2]->Rows() << " "  << std::endl;
                            break;
                        }

                        if(l-1 >= 0)
                        {
                            if(!(c_33 > -1 && r_33 > -1 && c_33 < flowstates[l-1][2]->Cols() && r_33 < flowstates[l-1][2]->Rows()))
                            {
                                //FF::cout << " aob " << l <<std::endl;
                                break;
                            }

                            if(flowstates[l-1][2]->data[r_33][c_33] > 0.5)
                            {
                                //FF::cout << "parent is subbed " <<std::endl;
                                 l_rs = l;
                                do_rs = true;
                                h_surf_av += frac * flowstates[l][4]->data[r_3][c_3];
                                n_pix += frac;
                                r_rs = r_3;
                                c_rs = c_3;
                                break;
                            }
                        }else
                        {
                            if(!pcr::isMV(flowstates[l][0]->data[r_3][c_3]))
                            {
                                //FF::cout << "no parent " <<std::endl;
                                l_rs = l;
                                do_rs = true;
                                h_surf_av += frac * flowstates[l][4]->data[r_3][c_3];
                                n_pix += frac;
                                r_rs = r_3;
                                c_rs = c_3;
                                break;
                            }
                        }

                        l--;
                    }
                }

                h_surf_av = h_surf_av/std::max(1e-12f,n_pix);

                //FF::cout << r << "  "  << c  << " " << do_r1 << " " << do_r2 << " " << do_rs << "    " << l_rs << " " <<  h_surf_av << " " << n_pix << std::endl;
                //CHAN_H->data[r][c] = h_surf_av;

                //FF::cout<<"channel check rs done"<< std::endl;
                if(n_pix > 1e-12)
                {
                    //get current volume in river
                    float vol_chan = CHAN_W->data[r][c] * CHAN_H->data[r][c] * dx_chan;

                    //get current space left over/short
                    float vol_excess = std::max(0.0f,CHAN_W->data[r][c] * (CHAN_H->data[r][c]-CHAN_D->data[r][c]) * dx_chan);
                    float vol_room = std::max(0.0f,CHAN_W->data[r][c] * (CHAN_D->data[r][c]-CHAN_H->data[r][c]) * dx_chan);


                    //FF::cout << " channel " << r << " " << c << " " << n_pix << " " << h_surf_av << " " << vol_chan << " " << vol_excess << " " << vol_room << std::endl;
                    if(vol_room > 1e-12f)
                    {
                        float vol = h_surf_av* n_pix * dx_h*dx_h;
                        //correct for travel distance
                        float vel = (std::pow(h_surf_av,2.0/3.0) * std::sqrt(std::max(0.001,CHAN_S->data[r][c]*0.5 + h_surf_av/dx_chan))/0.03);
                        float frac_max = vol_room/std::max(1e-12f,vol);
                        float frac_in = std::min(1.0f,std::max(0.0f,std::min(frac_max,(float)(vel*_dt/dx_h))));

                        //
                        //FF::cout << "flow in " << frac_in << " " << vel << "  " << CHAN_D->data[r][c] << " " << CHAN_S->data[r][c] << " " << do_r1 << " " << do_r2 << " " << do_rs << std::endl;

                        if(do_r1)
                        {
                            int l = flowstates.size()-1;
                            for(int r2 = 0; r2 < 3; r2++)
                            {
                                for(int c2 = 0; c2 < 3; c2++)
                                {
                                    int r1 = r *3 + r2;
                                    int c1 = c *3 + c2;
                                    if(r1 > -1 && c1 > -1 && r1 < flowstates[l][4]->Rows() && c1 < flowstates[l][4]->Cols())
                                    {
                                        if(!pcr::isMV(flowstates[l][0]->data[r1][c1]))
                                        {
                                            flowstates[l][4]->data[r1][c1] *= (1.0 - frac_in);
                                        }
                                    }
                                }
                            }
                        }else if(do_r2)
                        {
                            int l = flowstates.size()-2;
                            flowstates[l][4]->data[r][c] *= (1.0 - frac_in);
                        }else if(do_rs)
                        {
                            int l = l_rs;
                            flowstates[l][4]->data[r_rs][c_rs] *= ((frac*(1.0 - frac_in)) + (1.0-frac) *1.0);
                        }

                        CHAN_H->data[r][c] += (frac_in * vol)/std::max(1e-12f,CHAN_W->data[r][c] *dx_chan);
                        CHAN_H->data[r][c] = std::min(CHAN_H->data[r][c],CHAN_D->data[r][c]);

                    }
                    if(vol_excess > 1e-12f)
                    {
                        //transfer to cells
                        float height_add = vol_excess/std::max(1e-12f,n_pix * dx_h*dx_h);
                        if(do_r1)
                        {
                            int l = flowstates.size()-1;
                            for(int r2 = 0; r2 < 3; r2++)
                            {
                                for(int c2 = 0; c2 < 3; c2++)
                                {
                                    int r1 = r *3 + r2;
                                    int c1 = c *3 + c2;
                                    if(r1 > -1 && c1 > -1 && r1 < flowstates[l][4]->Rows() && c1 < flowstates[l][4]->Cols())
                                    {
                                        if(!pcr::isMV(flowstates[l][0]->data[r1][c1]))
                                        {
                                            flowstates[l][4]->data[r1][c1] += height_add;
                                        }
                                    }
                                }
                            }
                        }else if(do_r2)
                        {
                            int l = flowstates.size()-2;
                            flowstates[l][4]->data[r][c] += height_add;
                        }else if(do_rs)
                        {
                            int l = l_rs;
                            flowstates[l][4]->data[r_rs][c_rs] += frac * height_add;
                        }

                        //update channel height
                        CHAN_H->data[r][c] = std::min(CHAN_H->data[r][c],CHAN_D->data[r][c]);
                    }

                    //update channel height
                }
            }
        }
    }


    //FF::cout << "Flow Kinematic 1D2d done" << std::endl;
    //Channel Kinematic Flow
    /*writeRaster(*CHAN_H,"chha.tif");
    writeRaster(*CHAN_INFIL,"infil.tif");
    writeRaster(*CHAN_ALPHA,"alpha.tif");
    writeRaster(*CHAN_LDD,"chldd.tif");
    writeRaster(*CHAN_W,"chw.tif");
    FF::cout << "Flow Kinematic 1D2d done, timestepping channel with dt " << _dt << " dx " << dx_chan << std::endl;
    */

    FlowKinematic(CHAN_LDD, CHAN_S, CHAN_W, 0.03, CHAN_H, CHAN_HN, _dt,dx_chan,CHAN_ALPHA, CHAN_QN, CHAN_INFIL);

    //writeRaster(*CHAN_HN,"chhn.tif");
    //writeRaster(*CHAN_ALPHA,"alphan.tif");
    for(int r = 0; r < CHAN_LDD->nrRows(); r++)
    {
        for(int c = 0; c < CHAN_LDD->nrCols(); c++)
        {
            CHAN_H->data[r][c] = CHAN_HN->data[r][c];
        }
    }


}




#endif // FULLFLOOD_SVSOLVE_H
