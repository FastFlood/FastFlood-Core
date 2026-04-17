#ifndef FULLFLOOD_SIMPLEGRID_H
#define FULLFLOOD_SIMPLEGRID_H

static inline void flow_saintvenantdtandflow(double _dt,float dx, cTMapf * DEM,cTMapf * N,cTMapf * H,cTMapf * VX,cTMapf * VY, cTMapf * HN,cTMapf * VXN,cTMapf * VYN, cTMapf * QX1 = nullptr, cTMapf * QX2 = nullptr,cTMapf * QY1 = nullptr,cTMapf * QY2 = nullptr, float courant = 0.1, float mindt = 0.0, float courant_forced = 0.25f)
{



    int dim0 = DEM->nrCols();
    int dim1 = DEM->nrRows();


    double dt = 1e6;
    float dt_min = 1e6;
    float dt_store = _dt;
    bool stop = false;
    #pragma omp parallel shared(stop,dt_min,dt)
    {
        double t = 0;
        double t_end = _dt;

        int iter = 0;
        while(true)
        {

            #pragma omp single
            {
                dt_min = 1e6;
                dt = 1e6;
            }

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
             //get dt
            #pragma omp for collapse(2) reduction(min:dt_min,dt)
            for(int r = 0; r < DEM->nrRows(); r++)
            {
                for (int c = 0; c < DEM->nrCols(); c++)
                {

                    if(!pcr::isMV(DEM->data[r][c]))
                    {

                        int x = r;
                        int y = c;

                        const int gy = std::min(dim1-(int)(1),std::max((int)(0),(int)(x)));
                        const int gx = std::min(dim0-(int)(1),std::max((int)(0),(int)(y)));


                        const int gx_x1 = std::min(dim0-(int)(1),std::max((int)(0),(int)(gx - 1 )));
                        const int gy_x1 = gy;
                        const int gx_x2 = std::min(dim0-(int)(1),std::max((int)(0),(int)(gx + 1)));
                        const int gy_x2 = gy;
                        const int gx_y1 = gx;
                        const int gy_y1 = std::min(dim1-(int)(1),std::max((int)(0),(int)(gy - 1)));
                        const int gx_y2 = gx;
                        const int gy_y2 = std::min(dim1-(int)(1),std::max((int)(0),(int)(gy + 1)));


                        float z = DEM->data[gy][gx];

                        float z_x1 = DEM->data[gy][gx_x1];
                        float z_x2 = DEM->data[gy][gx_x2];
                        float z_y1 = DEM->data[gy_y1][gx];
                        float z_y2 = DEM->data[gy_y2][gx];


                        if(gx + 1 > dim0-1)
                        {
                            z_x2 = -99999;
                        }
                        if(gx - 1 < 0)
                        {
                            z_x1 = -99999;
                        }
                        if(gy + 1 > dim1-1)
                        {
                            z_y2 = -99999;
                        }
                        if(gy - 1 < 0)
                        {
                            z_y1 = -99999;
                        }
                        float zc_x1 = z_x1 < -99995? z : z_x1;
                        float zc_x2 = z_x2 < -99995? z : z_x2;
                        float zc_y1 = z_y1 < -99995? z : z_y1;
                        float zc_y2 = z_y2 < -99995? z : z_y2;

                        float h = std::max(0.0f,H->data[gy][gx]);
                        float h_x1 = std::max(0.0f,H->data[gy][gx_x1]);
                        float h_x2 = std::max(0.0f,H->data[gy][gx_x2]);
                        float h_y1 = std::max(0.0f,H->data[gy_y1][gx]);
                        float h_y2 = std::max(0.0f,H->data[gy_y2][gx]);

                        float hw_x1 = h > h_x1? std::max(0.0f,std::min(h,h+z - (h_x1 +z_x1))) : 0.0;
                        float hw_x2 = h > h_x2? std::max(0.0f,std::min(h,h+z - (h_x2 +z_x2))) : 0.0;
                        float hw_y1 = h > h_y1? std::max(0.0f,std::min(h,h+z - (h_y1 +z_y1))) : 0.0;
                        float hw_y2 = h > h_y2? std::max(0.0f,std::min(h,h+z - (h_y2 +z_y2))) : 0.0;

                        float vmax = 1000.0;

                        float vx = std::max(std::sqrt(9.81f * std::fabs(hw_x2)),std::max(std::sqrt(9.81f * std::fabs(hw_x1)),std::min(vmax,std::max(-vmax,VX->data[gy][gx]))));
                        float vy = std::max(std::sqrt(9.81f * std::fabs(hw_y2)),std::max(std::sqrt(9.81f * std::fabs(hw_y1)),std::min(vmax,std::max(-vmax,VY->data[gy][gx]))));

                        double dt_req = courant *dx/( std::min(1000.0f,std::max(0.01f,(sqrt(vx*vx + vy * vy)))));
                        dt_min = std::min(_dt-t,std::max((double)mindt,std::min((double)dt_min,dt_req)));
                        dt = std::min(_dt-t,std::max((double)mindt,std::min((double)dt_min,dt)));
                    }

                }
            }

            #pragma omp critical
            {
                //dt = dt_min;
                //dt_min = std::min(dt_min,(float)dt_min_this);
                //dt_store = std::min(t_end - t,std::max((double)dt_min,1e-6));


                //reset already for next time
                //dt_min = 1e6;
            }

            //wait
            #pragma omp barrier

            //reduce for minimum dt
            //#pragma omp single
            {
               //dt = dt_store;
            }


            //wait
            #pragma omp barrier


        #pragma omp single
        {
            //FF::cout << "dt  " << dt << std::endl;
        }

            //do flow calc


            float tx = dt/dx;

            #pragma omp for collapse(2)
            for(int r = 0; r < DEM->nrRows(); r++)
            {
                for (int c = 0; c < DEM->nrCols(); c++)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        int x = r;
                        int y = c;
                        const int gy = std::min(dim1-(int)(1),std::max((int)(0),(int)(x)));
                        const int gx = std::min(dim0-(int)(1),std::max((int)(0),(int)(y)));

                        const int gx_x1 = std::min(dim0-(int)(1),std::max((int)(0),(int)(gx - 1 )));
                        const int gy_x1 = gy;
                        const int gx_x2 = std::min(dim0-(int)(1),std::max((int)(0),(int)(gx + 1)));
                        const int gy_x2 = gy;
                        const int gx_y1 = gx;
                        const int gy_y1 = std::min(dim1-(int)(1),std::max((int)(0),(int)(gy - 1)));
                        const int gx_y2 = gx;
                        const int gy_y2 = std::min(dim1-(int)(1),std::max((int)(0),(int)(gy + 1)));

                        float z = DEM->data[gy][gx];
                        float n = std::max(0.00001f,N->data[gy][gx]);

                        float z_x1 = DEM->data[gy][gx_x1];
                        float z_x2 = DEM->data[gy][gx_x2];
                        float z_y1 = DEM->data[gy_y1][gx];
                        float z_y2 = DEM->data[gy_y2][gx];


                        if(gx + 1 > dim0-1)
                        {
                            z_x2 = -99999;
                        }
                        if(gx - 1 < 0)
                        {
                            z_x1 = -99999;
                        }
                        if(gy + 1 > dim1-1)
                        {
                            z_y2 = -99999;
                        }
                        if(gy - 1 < 0)
                        {
                            z_y1 = -99999;
                        }
                        float zc_x1 = z_x1 < -99995? z : z_x1;
                        float zc_x2 = z_x2 < -99995? z : z_x2;
                        float zc_y1 = z_y1 < -99995? z : z_y1;
                        float zc_y2 = z_y2 < -99995? z : z_y2;

                        float h = std::max(0.0f,H->data[gy][gx]);
                        float h_x1 = std::max(0.0f,H->data[gy][gx_x1]);
                        float h_x2 = std::max(0.0f,H->data[gy][gx_x2]);
                        float h_y1 = std::max(0.0f,H->data[gy_y1][gx]);
                        float h_y2 = std::max(0.0f,H->data[gy_y2][gx]);

                        float h_corr_x1 =  std::max(0.0f,(h_x1 - std::max(0.0f,z - zc_x1)));
                        float h_corr_x1b =  std::max(0.0f,(h - std::max(0.0f,zc_x1-z)));
                        float h_corr_x2 =  std::max(0.0f,(h - std::max(0.0f,zc_x2 - z)));
                        float h_corr_x2b =  std::max(0.0f,(h_x2 - std::max(0.0f,z - zc_x2)));
                        float h_corr_y1 =  std::max(0.0f,(h_y1 - std::max(0.0f,z - zc_y1)));
                        float h_corr_y1b =  std::max(0.0f,(h - std::max(0.0f,zc_y1 - z)));
                        float h_corr_y2 =  std::max(0.0f,(h - std::max(0.0f,zc_y2 - z)));
                        float h_corr_y2b =  std::max(0.0f,(h_y2 - std::max(0.0f,z - zc_y2)));

                        float vmax = 0.1 * dx/dt;

                        float vx = std::min(vmax,std::max(-vmax,VX->data[gy][gx]));
                        float vy = std::min(vmax,std::max(-vmax,VY->data[gy][gx]));

                        float vx_x1 = std::min(vmax,std::max(-vmax,VX->data[gy][gx_x1]));
                        float vy_x1 = std::min(vmax,std::max(-vmax,VY->data[gy][gx_x1]));
                        float vx_x2 = std::min(vmax,std::max(-vmax,VX->data[gy][gx_x2]));
                        float vy_x2 = std::min(vmax,std::max(-vmax,VY->data[gy][gx_x2]));
                        float vx_y1 = std::min(vmax,std::max(-vmax,VX->data[gy_y1][gx]));
                        float vy_y1 = std::min(vmax,std::max(-vmax,VY->data[gy_y1][gx]));
                        float vx_y2 = std::min(vmax,std::max(-vmax,VX->data[gy_y2][gx]));
                        float vy_y2 = std::min(vmax,std::max(-vmax,VY->data[gy_y2][gx]));

                        float zh = z + h;
                        float zh_x1 =std::min(1.0f,h_x1/std::max(0.0001f,h)) * (zc_x1 + h_x1) + (1.0f-std::min(1.0f,h_x1/std::max(0.0001f,h)))*std::min(zh,zc_x1 + h_x1);
                        float zh_x2 =std::min(1.0f,h_x2/std::max(0.0001f,h)) * (zc_x2 + h_x2) + (1.0f-std::min(1.0f,h_x2/std::max(0.0001f,h)))*std::min(zh,zc_x2 + h_x2);
                        float zh_y1 =std::min(1.0f,h_y1/std::max(0.0001f,h)) * (zc_y1 + h_y1) + (1.0f-std::min(1.0f,h_y1/std::max(0.0001f,h)))*std::min(zh,zc_y1 + h_y1);
                        float zh_y2 =std::min(1.0f,h_y2/std::max(0.0001f,h)) * (zc_y2 + h_y2) + (1.0f-std::min(1.0f,h_y2/std::max(0.0001f,h)))*std::min(zh,zc_y2 + h_y2);

                        float sx_zh_x2 = std::min((float)(0.5),std::max((float)(-0.5),(float)((zh_x2-zh)/dx)));
                        float sy_zh_y1 = std::min((float)(0.5),std::max((float)(-0.5),(float)((zh-zh_y1)/dx)));
                        float sx_zh_x1 = std::min((float)(0.5),std::max((float)(-0.5),(float)((zh-zh_x1)/dx)));
                        float sy_zh_y2 = std::min((float)(0.5),std::max((float)(-0.5),(float)((zh_y2-zh)/dx)));

                        float sx_zh = std::min(1.0f,std::max(-1.0f,0.5f * (sx_zh_x1 + sx_zh_x2)));
                        float sy_zh = std::min(1.0f,std::max(-1.0f,0.5f * (sy_zh_y1 + sy_zh_y2)));

                        LSMVector3 hll_x1 = z_x1 < -99995? LSMVector3(0.0,0.0,0.0):F_HLL2FL(h_corr_x1,vx_x1,vy_x1,h_corr_x1b,vx,vy);
                        LSMVector3 hll_x2 = z_x2 < -99995? LSMVector3(0.0,0.0,0.0):F_HLL2FL(h_corr_x2,vx,vy,h_corr_x2b,vx_x2,vy_x2);
                        LSMVector3 hll_y1 = z_y1 < -99995? LSMVector3(0.0,0.0,0.0):F_HLL2FL(h_corr_y1,vy_y1,vx_y1,h_corr_y1b,vy,vx);
                        LSMVector3 hll_y2 = z_y2 < -99995? LSMVector3(0.0,0.0,0.0):F_HLL2FL(h_corr_y2,vy,vx,h_corr_y2b,vy_y2,vx_y2);

                        float C = courant_forced;

                        float flux_x1 = z_x1 < -99995? std::max(-h * C,(float)(-h * sqrt(h) *dt*  sqrt(h)/(dx*(0.001+n)))):std::max(-h * C,std::min((float)(+tx*(hll_x1.x)),h_x1 * C));//std::max(-h * C,(float)(-h * sqrt(h) *0.0*dt*  sqrt(h)/(dx*(0.001+n))))
                        float flux_x2 = z_x2 < -99995? std::max(-h * C,(float)(-h * sqrt(h) *dt*  sqrt(h)/(dx*(0.001+n)))):std::max(-h * C,std::min((float)(-tx*(hll_x2.x)),h_x2 * C));//std::max(-h * C,(float)(-h * sqrt(h) *0.0*dt*  sqrt(h)/(dx*(0.001+n))))
                        float flux_y1 = z_y1 < -99995? std::max(-h * C,(float)(-h * sqrt(h) *dt*  sqrt(h)/(dx*(0.001+n)))):std::max(-h * C,std::min((float)(+tx*(hll_y1.x)),h_y1 * C));//std::max(-h * C,(float)(-h * sqrt(h) *0.0*dt*  sqrt(h)/(dx*(0.001+n))))
                        float flux_y2 = z_y2 < -99995? std::max(-h * C,(float)(-h * sqrt(h) *dt*  sqrt(h)/(dx*(0.001+n)))):std::max(-h * C,std::min((float)(-tx*(hll_y2.x)),h_y2 * C));//std::max(-h * C,(float)(-h * sqrt(h) *0.0*dt*  sqrt(h)/(dx*(0.001+n))))

                        int edges = ((z_x1 < -99995)?1:0) +((z_x2 < -99995)?1:0)+((z_y1 < -99995)?1:0)+((z_y2 < -99995)?1:0);

                        bool edge = (z_x1 < -99995 || z_x2 < -99995 || z_y1 < -99995 || z_y2 < -99995);
                        float hold = h;
                        float hn = ((std::max(0.00f,(float)(h + flux_x1 + flux_x2 + flux_y1 + flux_y2))));


                        float f_centre_x = 0.5 * GRAV*((z<zc_x2?std::min(zc_x2,z + h) : std::min(z,zc_x2 + h_x2)) - (z<zc_x1?std::min(zc_x1,z + h) : std::min(z,zc_x1 + h_x1)))*(h_corr_x1 + h_corr_x2b);
                        float f_centre_y = 0.5 * GRAV*((z<zc_y2?std::min(zc_y2,z + h) : std::min(z,zc_y2 + h_y2)) - (z<zc_y1?std::min(zc_y1,z + h) : std::min(z,zc_y1 + h_y1)))*(h_corr_y1 + h_corr_y2b);


                        if(h > 0.01f )
                        {
                            //FF::cout << "h: " << h << " - " << sx_zh << " " << sy_zh << " " << h_x1 << " " << h_x2 << " " << h_corr_x1 << " " << h_corr_x2 << " " << h_corr_x1b << " " << h_corr_x2b << " " << hll_x1.y << "  " << hll_x2.y  <<  " " << zc_x1 << " " << zc_x2 << std::endl;
                            //FF::cout << "f: " << tx*(hll_x2.y - hll_x1.y + f_centre_x) << " " << hll_x2.y - hll_x1.y << " " <<  f_centre_x << " " << hll_y2.z - hll_y1.z << std::endl;
                        }

                        float qxn = h * vx - tx*(hll_x2.y - hll_x1.y + f_centre_x) - tx*(hll_y2.z - hll_y1.z);//- 0.5 * GRAV *hn*sx_zh * dt;
                        float qyn = h * vy - tx*(hll_x2.z - hll_x1.z) - tx*(hll_y2.y - hll_y1.y + f_centre_y);//- 0.5 * GRAV *hn*sy_zh * dt;

                        float vsq = sqrt((float)(vx * vx + vy * vy));
                        float nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01f,pow((float)(hn),(float)(4.0/3.0)));
                        float nsq = nsq1*vsq*dt;

                        float vxn = (float)((qxn/(1.0f+nsq)))/std::max(0.01f,(float)(hn));
                        float vyn = (float)((qyn/(1.0f+nsq)))/std::max(0.01f,(float)(hn));

                        float threshold = std::min(0.1f,0.01f * dx);
                        if(hn < threshold)
                        {
                            float kinfac = std::max(0.0f,(threshold - hn) / (0.025f * dx));
                            float acc_eff = (vxn -vx)/std::max(0.0001,dt);

                            float v_kin = (sx_zh>0.0f?1.0f:-1.0f) * hn * sqrt(hn) * sqrt(std::max(0.0f,sx_zh>0.0f?sx_zh:-sx_zh))/(0.001f+n);

                            vxn = kinfac * v_kin + vxn*(1.0f-kinfac);

                        }

                        if(hn < threshold)
                        {
                            float kinfac = std::max(0.0f,(threshold - hn) / (0.025f * dx));
                            float acc_eff = (vyn -vy)/std::max(0.0001,dt);

                            float v_kin = (sy_zh>0.0f?1.0f:-1.0f) * hn * sqrt(hn) * sqrt(std::max(0.0f,sy_zh>0.0f?sy_zh:-sy_zh))/(0.001f+n);

                            vyn = kinfac * v_kin + vyn*(1.0f-kinfac);

                        }

                        hn = edges > 2? 0.0f:(isnan(hn)? 0.0:hn);
                        vxn = edges > 2? 0.0f:isnan(vxn)? 0.0:vxn;
                        vyn = edges > 2? 0.0f:isnan(vyn)? 0.0:vyn;


                        //write output
                        HN->data[gy][gx] = hn;
                        VXN->data[gy][gx] = vxn;
                        VYN->data[gy][gx] = vyn;

                        if(QX1 != nullptr && QX2 != nullptr && QY1 != nullptr && QY2 != nullptr)
                        {

                            QX1->data[gy][gx] = flux_x1 *(dx*dx);
                            QX2->data[gy][gx] = flux_x2 *(dx*dx);
                            QY1->data[gy][gx] = flux_y1 *(dx*dx);
                            QY2->data[gy][gx] = flux_y2 *(dx*dx);
                        }
                    }
                }
            }

            //wait
            #pragma omp barrier

            #pragma omp for collapse(2)
            for(int r = 0; r < H->nrRows(); r++)
            {
                for(int c = 0; c < H->nrCols(); c++)
                {
                    H->Drc = HN->Drc;
                    VX->Drc = VXN->Drc;
                    VY->Drc = VYN->Drc;


                }
            }

            //wait
            #pragma omp barrier

            t = t+dt;
            iter ++;


        }

        #pragma omp master
        {

        }

    }
    return;
}
#endif // FULLFLOOD_SIMPLEGRID_H
