#ifndef FULLFLOOD_MMSOLVE_H
#define FULLFLOOD_MMSOLVE_H

/*
 * @brief
 * This function does the actual HLL solver for momentum and mass transfer terms
 * but specifically only for one gridcell boundary. Give it the current cell properties, and
 * the neighbor to look at, and it will fill the hll vector thats appropriate (depending on ca and ca2, the directionality indicators)
 * these contain the x and y momentum transfers.
 * h_add will be incremented for the new height change rate.
 * it will also fill the center balancing term fx_center and fy_center,
 * and it will fill the slope, to track the average slope in a direction (sx_av and sy_av)
 * NOTE: this function is carefully crafted so that it gives identical results (but opposite),
 * when called from pixel r1,c1 to r,c as opposed to from r,c to r1,c1.
 * this means its strictly mass conserving, and we can parallalize our scheme, although it requires
 * we do all calculations for each cell edge double, this still wins out in the end
 *
 *
 *
 */

inline void MomentumFluxGet(double tglob, int l, int r, int c, float _dt, float dx, std::vector<std::vector<cTMapf *>> &flowstates,
                            bool do_ocean, float ocean_level,float h,float vx,float vy,float z,float n,
                            float &fx_center /*out*/,float &fy_center /*out*/,
                            float &h_add/*out*/, LSMVector3 & hll_x1 /*out*/, LSMVector3 & hll_x2 /*out*/, LSMVector3 & hll_y1 /*out*/, LSMVector3 & hll_y2 /*out*/,
                            float & sx_av /*out*/, float & sy_av /*out*/,

                            int r2, int c2, int l2, float dy, int ca, int ca2)
{



    //is the neighboring cell dormatn due to sdt, then use the

    float dxh = 0.5*(std::fabs(dy) + std::fabs(dx));
    float dx_w = std::min(std::fabs(dy),std::fabs(dx));
    //do mass and momentum calculations
    int dim0 = flowstates[l2][0]->Cols();
    int dim1 = flowstates[l2][0]->Rows();
    int x = c2;
    int y = r2;
    const int gy_x1 = std::min(dim1-(int)(1),std::max((int)(0),(int)(y)));
    const int gx_x1 = std::min(dim0-(int)(1),std::max((int)(0),(int)(x)));

    //does this neighbor have a timestep?
    if(flowstates[l][21] != nullptr)
    {
        if(!(r2 +1 >= dim1 || r2 < 0 || c2 +1 >= dim0 || c2 < 0))
        {
            //if not, use the cache here
            float dt_other = flowstates[l2][13]->data[r2][c2];
            int r_base =r2;
            int c_base =c2;
            if(flowstates[0][17] != nullptr)
            {
                int l_c = l2;
                while(l_c > 0)
                {
                    r_base /= 3;
                    c_base /= 3;
                    l_c -= 1;
                }
                dt_other = flowstates[0][13]->data[r_base][c_base];
            }
            if(dt_other < 1e-12f)
            {
                int lt=l2;
                int rt=r2;
                int ct =c2;
                float m = -1;
                float haddrep = 0.0;
                int ca2temp = ca2 *(-1);
                if( dx < dy || (dx == dy && (r < r2 || c < c2)))
                {
                    m =1;
                    lt = l;
                    rt = r;
                    ct =c;
                    ca2temp = ca2temp *(-1);
                }
                if(ca==1)
                {
                    double fmfy = 0.0;
                    double fmfz = 0.0;
                  if(ca2temp == 1)
                  {
                     haddrep = flowstates[lt][24]->data[rt][ct];
                    h_add += m*dx_w *flowstates[lt][24]->data[rt][ct]*_dt;
                    sx_av += m*dx_w*flowstates[lt][40]->data[rt][ct];
                    fmfy = flowstates[lt][25]->data[rt][ct]*_dt;
                    fmfz = flowstates[lt][26]->data[rt][ct]*_dt;
                  }else{
                    haddrep = flowstates[lt][21]->data[rt][ct];
                    h_add += m*1.0*dx_w *flowstates[lt][21]->data[rt][ct]*_dt;
                    sx_av += m*dx_w*flowstates[lt][38]->data[rt][ct];
                    fmfy = flowstates[lt][22]->data[rt][ct]*_dt;
                    fmfz = flowstates[lt][23]->data[rt][ct]*_dt;
                  }

                  if(ca2 == 1)
                  {
                      hll_x2.y = hll_x2.y+ m*dx_w * fmfy;
                      hll_x2.z = hll_x2.z+ m*dx_w * fmfz;
                  }else
                  {
                      hll_x1.y = hll_x1.y+ m*dx_w * fmfy;
                      hll_x1.z = hll_x1.z+ m*dx_w * fmfz;
                  }

                }else if(ca == 2)
                {
                    double fmfy = 0.0;
                    double fmfz = 0.0;
                  if(ca2temp == 1)
                  {
                      haddrep = flowstates[lt][30]->data[rt][ct];
                    h_add += m*1.0*dx_w *flowstates[lt][30]->data[rt][ct]*_dt;
                    sy_av += m*dx_w*flowstates[lt][45]->data[rt][ct];
                    fmfy = flowstates[lt][31]->data[rt][ct]*_dt;
                    fmfz = flowstates[lt][32]->data[rt][ct]*_dt;
                  }else{
                      haddrep =flowstates[lt][27]->data[rt][ct];
                    h_add += m*1.0*dx_w *flowstates[lt][27]->data[rt][ct]*_dt;
                    sy_av += m*dx_w*flowstates[lt][43]->data[rt][ct];
                    fmfy = flowstates[lt][28]->data[rt][ct]*_dt;
                    fmfz = flowstates[lt][29]->data[rt][ct]*_dt;
                  }

                  if(ca2 == 1)
                  {
                      hll_y2.y = hll_y2.y+ m*dx_w * fmfy;
                      hll_y2.z = hll_y2.z+ m*dx_w * fmfz;
                  }else
                  {
                      hll_y1.y = hll_y1.y+ m*dx_w * fmfy;
                      hll_y1.z = hll_y1.z+ m*dx_w * fmfz;
                  }
                }
                if((r == 688 && c== 591 && l ==3) || (r2 == 688 && c2== 591 && l2 ==3))
                {
                   // FF::cout << "fsolve_nodt " << r << " " << c << " " << l << " loc2: " << r2 << " " << c2 << " " << l2 << " hadd :" <<  haddrep<< " dx " << dx << " " << dy << " " << dx_w << std::endl;
                }
                return;
            }
        }
    }

    //is our neighbor a normal cell or a cache cell?


    //is the other gridcells base level superpixel a cached one?
    //get scaling and apply it

    int r_base =r2;
    int c_base =c2;
    if(flowstates[0][17] != nullptr)
    {
        int l_c = l2;
        while(l_c > 0)
        {
            r_base /= 3;
            c_base /= 3;
            l_c -= 1;
        }
    }

    float z_x1 = flowstates[l2][0]->data[gy_x1][gx_x1];
    float h_x1 = flowstates[l2][4]->data[gy_x1][gx_x1];
    float vx_x1 = flowstates[l2][5]->data[gy_x1][gx_x1];
    float vy_x1 = flowstates[l2][6]->data[gy_x1][gx_x1];

    float dtlo = _dt;
    if(flowstates[l2][13] != nullptr)
    {
        dtlo = std::max(dtlo,flowstates[l2][13]->data[gy_x1][gx_x1]);
    }

    bool is_cache_x1 = false;
    if(flowstates[0][17] != nullptr)
    {
        if(r_base < flowstates[0][17]->nrRows() && c_base < flowstates[0][17]->nrCols())
        {
            if(flowstates[0][17]->data[r_base][c_base] > 0.5f)
            {
                is_cache_x1 = true;
                //FF::cout << "neihbor is cached : "<< flowstates[0][18]->data[r_base][c_base]/(std::max(1e-12f,flowstates[0][4]->data[r_base][c_base])) << std::endl;
                float h_cache = flowstates[0][4]->data[r_base][c_base];
                float vx_cache = flowstates[0][5]->data[r_base][c_base];
                float vy_cache = flowstates[0][6]->data[r_base][c_base];

                if(flowstates[0][21] != nullptr)
                {
                    float dt_frac_r = std::min(1.0,std::max(0.0,(flowstates[0][12]->data[r_base][c_base]-tglob))/std::max(1e-12f,flowstates[0][34]->data[r_base][c_base]));
                    h_cache = std::max(0.0f,h_cache -flowstates[0][34]->data[r_base][c_base] * flowstates[0][14]->data[r_base][c_base] * dt_frac_r);
                    vx_cache -= flowstates[0][34]->data[r_base][c_base]*flowstates[0][15]->data[r_base][c_base] * dt_frac_r;
                    vy_cache -= flowstates[0][34]->data[r_base][c_base]*flowstates[0][16]->data[r_base][c_base] * dt_frac_r;
                }

                h_x1 = h_x1 * h_cache/(std::max(1e-12f,flowstates[0][18]->data[r_base][c_base]));
                vx_x1 = vx_x1 * flowstates[0][5]->data[r_base][c_base]/(flowstates[0][19]->data[r_base][c_base] > 0.0f? std::max(1e-12f,flowstates[0][19]->data[r_base][c_base]):std::min(-1e-12f,flowstates[0][19]->data[r_base][c_base]));
                vy_x1 = vy_x1 * flowstates[0][6]->data[r_base][c_base]/(flowstates[0][20]->data[r_base][c_base] > 0.0f? std::max(1e-12f,flowstates[0][20]->data[r_base][c_base]):std::min(-1e-12f,flowstates[0][20]->data[r_base][c_base]));
            }
        }
    }
    if(!is_cache_x1)
    {
        //potentiall go back in time a bit for the neigbor:
        if(flowstates[l][21] != nullptr)
        {
            if(r2 > -1 && r2 < flowstates[l2][0]->nrRows() && c2 > -1 && c2 < flowstates[l2][0]->nrCols())
            {
                float dt_frac_r = std::min(1.0,std::max(0.0,(flowstates[l2][12]->data[r2][c2]-tglob))/std::max(1e-12f,flowstates[l2][34]->data[r2][c2]));
                float temp = h_x1;
                h_x1 = std::max(0.0f,h_x1-flowstates[l2][34]->data[r2][c2] * flowstates[l2][14]->data[r2][c2] * dt_frac_r);
                vx_x1 -= flowstates[l2][34]->data[r2][c2]*flowstates[l2][15]->data[r2][c2] * dt_frac_r;
                vy_x1 -= flowstates[l2][34]->data[r2][c2]*flowstates[l2][16]->data[r2][c2] * dt_frac_r;
                if((r == 688 && c== 591 && l ==3) || (r2 == 688 && c2== 591 && l2 ==3))
                {
                    //FF::cout << "h corr x1 horg " << temp << " h_x1 " << h_x1 << " " << dt_frac_r << " dhdt " << flowstates[l2][14]->data[r2][c2] << " dt " <<  flowstates[l2][34]->data[r2][c2] << std::endl;
                }
            }

        }
    }


    if(pcr::isMV(z_x1) || r2 +1 >= dim1 || r2 < 1 || c2 +1 >= dim0 || c2 < 1)
    {
        if(r2 +1 >= dim1 || r2 < 1 || c2 +1 >= dim0 || c2 < 1)
        {

            z_x1 = -9999999999;
            h_x1 = 0.0;
            vx_x1 = 0.0;
            vy_x1 = 0.0;
        }else if(do_ocean)
        {
            z_x1 = z;
            h_x1 = std::max(0.0f,ocean_level - z);
            vx_x1 = 0.0;
            vy_x1 = 0.0;

        }else{

            z_x1 = -9999999999;
            h_x1 = 0.0;
            vx_x1 = 0.0;
            vy_x1 = 0.0;
        }
    }


    float zc_x1 = z_x1 < -99995? z : z_x1;

    float h_corr_x1 =  std::max(0.0f,(h_x1 - std::max(0.0f,z - zc_x1)));
    float h_corr_x1b =  std::max(0.0f,(h - std::max(0.0f,zc_x1-z)));

    float vmax = 0.25 * dxh/std::max(1e-6f,std::max(dtlo,_dt));

    vx_x1 = std::min(vmax,std::max(-vmax,vx_x1));
    vy_x1 = std::min(vmax,std::max(-vmax,vy_x1));
    float zh = z + h;
    float zh_x1 =std::min(1.0f,h_x1/std::max(0.0001f,h)) * (zc_x1 + h_x1) + (1.0f-std::min(1.0f,h_x1/std::max(0.0001f,h)))*std::min(zh,zc_x1 + h_x1);

    float sx_zh_x1 =((float)(ca2))*std::min((float)(0.5),std::max((float)(-0.5),(float)((zh_x1-zh)/dxh)));

    float sx_zh = std::min(1.0f,std::max(-1.0f,sx_zh_x1));

    float vxh = vx;
    float vyh = vy;

    vxh = std::min(vmax,std::max(-vmax,vx));
    vyh = std::min(vmax,std::max(-vmax,vy));


    if(ca == 2)
    {
        float temp = vy_x1;
        vy_x1 = vx_x1;
        vx_x1= temp;

        temp = vyh;
        vyh = vxh;
        vxh = temp;
    }

    LSMVector3 hll_x1h;

    if(ca2 == 1)
    {
        hll_x1h = z_x1 < -99995? LSMVector3(0.0,0.0,0.0):F_HLL2FL(h_corr_x1b,vxh,vyh,h_corr_x1,vx_x1,vy_x1);
    }else{
        hll_x1h = z_x1 < -99995? LSMVector3(0.0,0.0,0.0):F_HLL2FL(h_corr_x1,vx_x1,vy_x1,h_corr_x1b,vxh,vyh);
    }

    //fluxes are positive in direction of

    float C = 0.125;
    //this makes our scheme strictly mass conserving, but requires still a quite low timestep
    //what we can do though is make this bigger if we have homogeneous flow conditions
    //check the in and outflow of the neighbors
    //



    float tx = _dt/(dx);
    float flux_edge = std::max(-h * C,(float)(-h * sqrt(h) *_dt*  sqrt(h)/(dxh*(0.001+n))));
    float flux_x1 = (float)(tx*(- ((float)(ca2))*hll_x1h.x));

    //limiting for strict mass conservation
    float dtmax = std::max(1e-6f,std::max(dtlo,_dt));
    flux_x1 = flux_x1 * dtmax/_dt;
    flux_x1 = std::max(-h * C,std::min(h_x1 * C,flux_x1));
    flux_x1 = flux_x1 * _dt/dtmax;

    if(z_x1 < -99995)
    {
        flux_x1 = flux_edge;
    }
    int edges = ((z_x1 < -99995)?1:0);

    bool edge = (z_x1 < -99995);
    float hold = h;
    h_add += 1.0*dx_w *(flux_x1);

    if((r == 688 && c== 591 && l ==3) || (r2 == 688 && c2== 591 && l2 ==3))
    {
        //FF::cout << "fsolve " << r << " " << c << " " << l << " loc2: " << r2 << " " << c2 << " " << l2 << " h:" <<  h << " vx "<< vxh <<  " " << vyh << " hx1:" << h_x1 << " vx1: " << vx_x1 << " " << vy_x1 << " dt: " << _dt << " there:" << flowstates[l][12]->data[r][c] << " "  << tglob << "  dhdt"  <<  flowstates[l][14]->data[r][c] << "  frac:" << std::min(1.0,std::max(0.0,(flowstates[l][12]->data[r][c]-tglob))/std::max(1e-12f,flowstates[l][34]->data[r][c])) << "  ha_add " << flux_x1/_dt << " " << flux_nb/_dt <<  "  dx " << dx << " "  << dy << "  " << dx_w << std::endl;

    }
    float f_centre = 0.5 * GRAV*((z+h) - (z<zc_x1?std::min(zc_x1,z + h) : std::min(z,zc_x1 + h_x1)))*(h_corr_x1 + h_corr_x1b);
    if(ca==1)
    {
      fx_center += dx_w * f_centre;

      if(ca2 == 1)
      {
        sx_av += dx_w*sx_zh;
        hll_x2.y = hll_x2.y+ dx_w * hll_x1h.y;
        hll_x2.z = hll_x2.z+ dx_w * hll_x1h.z;
      }else{
        sx_av += dx_w*sx_zh;
        hll_x1.y = hll_x1.y+dx_w * hll_x1h.y;
        hll_x1.z = hll_x1.z+dx_w * hll_x1h.z;
      }
    }else if(ca == 2)
    {
      fy_center += dx_w * f_centre;
      if(ca2 == 1)
      {
        sy_av += dx_w*sx_zh;
        hll_y2.y = hll_y2.y+dx_w * hll_x1h.y;
        hll_y2.z = hll_y2.z+dx_w * hll_x1h.z;
      }else{
        sy_av += dx_w*sx_zh;
        hll_y1.y = hll_y1.y+dx_w * hll_x1h.y;
        hll_y1.z = hll_y1.z+dx_w * hll_x1h.z;
      }
    }

    //does this neighbor have a timestep?
    if(flowstates[l][21] != nullptr)
    {
        //if not, use the cache here
        if(flowstates[l][13]->data[r][c] > 1e-12f)
        {
            float m = -1.0;
            int lt=l2;
            int rt=r2;
            int ct =c2;
            ca2 = ca2 *(-1);
            if(( dx < dy || (dx == dy && (r < r2 || c < c2))) || (r2 +1 >= dim1 || r2 < 0 || c2 +1 >= dim0 || c2 < 0))
            {
                m = 1.0;
                lt = l;
                rt = r;
                ct =c;
                ca2 = ca2 *(-1);
            }
            if(ca==1)
            {
              if(ca2 == 1)
              {
                flowstates[lt][24]->data[rt][ct] = m*flux_x1/std::max(1e-12f,_dt);
                flowstates[lt][40]->data[rt][ct] = m*sx_zh;
                flowstates[lt][25]->data[rt][ct] = m*hll_x1h.y/std::max(1e-12f,_dt);
                flowstates[lt][26]->data[rt][ct] = m*hll_x1h.z/std::max(1e-12f,_dt);
              }else{
                flowstates[lt][21]->data[rt][ct] = m*flux_x1/std::max(1e-12f,_dt);
                flowstates[lt][38]->data[rt][ct] = m*sx_zh;
                flowstates[lt][22]->data[rt][ct] = m*hll_x1h.y/std::max(1e-12f,_dt);
                flowstates[lt][23]->data[rt][ct] = m*hll_x1h.z/std::max(1e-12f,_dt);
              }
            }else if(ca == 2)
            {
              if(ca2 == 1)
              {
                flowstates[lt][30]->data[rt][ct] = m*flux_x1/std::max(1e-12f,_dt);
                flowstates[lt][45]->data[rt][ct] = m*sx_zh;
                flowstates[lt][31]->data[rt][ct] = m*hll_x1h.y/std::max(1e-12f,_dt);
                flowstates[lt][32]->data[rt][ct] = m*hll_x1h.z/std::max(1e-12f,_dt);
              }else{
                flowstates[lt][27]->data[rt][ct] = m*flux_x1/std::max(1e-12f,_dt);
                flowstates[lt][43]->data[rt][ct] = m*sx_zh;
                flowstates[lt][28]->data[rt][ct] = m*hll_x1h.y/std::max(1e-12f,_dt);
                flowstates[lt][29]->data[rt][ct] = m*hll_x1h.z/std::max(1e-12f,_dt);
              }
            }
        }
    }


}

inline double GetDHDT(float rain,float h_rate)
{
    return rain + h_rate;
}
inline LSMVector2 GetDUDTAddTerms(float hn, float sx_av, float sy_av,float _dt)
{
    return LSMVector2(- 0.5 * GRAV *hn*sx_av * _dt,- 0.5 * GRAV *hn*sy_av * _dt);
}
inline LSMVector2 GetDUDTNScale(float n, float hn, float dx, double _dt, float vx, float vy, float sx_av, float sy_av)
{
    float qxn = vx*hn;
    float qyn = vy*hn;
    float vsq = sqrt((float)(vx * vx + vy * vy));
    float nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01f,pow((float)(hn),(float)(4.0/3.0)));
    float nsq = nsq1*vsq*_dt;

    float vxn = (float)((qxn/(1.0f+nsq)))/std::max(0.01f,(float)(hn));
    float vyn = (float)((qyn/(1.0f+nsq)))/std::max(0.01f,(float)(hn));

    float threshold = std::min(0.1f,0.01f * dx);
    if(hn < threshold)
    {
        float kinfac = std::max(0.0f,(threshold - hn) / (0.025f * dx));
        float acc_eff = (vxn -vx)/std::max(0.0001,_dt);

        float v_kin = (sx_av>0.0f?1.0f:-1.0f) * hn * sqrt(hn) * sqrt(std::max(0.0f,sx_av>0.0f?sx_av:-sx_av))/(0.001f+n);

        vxn = kinfac * v_kin + vxn*(1.0f-kinfac);
    }

    if(hn < threshold)
    {
        float kinfac = std::max(0.0f,(threshold - hn) / (0.025f * dx));
        float acc_eff = (vyn -vy)/std::max(0.0001,_dt);

        float v_kin = (sy_av>0.0f?1.0f:-1.0f) * hn * sqrt(hn) * sqrt(std::max(0.0f,sy_av>0.0f?sy_av:-sy_av))/(0.001f+n);

        vyn = kinfac * v_kin + vyn*(1.0f-kinfac);
    }

    return LSMVector2(vxn/vx, vyn/vy);
}
inline LSMVector2 GetDUDT(float n, float h, float hn, float dx, double _dt, float vx, float vy, LSMVector3 hll_x1, LSMVector3 hll_x2, LSMVector3 hll_y1, LSMVector3 hll_y2, float sx_av, float sy_av, float fx_center, float fy_center)
{
    float tx = _dt/dx;
    float qxn = h * vx - (tx/dx)*(hll_x2.y - hll_x1.y) - (tx/dx)*(hll_y2.z - hll_y1.z) - 0.5 * GRAV *hn*sx_av * _dt;//+ fx_center;
    float qyn = h * vy - (tx/dx)*(hll_x2.z - hll_x1.z) - (tx/dx)*(hll_y2.y - hll_y1.y) - 0.5 * GRAV *hn*sy_av * _dt;//+ fy_center;

    float vsq = sqrt((float)(vx * vx + vy * vy));
    float nsq1 = (0.001+n)*(0.001+n)*GRAV/std::max(0.01f,pow((float)(hn),(float)(4.0/3.0)));
    float nsq = nsq1*vsq*_dt;

    float vxn = (float)((qxn/(1.0f+nsq)))/std::max(0.01f,(float)(hn));
    float vyn = (float)((qyn/(1.0f+nsq)))/std::max(0.01f,(float)(hn));

    float threshold = std::min(0.1f,0.01f * dx);
    if(hn < threshold)
    {
        float kinfac = std::max(0.0f,(threshold - hn) / (0.025f * dx));
        float acc_eff = (vxn -vx)/std::max(0.0001,_dt);

        float v_kin = (sx_av>0.0f?1.0f:-1.0f) * hn * sqrt(hn) * sqrt(std::max(0.0f,sx_av>0.0f?sx_av:-sx_av))/(0.001f+n);

        vxn = kinfac * v_kin + vxn*(1.0f-kinfac);
    }

    if(hn < threshold)
    {
        float kinfac = std::max(0.0f,(threshold - hn) / (0.025f * dx));
        float acc_eff = (vyn -vy)/std::max(0.0001,_dt);

        float v_kin = (sy_av>0.0f?1.0f:-1.0f) * hn * sqrt(hn) * sqrt(std::max(0.0f,sy_av>0.0f?sy_av:-sy_av))/(0.001f+n);

        vyn = kinfac * v_kin + vyn*(1.0f-kinfac);
    }

    vxn = std::min(25.0f,std::max(-25.0f,vxn));
    vyn = std::min(25.0f,std::max(-25.0f,vyn));
    if(!std::isnormal(vxn))
    {
        vxn = 0.0;
    }
    if(!std::isnormal(vyn))
    {
        vyn = 0.0;
    }

    return LSMVector2((vxn-vx)/_dt, (vyn-vy)/_dt);

}


//only do one edge direction, for a cached supercell
inline void DoMassMomentum3CascCalcCached(double tglob, double courant, double _dt,float dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l,
                                   bool do_ocean, float ocean_level, float rain_fac,
                                   float dx_cache, float *cell_cache_h , float *cell_cache_vx, float * cell_cache_vy,float * area_total,
                                   float cell_cache_hscale, float cell_cache_vxscale, float cell_cache_vyscale,
                                   float *cell_cache_dhdt, float * cell_cache_dvxdt, float * cell_cache_dvydt,
                                   int cae, int ca2e)
{
    if(!(r < flowstates[l][0]->nrRows() && c < flowstates[l][0]->nrCols()))
    {
        return;
    }
    LSMVector3 hll_x1;
    LSMVector3 hll_y1;
    LSMVector3 hll_x2;
    LSMVector3 hll_y2;

    float h_add = 0.0;
    float fx_center = 0.0;
    float fy_center= 0.0;
    float sx_av = 0.0;
    float sy_av = 0.0;

    float z = flowstates[l][0]->data[r][c];
    float h = flowstates[l][4]->data[r][c];
    float vx = flowstates[l][5]->data[r][c];
    float vy = flowstates[l][6]->data[r][c];
    float n = flowstates[l][3]->data[r][c];

    float h_org= h;
    if(cell_cache_h != nullptr)
    {
        //FF::cout << " run cached gridcell " << r << " "  << c << " " << l << " " << cell_cache_hscale << " " << cell_cache_vxscale << "  "  << cell_cache_vyscale <<  "  " << dx <<  " "  << dx_cache << std::endl;
        h = h* cell_cache_hscale;
        vx = vx * cell_cache_vxscale;
        vy = vy * cell_cache_vyscale;
    }else if(flowstates[l][21] != nullptr)
    {
        //potentially go back in time?
        float dt_frac_r = std::min(1.0,std::max(0.0,(flowstates[l][12]->data[r][c]-tglob))/std::max(1e-12f,flowstates[l][34]->data[r][c]));
        float h_org = h;
        h = std::max(0.0f, h - flowstates[l][34]->data[r][c]* flowstates[l][14]->data[r][c] * dt_frac_r);
        vx -= flowstates[l][34]->data[r][c]*flowstates[l][15]->data[r][c] * dt_frac_r;
        vy -= flowstates[l][34]->data[r][c]*flowstates[l][16]->data[r][c] * dt_frac_r;


        if((r == 688 && c== 591 && l ==3))
        {
            //FF::cout << "!!!!!run our cell tglob: " << tglob << " there: " << flowstates[l][12]->data[r][c] << "dt:" << _dt << " frac: " << dt_frac_r << " h:" << h << " h_org: " << h_org << "  tpast:" << flowstates[l][34]->data[r][c] << " dhdt " << flowstates[l][14]->data[r][c] <<   std::endl;
        }

    }

    auto mmlambda = [tglob, l, r,c,_dt,dx,  &flowstates,&fx_center, do_ocean, ocean_level, &fy_center, h,vx,vy,z,n,&h_add,&hll_x1,&hll_y1, &hll_x2, &hll_y2, &sx_av, &sy_av](int r2,int c2,int l2,float dy, int ca, int ca2)
    {

        //call the actual function
        MomentumFluxGet(tglob, l, r, c, _dt, dx, flowstates,
                        do_ocean, ocean_level, h, vx, vy, z, n,
                        fx_center , fy_center ,
                        h_add, hll_x1 , hll_x2 , hll_y1 , hll_y2 ,
                        sx_av , sy_av ,
                        r2, c2, l2, dy, ca, ca2);
    };


    //pass function pointer instead of lamba for wider access
    if(cell_cache_h == nullptr)
    {
        ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r,c+1,l,1,0,true, mmlambda,1,1);
        ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r,c-1,l,-1,0,true, mmlambda,1,-1);
        ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r+1,c,l,0,1,true, mmlambda,2,1);
        ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r-1,c,l,0,-1,true, mmlambda,2,-1);
    }else
    {
        //ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r,c+1,l,1,0,true, mmlambda,cae,ca2e);

        if(cae == 1 && ca2e == 1)
        {
            ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r,c+1,l,1,0,true, mmlambda,1,1);
        }
        if(cae == 1 && ca2e == -1)
        {
            ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r,c-1,l,-1,0,true, mmlambda,1,-1);
        }
        if(cae == 2 && ca2e == 1)
        {
            ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r+1,c,l,0,1,true, mmlambda,2,1);
        }
        if(cae == 2 && ca2e == -1)
        {
            ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r-1,c,l,0,-1,true, mmlambda,2,-1);
        }
    }

    sx_av = sx_av/dx;
    sy_av = sy_av/dx;
    h_add = h_add/dx;

    //now finally get dh/dt and dux/dy and duy/dt
    double dhdt = 0.0;
    double raindhdt = 0.0;//
    if(cell_cache_h == nullptr)
    {
        raindhdt= (rain_fac * flowstates[l][10]->data[r][c]/(3600000.0));
        //FF::cout << " rain normal " << (rain_fac * flowstates[l][10]->data[r][c]/(3600000.0)) * _dt << std::endl;
        dhdt = GetDHDT(raindhdt, h_add/(_dt));
    }else
    {
        //*area_total += dx* dx;
        dhdt = GetDHDT(0.0, h_add/(_dt));
    }

    //and alter the values
    double hn = h + dhdt * _dt;
    hn= std::min(1e4,std::max(0.0,hn));
    if(!std::isfinite(hn))
    {
        hn =0.0;
    }


    //FF::cout << "rain add " << r << " " << c << " " << l << " rain: " << (rain_fac * flowstates[l][10]->data[r][c]/(3600000.0)) << " " << flowstates[l][10]->data[r][c] <<  "  dt: " << _dt << " h:" << h << " " << hn << "  dhdt: " << dhdt << std::endl;

    if( hn > 5.0)
    {
        //FF::cout << "msolve " << r << " " << c << " " << l << " h:" <<  h << " hn:" << hn << " dt: " << _dt << " dhdt " << dhdt << "  there:" << flowstates[l][12]->data[r][c] << " "  << tglob << "  horg:" << h_org << "  dhdt"  <<  flowstates[l][14]->data[r][c] << "  frac:" << std::min(1.0,std::max(0.0,(flowstates[l][12]->data[r][c]-tglob))/std::max(1e-12f,flowstates[l][34]->data[r][c])) << std::endl;
    }

    if(r == 15 && c== 15 && l ==0)
    {
        //FF::cout << "rain add " << r << " " << c << " " << l << " rain: " << (rain_fac * flowstates[l][10]->data[r][c]/(3600000.0)) << " " << flowstates[l][10]->data[r][c] <<  " tglob: " << tglob << " t:" << flowstates[l][12]->data[r][c] << " dt: " << _dt << " h:" << h << " " << hn << "  dhdt: " << dhdt << std::endl;

        //FF::cout << "msolve " << r << " " << c << " " << l << " h:" <<  h << " hn:" << hn << " dt: " << _dt << " dhdt " << dhdt << "  there:" << flowstates[l][12]->data[r][c] << " "  << tglob << "  horg:" << h_org << "  dhdt"  <<  flowstates[l][14]->data[r][c] << "  frac:" << std::min(1.0,std::max(0.0,(flowstates[l][12]->data[r][c]-tglob))/std::max(1e-12f,flowstates[l][34]->data[r][c])) << std::endl;

    }

    //convert to dhdt for the cache superpixel

    LSMVector2 dudt = GetDUDT(n,h,hn, dx,_dt,vx,vy,hll_x1,hll_x2,hll_y1,hll_y2,sx_av,sy_av,fx_center,fy_center);
    //convert to duxdt for cache superpixel

    float vxn = vx + _dt * dudt.x;
    float vyn = vy + _dt * dudt.y;

    //write output
    if(cell_cache_h == nullptr)
    {
        flowstates[l][7]->data[r][c] = hn;
        flowstates[l][8]->data[r][c] = vxn;
        flowstates[l][9]->data[r][c] = vyn;

        if(flowstates[l][21] != nullptr)
        {
            flowstates[l][35]->data[r][c] = dhdt;
            flowstates[l][36]->data[r][c] = dudt.x;
            flowstates[l][37]->data[r][c] = dudt.y;
        }

    }else
    {


        *cell_cache_h = std::min(1e4,std::max(0.0,*cell_cache_h + dhdt * _dt * dx*dx/(dx_cache*dx_cache)));
        //*cell_cache_vx = *cell_cache_vx + dudt.x * _dt * dx*dx/(dx_cache*dx_cache);
        //*cell_cache_vy = *cell_cache_vy + dudt.y * _dt * dx*dx/(dx_cache*dx_cache);

        if(flowstates[l][21] != nullptr)
        {
            *cell_cache_dhdt += raindhdt * dx*dx/(dx_cache*dx_cache) +  dhdt * dx*dx/(dx_cache*dx_cache);
            *cell_cache_dvxdt += dudt.x * dx*dx/(dx_cache*dx_cache);
            *cell_cache_dvydt += dudt.y * dx*dx/(dx_cache*dx_cache);
        }
    }

    //update cell time

}


inline void DoMassMomentum3CascCalc(double tglob,double courant, double _dt,float dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l,
                             bool do_ocean, float ocean_level, float rain_fac)
{
    LSMVector3 hll_x1;
    LSMVector3 hll_y1;
    LSMVector3 hll_x2;
    LSMVector3 hll_y2;

    float h_add = 0.0;
    float fx_center = 0.0;
    float fy_center= 0.0;
    float sx_av = 0.0;
    float sy_av = 0.0;

    float z = flowstates[l][0]->data[r][c];
    float h = flowstates[l][4]->data[r][c];
    float vx = flowstates[l][5]->data[r][c];
    float vy = flowstates[l][6]->data[r][c];
    float n = flowstates[l][3]->data[r][c];


    auto mmlambda = [tglob,l, r,c,_dt,dx,  &flowstates,&fx_center, do_ocean, ocean_level, &fy_center, h,vx,vy,z,n,&h_add,&hll_x1,&hll_y1, &hll_x2, &hll_y2, &sx_av, &sy_av](int r2,int c2,int l2,float dy, int ca, int ca2)
    {
        //call the actual function
        MomentumFluxGet(tglob, l, r, c, _dt, dx, flowstates,
                        do_ocean, ocean_level, h, vx, vy, z, n,
                        fx_center , fy_center ,
                        h_add, hll_x1 , hll_x2 , hll_y1 , hll_y2 ,
                        sx_av , sy_av ,
                        r2, c2, l2, dy, ca, ca2);
    };

    //pass function pointer instead of lamba for wider access

    ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r,c+1,l,1,0,true, mmlambda,1,1);
    ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r,c-1,l,-1,0,true, mmlambda,1,-1);
    ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r+1,c,l,0,1,true, mmlambda,2,1);
    ForEach3CascNeighbor(_dt,dx,flowstates, r, c, l, r-1,c,l,0,-1,true, mmlambda,2,-1);

    sx_av = sx_av/dx;
    sy_av = sy_av/dx;
    h_add = h_add/dx;


    //now finally get dh/dt and dux/dy and duy/dt
    float dhdt = GetDHDT((rain_fac * flowstates[l][10]->data[r][c]/(3600000.0)), h_add/(_dt));

    //and alter the values
    float hn = h + dhdt * _dt;
    hn= std::min(1e4f,std::max(0.0f,hn));
    if(!std::isfinite(hn))
    {
        hn =0.0;
    }

    LSMVector2 dudt = GetDUDT(n,h,hn, dx,_dt,vx,vy,hll_x1,hll_x2,hll_y1,hll_y2,sx_av,sy_av,fx_center,fy_center);

    float vxn = vx + _dt * dudt.x;
    float vyn = vy + _dt * dudt.y;

    //write output
    flowstates[l][7]->data[r][c] = hn;
    flowstates[l][8]->data[r][c] = vxn;
    flowstates[l][9]->data[r][c] = vyn;
}

/*
 * @brief
 * for the current pixel, solve the new height and velocities, by calling the flux solver function for each edge,
 * and each other gridcell along that edge (could be multiple due to the cascading grid)
 *
 */

inline void DoMassMomentum3Casc(double tglob, double courant, double _dt,float dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, bool do_ocean, float ocean_level, float rain_fac)
{

  if(!(r>-1 && c > -1 && r < flowstates[l][0]->Rows() && c < flowstates[l][0]->Cols()))
  {
    return;
  }


  //skip if the entire supercell is not required due to sdt




  if(flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
  {
      bool is_cache = false;
      if(l == 0)
      {
          if(flowstates[l][17] != nullptr)
          {
             if(flowstates[l][17]->data[r][c] > 0.5f)
             {
                 is_cache = true;
             }
          }
      }


      if(is_cache)
      {
          //check if we run a timestep here now
          //if current cell, or required by any neighbor


          if(flowstates[0][21] != nullptr)
          {
              bool do_thiscell = flowstates[l][13]->data[r][c] > 1e-12f;
              //any neighbors required?
              if(!do_thiscell)
              {
                  return;
              }
          }



          //set reference to these water properties
          float cell_cache_dx = dx;
          float *cell_cache_h = &(flowstates[l][7]->data[r][c]);
          float *cell_cache_vx = &(flowstates[l][8]->data[r][c]);
          float *cell_cache_vy = &(flowstates[l][9]->data[r][c]);

          float *cell_cache_dhdt = nullptr;
          float *cell_cache_dvxdt = nullptr;
          float *cell_cache_dvydt = nullptr;


          float h_org = flowstates[l][4]->data[r][c];
          float vx_org = flowstates[l][5]->data[r][c];
          float vy_org = flowstates[l][6]->data[r][c];

          //check if we must go back in time a bit for this cell
          if(flowstates[0][21] != nullptr)
          {
             float dt_frac_r = std::min(1.0,std::max(0.0,(flowstates[l][12]->data[r][c]-tglob))/std::max(1e-12f,flowstates[l][34]->data[r][c]));
             h_org -= flowstates[l][34]->data[r][c]* flowstates[l][14]->data[r][c] * dt_frac_r;
             vx_org -=flowstates[l][34]->data[r][c]* flowstates[l][15]->data[r][c] * dt_frac_r;
             vy_org -= flowstates[l][34]->data[r][c]*flowstates[l][16]->data[r][c] * dt_frac_r;
             flowstates[l][7]->data[r][c] = h_org;
             flowstates[l][8]->data[r][c] = vx_org;
             flowstates[l][9]->data[r][c] = vy_org;

             flowstates[l][35]->data[r][c] = 0.0;
             flowstates[l][36]->data[r][c] = 0.0;
             flowstates[l][37]->data[r][c] = 0.0;

             cell_cache_dhdt = &(flowstates[l][35]->data[r][c]);
             cell_cache_dvxdt = &(flowstates[l][36]->data[r][c]);
             cell_cache_dvydt = &(flowstates[l][37]->data[r][c]);

             _dt = flowstates[l][13]->data[r][c];

          }


          float area_total = 0.0;
          float cell_cache_hscale = h_org/(std::max(1e-12f,flowstates[0][18]->data[r][c]));
          float cell_cache_vxscale = vx_org/(flowstates[0][19]->data[r][c] > 0.0f? std::max(1e-12f,flowstates[0][19]->data[r][c]):std::min(-1e-12f,flowstates[0][19]->data[r][c]));
          float cell_cache_vyscale = vy_org/(flowstates[0][20]->data[r][c] > 0.0f? std::max(1e-12f,flowstates[0][20]->data[r][c]):std::min(-1e-12f,flowstates[0][20]->data[r][c]));

          //FF::cout << " cachce cell main " << r << "  " << c << " " << l << "  " << cell_cache_h << std::endl;
          //only do edges

          ForEach3CascEdgeP(_dt,dx,flowstates, r, c, l, r,c,l,true,[&](int r2e,int c2e,int l2e,float dye, int cae, int ca2e){

            DoMassMomentum3CascCalcCached(tglob,courant,_dt,dye,flowstates, r2e,c2e, l2e, do_ocean, ocean_level, rain_fac,
                                          cell_cache_dx, cell_cache_h, cell_cache_vx, cell_cache_vy, &area_total,
                                          cell_cache_hscale, cell_cache_vxscale, cell_cache_vyscale,
                                          cell_cache_dhdt, cell_cache_dvxdt, cell_cache_dvydt,
                                          cae, ca2e);
          });

          //collect changes to the values

          //FF::cout << "number of pixels in cached pixel: " << area_total << " " << cell_cache_dx*cell_cache_dx << std::endl;
          //compensate for missing precipitation
          float area_extra = std::max(0.0f,cell_cache_dx * cell_cache_dx - area_total);
          float h_add_rain = (rain_fac * flowstates[l][10]->data[r][c]/(3600000.0)) * _dt * area_extra/(cell_cache_dx*cell_cache_dx);

          //FF::cout << "h_add " << h_add_rain << std::endl;
          flowstates[l][7]->data[r][c] += h_add_rain;

          if(flowstates[0][21] != nullptr)
          {
              flowstates[l][35]->data[r][c] += h_add_rain/_dt;
          }

          //set TN and DTL
          if(flowstates[0][21] != nullptr)
          {
              flowstates[l][33]->data[r][c] = tglob + flowstates[l][13]->data[r][c];
              flowstates[l][46]->data[r][c] = flowstates[l][13]->data[r][c];
          }

          //FF::cout << "export ratio cache: "  << " " << r << "  " << c << " " << flowstates[l][7]->data[r][c] << " " << std::max(1e-12f,flowstates[l][18]->data[r][c]) << "  " << flowstates[l][7]->data[r][c]/std::max(1e-12f,flowstates[l][18]->data[r][c]) << std::endl;


      }else
      {
          //or do all
          for(int r2 = 0; r2 < 3; r2++)
          {
              for (int c2 = 0; c2 < 3; c2++)
              {
                  DoMassMomentum3Casc(tglob,courant,_dt,dx/ 3.0,flowstates, r * 3 + r2,c*3 + c2, l+1, do_ocean, ocean_level, rain_fac);
              }
          }
      }

  }else{

      //check if we do a timestep here now

      if(flowstates[0][21] != nullptr)
      {
          if(flowstates[0][21] != nullptr)
          {
              bool do_thiscell = flowstates[l][13]->data[r][c] > 1e-12f;
              //any neighbors required?
              if(!do_thiscell)
              {
                  //FF::cout << "dont model " << std::endl;
                  return;
              }else
              {
                  //FF::cout <<" !!!!!!!!!!  " << std::endl;
              }
          }
          _dt = flowstates[l][13]->data[r][c];
      }

      //FF::cout << "do model " << std::endl;

      DoMassMomentum3CascCalcCached(tglob,courant,_dt,dx,flowstates, r,c, l, do_ocean, ocean_level, rain_fac,
                                    dx,nullptr,nullptr,nullptr,nullptr,1.0,1.0,1.0,nullptr, nullptr, nullptr,0,0);

      //set TN and DTL
      if(flowstates[0][21] != nullptr)
      {
          flowstates[l][33]->data[r][c] = tglob + flowstates[l][13]->data[r][c];
          flowstates[l][46]->data[r][c] = flowstates[l][13]->data[r][c];

      }

  }
}


#endif
