#ifndef FLOWREFINE_H
#define FLOWREFINE_H

#include "map.h"
#include "upscale.h"


inline void UpscaleArtificialFlow(cTMapf * DEM,
                            cTMapf * H,
                            int iter,
                            float couranti, bool do_stable = false)
{


    //then the diffusive max wave
    cTMapf *map = H;
    cTMapf *map2 = DEM->GetCopy0();

    cTMapf *qxml = DEM->GetCopy0();
    cTMapf *qyml = DEM->GetCopy0();
    cTMapf *qxmr = DEM->GetCopy0();
    cTMapf *qymr = DEM->GetCopy0();


    const float coeff_boundary = 0.5;


    float dx = DEM->cellSizeX();
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sxl,syl,sxr,syr;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

    for(int r = 0; r < DEM->nrRows();r++)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                map->data[r][c] = (!pcr::isMV(H->data[r][c]))? H->Drc : 0.0;
                map2->data[r][c] = map->data[r][c];
            }
        }
    }

    for(int i = 0; i < iter; i++)
    {
        float courant = couranti;

        float fac_sourcelimit = 1.0;
        if(i > iter - 6)
        {
            fac_sourcelimit= 1.0;
            float in = i-(iter-6);
            courant = couranti * ((1.0- in/6.0)*1.0) + (in/6.0) * 0.25 * couranti;
            if(do_stable)
            {
                courant = couranti;
            }
        }
        {
            for(int r = 0; r < DEM->nrRows();r++)
            {

                for(int c = DEM->data.nr_cols() -1; c > -1 ;c--)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {

                        //reset variables
                        hn = 0.0;
                        qmax_xr = 0.0;
                        qmax_xl = 0.0;
                        qmax_yl = 0.0;
                        qmax_yr = 0.0;
                        signx = 0.0;
                        signy = 0.0;
                        vx = 0.0;
                        vy = 0.0;
                        qx = 0.0;
                        qy = 0.0;
                        sxl = 0.0;
                        syl = 0.0;
                        sxr = 0.0;
                        syr = 0.0;
                        mv_xl = false;
                        mv_xr = false;
                        mv_yl = false;
                        mv_yr = false;

                        dem =DEM->data[r][c];

                        //slope
                        if(c > 0)
                        {
                            dem_xl = pcr::isMV(DEM->data[r][c-1])?dem:DEM->data[r][c-1];
                        }else {
                            dem_xl = dem;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            dem_xr = pcr::isMV(DEM->data[r][c+1])?dem:DEM->data[r][c+1];
                        }else {
                            dem_xr = dem;
                        }
                        if(r > 0)
                        {
                            dem_yl = pcr::isMV(DEM->data[r-1][c])?dem:DEM->data[r-1][c];
                        }else {
                            dem_yl = dem;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            dem_yr = pcr::isMV(DEM->data[r+1][c])?dem:DEM->data[r+1][c];
                        }else {
                            dem_yr = dem;

                        }

                        h = map->Drc;

                        if(c > 0)
                        {
                            h_xl = pcr::isMV(map->data[r][c-1])? coeff_boundary*h:map->data[r][c-1];
                        }else {
                            h_xl = coeff_boundary*h;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            h_xr = pcr::isMV(map->data[r][c+1])?coeff_boundary*h:map->data[r][c+1];
                        }else {
                            h_xr = coeff_boundary*h;
                        }
                        if(r > 0)
                        {
                            h_yl = pcr::isMV(map->data[r-1][c])?coeff_boundary*h:map->data[r-1][c];
                        }else {
                            h_yl = coeff_boundary*h;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            h_yr = pcr::isMV(map->data[r+1][c])?coeff_boundary*h:map->data[r+1][c];
                        }else {
                            h_yr = coeff_boundary*h;
                        }


                        sxl = -((h_xl +dem_xl)-(h+dem));
                        syl = -((h_yl +dem_yl)-(h+dem));
                        sxr = -((h_xr +dem_xr)-(h+dem));
                        syr = -((h_yr +dem_yr)-(h+dem));


                        float sxz = (dem_xr)-(dem_xl);
                        float syz = (dem_yr)-(dem_yl);
                        //velocity
                        sl = std::max(0.0f,sxl) + std::max(0.0f,syl) + std::max(0.0f,sxr) + std::max(0.0f,syr);
                        float slz = std::sqrt(sxz*sxz + syz*syz);

                        //flux limiter
                        if(do_stable)
                        {
                            h_xl = (1.0-fac_sourcelimit)*h + fac_sourcelimit *std::min(0.15f*h,std::max(0.0f, (dem+ h - dem_xl - h_xl)*(0.5f*fac_sourcelimit)));
                            h_xr = (1.0-fac_sourcelimit)*h + fac_sourcelimit *std::min(0.15f*h,std::max(0.0f, (dem+ h - dem_xr - h_xr)*(0.5f*fac_sourcelimit)));
                            h_yl = (1.0-fac_sourcelimit)*h + fac_sourcelimit *std::min(0.15f*h,std::max(0.0f, (dem+ h - dem_yl - h_yl)*(0.5f*fac_sourcelimit)));
                            h_yr = (1.0-fac_sourcelimit)*h + fac_sourcelimit *std::min(0.15f*h,std::max(0.0f, (dem+ h - dem_yr - h_yr)*(0.5f*fac_sourcelimit)));
                        }else
                        {
                            h_xl = h;
                            h_xr = h;
                            h_yl = h;
                            h_yr = h;
                        }
                        //change height
                        qxml->data[r][c] = (courant * h_xl* std::max(0.0f,sxl)/(std::max(1e-12f,sl)));
                        qxmr->data[r][c] = (courant * h_xr* std::max(0.0f,sxr)/(std::max(1e-12f,sl)));
                        qyml->data[r][c] = (courant * h_yl* std::max(0.0f,syl)/(std::max(1e-12f,sl)));
                        qymr->data[r][c] = (courant * h_yr* std::max(0.0f,syr)/(std::max(1e-12f,sl)));

                    }
                }
            }

            for(int r = 0; r < DEM->nrRows();r++)
            {

                for(int c = DEM->data.nr_cols() -1; c > -1 ;c--)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        float qxl = qxml->data[r][c];
                        float qyl = qyml->data[r][c];
                        float qxr = qxmr->data[r][c];
                        float qyr = qymr->data[r][c];
                        map2->data[r][c] = std::max(0.0f,map->data[r][c] - qxl - qyl - qxr - qyr);

                        if(c > 0)
                        {
                            map2->data[r][c] += qxmr->data[r][c-1];
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            map2->data[r][c] += qxml->data[r][c+1];
                        }

                        if(r > 0)
                        {
                            map2->data[r][c] += qymr->data[r-1][c];
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            map2->data[r][c] += qyml->data[r+1][c];
                        }

                    }
                }
            }

        }

        cTMapf* temp = map2;
        map2 = map;
        map = temp;
    }

    if(map2 != H)
    {
        delete map2;
    }else
    {
        delete map;
    }
    delete qxml;
    delete qyml;
    delete qxmr;
    delete qymr;
}

#endif // CASCADEREFINE_H
