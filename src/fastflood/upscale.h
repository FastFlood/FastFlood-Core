#ifndef UPSCALE_H
#define UPSCALE_H
#include "../common/defines.h"

extern "C"
{
    
inline void EMSCRIPTEN_KEEPALIVE do_algo_upscale(
    float * ret,
    int m_w,
    int m_h,
    float * data_dem,
    float dem_ulx,
    float dem_uly,
    float dem_dx,
    float dem_dy,
    float * data_wh,
    int m_w2,
    int m_h2,
    float wh_ulx,
    float wh_uly,
    float wh_dx,
    float wh_dy,
    int iter,
    float couranti,
    int ocean,
    float ocean_level
)
{
    FF::cout << "upscale " << dem_ulx << " " << dem_uly << " " << dem_dx << " " << dem_dy << std::endl;
    FF::cout << "upscale " << wh_ulx << " " << wh_uly << " " << wh_dx << " " << wh_dy << std::endl;


    cTMapf * DEM = new cTMapf(data_dem,m_w,m_h,0.0,0.0, dem_dx,dem_dy);
	cTMapf * WH_ORG = new cTMapf(data_wh,m_w2,m_h2,0.0,0.0, wh_dx,wh_dy);

    cTMapf * H = DEM->GetCopy0();

    //interpolation first
    GetTPool().parallel_for(0,m_h,[&](int r)
    {
        for(int c = 0; c < m_w; c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                int r2 = ((dem_uly + r * dem_dy) - wh_uly)/wh_dy;
                int c2 = ((dem_ulx + c * dem_dx) - wh_ulx)/wh_dx;

                //FF::cout << "r/c " << r2 << " " << c2 << std::endl;
                if(r2 > -1 && r2 < m_h2)
                {
                    if(c2 > -1 && c2 < m_w2)
                    {
                        if(!pcr::isMV(WH_ORG->data[r2][c2]))
                        {
                            //FF::cout << WH_ORG->data[r2][c2] << std::endl;
                            H->data[r][c] = WH_ORG->data[r2][c2];
                        }
                    }
                }

            }
        }
    });

    //then the diffusive max wave
    cTMapf *map = DEM->GetCopy0();
    cTMapf *map2 = DEM->GetCopy0();


    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

    GetTPool().parallel_for(0,DEM->nrRows(),[&](int r)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                map->data[r][c] = (!pcr::isMV(H->data[r][c]))? H->Drc : 0.0;
                map2->data[r][c] = map->data[r][c];
            }
        }
    });
    GetTPool().waitext();

    for(int i = 0; i < iter; i++)
    {
        float courant = couranti;
        if(i > iter - 10)
        {
            float in = i-(iter-10);
            courant = couranti * ((1.0- in/10.0)*1.0) + (in/10.0) * 0.05 * couranti;
        FF::cout << "iter " << i << " " << courant << std::endl;
        }

        if(i%2 == 0)
        {
            GetTPool().parallel_for(0,DEM->nrRows(),[&](int r)
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
                        sx = 0.0;
                        sy = 0.0;
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
                            h_xl = pcr::isMV(map->data[r][c-1])?h:map->data[r][c-1];
                        }else {
                            h_xl = h;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            h_xr = pcr::isMV(map->data[r][c+1])?h:map->data[r][c+1];
                        }else {
                            h_xr = h;
                        }
                        if(r > 0)
                        {
                            h_yl = pcr::isMV(map->data[r-1][c])?h:map->data[r-1][c];
                        }else {
                            h_yl = h;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            h_yr = pcr::isMV(map->data[r+1][c])?h:map->data[r+1][c];
                        }else {
                            h_yr = h;
                        }


                        sx = (h_xr +dem_xr)-(h_xl+dem_xl);
                        sy = (h_yr +dem_yr)-(h_yl+dem_yl);


                        //velocity
                        sl = std::sqrt(sx*sx + sy*sy);
                        vx = sx/sl;
                        vy = sy/sl;
                        signx = vx > 0.0? -1.0:1.0;
                        signy = vy > 0.0? -1.0:1.0;

                        //flux limiter

                        qmax_xl = std::max(0.0f,dem+h - dem_xl - h_xl);
                        qmax_xr = std::max(0.0f,dem+h - dem_xr - h_xr);
                        qmax_yl = std::max(0.0f,dem+h - dem_yl - h_yl);
                        qmax_yr = std::max(0.0f,dem+h - dem_yr - h_yr);

                        //discharge

                        qmaxx = signx > 0.0? qmax_xr : qmax_xl;
                        qmaxy = signy > 0.0? qmax_yr : qmax_yl;

                        qx = std::min(courant*h,std::min(qmaxx,std::abs(courant*h * vx)));
                        qy = std::min(courant*h,std::min(qmaxy,std::abs(courant*h * vy)));
                        //change height

                        hn = h - qx - qy;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy);

                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                            }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                            }
                        }

                    }
                }
            });
            GetTPool().waitext();

        }else {
            GetTPool().parallel_for(0,DEM->nrRows(),[&](int r)
            {
                for(int c = 0; c < DEM->data.nr_cols();c++)
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
                        sx = 0.0;
                        sy = 0.0;
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
                            h_xl = pcr::isMV(map->data[r][c-1])?h:map->data[r][c-1];
                        }else {
                            h_xl = h;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            h_xr = pcr::isMV(map->data[r][c+1])?h:map->data[r][c+1];
                        }else {
                            h_xr = h;
                        }
                        if(r > 0)
                        {
                            h_yl = pcr::isMV(map->data[r-1][c])?h:map->data[r-1][c];
                        }else {
                            h_yl = h;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            h_yr = pcr::isMV(map->data[r+1][c])?h:map->data[r+1][c];
                        }else {
                            h_yr = h;
                        }


                        sx = (h_xr +dem_xr)-(h_xl+dem_xl);
                        sy = (h_yr +dem_yr)-(h_yl+dem_yl);


                        //velocity
                        sl = std::sqrt(sx*sx + sy*sy);
                        vx = sx/sl;
                        vy = sy/sl;
                        signx = vx > 0.0? -1.0:1.0;
                        signy = vy > 0.0? -1.0:1.0;

                        //flux limiter

                        qmax_xl = std::max(0.0f,dem+h - dem_xl - h_xl);
                        qmax_xr = std::max(0.0f,dem+h - dem_xr - h_xr);
                        qmax_yl = std::max(0.0f,dem+h - dem_yl - h_yl);
                        qmax_yr = std::max(0.0f,dem+h - dem_yr - h_yr);

                        //discharge

                        qmaxx = signx > 0.0? qmax_xr : qmax_xl;
                        qmaxy = signy > 0.0? qmax_yr : qmax_yl;

                        qx = std::min(courant*h,std::min(qmaxx,std::abs(courant*h * vx)));
                        qy = std::min(courant*h,std::min(qmaxy,std::abs(courant*h * vy)));
                        //change height

                        hn = h - qx - qy;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy);

                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                            }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                            }
                        }

                    }
                }
            });
            GetTPool().waitext();


        }

        GetTPool().parallel_for(0,DEM->nrRows(),[&](int r)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];
            }
        });
        GetTPool().waitext();

    }


    if(ocean > 0)
    {
    
        printf("ocean level %g", ocean_level);
        cTMapf * WHO = do_fill_ocean(DEM,ocean_level,40);

        GetTPool().parallel_for(0,DEM->nrRows(),[&](int r) {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        map->data[r][c] = std::max(map->data[r][c],WHO->data[r][c]);
                    }
                }
        });
        GetTPool().waitext();
        delete WHO;
    }

    delete map2;


        GetTPool().parallel_for(0,DEM->nrRows(),[&](int r) {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                    int index = r* DEM->Cols() + c;
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        ret[index] = map->data[r][c];
                    }
                }
        });
        GetTPool().waitext();


    //finally re-enforce potential coastal boundary condition

    delete H;
    delete map;
}


}

#endif
