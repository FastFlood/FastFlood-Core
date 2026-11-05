#ifndef MAXWAVE_H
#define MAXWAVE_H

const float CULVERT_NEIGHBOR_FACTOR = 0.85f;

extern "C"
{

    
inline cTMapf * DiffusiveMaxWavePrefDir(cTMapf * DEM,cTMapf * H, cTMapf * DirX, cTMapf * DirY, int iter, float couranti)
{

    ////FF::cout << "max wave prefdir " << DEM << " " << H << " " << DirX << " " << DirY << " dx" <<DEM->data.cell_size() << " dy" << DEM->data.cell_sizeY() << std::endl;

    cTMapf *map = DEM->GetCopy0();


    MaskedRaster<float> raster_data2(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
    cTMapf *map2 = DEM->GetCopy0();
    cTMapf *mapx = map->GetCopy0();
    cTMapf *mapx2 = map->GetCopy0();
    cTMapf *mapy = map->GetCopy0();
    cTMapf *mapy2 = map->GetCopy0();



    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

#pragma omp parallel for collapse(2)
    for(int r = 0; r < DEM->data.nr_rows();r++)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                map->data[r][c] = (!pcr::isMV(H->data[r][c]))? H->Drc : 0.0;
                map2->data[r][c] = map->data[r][c];
                mapx->data[r][c] =  map->data[r][c] * DirX->data[r][c];
                mapy->data[r][c] =  map->data[r][c] * DirY->data[r][c];

            }
        }
    }

    for(int i = 0; i < iter; i++)
    {
        float courant = couranti;
        if(i > iter - 10)
        {
            float in = i-(iter-10);
            courant = couranti * ((1.0- in/10.0)*1.0) + (in/10.0) * 0.05 * couranti;
            ////FF::cout << "iter " << i << " " << courant << std::endl;
        }

        if(i%2 == 0)
        {
            for(int r = DEM->data.nr_rows()-1; r >-1 ;r--)
            {
                for(int c = DEM->data.nr_cols() -1; c > -1 ;c--)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {


                        float dirxh = mapx->data[r][c]/ std::max(0.0001f,map->data[r][c]);
                        float diryh = mapy->data[r][c]/ std::max(0.0001f,map->data[r][c]);

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


                        //oh no, where to apply the preferential direction?

                        float dot = 1.0-((sx * dirxh + sy * diryh)/(std::sqrt(sx*sx + sy*sy)*std::sqrt(dirxh * dirxh + diryh *diryh)));
                        sx = sx;// * (1.0-std::fabs(dirxh));
                        sy = sy;// * (1.0-std::fabs(diryh));

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

                        qx = std::min(courant*h,std::min(qmaxx,std::abs(courant*h * vx  *dot)));
                        qy = std::min(courant*h,std::min(qmaxy,std::abs(courant*h * vy *dot )));

                        //
                        //change height

                        hn = h - qx - qy;
                        mapx2->data[r][c] -= (qx+qy) * dirxh;
                        mapy2->data[r][c] -= (qx+qy) * diryh;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy);

                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                                mapx2->data[r][c-1] = mapx->data[r][c-1] + qx * dirxh;
                                mapy2->data[r][c-1] = mapy->data[r][c-1] + qx * diryh;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                                mapx2->data[r][c+1] = mapx->data[r][c+1] + qx * dirxh;
                                mapy2->data[r][c+1] = mapy->data[r][c+1] + qx * diryh;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                                mapx2->data[r-1][c] = mapx->data[r-1][c] + qy * dirxh;
                                mapy2->data[r-1][c] = mapy->data[r-1][c] + qy * diryh;
                            }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                                mapx2->data[r+1][c] = mapx->data[r+1][c] + qy * dirxh;
                                mapy2->data[r+1][c] = mapy->data[r+1][c] + qy * diryh;
                            }
                        }

                    }
                }
            }

        }else {
            for(int r = 0; r < DEM->data.nr_rows();r++)
            {
                for(int c = 0; c < DEM->data.nr_cols();c++)
                {

                    if(!pcr::isMV(DEM->data[r][c]))
                    {

                        float dirxh = mapx->data[r][c]/ std::max(0.0001f,map->data[r][c]);
                        float diryh = mapy->data[r][c]/ std::max(0.0001f,map->data[r][c]);


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

                        float dot = 1.0- ((sx * dirxh + sy * diryh)/(std::sqrt(sx*sx + sy*sy)*std::sqrt(dirxh * dirxh + diryh *diryh)));
                        sx = sx;// * (1.0-std::fabs(dirxh));
                        sy = sy;// * (1.0-std::fabs(diryh));


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

                        qx = std::min(courant*h,std::min(qmaxx,std::abs(courant*h * vx*dot)));
                        qy = std::min(courant*h,std::min(qmaxy,std::abs(courant*h * vy*dot)));
                        //change height

                        hn = h - qx - qy;
                        mapx2->data[r][c] -= (qx+qy) * dirxh;
                        mapy2->data[r][c] -= (qx+qy) * diryh;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy);

                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                                mapx2->data[r][c-1] = mapx->data[r][c-1] + qx * dirxh;
                                mapy2->data[r][c-1] = mapy->data[r][c-1] + qx * diryh;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                                mapx2->data[r][c+1] = mapx->data[r][c+1] + qx * dirxh;
                                mapy2->data[r][c+1] = mapy->data[r][c+1] + qx * diryh;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                                mapx2->data[r-1][c] = mapx->data[r-1][c] + qy * dirxh;
                                mapy2->data[r-1][c] = mapy->data[r-1][c] + qy * diryh;
                            }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                                mapx2->data[r+1][c] = mapx->data[r+1][c] + qy * dirxh;
                                mapy2->data[r+1][c] = mapy->data[r+1][c] + qy * diryh;
                            }
                        }
                    }
                }
            }
        }

        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];
                mapx->data[r][c] = mapx2->data[r][c];
                mapy->data[r][c] = mapy2->data[r][c];
            }
        }

    }

    delete map2;
    delete mapx2;
    delete mapy2;
    delete mapx;
    delete mapy;

    return map;
}

inline cTMapf * DiffusiveMaxWavePrefDirDeposition(cTMapf * DEM,cTMapf * H, cTMapf * DirX, cTMapf * DirY, int iter, float couranti, cTMapf * Deposition_in, cTMapf * Deposition_out)
{

    ////FF::cout << "max wave prefdir " << DEM << " " << H << " " << DirX << " " << DirY << " dx" <<DEM->data.cell_size() << " dy" << DEM->data.cell_sizeY() << std::endl;

    cTMapf *map = DEM->GetCopy0();


    MaskedRaster<float> raster_data2(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
    cTMapf *map2 = DEM->GetCopy0();
    cTMapf *mapx = map->GetCopy0();
    cTMapf *mapx2 = map->GetCopy0();
    cTMapf *mapy = map->GetCopy0();
    cTMapf *mapy2 = map->GetCopy0();



    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

#pragma omp parallel for collapse(2)
    for(int r = 0; r < DEM->data.nr_rows();r++)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                map->data[r][c] = (!pcr::isMV(H->data[r][c]))? H->Drc : 0.0;
                map2->data[r][c] = map->data[r][c];
                mapx->data[r][c] =  map->data[r][c] * DirX->data[r][c];
                mapy->data[r][c] =  map->data[r][c] * DirY->data[r][c];
                Deposition_out->data[r][c] = Deposition_in->data[r][c];
            }
        }
    }

    for(int i = 0; i < iter; i++)
    {
        float courant = couranti;
        if(i > iter - 10)
        {
            float in = i-(iter-10);
            courant = couranti * ((1.0- in/10.0)*1.0) + (in/10.0) * 0.05 * couranti;
            ////FF::cout << "iter " << i << " " << courant << std::endl;
        }

        if(i%2 == 0)
        {
            for(int r = DEM->data.nr_rows()-1; r >-1 ;r--)
            {
                for(int c = DEM->data.nr_cols() -1; c > -1 ;c--)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {


                        float dirxh = mapx->data[r][c]/ std::max(0.0001f,map->data[r][c]);
                        float diryh = mapy->data[r][c]/ std::max(0.0001f,map->data[r][c]);

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


                        //oh no, where to apply the preferential direction?

                        float dot = 1.0-((sx * dirxh + sy * diryh)/(std::sqrt(sx*sx + sy*sy)*std::sqrt(dirxh * dirxh + diryh *diryh)));
                        sx = sx;// * (1.0-std::fabs(dirxh));
                        sy = sy;// * (1.0-std::fabs(diryh));

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

                        qx = std::min(courant*h,std::min(qmaxx,std::abs(courant*h * vx  *dot)));
                        qy = std::min(courant*h,std::min(qmaxy,std::abs(courant*h * vy *dot )));

                        //
                        //change height

                        hn = h - qx - qy;
                        mapx2->data[r][c] -= (qx+qy) * dirxh;
                        mapy2->data[r][c] -= (qx+qy) * diryh;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy);

                        float frach = std::max(0.0f,std::min(1.0f,Deposition_in->data[r][c]/map->data[r][c]));
                        float qdepx = qx * frach;
                        float qdepy = qy * frach;
                        Deposition_out->data[r][c] = std::max(0.0f,Deposition_out->data[r][c] - qdepx - qdepy); 

                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                                mapx2->data[r][c-1] = mapx->data[r][c-1] + qx * dirxh;
                                mapy2->data[r][c-1] = mapy->data[r][c-1] + qx * diryh;
                                Deposition_out->data[r][c-1] = Deposition_out->data[r][c-1] + qdepx;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                                mapx2->data[r][c+1] = mapx->data[r][c+1] + qx * dirxh;
                                mapy2->data[r][c+1] = mapy->data[r][c+1] + qx * diryh;
                                Deposition_out->data[r][c+1] = Deposition_out->data[r][c+1] + qdepx;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                                mapx2->data[r-1][c] = mapx->data[r-1][c] + qy * dirxh;
                                mapy2->data[r-1][c] = mapy->data[r-1][c] + qy * diryh;
                                Deposition_out->data[r-1][c] = Deposition_out->data[r-1][c] + qdepy;
                            }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                                mapx2->data[r+1][c] = mapx->data[r+1][c] + qy * dirxh;
                                mapy2->data[r+1][c] = mapy->data[r+1][c] + qy * diryh;
                                Deposition_out->data[r+1][c] = Deposition_out->data[r+1][c] + qdepy;
                            }
                        }

                    }
                }
            }

        }else {
            for(int r = 0; r < DEM->data.nr_rows();r++)
            {
                for(int c = 0; c < DEM->data.nr_cols();c++)
                {

                    if(!pcr::isMV(DEM->data[r][c]))
                    {

                        float dirxh = mapx->data[r][c]/ std::max(0.0001f,map->data[r][c]);
                        float diryh = mapy->data[r][c]/ std::max(0.0001f,map->data[r][c]);


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

                        float dot = 1.0- ((sx * dirxh + sy * diryh)/(std::sqrt(sx*sx + sy*sy)*std::sqrt(dirxh * dirxh + diryh *diryh)));
                        sx = sx;// * (1.0-std::fabs(dirxh));
                        sy = sy;// * (1.0-std::fabs(diryh));


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

                        qx = std::min(courant*h,std::min(qmaxx,std::abs(courant*h * vx*dot)));
                        qy = std::min(courant*h,std::min(qmaxy,std::abs(courant*h * vy*dot)));
                        //change height

                        hn = h - qx - qy;
                        mapx2->data[r][c] -= (qx+qy) * dirxh;
                        mapy2->data[r][c] -= (qx+qy) * diryh;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy);

                        float frach = std::max(0.0f,std::min(1.0f,Deposition_in->data[r][c]/map->data[r][c]));
                        float qdepx = qx * frach;
                        float qdepy = qy * frach;
                        Deposition_out->data[r][c] = std::max(0.0f,Deposition_out->data[r][c] - qdepx - qdepy); 

                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                                mapx2->data[r][c-1] = mapx->data[r][c-1] + qx * dirxh;
                                mapy2->data[r][c-1] = mapy->data[r][c-1] + qx * diryh;
                                Deposition_out->data[r][c-1] = Deposition_out->data[r][c-1] + qdepx;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                                mapx2->data[r][c+1] = mapx->data[r][c+1] + qx * dirxh;
                                mapy2->data[r][c+1] = mapy->data[r][c+1] + qx * diryh;
                                Deposition_out->data[r][c+1] = Deposition_out->data[r][c+1] + qdepx;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                                mapx2->data[r-1][c] = mapx->data[r-1][c] + qy * dirxh;
                                mapy2->data[r-1][c] = mapy->data[r-1][c] + qy * diryh;
                                Deposition_out->data[r-1][c] = Deposition_out->data[r-1][c] + qdepy;
                            }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                                mapx2->data[r+1][c] = mapx->data[r+1][c] + qy * dirxh;
                                mapy2->data[r+1][c] = mapy->data[r+1][c] + qy * diryh;
                                Deposition_out->data[r+1][c] = Deposition_out->data[r+1][c] + qdepy;
                            }
                        }
                    }
                }
            }
        }

        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];
                mapx->data[r][c] = mapx2->data[r][c];
                mapy->data[r][c] = mapy2->data[r][c];

                Deposition_in->data[r][c] = Deposition_out->data[r][c];
            }
        }

    }

    delete map2;
    delete mapx2;
    delete mapy2;
    delete mapx;
    delete mapy;

    return map;
}



inline cTMapf * DiffusiveMaxWaveCascade(cTMapf * DEM,cTMapf * DEMM,cTMapf * H, int iter, float couranti)
{


    //printf("start diffusive max wave\n");

    cTMapf *map = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
    ////FF::cout << "created map 1" << std::endl;
    cTMapf *map2 = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());

    

    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr,demm_xl,demm_xr,demm_yl,demm_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

    ////FF::cout << "get initial height " << std::endl;
    #pragma omp parallel for collapse(2)
    for(int r = 0; r < DEM->data.nr_rows();r++)
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
        int n_ease = 10;
        n_ease = std::min(10,(int)(iter/2.0));
        float n_easef = n_ease;
        if(i > iter - n_ease)
        {
            float in = i-(iter-n_easef);
            courant = couranti * ((1.0- in/n_easef)*1.0) + (in/n_easef) * 0.05 * couranti;
        }

        ////FF::cout << "diff max wave " << i << " " << courant << std::endl;
#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        if(i%2 == 0)
        {
            for(int r = DEM->data.nr_rows()-1; r >-1 ;r--)
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

                        if(c > 0)
                        {
                            demm_xl = pcr::isMV(DEMM->data[r][c-1])?dem:DEMM->data[r][c-1];
                        }else {
                            demm_xl = dem;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            demm_xr = pcr::isMV(DEMM->data[r][c+1])?dem:DEMM->data[r][c+1];
                        }else {
                            demm_xr = dem;
                        }
                        if(r > 0)
                        {
                            demm_yl = pcr::isMV(DEMM->data[r-1][c])?dem:DEMM->data[r-1][c];
                        }else {
                            demm_yl = dem;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            demm_yr = pcr::isMV(DEMM->data[r+1][c])?dem:DEMM->data[r+1][c];
                        }else {
                            demm_yr = dem;

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

                        qmax_xl = std::max(0.0f,dem+h - std::max(demm_xl, dem_xl +h_xl));
                        qmax_xr = std::max(0.0f,dem+h - std::max(demm_xr, dem_xr +h_xr));
                        qmax_yl = std::max(0.0f,dem+h - std::max(demm_yl, dem_yl +h_yl));
                        qmax_yr = std::max(0.0f,dem+h - std::max(demm_yr, dem_yr +h_yr));

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
            }

        }else {
            for(int r = 0; r < DEM->data.nr_rows();r++)
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

                         if(c > 0)
                        {
                            demm_xl = pcr::isMV(DEMM->data[r][c-1])?dem:DEMM->data[r][c-1];
                        }else {
                            demm_xl = dem;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            demm_xr = pcr::isMV(DEMM->data[r][c+1])?dem:DEMM->data[r][c+1];
                        }else {
                            demm_xr = dem;
                        }
                        if(r > 0)
                        {
                            demm_yl = pcr::isMV(DEMM->data[r-1][c])?dem:DEMM->data[r-1][c];
                        }else {
                            demm_yl = dem;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            demm_yr = pcr::isMV(DEMM->data[r+1][c])?dem:DEMM->data[r+1][c];
                        }else {
                            demm_yr = dem;

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

                        qmax_xl = std::max(0.0f,dem+h - std::max(demm_xl, dem_xl +h_xl));
                        qmax_xr = std::max(0.0f,dem+h - std::max(demm_xr, dem_xr +h_xr));
                        qmax_yl = std::max(0.0f,dem+h - std::max(demm_yl, dem_yl +h_yl));
                        qmax_yr = std::max(0.0f,dem+h - std::max(demm_yr, dem_yr +h_yr));

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
            }


        }

        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];
            }
        }

    }

    delete map2;

    return map;
}


inline cTMapf * DiffusiveMaxWave(cTMapf * DEM,cTMapf * H, int iter, float couranti)
{


    //printf("start diffusive max wave\n");

    //MaskedRaster<float> raster_data(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


    //MaskedRaster<float> raster_data2(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map2 = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

    const float coeff_boundary = 0.5;
#pragma omp parallel for collapse(2)
    for(int r = 0; r < DEM->data.nr_rows();r++)
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
        if(i > iter - 10)
        {
            float in = i-(iter-10);
            courant = couranti * ((1.0- in/10.0)*1.0) + (in/10.0) * 0.05 * couranti;
        }
#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        if(i%2 == 0)
        {
            for(int r = DEM->data.nr_rows()-1; r >-1 ;r--)
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
                            h_xl = pcr::isMV(map->data[r][c-1])? coeff_boundary* h:map->data[r][c-1];
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
            }

        }else {
            for(int r = 0; r < DEM->data.nr_rows();r++)
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
                            h_xl = pcr::isMV(map->data[r][c-1])?coeff_boundary*h:map->data[r][c-1];
                        }else {
                            h_xl = coeff_boundary*h;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            h_xr = pcr::isMV(map->data[r][c+1])?h:map->data[r][c+1];
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
            }


        }

        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];
            }
        }

    }

    delete map2;

    return map;
}


inline cTMapf * DiffusiveMaxWaveCulverts(cTMapf * DEM,cTMapf * H, int iter, float couranti,
                                         cTMapf * CULVERTID = nullptr,
                                         const std::vector<float>& CULVERT_Q = {},
                                         const std::vector<float>& CULVERT_FRAC = {},
                                         const std::vector<float>& CULVERT_QTOT = {},
                                         const std::vector<int>& CULVERT_RS = {},
                                         const std::vector<int>& CULVERT_CS = {},
                                         const std::vector<int>& CULVERT_RT = {},
                                         const std::vector<int>& CULVERT_CT = {},
                                         const std::vector<float>& culvert_ffactors = {})
{
    bool transfer_all = false;


    printf("start diffusive max wave with culverts\n");

    //MaskedRaster<float> raster_data(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


     //MaskedRaster<float> raster_data2(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map2 = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");

    std::vector<int> eff_culvert_id(DEM->nrRows() * DEM->nrCols(), -1);
    std::vector<float> eff_culvert_mult(DEM->nrRows() * DEM->nrCols(), 1.0f);

    for (int idc = 0; idc < CULVERT_RS.size(); idc++)
    {
        int r_src = CULVERT_RS[idc];
        int c_src = CULVERT_CS[idc];
        if (r_src >= 0 && r_src < DEM->nrRows() && c_src >= 0 && c_src < DEM->nrCols())
        {
            eff_culvert_id[r_src * DEM->nrCols() + c_src] = idc;
            eff_culvert_mult[r_src * DEM->nrCols() + c_src] = 1.0f;
        }
    }

    if (DEM->cellSizeX() < 50.0f)
    {
        int dr[4] = {-1, 1, 0, 0};
        int dc[4] = {0, 0, -1, 1};
        for (int idc = 0; idc < CULVERT_RS.size(); idc++)
        {
            int r_src = CULVERT_RS[idc];
            int c_src = CULVERT_CS[idc];
            if (r_src >= 0 && r_src < DEM->nrRows() && c_src >= 0 && c_src < DEM->nrCols())
            {
                for (int n = 0; n < 4; n++)
                {
                    int rn = r_src + dr[n];
                    int cn = c_src + dc[n];
                    if (rn >= 0 && rn < DEM->nrRows() && cn >= 0 && cn < DEM->nrCols())
                    {
                        // check if neighbor is target cell of any culvert
                        bool is_tgt = false;
                        for (int t = 0; t < CULVERT_RT.size(); t++)
                        {
                            if (rn == CULVERT_RT[t] && cn == CULVERT_CT[t])
                            {
                                is_tgt = true;
                                break;
                            }
                        }
                        if (is_tgt) continue;

                        int idx = rn * DEM->nrCols() + cn;
                        if (eff_culvert_id[idx] == -1)
                        {
                            eff_culvert_id[idx] = idc;
                            eff_culvert_mult[idx] = CULVERT_NEIGHBOR_FACTOR;
                        }
                    }
                }
            }
        }
    }


    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

    const float coeff_boundary = 0.5;
#pragma omp parallel for collapse(2)
    for(int r = 0; r < DEM->data.nr_rows();r++)
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
        if(i > iter - 10)
        {
            float in = i-(iter-10);
            courant = couranti * ((1.0- in/10.0)*1.0) + (in/10.0) * 0.05 * couranti;
        }
#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        if(i%2 == 0)
        {
            for(int r = DEM->data.nr_rows()-1; r >-1 ;r--)
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
                            h_xl = pcr::isMV(map->data[r][c-1])? coeff_boundary* h:map->data[r][c-1];
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

                        // Calculate culvert flow
                        float q_c = 0.0f;
                        int idc = eff_culvert_id[r * DEM->nrCols() + c];
                        int r_tgt = -1;
                        int c_tgt = -1;

                        if (idc >= 0)
                        {
                            r_tgt = CULVERT_RT[idc];
                            c_tgt = CULVERT_CT[idc];

                            if (r_tgt >= 0 && r_tgt < DEM->nrRows() && c_tgt >= 0 && c_tgt < DEM->nrCols())
                            {
                                float h_avail = map2->data[r][c];
                                float mult = eff_culvert_mult[r * DEM->nrCols() + c];
                                if (transfer_all)
                                {
                                    q_c = h_avail * mult;
                                }
                                else
                                {
                                    float actual_q = CULVERT_QTOT[idc];
                                    float max_q = CULVERT_Q[idc] * (idc < culvert_ffactors.size() ? culvert_ffactors[idc] : 1.0f);
                                    float frac_water = 1.0f;
                                    if (actual_q > max_q && actual_q > 1e-6f)
                                    {
                                        frac_water = max_q / actual_q;
                                    }
                                    
                                    q_c = h_avail * frac_water * mult;
                                }
                            }
                        }

                        // scale fluxes if they exceed available height in map2
                        float h_avail = map2->data[r][c];
                        float h_after_culv = std::max(0.0f, h_avail - q_c);

                        float q_diff_total = qx + qy;
                        if (q_diff_total > h_after_culv && q_diff_total > 1e-12f)
                        {
                            float factor = h_after_culv / q_diff_total;
                            qx *= factor;
                            qy *= factor;
                        }

                        //change height

                        hn = h - qx - qy;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy - q_c);

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

                        if (q_c > 0.0f)
                        {
                            map2->data[r_tgt][c_tgt] = map2->data[r_tgt][c_tgt] + q_c;
                        }

                    }
                }
            }

        }else {
            for(int r = 0; r < DEM->data.nr_rows();r++)
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
                            h_xl = pcr::isMV(map->data[r][c-1])?coeff_boundary*h:map->data[r][c-1];
                        }else {
                            h_xl = coeff_boundary*h;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            h_xr = pcr::isMV(map->data[r][c+1])?h:map->data[r][c+1];
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

                        // Calculate culvert flow
                        float q_c = 0.0f;
                        int idc = eff_culvert_id[r * DEM->nrCols() + c];
                        int r_tgt = -1;
                        int c_tgt = -1;

                        if (idc >= 0)
                        {
                            r_tgt = CULVERT_RT[idc];
                            c_tgt = CULVERT_CT[idc];

                            if (r_tgt >= 0 && r_tgt < DEM->nrRows() && c_tgt >= 0 && c_tgt < DEM->nrCols())
                            {
                                float h_avail = map2->data[r][c];
                                float mult = eff_culvert_mult[r * DEM->nrCols() + c];
                                if (transfer_all)
                                {
                                    q_c = h_avail * mult;
                                }
                                else
                                {
                                    float actual_q = CULVERT_QTOT[idc];
                                    float max_q = CULVERT_Q[idc] * (idc < culvert_ffactors.size() ? culvert_ffactors[idc] : 1.0f);
                                    float frac_water = 1.0f;
                                    if (actual_q > max_q && actual_q > 1e-6f)
                                    {
                                        frac_water = max_q / actual_q;
                                    }
                                    
                                    q_c = h_avail * frac_water * mult;
                                }
                            }
                        }

                        // scale fluxes if they exceed available height in map2
                        float h_avail = map2->data[r][c];
                        float h_after_culv = std::max(0.0f, h_avail - q_c);

                        float q_diff_total = qx + qy;
                        if (q_diff_total > h_after_culv && q_diff_total > 1e-12f)
                        {
                            float factor = h_after_culv / q_diff_total;
                            qx *= factor;
                            qy *= factor;
                        }

                        //change height

                        hn = h - qx - qy;
                        map2->data[r][c] = std::max(0.0f,map2->data[r][c] - qx - qy - q_c);

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

                        if (q_c > 0.0f)
                        {
                            map2->data[r_tgt][c_tgt] = map2->data[r_tgt][c_tgt] + q_c;
                        }

                    }
                }
            }


        }

        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];

                //if (CULVERTID != nullptr && (int)CULVERTID->data[r][c] >= 0)
                //{
                //    map->data[r][c] = 0.0;
                //}
            }
        }

    }

    delete map2;

    return map;
}

inline cTMapf * DiffusiveQMaxWave(cTMapf * DEM,cTMapf * H, int iter, float couranti, float power = 0.6)
{

	printf("start diffusive max q wave\n");


    float power2 = 1.0/std::max(1e-12f,power);

   

    //MaskedRaster<float> raster_data(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
    cTMapf *map = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


    //MaskedRaster<float> raster_data2(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map2 = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

    float total = 0.0;
    for(int r = 0; r < DEM->data.nr_rows();r++)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                map->data[r][c] = (!pcr::isMV(H->data[r][c]))? std::pow(H->Drc,power2) : 0.0;
                map2->data[r][c] = map->data[r][c];
                total += map->data[r][c];
            }
        }
    }
#ifdef FASTFLOOD_WASM
    emscripten_run_script("increment_progress_fastflood()");
#endif
    for(int i = 0; i < iter; i++)
    {
        float courant = couranti;
        if(i > iter - 10)
        {
            float in = i-(iter-10);
            courant = couranti * ((1.0- in/10.0)*1.0) + (in/10.0) * 0.05 * couranti;

        }
         //FF::cout << "iter " << i << " " << courant << std::endl;


        if(true)//i%2 == 0)
        {
            for(int r = DEM->data.nr_rows()-1; r >-1 ;r--)
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

                        float h_org = h;
                        float h_org_xl = h_xl;
                        float h_org_xr = h_xr;
                        float h_org_yl = h_yl;
                        float h_org_yr = h_yr;

                        h = std::pow(h,power);
                        h_xl = std::pow(h_xl,power);
                        h_xr = std::pow(h_xr,power);
                        h_yl = std::pow(h_yl,power);
                        h_yr = std::pow(h_yr,power);

                        sx = (h_xr +dem_xr)-(h_xl+dem_xl);
                        sy = (h_yr +dem_yr)-(h_yl+dem_yl);


                        //velocity
                        sl = std::sqrt(sx*sx + sy*sy);
                        vx = sx/sl;
                        vy = sy/sl;
                        signx = vx > 0.0? -1.0:1.0;
                        signy = vy > 0.0? -1.0:1.0;

                        //how much can we transfer so that the height in the neighbor is the same as ours?


                        //flux limiter

                        qmax_xl = std::max(0.0f,std::pow(std::max(0.0f,dem+h - dem_xl),power2) -h_org_xl);
                        qmax_xr = std::max(0.0f,std::pow(std::max(0.0f,dem+h - dem_xr),power2) -h_org_xr);
                        qmax_yl = std::max(0.0f,std::pow(std::max(0.0f,dem+h - dem_yl),power2) -h_org_yl);
                        qmax_yr = std::max(0.0f,std::pow(std::max(0.0f,dem+h - dem_yr),power2) -h_org_yr);

                        //discharge

                        qmaxx = signx > 0.0? qmax_xr : qmax_xl;
                        qmaxy = signy > 0.0? qmax_yr : qmax_yl;

                        qx = std::min(courant*h_org,std::min(0.5f*qmaxx,std::abs(courant*h_org * vx)));
                        qy = std::min(courant*h_org,std::min(0.5f*qmaxy,std::abs(courant*h_org * vy)));
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
            }

        }else {
            for(int r = 0; r < DEM->data.nr_rows();r++)
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

                        h = std::pow(h,power);
                        h_xl = std::pow(h_xl,power);
                        h_xr = std::pow(h_xr,power);
                        h_yl = std::pow(h_yl,power);
                        h_yr = std::pow(h_yr,power);




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
            }


        }

        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];
            }
        }

    }

    delete map2;

    total = 0.0;
    for(int r = 0; r < DEM->data.nr_rows();r++)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(map->data[r][c]))
            {
                total += map->data[r][c];
                map->data[r][c] = std::pow(map->data[r][c],power);
            }

        }
    }

    //FF::cout << "total2 " << total <<  std::endl;


    return map;
}




inline cTMapf * DiffusiveMaxWaveDeposition(cTMapf * DEM,cTMapf * H, int iter, float couranti, cTMapf * Deposition_in, cTMapf * Deposition_out)
{


	printf("start diffusive max wave\n");

    //MaskedRaster<float> raster_data(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


    //MaskedRaster<float> raster_data2(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map2 = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


    float dx = 1.0;
    float sl;
    float dem,dem_xl,dem_xr,dem_yl,dem_yr;
    float h, h_xl, h_xr, h_yl, h_yr;
    float hn, qmax_xr,qmax_xl,qmax_yl,qmax_yr,signx,signy,vx,vy,qx,qy,sx,sy;
    bool mv_xl,mv_xr,mv_yl,mv_yr;
    float qmaxx, qmaxy;

#pragma omp parallel for collapse(2)
    for(int r = 0; r < DEM->data.nr_rows();r++)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                map->data[r][c] = (!pcr::isMV(H->data[r][c]))? H->Drc : 0.0;
                map2->data[r][c] = map->data[r][c];
                Deposition_out->data[r][c] = Deposition_in->data[r][c];
            }
        }
    }

    for(int i = 0; i < iter; i++)
    {
        float courant = couranti;
        if(i > iter - 10)
        {
            float in = i-(iter-10);
            courant = couranti * ((1.0- in/10.0)*1.0) + (in/10.0) * 0.05 * couranti;
        }
#ifdef FASTFLOOD_WASM
		emscripten_run_script("increment_progress_fastflood()");
#endif
        if(i%2 == 0)
        {
            for(int r = DEM->data.nr_rows()-1; r >-1 ;r--)
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

                         float frach = std::max(0.0f,std::min(1.0f,Deposition_in->data[r][c]/map->data[r][c]));
                        float qdepx = qx * frach;
                        float qdepy = qy * frach;
                        Deposition_out->data[r][c] = std::max(0.0f,Deposition_out->data[r][c] - qdepx - qdepy); 


                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                                Deposition_out->data[r][c-1] = Deposition_out->data[r][c-1] + qdepx;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                             Deposition_out->data[r][c+1] = Deposition_out->data[r][c+1] + qdepx;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                            Deposition_out->data[r-1][c] = Deposition_out->data[r-1][c] + qdepy;
                           }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                            Deposition_out->data[r+1][c] = Deposition_out->data[r+1][c] + qdepy;
                            }
                        }

                    }
                }
            }

        }else {
            for(int r = 0; r < DEM->data.nr_rows();r++)
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

                         float frach = std::max(0.0f,std::min(1.0f,Deposition_in->data[r][c]/map->data[r][c]));
                        float qdepx = qx * frach;
                        float qdepy = qy * frach;
                        Deposition_out->data[r][c] = std::max(0.0f,Deposition_out->data[r][c] - qdepx - qdepy); 


                        if(c > 0)
                        {
                            if(signx< 0 && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                map2->data[r][c-1] = map2->data[r][c-1]+qx;
                                Deposition_out->data[r][c-1] = Deposition_out->data[r][c-1] + qdepx;
                            }
                        }
                        if(c < DEM->nrCols() -1)
                        {
                            if(signx > 0 && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                map2->data[r][c+1] = map2->data[r][c+1]+qx;
                             Deposition_out->data[r][c+1] = Deposition_out->data[r][c+1] + qdepx;
                            }
                        }

                        if(r > 0)
                        {
                            if(signy < 0 && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                map2->data[r-1][c] = map2->data[r-1][c]+qy;
                            Deposition_out->data[r-1][c] = Deposition_out->data[r-1][c] + qdepy;
                           }
                        }
                        if(r < DEM->nrRows() -1)
                        {
                            if(signy > 0 && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                map2->data[r+1][c] = map2->data[r+1][c]+qy;
                            Deposition_out->data[r+1][c] = Deposition_out->data[r+1][c] + qdepy;
                            }
                        }

                    }
                }
            }


        }

        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                map->data[r][c] = map2->data[r][c];

                Deposition_in->data[r][c] = Deposition_out->data[r][c];
            }
        }

    }

    delete map2;

    return map;
}





}

#endif
