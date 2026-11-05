#ifndef VELCOR_H
#define VELCOR_H

extern "C"
{
inline cTMapf * WindowDistributeWRDir(cTMapf  *Other, cTMapf  *sx,cTMapf  *sy,cTMapf * WindowSizeL,cTMapf * WindowSizeR, cTMapf * WindowW ,float dx)
{

    cTMapf * map = Other->GetCopy0();
    cTMapf * mapw = Other->GetCopy0();

   for(int r = 0; r < Other->data.nr_rows();r++)
    {
        for(int c = 0; c < Other->data.nr_cols();c++)
        {
            if(pcr::isMV(Other->data[r][c]) )
            {
                pcr::setMV(map->data[r][c]);
            }else
            {

                float sizel = std::max(1e-12f,WindowSizeL->data[r][c]);
                float sizer = std::max(1e-12f,WindowSizeR->data[r][c]);
                
                int cells = std::min(100.f,std::max(1.0f,std::max(sizel,sizer)/dx));

                int rbegin = std::max(r-cells,0);
                int rend = std::min(r+cells,Other->nrRows()-1);
                int cbegin = std::max(c-cells,0);
                int cend = std::min(c+cells,Other->nrCols()-1);

                float val = Other->data[r][c];
                float w = WindowW->data[r][c];

                float sxh = sx->data[r][c];
                float syh = sy->data[r][c];

                float st = std::max(1e-12f,std::sqrt(sxh*sxh+syh*syh));
                sxh = sxh/st;
                syh = syh/st;


                for(int rn = rbegin; rn < rend + 1; rn++)
                {
                    for(int cn = cbegin; cn < cend + 1; cn++)
                    {
                        if(!pcr::isMV(Other->data[rn][cn]))
                        {
                           //directional weight based on dot product
                           float distc = std::sqrt(((float)(rn-r))*((float)(rn-r))+((float)(cn -c))*((float)(cn-c)));
                           float dirwl = std::max(0.0f, 0.5f +  (-syh * ((float)(cn -c)) + sxh * ((float)(rn -r)))/distc);
                           float dirwr = std::max(0.0f, 0.5f +  (syh * ((float)(cn -c)) - sxh * ((float)(rn -r)))/distc);
                           float dist = distc*dx;
                           float distf_l = std::max(0.0f,1.0f-dist/sizel);
                           float distf_r = std::max(0.0f,1.0f-dist/sizer);
                           float dirw = dirwl * distf_l + dirwr * distf_r;
                           map->data[rn][cn] +=val*w*dirw;
                           mapw->data[rn][cn] += w*dirw;
                        }

                    }
                }
            }
        }
    }

    for(int r = 0; r < Other->data.nr_rows();r++)
    {
        for(int c = 0; c < Other->data.nr_cols();c++)
        {
            if(!pcr::isMV(Other->data[r][c]) )
            {
                map->data[r][c] = map->data[r][c]/std::max(1e-12f,mapw->data[r][c]);
            }
        }
    }
        
    delete mapw;
    return map;

}
inline cTMapf * WindowDistributeWR(cTMapf  *Other, cTMapf * WindowSize, cTMapf * WindowW ,float dx)
{

    cTMapf * map = Other->GetCopy0();
    cTMapf * mapw = Other->GetCopy0();

   for(int r = 0; r < Other->data.nr_rows();r++)
    {
        for(int c = 0; c < Other->data.nr_cols();c++)
        {
            if(pcr::isMV(Other->data[r][c]) )
            {
                pcr::setMV(map->data[r][c]);
            }else
            {
                int cells = std::min(100.f,std::max(1.0f, WindowSize->data[r][c]/dx));

                int rbegin = std::max(r-cells,0);
                int rend = std::min(r+cells,Other->nrRows()-1);
                int cbegin = std::max(c-cells,0);
                int cend = std::min(c+cells,Other->nrCols()-1);

                float val = Other->data[r][c];
                float w = WindowW->data[r][c];

                for(int rn = rbegin; rn < rend + 1; rn++)
                {
                    for(int cn = cbegin; cn < cend + 1; cn++)
                    {
                        if(!pcr::isMV(Other->data[rn][cn]))
                        {

                           map->data[rn][cn] +=val*w;
                           mapw->data[rn][cn] += w;
                        }

                    }
                }
            }
        }
    }

    for(int r = 0; r < Other->data.nr_rows();r++)
    {
        for(int c = 0; c < Other->data.nr_cols();c++)
        {
            if(!pcr::isMV(Other->data[r][c]) )
            {
                map->data[r][c] = map->data[r][c]/std::max(1e-12f,mapw->data[r][c]);
            }
        }
    }
        
    delete mapw;
    return map;

}
inline cTMapf * WindowAverageWR(cTMapf * Other,cTMapf * WindowSizepx, float dx)
{

    cTMapf * map = Other->GetCopy0();

    for(int r = 0; r < Other->data.nr_rows();r++)
    {
        for(int c = 0; c < Other->data.nr_cols();c++)
        {
            if(pcr::isMV(Other->data[r][c]) )
            {
                pcr::setMV(map->data[r][c]);
            }else
            {
                int cells = std::min(25.f,std::max(4.0f,WindowSizepx->data[r][c]/dx));

                int rbegin = std::max(r-cells,0);
                int rend = std::min(r+cells,Other->nrRows()-1);
                int cbegin = std::max(c-cells,0);
                int cend = std::min(c+cells,Other->nrCols()-1);

                float val = 0.0f;
                float w = 0.0f;

                for(int rn = rbegin; rn < rend + 1; rn++)
                {
                    for(int cn = cbegin; cn < cend + 1; cn++)
                    {
                        if(!pcr::isMV(Other->data[rn][cn]))
                        {

                            if(rn == r && cn == c)
                            {
                                w += WindowSizepx->data[rn][cn];
                                val +=WindowSizepx->data[rn][cn] * Other->data[rn][cn];
                            }else
                            {
                                float dist = (float)(std::min(std::abs(c - cn),std::abs(r - rn))) * Other->cellSize();
                                float wrncn = WindowSizepx->data[rn][cn];
                                w += wrncn;
                                val += wrncn * Other->data[rn][cn];
                            }
                        }

                    }
                }

                if(w > 0.0f)
                {
                    map->data[r][c] = val/w;
                }else {
                    pcr::setMV(map->data[r][c]);
                }
            }
        }
    }

    return map;


}

inline cTMapf * WindowAverageWR2(cTMapf * Other,cTMapf * WindowSizepx,cTMapf * WindowW, float dx)
{

    cTMapf * map = Other->GetCopy0();

    for(int r = 0; r < Other->data.nr_rows();r++)
    {
        for(int c = 0; c < Other->data.nr_cols();c++)
        {
            if(pcr::isMV(Other->data[r][c]) )
            {
                pcr::setMV(map->data[r][c]);
            }else
            {
                int cells = std::min(25.f,std::max(1.0f,WindowSizepx->data[r][c]/dx));

                int rbegin = std::max(r-cells,0);
                int rend = std::min(r+cells,Other->nrRows()-1);
                int cbegin = std::max(c-cells,0);
                int cend = std::min(c+cells,Other->nrCols()-1);

                float val = 0.0f;
                float w = 0.0f;

                for(int rn = rbegin; rn < rend + 1; rn++)
                {
                    for(int cn = cbegin; cn < cend + 1; cn++)
                    {
                        if(!pcr::isMV(Other->data[rn][cn]))
                        {

                            if(rn == r && cn == c)
                            {
                                w += WindowW->data[rn][cn];
                                val +=WindowW->data[rn][cn] * Other->data[rn][cn];
                            }else
                            {
                                float dist = (float)(std::min(std::abs(c - cn),std::abs(r - rn))) * Other->cellSize();
                                float wrncn = WindowW->data[rn][cn];
                                w += wrncn;
                                val += wrncn * Other->data[rn][cn];
                            }
                        }

                    }
                }

                if(w > 0.0f)
                {
                    map->data[r][c] = val/w;
                }else {
                    pcr::setMV(map->data[r][c]);
                }
            }
        }
    }

    return map;


}
inline cTMapf* EMSCRIPTEN_KEEPALIVE SmoothWHOscillations(cTMapf * DEM, cTMapf *wh)
{
    cTMapf * whn= wh->GetCopy();

    for(int i = 0; i < 50; i++)
    {
         
        //#pragma omp parallel for collapse(2)
            for (int r = 0; r < DEM->Rows(); r++) {
                for (int c = 0; c < DEM->Cols(); c++) {
                    if (!pcr::isMV(DEM->data[r][c])) {

                        float whh = wh->data[r][c];
                        float dem = DEM->data[r][c];
                        float wh_x1 = whh;
                        float wh_x2 = whh;
                        float wh_y1 = whh;
                        float wh_y2 = whh;

                        float dem_x1 = dem;
                        float dem_y1 = dem;
                        float dem_x2 = dem;
                        float dem_y2 = dem;

                        float qx1 = 0.0;
                        float qx2 = 0.0;
                        float qy1 = 0.0;
                        float qy2 = 0.0;

                        float courant = 0.01;
                        if(c > 0)
                        {
                            if(!pcr::isMV(DEM->data[r][c-1]))
                            {
                                wh_x1 = wh->data[r][c-1];
                                dem_x1 = DEM->data[r][c-1];
                                qx1 = std::min(courant*wh_x1,courant*std::max(dem_x1 + wh_x1 - dem - whh,0.0f));
                                whn->data[r][c-1] -= qx1;
                            }
                        }
                        if(c <  DEM->Cols()-1)
                        {
                            if(!pcr::isMV(DEM->data[r][c+1]))
                            {
                                wh_x2 = wh->data[r][c+1];
                                dem_x2 = DEM->data[r][c+1];
                                qx2 = std::min(courant*wh_x2,courant*std::max(dem_x2 + wh_x2 - dem - whh,0.0f));
                                whn->data[r][c+1] -= qx2;
                            }
                        }
                        if(r > 0)
                        {
                            if(!pcr::isMV(DEM->data[r-1][c]))
                            {
                                wh_y1 = wh->data[r-1][c];
                                dem_y1 = DEM->data[r-1][c];
                                qy1 = std::min(courant*wh_y1,courant*std::max(dem_y1 + wh_y1 - dem - whh,0.0f));
                                whn->data[r-1][c] -= qy1;
                            }
                        }
                        if(r < DEM->Rows()-1)
                        {
                            if(!pcr::isMV(DEM->data[r+1][c]))
                            {
                                wh_y2 = wh->data[r+1][c];
                                dem_y2 = DEM->data[r+1][c];
                                qy2 = std::min(courant*wh_y2,courant*std::max(dem_y2 + wh_y2 - dem - whh,0.0f));
                                whn->data[r+1][c] -= qy2;
                            }
                        }

                        whn->data[r][c] = whn->data[r][c] +qx1 + qx2 + qy1 + qy2;

                    }
                }
            }

            for (int r = 0; r < DEM->Rows(); r++) {
                for (int c = 0; c < DEM->Cols(); c++) {
                    wh->data[r][c] = whn->data[r][c];
                }
            }

    }

    return whn;

}
}

#endif
