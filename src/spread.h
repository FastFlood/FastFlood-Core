#ifndef SPREAD_H
#define SPREAD_H

#include "map.h"

#define FF_MULT_DELT 0.0000002
extern "C"
{


inline cTMapf * SpreadMonotonicReconstruct(cTMapf * DEM, float delta,int iter_max = 0)
{



    //MaskedRaster<float> raster_data(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
    cTMapf *map = DEM->GetCopy0();

    //initialize the map with friction values
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(DEM->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 1e31;
            }
        }
    }


    float dx = map->cellSize();

    //we keep iterating through this algorithm untill there is no change left to make
    bool change = true;
    bool first = true;

    int iter = 0;

    while(change && ((iter_max <= 0) || (iter_max > 0 && iter < iter_max)))
    {
        iter ++;
        change = false;

        if(iter%2 == 0)
        {
            //first we move in right-lower direction
            for(int r = 0; r < map->data.nr_rows();r++)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = DEM->data[r][c];

                    if(!pcr::isMV(v_points))
                    {
                        if((r-1 > -1))
                        {
                            float vn_points = DEM->data[r-1][c];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r-1][c] + (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]) ,DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }

                        if((c-1 > -1))
                        {
                            float vn_points = DEM->data[r][c-1];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r][c-1]+ (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]),DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }
                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = DEM->data[r][c];

                    if(!pcr::isMV(v_points))
                    {
                        if((r+1 < map->data.nr_rows()))
                        {
                            float vn_points = DEM->data[r+1][c];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r+1][c]+ (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]),DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }

                        if((c+1 < map->data.nr_cols()))
                        {
                            float vn_points = DEM->data[r][c+1];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r][c+1]+ (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]),DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }
                    }

                }

            }
        }else
            {
                //first we move in right-lower direction
                for(int r = 0; r < map->data.nr_rows();r++)
                {
                    for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                    {
                        float v_points = DEM->data[r][c];

                        if(!pcr::isMV(v_points))
                        {
                            if((r-1 > -1))
                            {
                                float vn_points = DEM->data[r-1][c];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;

                                        demh = std::max(map->data[r-1][c]+ (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]),DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }


                            if((c+1 < map->data.nr_cols()))
                            {
                                float vn_points = DEM->data[r][c+1];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;

                                        demh = std::max(map->data[r][c+1]+ (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]),DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }
                        }
                    }
                }


                //then we move in left-upper direction
                for(int r = map->data.nr_rows()-1; r > -1 ;r--)
                {
                    for(int c = 0; c < map->data.nr_cols();c++)
                    {
                        float v_points = DEM->data[r][c];

                        if(!pcr::isMV(v_points))
                        {
                            if((r+1 < map->data.nr_rows()))
                            {
                                float vn_points = DEM->data[r+1][c];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;

                                        demh = std::max(map->data[r+1][c]+ (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]),DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }

                            if((c-1 > -1))
                            {
                                float vn_points = DEM->data[r][c-1];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;


                                        demh = std::max(map->data[r][c-1]+ (float)std::max(0.0001,FF_MULT_DELT*(double)map->data[r][c]),DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }
                        }

                    }
            }


        }


        first = false;

    }

    return map;

}


inline cTMapf * SpreadMonotonicReconstructIncr(cTMapf * DEM, float delta, int iter_max = 0)
{
    cTMapf *map = DEM->GetCopy0();

    // initialize the map with friction values
    for(int r = 0; r < map->data.nr_rows(); r++)
    {
        for(int c = 0; c < map->data.nr_cols(); c++)
        {
            if(pcr::isMV(DEM->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            } else {
                map->data[r][c] = 1e31;
            }
        }
    }

    bool change = true;
    int iter = 0;

    while(change && ((iter_max <= 0) || (iter_max > 0 && iter < iter_max)))
    {
        iter++;
        change = false;

        if(iter % 2 == 0)
        {
            // --- First Pass: Right-Lower Direction (r++, c++) ---
            for(int r = 0; r < map->data.nr_rows(); r++)
            {
                for(int c = 0; c < map->data.nr_cols(); c++)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        // Check North neighbor (r-1)
                        if(r - 1 > -1)
                        {
                            if(pcr::isMV(DEM->data[r-1][c])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r-1][c];
                                float z_x1 = OUTORMV(DEM, r-1, c-1) ? z_h : DEM->data[r-1][c-1];
                                float z_x2 = OUTORMV(DEM, r-1, c+1) ? z_h : DEM->data[r-1][c+1];
                                float z_y1 = OUTORMV(DEM, r-2, c) ? z_h : DEM->data[r-2][c];
                                float z_y2 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sy1 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r-1][c] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT* (double)map->data[r][c]) + 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }

                        // Check West neighbor (c-1)
                        if(c - 1 > -1)
                        {
                            if(pcr::isMV(DEM->data[r][c-1])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r][c-1];
                                float z_x1 = OUTORMV(DEM, r, c-2) ? z_h : DEM->data[r][c-2];
                                float z_x2 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float z_y1 = OUTORMV(DEM, r-1, c-1) ? z_h : DEM->data[r-1][c-1];
                                float z_y2 = OUTORMV(DEM, r+1, c-1) ? z_h : DEM->data[r+1][c-1];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sx1 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r][c-1] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT * (double)map->data[r][c]) + 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }
                    }
                }
            }

            // --- Second Pass: Left-Upper Direction (r--, c--) ---
            for(int r = map->data.nr_rows() - 1; r > -1; r--)
            {
                for(int c = map->data.nr_cols() - 1; c > -1; c--)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        // Check South neighbor (r+1)
                        if(r + 1 < map->data.nr_rows())
                        {
                            if(pcr::isMV(DEM->data[r+1][c])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r+1][c];
                                float z_x1 = OUTORMV(DEM, r+1, c-1) ? z_h : DEM->data[r+1][c-1];
                                float z_x2 = OUTORMV(DEM, r+1, c+1) ? z_h : DEM->data[r+1][c+1];
                                float z_y1 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float z_y2 = OUTORMV(DEM, r+2, c) ? z_h : DEM->data[r+2][c];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sy2 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r+1][c] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT * (double)map->data[r][c]) + + 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }

                        // Check East neighbor (c+1)
                        if(c + 1 < map->data.nr_cols())
                        {
                            if(pcr::isMV(DEM->data[r][c+1])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r][c+1];
                                float z_x1 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float z_x2 = OUTORMV(DEM, r, c+2) ? z_h : DEM->data[r][c+2];
                                float z_y1 = OUTORMV(DEM, r-1, c+1) ? z_h : DEM->data[r-1][c+1];
                                float z_y2 = OUTORMV(DEM, r+1, c+1) ? z_h : DEM->data[r+1][c+1];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sx2 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r][c+1] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT * (double)map->data[r][c]) + 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }
                    }
                }
            }
        } else {
            // --- First Pass: Right-Upper Direction (r++, c--) ---
            for(int r = 0; r < map->data.nr_rows(); r++)
            {
                for(int c = map->data.nr_cols() - 1; c > -1; c--)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        // Check North neighbor (r-1)
                        if(r - 1 > -1)
                        {
                            if(pcr::isMV(DEM->data[r-1][c])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r-1][c];
                                float z_x1 = OUTORMV(DEM, r-1, c-1) ? z_h : DEM->data[r-1][c-1];
                                float z_x2 = OUTORMV(DEM, r-1, c+1) ? z_h : DEM->data[r-1][c+1];
                                float z_y1 = OUTORMV(DEM, r-2, c) ? z_h : DEM->data[r-2][c];
                                float z_y2 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sy1 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r-1][c] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT * (double)map->data[r][c]) + 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }

                        // Check East neighbor (c+1)
                        if(c + 1 < map->data.nr_cols())
                        {
                            if(pcr::isMV(DEM->data[r][c+1])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r][c+1];
                                float z_x1 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float z_x2 = OUTORMV(DEM, r, c+2) ? z_h : DEM->data[r][c+2];
                                float z_y1 = OUTORMV(DEM, r-1, c+1) ? z_h : DEM->data[r-1][c+1];
                                float z_y2 = OUTORMV(DEM, r+1, c+1) ? z_h : DEM->data[r+1][c+1];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sx2 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r][c+1] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT * (double)map->data[r][c]) + 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }
                    }
                }
            }

            // --- Second Pass: Left-Lower Direction (r--, c++) ---
            for(int r = map->data.nr_rows() - 1; r > -1; r--)
            {
                for(int c = 0; c < map->data.nr_cols(); c++)
                {
                    if(!pcr::isMV(DEM->data[r][c]))
                    {
                        // Check South neighbor (r+1)
                        if(r + 1 < map->data.nr_rows())
                        {
                            if(pcr::isMV(DEM->data[r+1][c])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r+1][c];
                                float z_x1 = OUTORMV(DEM, r+1, c-1) ? z_h : DEM->data[r+1][c-1];
                                float z_x2 = OUTORMV(DEM, r+1, c+1) ? z_h : DEM->data[r+1][c+1];
                                float z_y1 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float z_y2 = OUTORMV(DEM, r+2, c) ? z_h : DEM->data[r+2][c];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sy2 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r+1][c] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT * (double)map->data[r][c])+ 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }

                        // Check West neighbor (c-1)
                        if(c - 1 > -1)
                        {
                            if(pcr::isMV(DEM->data[r][c-1])) { map->data[r][c] = DEM->data[r][c]; }
                            else {
                                float z_h = DEM->data[r][c-1];
                                float z_x1 = OUTORMV(DEM, r, c-2) ? z_h : DEM->data[r][c-2];
                                float z_x2 = OUTORMV(DEM, r, c) ? z_h : DEM->data[r][c];
                                float z_y1 = OUTORMV(DEM, r-1, c-1) ? z_h : DEM->data[r-1][c-1];
                                float z_y2 = OUTORMV(DEM, r+1, c-1) ? z_h : DEM->data[r+1][c-1];
                                float sx1 = std::max(0.0f, z_x1 - z_h); float sx2 = std::max(0.0f, z_x2 - z_h);
                                float sy1 = std::max(0.0f, z_y1 - z_h); float sy2 = std::max(0.0f, z_y2 - z_h);
                                double sh = sx1 / std::max(0.00001f, sx1 + sx2 + sy1 + sy2);

                                float demh = std::max(map->data[r][c-1] + (float)std::max(0.0001+ 0.0001f * (float)sh, FF_MULT_DELT * (double)map->data[r][c]) + 0.0001f * (float)sh, DEM->data[r][c]);
                                float mapn = std::min(demh, map->data[r][c]);
                                if(mapn < map->data[r][c]) { change = true; map->data[r][c] = mapn; }
                            }
                        } else { map->data[r][c] = DEM->data[r][c]; }
                    }
                }
            }
        }
    }
    return map;
}



inline cTMapf * SpreadDepressionFind(cTMapf * DEM, int iter_max = 0)
{

    //MaskedRaster<float> raster_data(DEM->data.nr_rows(), DEM->data.nr_cols(), DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());
     cTMapf *map = new cTMapf(nullptr,DEM->data.nr_cols(),DEM->data.nr_rows(),  DEM->data.north(), DEM->data.west(), DEM->data.cell_size(),DEM->data.cell_sizeY());//std::move(raster_data),DEM->projection(),"");


    //initialize the map with friction values
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(DEM->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 1e31;
            }
        }
    }


    float dx = map->cellSize();

    //we keep iterating through this algorithm untill there is no change left to make
    bool change = true;
    bool first = true;

    int iter = 0;

    while(change && ((iter_max <= 0) || (iter_max > 0 && iter < iter_max)))
    {
        iter ++;
        change = false;

        if(iter%2 == 0)
        {
            //first we move in right-lower direction
            for(int r = 0; r < map->data.nr_rows();r++)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = DEM->data[r][c];

                    if(!pcr::isMV(v_points))
                    {
                        if((r-1 > -1))
                        {
                            float vn_points = DEM->data[r-1][c];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r-1][c],DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }

                        if((c-1 > -1))
                        {
                            float vn_points = DEM->data[r][c-1];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r][c-1],DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }
                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = DEM->data[r][c];

                    if(!pcr::isMV(v_points))
                    {
                        if((r+1 < map->data.nr_rows()))
                        {
                            float vn_points = DEM->data[r+1][c];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r+1][c],DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }

                        if((c+1 < map->data.nr_cols()))
                        {
                            float vn_points = DEM->data[r][c+1];
                            if(pcr::isMV(vn_points))
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }else
                            {
                                float demh;

                                    demh = std::max(map->data[r][c+1],DEM->data[r][c]);

                                float mapn = std::min(demh,map->data[r][c]);
                                if(mapn < map->data[r][c])
                                {
                                    change = true;
                                    map->data[r][c] = mapn;
                                }
                            }
                        }else
                        {
                            map->data[r][c] = DEM->data[r][c];
                        }
                    }

                }

            }
        }else
            {
                //first we move in right-lower direction
                for(int r = 0; r < map->data.nr_rows();r++)
                {
                    for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                    {
                        float v_points = DEM->data[r][c];

                        if(!pcr::isMV(v_points))
                        {
                            if((r-1 > -1))
                            {
                                float vn_points = DEM->data[r-1][c];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;

                                        demh = std::max(map->data[r-1][c],DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }


                            if((c+1 < map->data.nr_cols()))
                            {
                                float vn_points = DEM->data[r][c+1];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;

                                        demh = std::max(map->data[r][c+1],DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }
                        }
                    }
                }


                //then we move in left-upper direction
                for(int r = map->data.nr_rows()-1; r > -1 ;r--)
                {
                    for(int c = 0; c < map->data.nr_cols();c++)
                    {
                        float v_points = DEM->data[r][c];

                        if(!pcr::isMV(v_points))
                        {
                            if((r+1 < map->data.nr_rows()))
                            {
                                float vn_points = DEM->data[r+1][c];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;

                                        demh = std::max(map->data[r+1][c],DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }

                            if((c-1 > -1))
                            {
                                float vn_points = DEM->data[r][c-1];
                                if(pcr::isMV(vn_points))
                                {
                                    map->data[r][c] = DEM->data[r][c];
                                }else
                                {
                                    float demh;


                                        demh = std::max(map->data[r][c-1],DEM->data[r][c]);

                                    float mapn = std::min(demh,map->data[r][c]);
                                    if(mapn < map->data[r][c])
                                    {
                                        change = true;
                                        map->data[r][c] = mapn;
                                    }
                                }
                            }else
                            {
                                map->data[r][c] = DEM->data[r][c];
                            }
                        }

                    }
            }


        }


        first = false;

    }

    return map;

}



}

#endif
