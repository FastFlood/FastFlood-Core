#ifndef OCEANFILL_H
#define OCEANFILL_H

#include "../geospatial/map/map.h"
#include "../math/fastmath.h"



extern "C"
{
    
inline cTMapf *  do_fill_ocean_wave_decay(cTMapf * map /*dem*/,  cTMapf * N, float wavel, float freq, float level, int iter_max, float wave_decay_cal=1.35)
{

    auto wh_newf = [wave_decay_cal,level](float wh_1, float dem_1, float wh_2, float dem_2, float N_1, float N_2, float dx,float freqh, float wavelh) {

        float wh_nb =std::max(0.0f,wh_1) ;
        float wh_nb_wavel =std::min(wh_nb,std::max(0.0f,wh_1 + dem_1 - level));
        float wh_nb_nwavel = std::max(0.0f,wh_nb - wh_nb_wavel);
        float wh_nb_decay = wh_nb_nwavel  + wh_nb_wavel * fastexp(- dx * (wave_decay_cal * (0.5*(N_1+N_2)) * freqh / std::max(0.01f,std::max( std::max(wh_1,wh_2),wavelh))));
        //should probably be:
        //float wh_nb_decay = wh_nb_nwavel  + wh_nb_wavel * fastexp(- dx * (wave_decay_cal * (0.5*(N_1+N_2)) * freqh / std::max(0.01f,std::max(0.1f* std::max(wh_1,wh_2),wh_nb_wavel))));
        float available_here = std::max(0.0f,dem_1 + wh_nb_decay - dem_2);
        return available_here;
    };

    FF_INFO("Ocean fill with decay " );
    FF_INFO("Start ocean fill wave decay " + to_string<int>(iter_max) + " " + to_string<float>(level)+" " + to_string<float>(wavel) +" ");

    cTMapf * WH = map->GetCopy0();

    //sweeping ocean fill

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
            for(int r = 1; r < map->data.nr_rows();r++)
            {
                for(int c = 1; c < map->data.nr_cols();c++)
                {
                     if(!pcr::isMV(map->data[r][c]))
                    {
                       {
                            float wh = WH->data[r-1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r-1][c];
                            if(pcr::isMV(map->data[r-1][c]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r-1][c],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c-1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c-1];
                            if(pcr::isMV(map->data[r][c-1]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r][c-1],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                    }

                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-2; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-2; c > -1 ;c--)
                {
                     if(!pcr::isMV(map->data[r][c]))
                    {
                        {
                            float wh = WH->data[r+1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r+1][c];
                            if(pcr::isMV(map->data[r+1][c]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r+1][c],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c+1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c+1];
                            if(pcr::isMV(map->data[r][c+1]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r][c+1],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                    }


                }
            }
        }else
        {

            //first we move in right-lower direction
            for(int r = 1; r < map->data.nr_rows();r++)
            {
                for(int c = map->data.nr_cols()-2; c > -1 ;c--)
                {
                     if(!pcr::isMV(map->data[r][c]))
                    {
                        {
                            float wh = WH->data[r-1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r-1][c];
                            if(pcr::isMV(map->data[r-1][c]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r-1][c],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c+1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c+1];
                            if(pcr::isMV(map->data[r][c+1]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r][c+1],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-2; r > -1 ;r--)
            {
                for(int c = 1; c < map->data.nr_cols();c++)
                {

                     if(!pcr::isMV(map->data[r][c]))
                    {
                        {
                            float wh = WH->data[r+1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r+1][c];
                            if(pcr::isMV(map->data[r+1][c]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r+1][c],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c-1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c-1];
                            if(pcr::isMV(map->data[r][c-1]))
                            {
                                wh = level +wavel + 10;
                                dem = -10.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = wh_newf(wh,dem,whh,demh,N->data[r][c-1],N->data[r][c],dx,freq,wavel);
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = std::min(std::max(0.0f,wavel + level -demh),std::max(WH->data[r][c],wh_new));
                                    change = true;
                                }
                            }
                        }
                    }

                }
            }

        }
        first = false;

    }

    return WH;

}
inline cTMapf *  do_fill_ocean(cTMapf * map /*dem*/, float level, int iter_max)
{

    cTMapf * WH = map->GetCopy0();

    //sweeping ocean fill

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
            for(int r = 1; r < map->data.nr_rows();r++)
            {
                for(int c = 1; c < map->data.nr_cols();c++)
                {
                     if(!pcr::isMV(map->data[r][c]))
                    {
                       {
                            float wh = WH->data[r-1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r-1][c];
                            if(pcr::isMV(map->data[r-1][c]))
                            {
                                wh = 1.0;
                                dem = level - 1.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c-1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c-1];
                            if(pcr::isMV(map->data[r][c-1]))
                            {
                                wh = 1.0;
                                dem = level - 1.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        } 
                    }
                    
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-2; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-2; c > -1 ;c--)
                {
                     if(!pcr::isMV(map->data[r][c]))
                    {
                        {
                            float wh = WH->data[r+1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r+1][c];
                            if(pcr::isMV(map->data[r+1][c]))
                            {
                                wh = 1.0;
                                dem = level - 1.0;
                            }
                            float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                            if(wh > 0.0000000001f)
                            {
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c+1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c+1];
                            if(pcr::isMV(map->data[r][c+1]))
                            {
                                wh = 1.0;
                                dem = level - 1.0;
                            }
                            if(wh > 0.0000000001f)
                            {
                                float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        }
                    }


                }
            }
        }else
        {

            //first we move in right-lower direction
            for(int r = 1; r < map->data.nr_rows();r++)
            {
                for(int c = map->data.nr_cols()-2; c > -1 ;c--)
                {
                     if(!pcr::isMV(map->data[r][c]))
                    {
                        {
                            float wh = WH->data[r-1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r-1][c];
                            if(pcr::isMV(map->data[r-1][c]))
                            {
                                wh = level;
                                dem = demh;
                            }
                            float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                            if(wh > 0.0000000001f)
                            {
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c+1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c+1];
                            if(pcr::isMV(map->data[r][c+1]))
                            {
                                wh = level;
                                dem = demh;
                            }
                            float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                            if(wh > 0.0000000001f)
                            {
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        }
                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-2; r > -1 ;r--)
            {
                for(int c = 1; c < map->data.nr_cols();c++)
                {
                    
                     if(!pcr::isMV(map->data[r][c]))
                    {
                        {
                            float wh = WH->data[r+1][c];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r+1][c];
                            if(pcr::isMV(map->data[r+1][c]))
                            {
                                wh = level;
                                dem = demh;
                            }
                            float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                            if(wh > 0.0000000001f)
                            {
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        }
                        {
                            float wh = WH->data[r][c-1];
                            float whh= WH->data[r][c];
                            float demh= map->data[r][c];
                            float dem = map->data[r][c-1];
                            if(pcr::isMV(map->data[r][c-1]))
                            {
                                wh = level;
                                dem = demh;
                            }
                            float wh_new = std::min(std::max(0.0f,level-demh),std::max(0.0f,dem + wh - demh));
                            if(wh > 0.0000000001f)
                            {
                                if(wh_new != whh)
                                {
                                    WH->data[r][c] = wh_new;
                                    change = true;
                                }
                            }
                        }
                    }

                }
            }

        }
        first = false;

    }

    return WH;

}

}

#endif
