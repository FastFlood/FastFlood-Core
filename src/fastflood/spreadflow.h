#ifndef FF_SPREADFLOWH
#define FF_SPREADFLOWH

#include "map.h"
extern "C"
{



inline double Combine_M2(double n_a,double avg_a, double M2_a, double n_b, double avg_b, double M2_b)
{
    double n = n_a + n_b;
    double delta = std::fabs(avg_b - avg_a);
    double M2 = M2_a + M2_b + std::fabs((delta * delta)/n) * std::fabs( n_a *n_b);

    return M2;
}


inline cTMapf * SpreadFlowMD2Variance(cTMapf * source ,cTMapf * source_acc, cTMapf * prop_avg, cTMapf * fracx1, cTMapf * fracx2, cTMapf * fracy1,cTMapf * fracy2, int iter_max = 0)
{
    //source is the source material,
    //prop are the properties of which we want to get the variance
    //
    /*MaskedRaster<float> raster_data(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *mx1 = new cTMapf(std::move(raster_data),source->projection(),"");
    MaskedRaster<float> raster_data2(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my1 = new cTMapf(std::move(raster_data2),source->projection(),"");
    MaskedRaster<float> raster_data3(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *mx2 = new cTMapf(std::move(raster_data3),source->projection(),"");
    MaskedRaster<float> raster_data4(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my2 = new cTMapf(std::move(raster_data4),source->projection(),"");
    */

    cTMapf *map = source->GetCopy0();


    cTMapf * M2 = source->GetCopy0();

    //initialize the map with friction values

   /* #pragma omp parallel for collapse(2)
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(source->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 0.0;
                mx1->data[r][c] = fracx1->data[r][c] * source->data[r][c];
                my1->data[r][c] = fracy1->data[r][c] * source->data[r][c];
                mx2->data[r][c] = fracx2->data[r][c] * source->data[r][c];
                my2->data[r][c] = fracy2->data[r][c] * source->data[r][c];

            }
        }
    }
*/
    float dx = map->cellSize();

    //we keep iterating through this algorithm untill there is no change left to make
    bool change = true;
    bool first = true;

    int iter = 0;

    while(change && ((iter_max <= 0) || (iter_max > 0 && iter < iter_max)))
    {
        iter ++;
        //change = false;

        if(iter%2 == 0)
        {

            //first we move in right-lower direction
            for(int r = 0; r < map->data.nr_rows();r++)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        //get all incoming weights

                        //get all incoming averages

                        //get all incoming M2 values

                        double w_x1 = 0.0;
                        double av_x1 = 0.0;
                        double m2_x1 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(source->data[r][c-1]))
                            {
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * source_acc->data[r][c-1];
                                av_x1 = prop_avg->data[r][c-1];
                                m2_x1 = std::fabs(fracx2->data[r][c-1]) * M2->data[r][c-1];
                            }
                        }

                        double w_x2 = 0.0;
                        double av_x2 = 0.0;
                        double m2_x2 = 0.0;
                        if(c < map->data.nr_cols() - 1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(source->data[r][c+1]))
                            {
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * source_acc->data[r][c+1];
                                av_x2 = prop_avg->data[r][c+1];
                                m2_x2 = std::fabs(fracx1->data[r][c+1]) * M2->data[r][c+1];
                            }
                        }

                        double w_y1 = 0.0;
                        double av_y1 = 0.0;
                        double m2_y1 = 0.0;
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10)  && !pcr::isMV(source->data[r-1][c]))
                            {
                                w_y1 = std::fabs(fracy2->data[r-1][c]) * source_acc->data[r-1][c];
                                av_y1 = prop_avg->data[r-1][c];
                                m2_y1 = std::fabs(fracy2->data[r-1][c]) * M2->data[r-1][c];
                            }
                        }

                        double w_y2 = 0.0;
                        double av_y2 = 0.0;
                        double m2_y2 = 0.0;
                        if(r < map->data.nr_rows() - 1)
                        {
                            if((std::fabs(fracy1->data[r+1][c])  > 1e-10) && !pcr::isMV(source->data[r+1][c]))
                            {
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * source_acc->data[r+1][c];
                                av_y2 = prop_avg->data[r+1][c];
                                m2_y2 = std::fabs(fracy1->data[r+1][c]) *M2->data[r+1][c];
                            }
                        }

                        double w_x = 0.0;
                        double av_x= 0.0;
                        double m2_x = 0.0;
                        double w_y = 0.0;
                        double av_y = 0.0;
                        double m2_y = 0.0;

                        //Combine_M2(double n_a,double avg_a, double M2_a, double n_b, double avg_b, double M2_b)

                        //M2->data[r][c] = m2_x1 + m2_x2 + m2_y1 + m2_y2 + 1.0;



                        if((w_x1 > 1e-10) && (w_x2 >  1e-10))
                        {
                            w_x = w_x1 + w_x2;
                            av_x = (w_x1 * av_x1 + w_x2 * av_x2)/std::max(1e-12,w_x);
                            m2_x = Combine_M2(w_x1,av_x1, m2_x1, w_x2,av_x2,m2_x2);

                        }else if (w_x2 >  1e-10)
                        {
                            w_x = w_x2;
                            av_x = av_x2;
                            m2_x = m2_x2;
                        }else if (w_x1 >  1e-10)
                        {
                            w_x = w_x1;
                            av_x = av_x1;
                            m2_x = m2_x1;
                        }

                        if((w_y2 >  1e-10) && (w_y1 >  1e-10))
                        {
                            w_y = w_y1 + w_y2;
                            av_y = (w_y1 * av_y1 + w_y2 * av_y2)/std::max(1e-12,w_y);
                            m2_y = Combine_M2(w_y1,av_y1, m2_y1, w_y2,av_y2,m2_y2);

                        }else if (w_y2 >  1e-10)
                        {
                            w_y = w_y2;
                            av_y = av_y2;
                            m2_y = m2_y2;
                        }else if (w_y1 >  1e-10)
                        {
                            w_y = w_y1;
                            av_y = av_y1;
                            m2_y = m2_y1;
                        }

                        if((w_x >  1e-10) && (w_y >  1e-10))
                        {
                            M2->data[r][c] = Combine_M2(w_y,av_y, m2_y, w_x,av_x,m2_x);

                        }else if(w_x > 1e-10 )
                        {
                            M2->data[r][c] = m2_x + 1.0;

                        }else if(w_y >  1e-10)
                        {
                            M2->data[r][c] = m2_y + 1.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {
                        //get all incoming weights

                        //get all incoming averages

                        //get all incoming M2 values

                        double w_x1 = 0.0;
                        double av_x1 = 0.0;
                        double m2_x1 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(source->data[r][c-1]))
                            {
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * source_acc->data[r][c-1];
                                av_x1 = prop_avg->data[r][c-1];
                                m2_x1 = std::fabs(fracx2->data[r][c-1]) * M2->data[r][c-1];
                            }
                        }

                        double w_x2 = 0.0;
                        double av_x2 = 0.0;
                        double m2_x2 = 0.0;
                        if(c < map->data.nr_cols() - 1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(source->data[r][c+1]))
                            {
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * source_acc->data[r][c+1];
                                av_x2 = prop_avg->data[r][c+1];
                                m2_x2 = std::fabs(fracx1->data[r][c+1]) * M2->data[r][c+1];
                            }
                        }

                        double w_y1 = 0.0;
                        double av_y1 = 0.0;
                        double m2_y1 = 0.0;
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10)  && !pcr::isMV(source->data[r-1][c]))
                            {
                                w_y1 = std::fabs(fracy2->data[r-1][c]) * source_acc->data[r-1][c];
                                av_y1 = prop_avg->data[r-1][c];
                                m2_y1 = std::fabs(fracy2->data[r-1][c]) * M2->data[r-1][c];
                            }
                        }

                        double w_y2 = 0.0;
                        double av_y2 = 0.0;
                        double m2_y2 = 0.0;
                        if(r < map->data.nr_rows() - 1)
                        {
                            if((std::fabs(fracy1->data[r+1][c])  > 1e-10) && !pcr::isMV(source->data[r+1][c]))
                            {
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * source_acc->data[r+1][c];
                                av_y2 = prop_avg->data[r+1][c];
                                m2_y2 = std::fabs(fracy1->data[r+1][c]) *M2->data[r+1][c];
                            }
                        }

                        double w_x = 0.0;
                        double av_x= 0.0;
                        double m2_x = 0.0;
                        double w_y = 0.0;
                        double av_y = 0.0;
                        double m2_y = 0.0;

                        //Combine_M2(double n_a,double avg_a, double M2_a, double n_b, double avg_b, double M2_b)

                        //M2->data[r][c] = m2_x1 + m2_x2 + m2_y1 + m2_y2 + 1.0;



                        if((w_x1 > 1e-10) && (w_x2 >  1e-10))
                        {
                            w_x = w_x1 + w_x2;
                            av_x = (w_x1 * av_x1 + w_x2 * av_x2)/std::max(1e-12,w_x);
                            m2_x = Combine_M2(w_x1,av_x1, m2_x1, w_x2,av_x2,m2_x2);

                        }else if (w_x2 >  1e-10)
                        {
                            w_x = w_x2;
                            av_x = av_x2;
                            m2_x = m2_x2;
                        }else if (w_x1 >  1e-10)
                        {
                            w_x = w_x1;
                            av_x = av_x1;
                            m2_x = m2_x1;
                        }

                        if((w_y2 >  1e-10) && (w_y1 >  1e-10))
                        {
                            w_y = w_y1 + w_y2;
                            av_y = (w_y1 * av_y1 + w_y2 * av_y2)/std::max(1e-12,w_y);
                            m2_y = Combine_M2(w_y1,av_y1, m2_y1, w_y2,av_y2,m2_y2);

                        }else if (w_y2 >  1e-10)
                        {
                            w_y = w_y2;
                            av_y = av_y2;
                            m2_y = m2_y2;
                        }else if (w_y1 >  1e-10)
                        {
                            w_y = w_y1;
                            av_y = av_y1;
                            m2_y = m2_y1;
                        }

                        if((w_x >  1e-10) && (w_y >  1e-10))
                        {
                            M2->data[r][c] = Combine_M2(w_y,av_y, m2_y, w_x,av_x,m2_x);

                        }else if(w_x > 1e-10 )
                        {
                            M2->data[r][c] = m2_x + 1.0;

                        }else if(w_y >  1e-10)
                        {
                            M2->data[r][c] = m2_y + 1.0;
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
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {
                        //get all incoming weights

                        //get all incoming averages

                        //get all incoming M2 values

                        double w_x1 = 0.0;
                        double av_x1 = 0.0;
                        double m2_x1 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(source->data[r][c-1]))
                            {
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * source_acc->data[r][c-1];
                                av_x1 = prop_avg->data[r][c-1];
                                m2_x1 = std::fabs(fracx2->data[r][c-1]) * M2->data[r][c-1];
                            }
                        }

                        double w_x2 = 0.0;
                        double av_x2 = 0.0;
                        double m2_x2 = 0.0;
                        if(c < map->data.nr_cols() - 1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(source->data[r][c+1]))
                            {
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * source_acc->data[r][c+1];
                                av_x2 = prop_avg->data[r][c+1];
                                m2_x2 = std::fabs(fracx1->data[r][c+1]) * M2->data[r][c+1];
                            }
                        }

                        double w_y1 = 0.0;
                        double av_y1 = 0.0;
                        double m2_y1 = 0.0;
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10)  && !pcr::isMV(source->data[r-1][c]))
                            {
                                w_y1 = std::fabs(fracy2->data[r-1][c]) * source_acc->data[r-1][c];
                                av_y1 = prop_avg->data[r-1][c];
                                m2_y1 = std::fabs(fracy2->data[r-1][c]) * M2->data[r-1][c];
                            }
                        }

                        double w_y2 = 0.0;
                        double av_y2 = 0.0;
                        double m2_y2 = 0.0;
                        if(r < map->data.nr_rows() - 1)
                        {
                            if((std::fabs(fracy1->data[r+1][c])  > 1e-10) && !pcr::isMV(source->data[r+1][c]))
                            {
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * source_acc->data[r+1][c];
                                av_y2 = prop_avg->data[r+1][c];
                                m2_y2 = std::fabs(fracy1->data[r+1][c]) *M2->data[r+1][c];
                            }
                        }

                        double w_x = 0.0;
                        double av_x= 0.0;
                        double m2_x = 0.0;
                        double w_y = 0.0;
                        double av_y = 0.0;
                        double m2_y = 0.0;

                        //Combine_M2(double n_a,double avg_a, double M2_a, double n_b, double avg_b, double M2_b)

                        //M2->data[r][c] = m2_x1 + m2_x2 + m2_y1 + m2_y2 + 1.0;



                        if((w_x1 > 1e-10) && (w_x2 >  1e-10))
                        {
                            w_x = w_x1 + w_x2;
                            av_x = (w_x1 * av_x1 + w_x2 * av_x2)/std::max(1e-12,w_x);
                            m2_x = Combine_M2(w_x1,av_x1, m2_x1, w_x2,av_x2,m2_x2);

                        }else if (w_x2 >  1e-10)
                        {
                            w_x = w_x2;
                            av_x = av_x2;
                            m2_x = m2_x2;
                        }else if (w_x1 >  1e-10)
                        {
                            w_x = w_x1;
                            av_x = av_x1;
                            m2_x = m2_x1;
                        }

                        if((w_y2 >  1e-10) && (w_y1 >  1e-10))
                        {
                            w_y = w_y1 + w_y2;
                            av_y = (w_y1 * av_y1 + w_y2 * av_y2)/std::max(1e-12,w_y);
                            m2_y = Combine_M2(w_y1,av_y1, m2_y1, w_y2,av_y2,m2_y2);

                        }else if (w_y2 >  1e-10)
                        {
                            w_y = w_y2;
                            av_y = av_y2;
                            m2_y = m2_y2;
                        }else if (w_y1 >  1e-10)
                        {
                            w_y = w_y1;
                            av_y = av_y1;
                            m2_y = m2_y1;
                        }

                        if((w_x >  1e-10) && (w_y >  1e-10))
                        {
                            M2->data[r][c] = Combine_M2(w_y,av_y, m2_y, w_x,av_x,m2_x);

                        }else if(w_x > 1e-10 )
                        {
                            M2->data[r][c] = m2_x + 1.0;

                        }else if(w_y >  1e-10)
                        {
                            M2->data[r][c] = m2_y + 1.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {
                        //get all incoming weights

                        //get all incoming averages

                        //get all incoming M2 values

                        double w_x1 = 0.0;
                        double av_x1 = 0.0;
                        double m2_x1 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(source->data[r][c-1]))
                            {
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * source_acc->data[r][c-1];
                                av_x1 = prop_avg->data[r][c-1];
                                m2_x1 = std::fabs(fracx2->data[r][c-1]) * M2->data[r][c-1];
                            }
                        }

                        double w_x2 = 0.0;
                        double av_x2 = 0.0;
                        double m2_x2 = 0.0;
                        if(c < map->data.nr_cols() - 1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(source->data[r][c+1]))
                            {
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * source_acc->data[r][c+1];
                                av_x2 = prop_avg->data[r][c+1];
                                m2_x2 = std::fabs(fracx1->data[r][c+1]) * M2->data[r][c+1];
                            }
                        }

                        double w_y1 = 0.0;
                        double av_y1 = 0.0;
                        double m2_y1 = 0.0;
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10)  && !pcr::isMV(source->data[r-1][c]))
                            {
                                w_y1 = std::fabs(fracy2->data[r-1][c]) * source_acc->data[r-1][c];
                                av_y1 = prop_avg->data[r-1][c];
                                m2_y1 = std::fabs(fracy2->data[r-1][c]) * M2->data[r-1][c];
                            }
                        }

                        double w_y2 = 0.0;
                        double av_y2 = 0.0;
                        double m2_y2 = 0.0;
                        if(r < map->data.nr_rows() - 1)
                        {
                            if((std::fabs(fracy1->data[r+1][c])  > 1e-10) && !pcr::isMV(source->data[r+1][c]))
                            {
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * source_acc->data[r+1][c];
                                av_y2 = prop_avg->data[r+1][c];
                                m2_y2 = std::fabs(fracy1->data[r+1][c]) *M2->data[r+1][c];
                            }
                        }

                        double w_x = 0.0;
                        double av_x= 0.0;
                        double m2_x = 0.0;
                        double w_y = 0.0;
                        double av_y = 0.0;
                        double m2_y = 0.0;

                        //Combine_M2(double n_a,double avg_a, double M2_a, double n_b, double avg_b, double M2_b)

                        //M2->data[r][c] = m2_x1 + m2_x2 + m2_y1 + m2_y2 + 1.0;



                        if((w_x1 > 1e-10) && (w_x2 >  1e-10))
                        {
                            w_x = w_x1 + w_x2;
                            av_x = (w_x1 * av_x1 + w_x2 * av_x2)/std::max(1e-12,w_x);
                            m2_x = Combine_M2(w_x1,av_x1, m2_x1, w_x2,av_x2,m2_x2);

                        }else if (w_x2 >  1e-10)
                        {
                            w_x = w_x2;
                            av_x = av_x2;
                            m2_x = m2_x2;
                        }else if (w_x1 >  1e-10)
                        {
                            w_x = w_x1;
                            av_x = av_x1;
                            m2_x = m2_x1;
                        }

                        if((w_y2 >  1e-10) && (w_y1 >  1e-10))
                        {
                            w_y = w_y1 + w_y2;
                            av_y = (w_y1 * av_y1 + w_y2 * av_y2)/std::max(1e-12,w_y);
                            m2_y = Combine_M2(w_y1,av_y1, m2_y1, w_y2,av_y2,m2_y2);

                        }else if (w_y2 >  1e-10)
                        {
                            w_y = w_y2;
                            av_y = av_y2;
                            m2_y = m2_y2;
                        }else if (w_y1 >  1e-10)
                        {
                            w_y = w_y1;
                            av_y = av_y1;
                            m2_y = m2_y1;
                        }

                        if((w_x >  1e-10) && (w_y >  1e-10))
                        {
                            M2->data[r][c] = Combine_M2(w_y,av_y, m2_y, w_x,av_x,m2_x);

                        }else if(w_x > 1e-10 )
                        {
                            M2->data[r][c] = m2_x + 1.0;

                        }else if(w_y >  1e-10)
                        {
                            M2->data[r][c] = m2_y + 1.0;
                        }

                    }


                }
            }

        }
        first = false;

    }

    #pragma omp parallel for collapse(2)
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(source->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = M2->data[r][c]/std::max(1e-12f,source_acc->data[r][c]);
            }
        }
    }

    delete M2;

    return map;

}


inline cTMapf * SpreadFlowMD2Culvert(cTMapf * source ,cTMapf * fracx1, cTMapf * fracx2, cTMapf * fracy1,cTMapf * fracy2, int iter_max ,cTMapf * CULVERTID, std::vector<float>  &CULVERTFRAC, std::vector<int>  CULVERTRS, std::vector<int>  CULVERTCS, std::vector<int>  CULVERTRT, std::vector<int> CULVERTCT)
{
    //printf("start spread flow md2 %f %f %f \n", source->cellSizeX(), source->cellSizeY(), source->cellSize());
    //MaskedRaster<float> raster_data(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    //printf("start spread flow md22 \n");

	cTMapf *mx1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data2();
    cTMapf *mx2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data3(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data4(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");

    //MaskedRaster<float> raster_data5(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *map = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");

    //printf("start spread flow md23\n");


    //initialize the map with friction values

    //#pragma omp parallel for collapse(2)
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(source->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 0.0;
                mx1->data[r][c] = fracx1->data[r][c] * source->data[r][c];
                my1->data[r][c] = fracy1->data[r][c] * source->data[r][c];
                mx2->data[r][c] = fracx2->data[r][c] * source->data[r][c];
                my2->data[r][c] = fracy2->data[r][c] * source->data[r][c];

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
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                           //FF::cout << "culvert cell found " << r << " " << c <<  std::endl;

                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {

                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] * 1.0;

                                add += qtoth * frac_culv;

                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                            }

                        }




                        if((r+1 < map->data.nr_rows()) && myf2 > 1e-10)
                        {
                            {
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;
                                map->data[r][c] += myf2;

                            }

                            my2->data[r][c] = 0.0;
                        }

                        if((c+1 < map->data.nr_cols()) && mxf2 > 1e-10)
                        {
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;
                                map->data[r][c] += mxf2;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                           //FF::cout << "culvert cell found " << r << " " << c <<  std::endl;

                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {

                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] * 1.0;

                                add += qtoth * frac_culv;

                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                            }

                        }

                        if((r-1 > -1) && myf1 < -1e-10)
                        {
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;
                                map->data[r][c] -= myf1;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if((c-1 > -1) && mxf1 < -1e-10)
                        {
                            {
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;
                                map->data[r][c] -= mxf1;

                            }

                            mx1->data[r][c] = 0.0;
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
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                           //FF::cout << "culvert cell found " << r << " " << c <<  std::endl;

                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {

                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] * 1.0;

                                add += qtoth * frac_culv;

                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                            }

                        }

                        if((r+1 < map->data.nr_rows()) && myf2 > 1e-10)
                        {
                            {
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;
                                map->data[r][c] += myf2;

                            }

                            my2->data[r][c] = 0.0;
                        }
                        if((c-1 > -1) && mxf1 < -1e-10)
                        {
                            {
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;
                                map->data[r][c] -= mxf1;

                            }

                            mx1->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                           //FF::cout << "culvert cell found " << r << " " << c <<  std::endl;

                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {

                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] * 1.0;

                                add += qtoth * frac_culv;

                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                            }

                        }

                        if((r-1 > -1) && myf1 < -1e-10)
                        {
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;
                                map->data[r][c] -= myf1;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if((c+1 < map->data.nr_cols()) && mxf2 > 1e-10)
                        {
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;
                                map->data[r][c] += mxf2;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }


                }
            }

        }
        first = false;

    }

    delete mx1;
    delete my1;
    delete mx2;
    delete my2;

    //printf("start spread flow md24\n");


    return map;

}


inline cTMapf * SpreadFlowMD2(cTMapf * source ,cTMapf * fracx1, cTMapf * fracx2, cTMapf * fracy1,cTMapf * fracy2, int iter_max = 0)
{
    //printf("start spread flow md2 %f %f %f \n", source->cellSizeX(), source->cellSizeY(), source->cellSize());
    //MaskedRaster<float> raster_data(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    //printf("start spread flow md22 \n");

	cTMapf *mx1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data2();
    cTMapf *mx2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data3(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data4(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");

    //MaskedRaster<float> raster_data5(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *map = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");

    //printf("start spread flow md23\n");


    //initialize the map with friction values

    //#pragma omp parallel for collapse(2)
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(source->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 0.0;
                mx1->data[r][c] = fracx1->data[r][c] * source->data[r][c];
                my1->data[r][c] = fracy1->data[r][c] * source->data[r][c];
                mx2->data[r][c] = fracx2->data[r][c] * source->data[r][c];
                my2->data[r][c] = fracy2->data[r][c] * source->data[r][c];

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
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];
                        map->data[r][c] += std::max(0.0f,myf2);
                        map->data[r][c] += std::max(0.0f,mxf2);
                        if((r+1 < map->data.nr_rows()) && myf2 > 1e-10)
                        {
                            {
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;
                            }
                        }
                        if((c+1 < map->data.nr_cols()) && mxf2 > 1e-10)
                        {
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;
                            }
                        }

                        my2->data[r][c] = 0.0;
                        mx2->data[r][c] = 0.0;
                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];
                        map->data[r][c] -= std::min(0.0f,myf1);
                        map->data[r][c] -= std::min(0.0f,mxf1);

                        if((r-1 > -1) && myf1 < -1e-10)
                        {
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                        }

                        if((c-1 > -1) && mxf1 < -1e-10)
                        {
                            {
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                        }
                        my1->data[r][c] = 0.0;
                        mx1->data[r][c] = 0.0;
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
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];
                        map->data[r][c] += std::max(0.0f,myf2);
                        map->data[r][c] -= std::min(0.0f,mxf1);
                        if((r+1 < map->data.nr_rows()) && myf2 > 1e-10)
                        {
                            {
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;

                            }

                        }
                        if((c-1 > -1) && mxf1 < -1e-10)
                        {
                            {
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                        }
                        my2->data[r][c] = 0.0;
                        mx1->data[r][c] = 0.0;

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];
                        map->data[r][c] += std::max(0.0f,mxf2);
                        map->data[r][c] -= std::min(0.0f,myf1);
                        if((r-1 > -1) && myf1 < -1e-10)
                        {
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                        }

                        if((c+1 < map->data.nr_cols()) && mxf2 > 1e-10)
                        {
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;

                            }

                        }

                        my1->data[r][c] = 0.0;
                        mx2->data[r][c] = 0.0;
                    }


                }
            }

        }
        first = false;

    }

    delete mx1;
    delete my1;
    delete mx2;
    delete my2;

    //printf("start spread flow md24\n");


    return map;

}


inline cTMapf * SpreadFlowMD2Depression(cTMapf * source ,cTMapf * fracx1, cTMapf * fracx2, cTMapf * fracy1,cTMapf * fracy2, int iter_max, cTMapf * DEPRESSIONID, std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, bool store_depressionvol, float scale_in, cTMapf* scale_adjust_mask, float scale_adjust_scale)
{

    //FF::cout << "spread with depression "<< DEPRESSIONVOLS.size() << " " << DEPRESSIONCVOLS.size() << " " << std::endl;

    std::vector<float> temp;
    std::vector<float> *DEPRESSIONCVOLSH;
    if(!store_depressionvol)
    {
        for(int i = 0; i < DEPRESSIONCVOLS.size(); i++)
        {
            temp.push_back(0.0);
        }
        DEPRESSIONCVOLSH = &temp;
    }else
    {
        DEPRESSIONCVOLSH = &DEPRESSIONCVOLS;
    }
    //FF::cout << "spread with depression "<< DEPRESSIONVOLS.size() << " " << DEPRESSIONCVOLS.size() << " " << DEPRESSIONCVOLSH->size() << std::endl;


    cTMapf *mx1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data2();
    cTMapf *mx2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data3(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data4(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");

    //MaskedRaster<float> raster_data5(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *map = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");


    //initialize the map with friction values

    //#pragma omp parallel for collapse(2)
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(source->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 0.0;
                mx1->data[r][c] = fracx1->data[r][c] * source->data[r][c];
                my1->data[r][c] = fracy1->data[r][c] * source->data[r][c];
                mx2->data[r][c] = fracx2->data[r][c] * source->data[r][c];
                my2->data[r][c] = fracy2->data[r][c] * source->data[r][c];

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
                    float scale = scale_in;
                    if(scale_adjust_mask)
                    {
                        float area = scale_adjust_mask->data[r][c]*scale_adjust_scale*scale_adjust_scale;
                        scale = std::max(scale,sqrt(area));
                    }
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];


                        if(myf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,myf2)*scale < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (frac_depr);

                                        myf2 = myf2 * (1.0-frac_depr);
                                    }
                                }else{
                                    //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                }

                            }

                            map->data[r][c] += myf2;
                            if((r+1 < map->data.nr_rows()))
                            {
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;

                            }

                            my2->data[r][c] = 0.0;
                        }

                        if( mxf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {

                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,mxf2*scale) *  (frac_depr);

                                        mxf2 = mxf2 * (1.0-frac_depr);
                                    }
                                }else{
                                    //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                    }
                            }

                            map->data[r][c] += mxf2;
                            if((c+1 < map->data.nr_cols()))
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float scale = scale_in;
                    if(scale_adjust_mask)
                    {
                        float area = scale_adjust_mask->data[r][c]*scale_adjust_scale*scale_adjust_scale;
                        scale = std::max(scale,sqrt(area));
                    }
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];



                        if( myf1 < -1e-10)
                        {

                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                 if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {

                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,-myf1)*scale < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-myf1*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-myf1*scale) *  (frac_depr);

                                        myf1 = myf1 * (1.0-frac_depr);
                                    }

                                }else{
                                    //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                }
                            }
                            map->data[r][c] -= myf1;
                            if((r-1 > -1))
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if(mxf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                 if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {

                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,-mxf1)*scale < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-mxf1*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-mxf1*scale) *  (frac_depr);

                                        mxf1 = mxf1 * (1.0-frac_depr);
                                    }
                                }else{
                                        //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                    }
                            }
                            map->data[r][c] -= mxf1;
                            if((c-1 > -1))
                            {

                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                            mx1->data[r][c] = 0.0;
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
                    float scale = scale_in;
                    if(scale_adjust_mask)
                    {
                        float area = scale_adjust_mask->data[r][c]*scale_adjust_scale*scale_adjust_scale;
                        scale = std::max(scale,sqrt(area));
                    }
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];


                        if(myf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];

                                 if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,myf2*scale) < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,myf2*scale) *  (frac_depr);

                                        myf2 = myf2 * (1.0-frac_depr);
                                    }
                                }else{
                                    //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                }
                            }
                            map->data[r][c] += myf2;
                            if((r+1 < map->data.nr_rows())){
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;

                            }

                            my2->data[r][c] = 0.0;
                        }
                        if(mxf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];

                                 if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                     float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,-mxf1*scale) < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-mxf1*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-mxf1*scale) *  (frac_depr);

                                        mxf1 = mxf1 * (1.0-frac_depr);
                                    }
                                }else{
                                    //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                }
                            }
                            map->data[r][c] -= mxf1;
                            if((c-1 > -1)){
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                            mx1->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float scale = scale_in;
                    if(scale_adjust_mask)
                    {
                        float area = scale_adjust_mask->data[r][c]*scale_adjust_scale*scale_adjust_scale;
                        scale = std::max(scale,sqrt(area));
                    }
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];



                        if(myf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];

                                 if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,-myf1*scale) < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-myf1*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-myf1*scale) *  (frac_depr);

                                        myf1 = myf1 * (1.0-frac_depr);
                                    }
                                }else{
                                    //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                }
                            }
                            map->data[r][c] -= myf1;
                            if((r-1 > -1))
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if(mxf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];

                                 if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);
                                    if(cvol < vol - 0.0001)
                                    {
                                        if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                        {
                                            frac_depr = 1.0;
                                        }else
                                        {
                                            frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                        }

                                        DEPRESSIONCVOLSH->at(id) += std::max(0.0f,mxf2*scale) *  (frac_depr);

                                        mxf2 = mxf2 * (1.0-frac_depr);
                                    }
                                }else{
                                    //FF::cout << "depression id not found " << id << "  " << DEPRESSIONVOLS.size() << std::endl;
                                }
                            }

                            map->data[r][c] += mxf2;
                            if((c+1 < map->data.nr_cols()))
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }


                }
            }

        }
        first = false;

    }

    delete mx1;
    delete my1;
    delete mx2;
    delete my2;

    return map;

}




inline cTMapf * SpreadFlowMD2DepressionReservoir(cTMapf * source ,cTMapf * fracx1, cTMapf * fracx2, cTMapf * fracy1,cTMapf * fracy2, int iter_max, cTMapf * DEPRESSIONID, std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, bool store_depressionvol, cTMapf * RESERVOIRID, std::vector<float> & RESERVOIRCVOLS, std::vector<float> & RESERVOIRCAPPASS, std::vector<float> &RESERVOIRCAP, float scale)
{

    //FF::cout << "reservoirID map " << RESERVOIRID << " " <<  DEPRESSIONID <<  std::endl;
    //FF::cout << "spread with depression and reservoir"<< DEPRESSIONVOLS.size() << " " << DEPRESSIONCVOLS.size() << " " <<RESERVOIRCVOLS.size() << " " << RESERVOIRCAP.size() <<  std::endl;

    std::vector<float> temp;
    std::vector<float> temp2;
    std::vector<float> *DEPRESSIONCVOLSH;
    std::vector<float> *RESERVOIRCVOLSH;
    if(!store_depressionvol)
    {
        for(int i = 0; i < DEPRESSIONCVOLS.size(); i++)
        {
            temp.push_back(0.0);

        }
        for(int i = 0; i < RESERVOIRCVOLS.size(); i++)
        {
            //FF::cout << "capacity " << i << " " << RESERVOIRCAP.at(i) << " " << RESERVOIRCAPPASS.at(i) << std::endl;
            temp2.push_back(RESERVOIRCVOLS.at(i));
        }
        DEPRESSIONCVOLSH = &temp;
        RESERVOIRCVOLSH = &temp2;
    }else
    {
        DEPRESSIONCVOLSH = &DEPRESSIONCVOLS;
        RESERVOIRCVOLSH = &RESERVOIRCVOLS;
    }
    cTMapf *mx1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data2();
    cTMapf *mx2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data3(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data4(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");

    //MaskedRaster<float> raster_data5(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *map = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");


    //initialize the map with friction values
    //#pragma omp parallel for collapse(2)
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(source->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 0.0;
                mx1->data[r][c] = fracx1->data[r][c] * source->data[r][c];
                my1->data[r][c] = fracy1->data[r][c] * source->data[r][c];
                mx2->data[r][c] = fracx2->data[r][c] * source->data[r][c];
                my2->data[r][c] = fracy2->data[r][c] * source->data[r][c];

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

                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];


                        if( myf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (frac_depr);

                                    myf2 = myf2 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);


                                //FF::cout << "encountered reservoir at " << r << " "<< c<< " " << id << " " << vol << " " << vol2 << " " << cvol << std::endl;

                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (1.0-frac_depr);

                                    myf2 = myf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (frac_depr2);

                                    myf2 = myf2 * (1.0-frac_depr2);

                                }
                            }

                            map->data[r][c] += myf2;
                            if((r+1 < map->data.nr_rows())){
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;

                            }

                            my2->data[r][c] = 0.0;
                        }

                        if(mxf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,mxf2*scale) *  (frac_depr);

                                    mxf2 = mxf2 * (1.0-frac_depr);
                                }
                            }
                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (1.0-frac_depr);

                                    mxf2 = mxf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (frac_depr2);

                                    mxf2 = mxf2 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] += mxf2;
                            if((c+1 < map->data.nr_cols())){
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];



                        if(myf1 < -1e-10)
                        {

                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-myf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-myf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-myf1*scale) *  (frac_depr);

                                    myf1 = myf1 * (1.0-frac_depr);
                                }
                            }
                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (1.0-frac_depr);

                                    myf1 = myf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (frac_depr2);

                                    myf1 = myf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= myf1;
                            if((r-1 > -1))
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if(mxf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-mxf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-mxf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-mxf1*scale) *  (frac_depr);

                                    mxf1 = mxf1 * (1.0-frac_depr);
                                }

                            }
                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (1.0-frac_depr);

                                    mxf1 = mxf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (frac_depr2);

                                    mxf1 = mxf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= mxf1;
                            if((c-1 > -1)){
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                            mx1->data[r][c] = 0.0;
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

                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];


                        if(myf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];

                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2*scale) < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,myf2*scale) *  (frac_depr);

                                    myf2 = myf2 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];

                               float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);

                                //FF::cout << "encountered reservoir at " << r << " "<< c<< " " << id << " " << vol << " " << vol2 << " " << cvol << std::endl;


                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (1.0-frac_depr);

                                    myf2 = myf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (frac_depr2);

                                    myf2 = myf2 * (1.0-frac_depr2);

                                }
                            }

                            map->data[r][c] += myf2;
                            if((r+1 < map->data.nr_rows()))
                            {
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;

                            }

                            my2->data[r][c] = 0.0;
                        }
                        if(mxf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-mxf1*scale) < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-mxf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-mxf1*scale) *  (frac_depr);

                                    mxf1 = mxf1 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (1.0-frac_depr);

                                    mxf1 = mxf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (frac_depr2);

                                    mxf1 = mxf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= mxf1;
                            if((c-1 > -1))
                            {
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                            mx1->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];



                        if(myf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-myf1*scale) < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-myf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-myf1*scale) *  (frac_depr);

                                    myf1 = myf1 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (1.0-frac_depr);

                                    myf1 = myf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (frac_depr2);

                                    myf1 = myf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= myf1;
                            if(r-1 > -1)
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if( mxf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,mxf2*scale) *  (frac_depr);

                                    mxf2 = mxf2 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (1.0-frac_depr);

                                    mxf2 = mxf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (frac_depr2);

                                    mxf2 = mxf2 * (1.0-frac_depr2);

                                }
                            }

                            map->data[r][c] += mxf2;
                            if((c+1 < map->data.nr_cols()))
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }


                }
            }

        }
        first = false;

    }

    delete mx1;
    delete my1;
    delete mx2;
    delete my2;

    ////printf("spread with depression and reservoir done");
    return map;

}



inline cTMapf * SpreadFlowMD2DepressionReservoirCulvert(cTMapf * source ,cTMapf * fracx1, cTMapf * fracx2, cTMapf * fracy1,cTMapf * fracy2, int iter_max, cTMapf * DEPRESSIONID, std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, bool store_depressionvol, cTMapf * RESERVOIRID, std::vector<float> & RESERVOIRCVOLS, std::vector<float> & RESERVOIRCAPPASS, std::vector<float> &RESERVOIRCAP, cTMapf * CULVERTID, std::vector<float>  &CULVERTQ, std::vector<float>  &CULVERTFRAC,std::vector<float>  &CULVERTQTOT, std::vector<int>  CULVERTRS, std::vector<int>  CULVERTCS, std::vector<int>  CULVERTRT, std::vector<int> CULVERTCT,float scale, cTMapf * drainage, cTMapf * drainage_tot)
{

    //FF::cout << "CukvertID map " << RESERVOIRID << " " <<  DEPRESSIONID <<  std::endl;
    //FF::cout << "spread with depression and reservoir"<< DEPRESSIONVOLS.size() << " " << DEPRESSIONCVOLS.size() << " " <<RESERVOIRCVOLS.size() << " " << RESERVOIRCAP.size() <<  std::endl;

    std::vector<float> temp;
    std::vector<float> temp2;
    std::vector<float> *DEPRESSIONCVOLSH;
    std::vector<float> *RESERVOIRCVOLSH;

    //std::vector<float> CULVERTQTOT;

    //for(int i = 0; i < CULVERTQ.size(); i++)
    ///{
    //    CULVERTQTOT.push(0.0);
    //}

    if(!store_depressionvol)
    {
        for(int i = 0; i < DEPRESSIONCVOLS.size(); i++)
        {
            temp.push_back(0.0);

        }
        for(int i = 0; i < RESERVOIRCVOLS.size(); i++)
        {
            //FF::cout << "capacity " << i << " " << RESERVOIRCAP.at(i) << " " << RESERVOIRCAPPASS.at(i) << std::endl;
            temp2.push_back(RESERVOIRCVOLS.at(i));
        }
        DEPRESSIONCVOLSH = &temp;
        RESERVOIRCVOLSH = &temp2;
    }else
    {
        DEPRESSIONCVOLSH = &DEPRESSIONCVOLS;
        RESERVOIRCVOLSH = &RESERVOIRCVOLS;
    }
    cTMapf *mx1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data2();
    cTMapf *mx2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data3(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my1 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");
    //MaskedRaster<float> raster_data4(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *my2 = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");

    //MaskedRaster<float> raster_data5(source->data.nr_rows(), source->data.nr_cols(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());
    cTMapf *map = new cTMapf(nullptr, source->data.nr_cols(),source->data.nr_rows(), source->data.north(), source->data.west(), source->data.cell_size(),source->data.cell_sizeY());//std::move(raster_data),source->projection(),"");


    //initialize the map with friction values
    //#pragma omp parallel for collapse(2)
    for(int r = 0; r < map->data.nr_rows();r++)
    {
        for(int c = 0; c < map->data.nr_cols();c++)
        {
            if(pcr::isMV(source->data[r][c]))
            {
                pcr::setMV(map->data[r][c]);
            }else {
                map->data[r][c] = 0.0;
                mx1->data[r][c] = fracx1->data[r][c] * source->data[r][c];
                my1->data[r][c] = fracy1->data[r][c] * source->data[r][c];
                mx2->data[r][c] = fracx2->data[r][c] * source->data[r][c];
                my2->data[r][c] = fracy2->data[r][c] * source->data[r][c];

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

                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];

                        if(drainage != nullptr)
                        {
                            if(drainage->data[r][c] > 1e-11)
                            {
                                float room = std::max(0.0f, drainage->data[r][c] - drainage_tot->data[r][c]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_sd = std::max(0.0f,std::min(1.0f,room/std::max(0.0001f,qtoth)));

                                //FF::cout << " sd: "  << drainage->data[r][c] << " "  << qtoth << " "  << room << "  " << frac_sd <<std::endl;

                                mx1->data[r][c] = (1.0-frac_sd) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_sd) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_sd) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_sd) * my2->data[r][c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                                drainage_tot->data[r][c] += qtoth * frac_sd;

                            }

                        }

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                           //FF::cout << "culvert cell found " << r << " " << c <<  std::endl;

                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            //FF::cout << " source?  " << is_source << "  " << target_r << " " << target_c << std::endl;

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {


                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] * std::min(1.0f,room/std::max(0.0001f,qtoth));

                                //FF::cout << "culvert cell found and valid source " << CULVERTFRAC[idc] << " " <<   room << " " << frac_culv << " " << CULVERTQ[idc] << " " << CULVERTQTOT[idc] << " " << map->data[r][c] << std::endl;


                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                add += qtoth * frac_culv;
                                CULVERTQTOT[idc] += qtoth * frac_culv;


                                //connect this to any depressions or reservoirs, that can also outflow through the culvert
                                if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                                {
                                    int id = DEPRESSIONID->data[r][c];
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);

                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol/ scale)));

                                    add += (vol/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol/ scale) * frac_culv;
                                    DEPRESSIONCVOLSH->at(id) = DEPRESSIONCVOLSH->at(id) * (1.0-frac_culv);
                                }

                                if(!pcr::isMV(RESERVOIRID->data[r][c]))
                                {
                                    int id = RESERVOIRID->data[r][c];
                                    float vol = RESERVOIRCAPPASS.at(id);
                                    float vol2 = RESERVOIRCAP.at(id);
                                    float cvol = RESERVOIRCVOLSH->at(id);

                                    //room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    //frac_culv = std::min(1.0f,(vol/ scale)/std::max(0.0001f,room));

                                    //add += (vol/ scale) * frac_culv;
                                    //CULVERTQTOT[idc] += (vol/ scale) * frac_culv;


                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol2/ scale)));

                                    add +=  (vol2/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol2/ scale) * frac_culv;
                                    RESERVOIRCAP.at(id) = RESERVOIRCAP.at(id) * (1.0-frac_culv);
                                }

                                //redistribute the added discharge into the culvert to the target cell
                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];

                            }

                        }

                        if(myf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (frac_depr);

                                    myf2 = myf2 * (1.0-frac_depr);
                                }

                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);

                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (1.0-frac_depr);

                                    myf2 = myf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (frac_depr2);

                                    myf2 = myf2 * (1.0-frac_depr2);

                                }
                            }

                            map->data[r][c] += myf2;
                            if((r+1 < map->data.nr_rows())){
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;

                            }

                            my2->data[r][c] = 0.0;
                        }

                        if(mxf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,mxf2*scale) *  (frac_depr);

                                    mxf2 = mxf2 * (1.0-frac_depr);
                                }
                            }
                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (1.0-frac_depr);

                                    mxf2 = mxf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (frac_depr2);

                                    mxf2 = mxf2 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] += mxf2;
                            if((c+1 < map->data.nr_cols()))
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];

                        if(drainage != nullptr)
                        {
                            if(drainage->data[r][c] > 1e-11)
                            {
                                float room = std::max(0.0f, drainage->data[r][c] - drainage_tot->data[r][c]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_sd = std::max(0.0f,std::min(1.0f,room/std::max(0.0001f,qtoth)));

                                mx1->data[r][c] = (1.0-frac_sd) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_sd) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_sd) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_sd) * my2->data[r][c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                                drainage_tot->data[r][c] += qtoth * frac_sd;

                            }

                        }

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {

                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] *std::min(1.0f,room/std::max(0.0001f,qtoth));

                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                add += qtoth * frac_culv;
                                CULVERTQTOT[idc] += qtoth * frac_culv;


                                //connect this to any depressions or reservoirs, that can also outflow through the culvert
                                if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                                {
                                    int id = DEPRESSIONID->data[r][c];
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);

                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol/ scale)));

                                    add += (vol/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol/ scale) * frac_culv;
                                    DEPRESSIONCVOLSH->at(id) = DEPRESSIONCVOLSH->at(id) * (1.0-frac_culv);
                                }

                                if(!pcr::isMV(RESERVOIRID->data[r][c]))
                                {
                                    int id = RESERVOIRID->data[r][c];
                                    float vol = RESERVOIRCAPPASS.at(id);
                                    float vol2 = RESERVOIRCAP.at(id);
                                    float cvol = RESERVOIRCVOLSH->at(id);

                                    //room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    //frac_culv = std::min(1.0f,(vol/ scale)/std::max(0.0001f,room));

                                    //add += (vol/ scale) * frac_culv;
                                    //CULVERTQTOT[idc] += (vol/ scale) * frac_culv;


                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol2/ scale)));

                                    add +=  (vol2/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol2/ scale) * frac_culv;
                                    RESERVOIRCAP.at(id) = RESERVOIRCAP.at(id) * (1.0-frac_culv);
                                }

                                //redistribute the added discharge into the culvert to the target cell
                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];

                            }

                        }

                        if(myf1 < -1e-10)
                        {

                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-myf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-myf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-myf1*scale) *  (frac_depr);

                                    myf1 = myf1 * (1.0-frac_depr);
                                }
                            }
                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (1.0-frac_depr);

                                    myf1 = myf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (frac_depr2);

                                    myf1 = myf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= myf1;
                            if((r-1 > -1))
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if(mxf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-mxf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-mxf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-mxf1*scale) *  (frac_depr);

                                    mxf1 = mxf1 * (1.0-frac_depr);
                                }

                            }
                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (1.0-frac_depr);

                                    mxf1 = mxf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (frac_depr2);

                                    mxf1 = mxf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= mxf1;
                            if((c-1 > -1))
                            {
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                            mx1->data[r][c] = 0.0;
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

                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];
                        if(drainage != nullptr)
                        {
                            if(drainage->data[r][c] > 1e-11)
                            {
                                float room = std::max(0.0f, drainage->data[r][c] - drainage_tot->data[r][c]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_sd = std::max(0.0f,std::min(1.0f,room/std::max(0.0001f,qtoth)));

                                mx1->data[r][c] = (1.0-frac_sd) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_sd) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_sd) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_sd) * my2->data[r][c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                                drainage_tot->data[r][c] += qtoth * frac_sd;

                            }

                        }

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {

                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] *std::min(1.0f,room/std::max(0.0001f,qtoth));

                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                add += qtoth * frac_culv;
                                CULVERTQTOT[idc] += qtoth * frac_culv;


                                //connect this to any depressions or reservoirs, that can also outflow through the culvert
                                if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                                {
                                    int id = DEPRESSIONID->data[r][c];
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);

                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol/ scale)));

                                    add += (vol/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol/ scale) * frac_culv;
                                    DEPRESSIONCVOLSH->at(id) = DEPRESSIONCVOLSH->at(id) * (1.0-frac_culv);
                                }

                                if(!pcr::isMV(RESERVOIRID->data[r][c]))
                                {
                                    int id = RESERVOIRID->data[r][c];
                                    float vol = RESERVOIRCAPPASS.at(id);
                                    float vol2 = RESERVOIRCAP.at(id);
                                    float cvol = RESERVOIRCVOLSH->at(id);

                                    //room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    //frac_culv = std::min(1.0f,(vol/ scale)/std::max(0.0001f,room));

                                    //add += (vol/ scale) * frac_culv;
                                    //CULVERTQTOT[idc] += (vol/ scale) * frac_culv;


                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol2/ scale)));

                                    add +=  (vol2/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol2/ scale) * frac_culv;
                                    RESERVOIRCAP.at(id) = RESERVOIRCAP.at(id) * (1.0-frac_culv);
                                }

                                //redistribute the added discharge into the culvert to the target cell
                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];

                            }

                        }

                        if( myf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];

                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2*scale) < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,myf2*scale) *  (frac_depr);

                                    myf2 = myf2 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];

                               float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);

                                //FF::cout << "encountered reservoir at " << r << " "<< c<< " " << id << " " << vol << " " << vol2 << " " << cvol << std::endl;


                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (1.0-frac_depr);

                                    myf2 = myf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf2*scale) * (frac_depr2);

                                    myf2 = myf2 * (1.0-frac_depr2);

                                }
                            }

                            map->data[r][c] += myf2;
                            if((r+1 < map->data.nr_rows()))
                            {
                                mx1->data[r+1][c] += myf2 * fracx1->data[r+1][c];
                                my1->data[r+1][c] += myf2 * fracy1->data[r+1][c];
                                mx2->data[r+1][c] += myf2 * fracx2->data[r+1][c];
                                my2->data[r+1][c] += myf2 * fracy2->data[r+1][c];

                                change  = true;

                            }

                            my2->data[r][c] = 0.0;
                        }
                        if(mxf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-mxf1*scale) < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-mxf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-mxf1*scale) *  (frac_depr);

                                    mxf1 = mxf1 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (1.0-frac_depr);

                                    mxf1 = mxf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf1*scale) * (frac_depr2);

                                    mxf1 = mxf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= mxf1;
                            if((c-1 > -1))
                            {
                                mx1->data[r][c-1] -= mxf1 * fracx1->data[r][c-1];
                                my1->data[r][c-1] -= mxf1 * fracy1->data[r][c-1];
                                mx2->data[r][c-1] -= mxf1 * fracx2->data[r][c-1];
                                my2->data[r][c-1] -= mxf1 * fracy2->data[r][c-1];
                                change  = true;

                            }

                            mx1->data[r][c] = 0.0;
                        }

                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = source->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        float mxf1 = mx1->data[r][c];
                        float myf1 = my1->data[r][c];
                        float mxf2 = mx2->data[r][c];
                        float myf2 = my2->data[r][c];
                        if(drainage != nullptr)
                        {
                            if(drainage->data[r][c] > 1e-11)
                            {
                                float room = std::max(0.0f, drainage->data[r][c] - drainage_tot->data[r][c]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_sd = std::max(0.0f,std::min(1.0f,room/std::max(0.0001f,qtoth)));

                                mx1->data[r][c] = (1.0-frac_sd) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_sd) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_sd) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_sd) * my2->data[r][c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];
                                drainage_tot->data[r][c] += qtoth * frac_sd;

                            }

                        }

                        if(CULVERTID->data[r][c] >= 0.0)
                        {
                            bool is_source = false;
                            int target_r = 0;
                            int target_c = 0;
                            for(int i = 0; i < CULVERTRS.size(); i++)
                            {
                                if(r ==  CULVERTRS[i] && c == CULVERTCS[i])
                                {
                                    is_source = true;
                                    target_r = CULVERTRT[i];
                                    target_c = CULVERTCT[i];
                                    break;
                                }
                            }

                            if(is_source && target_r > -1 && target_c > -1 && target_r < map->data.nr_rows() && target_c < map->data.nr_cols())
                            {

                                float add = 0.0;
                                int idc = CULVERTID->data[r][c];
                                float room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                float qtoth = std::fabs(mxf1)+std::fabs(mxf2)+std::fabs(myf1)+std::fabs(myf2);
                                float frac_culv = CULVERTFRAC[idc] *std::min(1.0f,room/std::max(0.0001f,qtoth));

                                mx1->data[r][c] = (1.0-frac_culv) * mx1->data[r][c];
                                my1->data[r][c] = (1.0-frac_culv) * my1->data[r][c];
                                mx2->data[r][c] = (1.0-frac_culv) * mx2->data[r][c];
                                my2->data[r][c] = (1.0-frac_culv) * my2->data[r][c];

                                add += qtoth * frac_culv;
                                CULVERTQTOT[idc] += qtoth * frac_culv;


                                //connect this to any depressions or reservoirs, that can also outflow through the culvert
                                if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                                {
                                    int id = DEPRESSIONID->data[r][c];
                                    float vol = DEPRESSIONVOLS.at(id);
                                    float cvol = DEPRESSIONCVOLSH->at(id);

                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol/ scale)));

                                    add += (vol/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol/ scale) * frac_culv;
                                    DEPRESSIONCVOLSH->at(id) = DEPRESSIONCVOLSH->at(id) * (1.0-frac_culv);
                                }

                                if(!pcr::isMV(RESERVOIRID->data[r][c]))
                                {
                                    int id = RESERVOIRID->data[r][c];
                                    float vol = RESERVOIRCAPPASS.at(id);
                                    float vol2 = RESERVOIRCAP.at(id);
                                    float cvol = RESERVOIRCVOLSH->at(id);

                                    //room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    //frac_culv = std::min(1.0f,(vol/ scale)/std::max(0.0001f,room));

                                    //add += (vol/ scale) * frac_culv;
                                    //CULVERTQTOT[idc] += (vol/ scale) * frac_culv;


                                    room = std::max(0.0f,CULVERTQ[idc] - CULVERTQTOT[idc]);
                                    frac_culv = std::min(1.0f,room/std::max(0.0001f,(vol2/ scale)));

                                    add +=  (vol2/ scale) * frac_culv;
                                    CULVERTQTOT[idc] += (vol2/ scale) * frac_culv;
                                    RESERVOIRCAP.at(id) = RESERVOIRCAP.at(id) * (1.0-frac_culv);
                                }

                                //redistribute the added discharge into the culvert to the target cell
                                mx1->data[target_r][target_c] += add * fracx1->data[target_r][target_c];
                                my1->data[target_r][target_c] += add * fracy1->data[target_r][target_c];
                                mx2->data[target_r][target_c] += add * fracx2->data[target_r][target_c];
                                my2->data[target_r][target_c] += add * fracy2->data[target_r][target_c];

                                mxf1 = mx1->data[r][c];
                                myf1 = my1->data[r][c];
                                mxf2 = mx2->data[r][c];
                                myf2 = my2->data[r][c];

                            }

                        }


                        if(myf1 < -1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,-myf1*scale) < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,-myf1*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,-myf1*scale) *  (frac_depr);

                                    myf1 = myf1 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (1.0-frac_depr);

                                    myf1 = myf1 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,myf1)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,myf1*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,myf1*scale) * (frac_depr2);

                                    myf1 = myf1 * (1.0-frac_depr2);

                                }
                            }
                            map->data[r][c] -= myf1;
                            if((r-1 > -1))
                            {
                                mx1->data[r-1][c] -= myf1 * fracx1->data[r-1][c];
                                my1->data[r-1][c] -= myf1 * fracy1->data[r-1][c];
                                mx2->data[r-1][c] -= myf1 * fracx2->data[r-1][c];
                                my2->data[r-1][c] -= myf1 * fracy2->data[r-1][c];
                                change  = true;

                            }

                            my1->data[r][c] = 0.0;
                        }

                        if( mxf2 > 1e-10)
                        {
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                float vol = DEPRESSIONVOLS.at(id);
                                float cvol = DEPRESSIONCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    DEPRESSIONCVOLSH->at(id) += std::max(0.0f,mxf2*scale) *  (frac_depr);

                                    mxf2 = mxf2 * (1.0-frac_depr);
                                }
                            }

                            if(!pcr::isMV(RESERVOIRID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                float frac_depr2 = 0.0;
                                int id = RESERVOIRID->data[r][c];
                                float vol = RESERVOIRCAPPASS.at(id);
                                float vol2 = RESERVOIRCAP.at(id);
                                float cvol = RESERVOIRCVOLSH->at(id);
                                if(cvol < vol - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol - cvol)
                                    {
                                        frac_depr = 1.0;
                                    }else
                                    {
                                        frac_depr = 1.0 - std::max(0.0001f,vol-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (1.0-frac_depr);

                                    mxf2 = mxf2 * (frac_depr);
                                }
                                if(cvol < vol2 - 0.0001)
                                {
                                    if(std::max(0.0f,mxf2)*scale < vol2 - cvol)
                                    {
                                        frac_depr2 = 1.0;
                                    }else
                                    {
                                        frac_depr2 = std::max(0.0001f,vol2-cvol)/(std::max(0.001f,mxf2*scale));
                                    }

                                    RESERVOIRCVOLSH->at(id) += std::max(0.0f,mxf2*scale) * (frac_depr2);

                                    mxf2 = mxf2 * (1.0-frac_depr2);

                                }
                            }

                            map->data[r][c] += mxf2;
                            if((c+1 < map->data.nr_cols()))
                            {
                                mx1->data[r][c+1] += mxf2 * fracx1->data[r][c+1];
                                my1->data[r][c+1] += mxf2 * fracy1->data[r][c+1];
                                mx2->data[r][c+1] += mxf2 * fracx2->data[r][c+1];
                                my2->data[r][c+1] += mxf2 * fracy2->data[r][c+1];
                                change  = true;

                            }

                            mx2->data[r][c] = 0.0;
                        }

                    }


                }
            }

        }
        first = false;

    }

    delete mx1;
    delete my1;
    delete mx2;
    delete my2;

    ////printf("spread with depression and reservoir done");
    return map;

}



}

#endif
