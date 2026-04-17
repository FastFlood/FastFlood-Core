#ifndef FULLFLOOD_CACHE_H
#define FULLFLOOD_CACHE_H


inline int GetMaxSubGridLevelGT(std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, int gt)
{
    if(r < 0 || c < 0 || r > flowstates[l][0]->Rows()-1 || c > flowstates[l][0]->Cols()-1)
    {
      return l;
    }
    if(l >= gt)
    {
        return l;
    }
    if(flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
    {
        int maxl = l;
        for(int r2 = 0; r2 < 3; r2++)
        {
            for (int c2 = 0; c2 < 3; c2++)
            {
                maxl = std::max(maxl,GetMaxSubGridLevelGT(flowstates,r * 3 + r2,c*3 + c2,l+1, gt));
            }
        }
        return maxl;
    }else{
      return l;
    }
}



inline std::tuple<float,float,float,LSMVector2,LSMVector2> GetFullNewFluxAlongEdge(double tglob,int r, int c, int l, int r2, int c2, int l2, float hin, float dem,float ux, float uy,
                                          double courant, double _dt,float dx,std::vector<std::vector<cTMapf *>> &flowstates, bool do_ocean, float ocean_level, float rain_fac)
{

    if(!(r>-1 && c > -1 && r < flowstates[l][0]->Rows() && c < flowstates[l][0]->Cols()))
    {
      return std::make_tuple(0.0,0.0,0.0,LSMVector2(0.0,0.0),LSMVector2(1.0,1.0));
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

    float z = dem;//flowstates[l][0]->data[r][c];
    float h = hin;//flowstates[l][4]->data[r][c];
    float vx = ux;//flowstates[l][5]->data[r][c];
    float vy = uy;//flowstates[l][6]->data[r][c];
    float n = flowstates[l][3]->data[r][c];

    auto mmlambda = [tglob,l, r,c,_dt,dx,  &flowstates,&fx_center, do_ocean, ocean_level, &fy_center, h,vx,vy,z,n,&h_add,&hll_x1,&hll_y1, &hll_x2, &hll_y2, &sx_av, &sy_av](int r2,int c2,int l2,float dy, int ca, int ca2)
    {
        //call the actual function
        MomentumFluxGet(tglob,l, r, c, _dt, dx, flowstates,
                        do_ocean, ocean_level, h, vx, vy, z, n,
                        fx_center /*out*/, fy_center /*out*/,
                        h_add/*out*/, hll_x1 /*out*/, hll_x2 /*out*/, hll_y1 /*out*/, hll_y2 /*out*/,
                        sx_av /*out*/, sy_av /*out*/,
                        r2, c2, l2, dy, ca, ca2);
    };

    //pass function pointer instead of lamba for wider access
    ForEach3CascNeighbor(_dt,dx,flowstates, r2, c2, l2, r,c,l,c-c2,r-r2,true, mmlambda,c2 != c? 1: 2,c2 -c > 0?1 : -1);

    //get the darcy-weisbach scaling
    sx_av = sx_av/dx;
    sy_av = sy_av/dx;
    h_add = h_add/dx;

    //now finally get dh/dt and dux/dy and duy/dt
    float dhdt = GetDHDT((rain_fac * flowstates[l][10]->data[r][c]/(3600000.0)), h_add/(_dt));
    //and alter the values
    float hn = h + dhdt * _dt;

    LSMVector2 dudtadd = GetDUDTAddTerms(hn,sx_av, sy_av,_dt);
    LSMVector2 dudyn = GetDUDTNScale(n,hn,dx,_dt,vx,vy, sx_av, sy_av);
    //LSMVector2 dudt = GetDUDT(n,h,hn, dx,_dt,vx,vy,hll_x1,hll_x2,hll_y1,hll_y2,sx_av,sy_av,fx_center,fy_center);

    if(c2 == c)
    {
        if(r2-r > 0)
        {
            return std::make_tuple(h_add,hll_x2.y,hll_x2.z,dudtadd,dudyn);
        }else
        {
            return std::make_tuple(h_add,hll_x1.y,hll_x1.z,dudtadd,dudyn);
        }
    }else
    {
        if(c2-c > 0)
        {
            return std::make_tuple(h_add,hll_y2.y,hll_y2.z,dudtadd,dudyn);
        }else
        {
            return std::make_tuple(h_add,hll_y1.y,hll_y1.z,dudtadd,dudyn);
        }
    }

}

inline std::tuple<float,float,float,LSMVector2,LSMVector2> GetFullCurrentFluxAlongEdge(double tglob, int r, int c, int l, int r2, int c2, int l2,
                                              double courant, double _dt,float dx,std::vector<std::vector<cTMapf *>> &flowstates, bool do_ocean, float ocean_level, float rain_fac)
{
    if(!(r>-1 && c > -1 && r < flowstates[l][0]->Rows() && c < flowstates[l][0]->Cols()))
    {
      return std::make_tuple(0.0,0.0,0.0,LSMVector2(0.0,0.0),LSMVector2(1.0,1.0));
    }

    float z = flowstates[l][0]->data[r][c];
    float h = flowstates[l][4]->data[r][c];
    float vx = flowstates[l][5]->data[r][c];
    float vy = flowstates[l][6]->data[r][c];
    float n = flowstates[l][3]->data[r][c];

    return GetFullNewFluxAlongEdge(tglob, r, c, l, r2, c2, l2, h, z, vx, vy,
            courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
}

inline std::tuple<float, float, float> GetCellAveraged(int r, int c, int l, std::vector<std::vector<cTMapf *>> &flowstates,float dx)
{

    float h_sum = 0.0;
    float ux_sum = 0.0;
    float uy_sum = 0.0;

    float w_sum = 0.0;

    ForEachSubCell(r,c,l,flowstates,dx,[&h_sum, &ux_sum, &uy_sum, &w_sum](int r2, int c2, int l2, std::vector<std::vector<cTMapf *>> &flowstates,float dy){

        //FF::cout << " subcell for cache mv " << r2 <<  " " << c2 << " " << l2 << std::endl;

        if(!pcr::isMV(flowstates[l2][0]->data[r2][c2]))
        {
            //FF::cout << " subcell for cache " << flowstates[l2][4]->data[r2][c2] <<  "  "  << r2 <<  " " << c2 << " " << l2 << std::endl;
            w_sum += dy*dy;
            h_sum += dy*dy * flowstates[l2][4]->data[r2][c2];
            ux_sum += dy*dy * flowstates[l2][5]->data[r2][c2];
            uy_sum += dy*dy * flowstates[l2][6]->data[r2][c2];
        }

    });

    h_sum = h_sum/std::max(1e-12f,w_sum);
    ux_sum = ux_sum/std::max(1e-12f,w_sum);
    uy_sum = uy_sum/std::max(1e-12f,w_sum);

    return std::make_tuple(h_sum,ux_sum,uy_sum);
}

inline void ApplyRatioToGrid(int r, int c, int l, float dx, std::vector<std::vector<cTMapf *>> &flowstates, float r_h, float r_ux, float r_uy)
{
    ForEachSubCell(r,c,l,flowstates,dx,[r_h,r_ux,r_uy](int r2, int c2, int l2, std::vector<std::vector<cTMapf *>> &flowstates,float dy){
        flowstates[l2][4]->data[r2][c2] *= r_h;
        flowstates[l2][5]->data[r2][c2] *= r_ux;
        flowstates[l2][6]->data[r2][c2] *= r_uy;
        flowstates[l2][7]->data[r2][c2] = flowstates[l2][4]->data[r2][c2];
        flowstates[l2][8]->data[r2][c2] = flowstates[l2][5]->data[r2][c2];
        flowstates[l2][9]->data[r2][c2] = flowstates[l2][6]->data[r2][c2];
    });

}

inline std::tuple<float,float,float,float> GetCellEdgeWHAveraged(float _dt, int r, int c, int l, int r2, int c2, int l2,std::vector<std::vector<cTMapf *>> &flowstates,float dx)
{
    float hwsum = 0.0;
    float h = 0;
    float ux = 0.0;
    float uy = 0.0;
    float dem = 0.0;
    if(!(r>-1 && c > -1 && r < flowstates[l][0]->Rows() && c < flowstates[l][0]->Cols()))
    {
      return std::make_tuple(0.0,0.0,0.0,0.0);
    }

    auto mmlambda = [l, r,c, dx, &h, &ux, &uy, &dem, &hwsum, &flowstates](int r2,int c2,int l2,float dy, int ca, int ca2)
    {
        float dxh = 0.5*(std::fabs(dy) + std::fabs(dx));
        float dx_w = std::min(std::fabs(dy),std::fabs(dx));

        int dim0 = flowstates[l2][0]->Cols();
        int dim1 = flowstates[l2][0]->Rows();
        int x = c2;
        int y = r2;
        const int gy_x1 = std::min(dim1-(int)(1),std::max((int)(0),(int)(y)));
        const int gx_x1 = std::min(dim0-(int)(1),std::max((int)(0),(int)(x)));

        if(h > 0.1)
        {
            //FF::cout << "mlamb " << r << " " << c << " " << gy_x1 << " " << gx_x1 << std::endl;

        }
        float z_x1 = flowstates[l2][0]->data[gy_x1][gx_x1];
        float h_x1 = flowstates[l2][4]->data[gy_x1][gx_x1];
        float vx_x1 = flowstates[l2][5]->data[gy_x1][gx_x1];
        float vy_x1 = flowstates[l2][6]->data[gy_x1][gx_x1];

        dem = dem + z_x1 * h_x1 * dx_w;
        h = h + h_x1;
        hwsum = hwsum + h_x1 * dx_w;
        ux += vx_x1 * h_x1 * dx_w;
        uy += vy_x1 * h_x1 * dx_w;
    };

    //pass function pointer instead of lamba for wider access

    //apply the helper function, but do it reversed, we want to search all the edge cells along [r,c,l] -> [r2,c2,l2]
    ForEach3CascNeighbor(_dt,dx,flowstates, r2, c2, l2, r,c,l,c-c2,r-r2,true, mmlambda,1,1);

    dem = dem/std::max(1e-12f,hwsum);
    ux = ux/std::max(1e-12f,hwsum);
    uy = uy/std::max(1e-12f,hwsum);
    h = h/dx;
    return std::make_tuple(h,ux,uy,dem);
}
inline bool rcf(float a, float b, float ratio, float minv)
{
    if(std::fabs(a-b) < minv)
    {

        return true;
    }else
    {
        if( (a < 0.0?1:0) != (b < 0.0?1:0))
        {
            return false;
        }
        float r = std::fabs(a)/std::fabs(b);
        r = std::max(r,1.0f/std::max(1e-12f,r));
        if(r > (1.0 + ratio))
        {
            return false;
        }
        return true;
    }
}
inline void CheckCacheifyDo(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates,cTMapf * HCache, cTMapf * UXCache, cTMapf * UYCache,float courant, int r, int c, int l,float scale_req, float cache_threshold_rel, bool do_ocean, float ocean_level, float rain_fac)
{
    //if we dont have cache here currently
    if(flowstates[0][17]->data[r][c] < 0.5f)
    {

        //check if we have a sub-level of at least R3 or R1
        int ml = GetMaxSubGridLevelGT(flowstates,r,c,0,2);
        bool hasmaxlevel = ml >= 2;

        //check if all the pixels have been at constant value
        if(hasmaxlevel)
        {
            bool allsame = true;
            int n_pix = 0;

            //what is the pixel count within the first level (R27 probably?)
            int flsize = std::pow(3,flowstates.size()-1);

            //for each pixel within this superpixel
            for(int r2 = r * flsize; r2 < r*flsize + flsize; r2++)
            {
                for(int c2 = c * flsize; c2 < c*flsize + flsize; c2++)
                {
                    if(r2 < HCache->nrRows() && c2 < HCache->nrCols())
                    {
                        n_pix ++ ;
                        int r3 = r2*3;
                        int c3 = c2*3;

                        float h = 0.0;
                        float ux = 0.0;
                        float uy = 0.0;

                        for(int i = flowstates.size() -1; i > -1 ; i--)
                        {
                            r3 = r3/3;
                            c3 = c3/3;

                            //get height and velocity for this layer
                            if(r3 < flowstates[i][0]->nrRows() && c3 < flowstates[i][0]->nrCols() )
                            {

                                h += flowstates[i][4]->data[r3][c3];
                                ux += flowstates[i][5]->data[r3][c3];
                                uy += flowstates[i][6]->data[r3][c3];
                            }
                        }

                        //FF::cout << "compare pix: " <<h << " " << HCache->data[r2][c2] << " "<< ux << " " << UXCache->data[r2][c2] <<std::endl;

                        if(!rcf(h,HCache->data[r2][c2],cache_threshold_rel,0.0001) || !rcf(ux,UXCache->data[r2][c2],cache_threshold_rel,0.0005) || !rcf(uy,UYCache->data[r2][c2],cache_threshold_rel,0.0005))
                        {
                            //FF::cout << "not, stop"<<std::endl;
                            allsame = false;
                            break;
                        }
                    }
                }
                if(!allsame){
                    break;
                }
            }

            //if the water heights are virtually the same, consider it the same
            if(allsame)
            {
                //Activate the Cache by solving the equations around the edge of the superpixel and then compensating

                //run the saint venant solution to calculate
                //weight-averaged DEM
                //weight-averaged water height

                //weight averaged slope term
                //weight average friction term
                /*std::tuple<float,float,float,float> props_x1 = GetCellEdgeWHAveraged(_dt, r,c,0,r,c-1,0,flowstates,dx);
                std::tuple<float,float,float,float> props_x2 = GetCellEdgeWHAveraged(_dt,r,c,0,r,c+1,0,flowstates,dx);
                std::tuple<float,float,float,float> props_y1 = GetCellEdgeWHAveraged(_dt,r,c,0,r-1,c,0,flowstates,dx);
                std::tuple<float,float,float,float> props_y2 = GetCellEdgeWHAveraged(_dt,r,c,0,r+1,c,0,flowstates,dx);

                //for each side
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_old_x1 =  GetFullCurrentFluxAlongEdge(r,c,0,r,c-1,0,courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_new_x1 =  GetFullNewFluxAlongEdge(r,c,0,r,c-1,0,std::get<0>(props_x1),std::get<1>(props_x1),std::get<2>(props_x1),std::get<3>(props_x1),courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_old_x2 =  GetFullCurrentFluxAlongEdge(r,c,0,r,c+1,0,courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_new_x2 =  GetFullNewFluxAlongEdge(r,c,0,r,c+1,0,std::get<0>(props_x2),std::get<1>(props_x2),std::get<2>(props_x2),std::get<3>(props_x2),courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_old_y1 =  GetFullCurrentFluxAlongEdge(r,c,0,r-1,c,0,courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_new_y1 =  GetFullNewFluxAlongEdge(r,c,0,r-1,c,0,std::get<0>(props_y1),std::get<1>(props_y1),std::get<2>(props_y1),std::get<3>(props_y1),courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_old_y2 =  GetFullCurrentFluxAlongEdge(r,c,0,r+1,c,0,courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
                std::tuple<float,float,float,LSMVector2,LSMVector2> flux_new_y2 =  GetFullNewFluxAlongEdge(r,c,0,r+1,c,0,std::get<0>(props_y2),std::get<1>(props_y2),std::get<2>(props_y2),std::get<3>(props_y2),courant, _dt, dx,flowstates, do_ocean, ocean_level, rain_fac);
*/
                std::tuple<float,float,float> fprop = GetCellAveraged(r,c,0,flowstates,dx);

                if( std::get<0>(fprop) > 0.0001)
                {
                    //set all the correcting ratio's
                    flowstates[0][17]->data[r][c] = 1.0;

                    //FF::cout << "make cache " << n_pix <<  " " << ml << " " << cache_threshold_rel << " " << std::get<0>(fprop) << " "<< std::get<1>(fprop) << " " << std::get<2>(fprop) << std::endl;

                    flowstates[0][4]->data[r][c] = std::get<0>(fprop);
                    flowstates[0][5]->data[r][c] = std::get<1>(fprop);
                    flowstates[0][6]->data[r][c] = std::get<2>(fprop);

                    flowstates[0][7]->data[r][c] = std::get<0>(fprop);
                    flowstates[0][8]->data[r][c] = std::get<1>(fprop);
                    flowstates[0][9]->data[r][c] = std::get<2>(fprop);

                    flowstates[0][18]->data[r][c] = std::get<0>(fprop);
                    flowstates[0][19]->data[r][c] = std::get<1>(fprop);
                    flowstates[0][20]->data[r][c] = std::get<2>(fprop);
                }

            }
        }
    }
}

inline void CheckCacheifyUnDo(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates,cTMapf * HCache, cTMapf * UXCache, cTMapf * UYCache,float courant, int r, int c, int l,float scale_req, float cache_threshold_rel,bool do_ocean, float ocean_level, float rain_fac)
{
    //if we have cache here currently
    if(flowstates[0][17]->data[r][c] > 0.5f)
    {

        //get the ratio's in h, ux, uy
        float ho = flowstates[0][18]->data[r][c];
        float uxo = flowstates[0][19]->data[r][c];
        float uyo = flowstates[0][20]->data[r][c];
        float h = flowstates[0][4]->data[r][c];
        float ux = flowstates[0][5]->data[r][c];
        float uy = flowstates[0][6]->data[r][c];

        //check if the relative change of the values surpasses the threshold
        bool allsame = true;
        if(!rcf(h,ho,cache_threshold_rel,0.0001) || !rcf(ux,uxo,cache_threshold_rel,0.0005) || !rcf(uy,uyo,cache_threshold_rel,0.0005))
        {
            allsame = false;
        }

        //FF::cout << "check uncache " << h << " "<< ho << std::endl;
        //FF::cout << "check check "<<  FF_FULLFLOOD.R27FlowH->data[r][c]/std::max(1e-12f,FF_FULLFLOOD.R27SGHO->data[r][c]) << " " << FF_FULLFLOOD.R27FlowH->data[r][c] << " " << std::max(1e-12f,FF_FULLFLOOD.R27SGHO->data[r][c]) << std::endl;
        if(!allsame)
        {
            //De-activate the Cache

            //FF::cout << "un make cache " << cache_threshold_rel << " " << h/std::max(ho,1e-12f) << " " << ux/((uxo < 0.0?-1.0:1.0)*std::max(std::fabs(uxo),1e-12f)) << " " << uy/((uyo < 0.0?-1.0:1.0)*std::max(std::fabs(uyo),1e-12f)) << std::endl;

            ApplyRatioToGrid(r,c,0,dx,flowstates,h/std::max(ho,1e-12f),ux/((uxo < 0.0?-1.0:1.0)*std::max(std::fabs(uxo),1e-12f)),uy/((uyo < 0.0?-1.0:1.0)*std::max(std::fabs(uyo),1e-12f)));

            flowstates[0][17]->data[r][c] = 0.0;

            flowstates[0][18]->data[r][c] = 0.0;
            flowstates[0][19]->data[r][c] = 0.0;
            flowstates[0][20]->data[r][c] = 0.0;

            //clear old values that we put:
            if(flowstates[0][2]->data[r][c] > 0.5)
            {
                flowstates[0][4]->data[r][c] = 0.0;
                flowstates[0][5]->data[r][c] = 0.0;
                flowstates[0][6]->data[r][c] = 0.0;
            }


        }
    }
}

static inline void FlowGridCachify(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates, cTMapf * HCache, cTMapf * UXCache, cTMapf * UYCache,float courant, float scale_req, double threshold_cache_rel,bool do_ocean, float ocean_level, float rain_fac)
{
    //detect when we want to do subgridding
    //strong gradient in flow surface?
    //strong discharge?


    cTMapf * DEM = flowstates[0][0];
    for(int r = 0; r < DEM->nrRows(); r++)
    {
        for (int c = 0; c < DEM->nrCols(); c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                //for all pixels in this supergrid pixel
                CheckCacheifyDo(_dt,dx,flowstates,HCache, UXCache, UYCache,courant,r,c,0, scale_req, threshold_cache_rel, do_ocean, ocean_level, rain_fac);
            }
        }
    }

    //set the Cache for the next check
    //what is the pixel count within the first level (R27 probably?)
    int flsize = std::pow(3,flowstates.size()-1);

    //for each pixel within this superpixel
    for(int r2 = 0; r2 < HCache->nrRows()*flsize; r2++)
    {
        for(int c2 = 0; c2 < HCache->nrCols()*flsize; c2++)
        {
            if(r2 < HCache->nrRows() && c2 < HCache->nrCols())
            {
                int r3 = r2*3;
                int c3 = c2*3;

                float h = 0.0;
                float ux = 0.0;
                float uy = 0.0;

                for(int i = flowstates.size() -1; i > -1 ; i--)
                {
                    r3 = r3/3;
                    c3 = c3/3;

                    if(i == 0 && flowstates[0][17]->data[r3][c3] > 0.5f)
                    {
                        h *= flowstates[0][4]->data[r3][c3]/(std::max(1e-12f,flowstates[0][18]->data[r3][c3]));
                        ux *= flowstates[0][5]->data[r3][c3]/(flowstates[0][19]->data[r3][c3] > 0.0f? std::max(1e-12f,flowstates[0][19]->data[r3][c3]):std::min(-1e-12f,flowstates[0][19]->data[r3][c3]));
                        uy *= flowstates[0][6]->data[r3][c3]/(flowstates[0][20]->data[r3][c3] > 0.0f? std::max(1e-12f,flowstates[0][20]->data[r3][c3]):std::min(-1e-12f,flowstates[0][20]->data[r3][c3]));
                    }else
                    {
                        //get height and velocity for this layer
                        if(r3 < flowstates[i][0]->nrRows() && c3 < flowstates[i][0]->nrCols() )
                        {
                            h += flowstates[i][4]->data[r3][c3];
                            ux += flowstates[i][5]->data[r3][c3];
                            uy += flowstates[i][6]->data[r3][c3];
                        }
                    }
                }

                HCache->data[r2][c2] = h;
                UXCache->data[r2][c2] = ux;
                UYCache->data[r2][c2] = uy;
            }
        }
    }



}


static inline void FlowGridUnCachify(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates, cTMapf * HCache, cTMapf * UXCache, cTMapf * UYCache,float courant , float scale_req, double threshold_cache_rel,bool do_ocean, float ocean_level, float rain_fac)
{
    //detect when we want to do subgridding
    //strong gradient in flow surface?
    //strong discharge?

    cTMapf * DEM = flowstates[0][0];
    for(int r = 0; r < DEM->nrRows(); r++)
    {
        for (int c = 0; c < DEM->nrCols(); c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                //for all pixels in this supergrid pixel
                CheckCacheifyUnDo(_dt,dx,flowstates, HCache, UXCache, UYCache,courant,r,c,0, scale_req, threshold_cache_rel, do_ocean, ocean_level, rain_fac);
            }
        }
    }

    //detect when we dont want subgridding any longer
}



#endif // FULLFLOOD_CACHE_H
