
//10 is seemingly sufficient as a maximum, but could potentially be less
#define FF_QH_ITERS 10
#define MANNING_MULT_CSCPROF 2.0

inline void SolveManningRectangular(float height,float s, float width, float man, float &Q, float &A_start, float &P_start, int iters = 5)
{

    float s_sqrt = std::sqrt(s);
    float d2_3 = 2.0/3.0;

    //calculate minimal and maximal height/discharge
    float h_min = 0;
    float h_max = height;

    float P_max = width + height * 2.0 + P_start;
    float P_min = width + P_start;
    float A_min = 0.0 + A_start;
    float A_max = width * height + A_start;

    float q_min = A_min * std::pow(A_min/P_min,d2_3)*s_sqrt/man;
    float q_max = A_max * std::pow(A_max/P_max,d2_3)*s_sqrt/man;

    
    if(q_max < Q)
    {
        Q = q_max;
        A_start = A_max;
        P_start = P_max;
        return;
    }
    if(q_min > Q)
    {
        Q = q_min;
        A_start =A_min;
        P_start =P_min;
        return;
    }

    int iter = 0;
    while(iter < iters)
    {
        //best guess for good Q, since Q is proportional to H^5/3
         float h_middle = h_min + (Q-q_min)/std::max(0.000001f,q_max-q_min) * (h_max - h_min);

        float P_middle = width + h_middle * 2.0 + P_start;
        float A_middle = width * h_middle + A_start;

        float q_middle = A_middle * std::pow(A_middle/P_middle,d2_3)*s_sqrt/man;
        
        if(std::fabs(q_middle - Q) < 1e-6)
        {
            Q= Q;
            A_start = A_middle;
            P_start = P_middle;

            return;
        }
        
        //based on which position the value lies, replace either the min or max estimate
        if(q_middle < Q)
        {
            h_min = h_middle;
        }else
        {
            h_max = h_middle;
        }

        if(iter == iters - 1)
        {
            //return best estimate
            Q= Q;
            A_start = A_middle;
            P_start = P_middle;
            return;
        }

        iter ++;
    } 
    
    return;

}

inline void SolveManningTrapezoid(float s_lat,float s, float width, float man, float &Q, float &A_start, float &P_start, int iters = 5)
{
    float slat_inv = 1.0/std::max(0.0001f,s_lat);
    float s_sqrt = std::sqrt(s);
    float d2_3 = 2.0/3.0;

    //calculate minimal and maximal height/discharge
    float h_min = 0;
    float h_max = s_lat * width;

    float P_max = (P_start == 0.0? width:0.0) + 2.0 *std::sqrt(h_max*h_max + width *width) +P_start;
    float P_min = (P_start == 0.0? width:0.0) + P_start;
    float A_min = 0.0 + A_start;
    float A_max = width * h_max + h_max * width +A_start;

    float q_min = A_min * std::pow(A_min/P_min,d2_3)*s_sqrt/man;
    float q_max = A_max * std::pow(A_max/P_max,d2_3)*s_sqrt/man;

    if(q_max < Q)
    {
        Q = q_max;
        A_start = A_max;
        P_start = P_max;
        return;
    }
    if(q_min > Q)
    {
        Q = q_min;
        A_start =A_min;
        P_start =P_min;
        return;
    }

    int iter = 0;
    while(iter < iters)
    {
        //best guess for good Q, since Q is proportional to H^5/3
        //we dont actually do the power law here, its too expensive
        float h_middle = h_min + (Q-q_min)/std::max(0.000001f,q_max-q_min) * (h_max - h_min);

        float w = (slat_inv * h_middle);
        float P_middle =  (P_start == 0.0? width:0.0) + 2.0 *std::sqrt(h_middle*h_middle + w *w ) +P_start;
        float A_middle = (width * h_middle) + h_middle *w   +A_start;

        float q_middle = A_middle * std::pow(A_middle/P_middle,d2_3)*s_sqrt/man;
        
        if(std::fabs(q_middle - Q) < 1e-6)
        {
            Q= Q;
            A_start = A_middle;
            P_start = P_middle;

            return;
        }
        
        //based on which position the value lies, replace either the min or max estimate
        if(q_middle < Q)
        {
            h_min = h_middle;
            q_min = q_middle;
        }else
        {
            h_max = h_middle;
            q_max = q_middle;
        }

        if(iter == iters - 1)
        {
            //return best estimate
            Q= Q;
            A_start = A_middle;
            P_start = P_middle;
            return;
        }

        iter ++;
    } 
    
    return;
}

inline void SolveManningTrapezoidOS(float s_lat,float s, float width, float width_middle, float man, float &Q, float &A_start, float &P_start, int iters = 5, float sides = 2.0)
{
    float slat_inv = 1.0/std::max(0.0001f,s_lat);
    float s_sqrt = std::sqrt(s);
    float d2_3 = 2.0/3.0;

    //calculate minimal and maximal height/discharge
    float h_min = 0;
    float h_max = s_lat * width;

    float P_max = (P_start == 0.0? width:0.0) + sides *std::sqrt(h_max*h_max + width *width) +P_start;
    float P_min = (P_start == 0.0? width:0.0) + P_start;
    float A_min = 0.0 + A_start;
    float A_max = width_middle * h_max + h_max * width * 0.5 * sides +A_start;

    float q_min = A_min * std::pow(A_min/P_min,d2_3)*s_sqrt/man;
    float q_max = A_max * std::pow(A_max/P_max,d2_3)*s_sqrt/man;

    if(Q > 100)
    {
        //FF::cout << "test " << q_min << " " << q_max << " " << h_min << " " << h_max << std::endl;
    }

    if(q_max < Q)
    {
        Q = q_max;
        A_start = A_max;
        P_start = P_max;
        return;
    }
    if(q_min > Q)
    {
        Q = q_min;
        A_start =A_min;
        P_start =P_min;
        return;
    }

    int iter = 0;
    while(iter < iters)
    {
        //best guess for good Q, since Q is proportional to H^5/3
        float h_middle = h_min + (Q-q_min)/std::max(0.000001f,q_max-q_min) * (h_max - h_min);

        float w = (slat_inv * h_middle);
        float P_middle =  (P_start == 0.0? width:0.0) + sides *std::sqrt(h_middle*h_middle + w *w ) +P_start;
        float A_middle = (width_middle * h_middle) + h_middle *w * 0.5 * sides   +A_start;

        float q_middle = A_middle * std::pow(A_middle/P_middle,d2_3)*s_sqrt/man;
        
        if(std::fabs(q_middle - Q) < 1e-6)
        {
            Q= Q;
            A_start = A_middle;
            P_start = P_middle;

            return;
        }
        
        //based on which position the value lies, replace either the min or max estimate
        if(q_middle < Q)
        {
            h_min = h_middle;
            q_min = q_middle;
        }else
        {
            h_max = h_middle;
            q_max = q_middle;
        }

        if(iter == iters - 1)
        {
            //return best estimate
            Q= Q;
            A_start = A_middle;
            P_start = P_middle;
            return;
        }

        iter ++;
    } 
    
    return;
}


inline cTMapf * GetFlowCrossSectionalVolumeFromRiverDischarge(cTMapf * DEM,cTMapf * q,  cTMapf * N, cTMapf * CHANW, cTMapf * CHANH,cTMapf * SLOPE, int iters)
{

    cTMapf * DEMAVI = DEM->GetCopy();
    cTMapf * DEMAVI2 = DEM->GetCopy();
    //cTMapf * WHAVI = DEM->GetCopy0();
    /*cTMapf * SlopeX1I = DEM->GetCopy0();
    cTMapf * SlopeX1IP = DEM->GetCopy0();
    cTMapf * SlopeX2I = DEM->GetCopy0();
    cTMapf * SlopeX2IP = DEM->GetCopy0();
    cTMapf * SlopeY1I = DEM->GetCopy0();
    cTMapf * SlopeY1IP = DEM->GetCopy0();
    cTMapf * SlopeY2I = DEM->GetCopy0();
    cTMapf * SlopeY2IP = DEM->GetCopy0();

    cTMapf * VMAXX1I = DEM->GetCopy0();
    cTMapf * VMAXX2I = DEM->GetCopy0();
    cTMapf * VMAXY1I = DEM->GetCopy0();
    cTMapf * VMAXY2I = DEM->GetCopy0();*/
    cTMapf * WidthI = DEM->GetCopy0();
    cTMapf * q_estI = q->GetCopy0();
    cTMapf * A_estI = q->GetCopy0();
    cTMapf * P_estI = q->GetCopy0();


    //cTMapf * DIRX = DEM->GetCopy0();
    //cTMapf * DIRY = DEM->GetCopy0();
    //cTMapf * DIRW = DEM->GetCopy0();

    float dx = std::fabs(DEM->cellSizeX());

   
    for(int i = 0; i < iters; i++)
    {
        //we look 2^i pixels away for slope
        float dxi = dx;

        //current with in logaritmic search
        dxi = dx * (1.0);
        if(i > 0)
        {
            dxi =dx * ((float)(std::pow(3,i)));
        }
        int idist = 1;
        if(i>0)
        {
            idist = std::pow(3,i);
        }


        #pragma omp parallel for collapse(2)
        for (int r = 0; r < DEM->Rows(); r++)
        {
            for(int c = 0; c < DEM->Cols(); c++)
            {
                if(!pcr::isMV(DEM->data[r][c]) && (!(std::fabs(q_estI->data[r][c] - q->data[r][c]) < 1e-4)))
                {
                    float dem = DEM->data[r][c];
                    float z= dem;
                    float zx1;
                    float zx2;
                    float zy1;
                    float zy2;



                    if(i == 0)
                    {
                        if(c > 0)
                        {
                            zx1 = pcr::isMV(DEM->data[r][c-1])?dem:DEM->data[r][c-1];
                        }else {
                            zx1 = dem;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            zx2 = pcr::isMV(DEM->data[r][c+1])?dem:DEM->data[r][c+1];
                        }else {
                            zx2 = dem;
                        }
                        if(r > 0)
                        {
                            zy1 = pcr::isMV(DEM->data[r-1][c])?dem:DEM->data[r-1][c];
                        }else {
                            zy1 = dem;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            zy2 = pcr::isMV(DEM->data[r+1][c])?dem:DEM->data[r+1][c];
                        }else {
                            zy2 = dem;

                        }
                    }else
                    {
                        dem = DEMAVI->data[r][c];

                        if(c > idist)
                        {
                            zx1 = pcr::isMV(DEMAVI->data[r][c-idist])?dem:DEMAVI->data[r][c-idist];
                        }else {
                            zx1 = dem;
                        }
                        if(c < DEM->nrCols()-1-idist)
                        {
                            zx2 = pcr::isMV(DEMAVI->data[r][c+idist])?dem:DEMAVI->data[r][c+idist];
                        }else {
                            zx2 = dem;
                        }
                        if(r > idist)
                        {
                            zy1 = pcr::isMV(DEMAVI->data[r-idist][c])?dem:DEMAVI->data[r-idist][c];
                        }else {
                            zy1 = dem;
                        }
                        if(r < DEM->nrRows()-1-idist)
                        {
                            zy2 = pcr::isMV(DEMAVI->data[r+idist][c])?dem:DEMAVI->data[r+idist][c];
                        }else {
                            zy2 = dem;

                        }
                    }

                    float sx1 =(zx1 -z)/dxi;
                    float sy1 =(zy1 -z)/dxi;
                    float sx2 =(zx2 -z)/dxi;
                    float sy2 =(zy2 -z)/dxi;

                    //get direction of steepest descent
                    float ux = (sx1>sx2?-1.0 * sx1:1.0*sx2);
                    float uy = (sy1>sy2?-1.0 * sy1:1.0*sy2);
                    float norm = std::sqrt(ux*ux+uy*uy);
                    ux= ux/std::max(1e-12f,norm);
                    uy= uy/std::max(1e-12f,norm);

                    //get left and right slope values (multiply vectors by 2d rotation matrix of 90 and -90 degrees)
                    float ulx = -uy;
                    float uly = ux;
                    float urx = uy;
                    float ury = -ux;
                    float sl= ((ulx > 0.0?std::fabs(ulx)*sx2: std::fabs(ulx)*sx1)) + std::fabs((uly > 0.0?std::fabs(uly)*sy2: std::fabs(uly)*sy1));
                    float sr= ((urx > 0.0?std::fabs(urx)*sx2: std::fabs(urx)*sx1)) + std::fabs((ury > 0.0?std::fabs(ury)*sy2: std::fabs(ury)*sy1));
                    float zl = z + dxi * sl;
                    float zr = z + dxi * sr;

                    //did we already meet a volume larger then ours there?

                    //redistribute the volume and see if we have left-over
                

                    sl = sl * 4.0/3.0;
                    sr = sr * 4.0/3.0;
                    float s = std::max(1e-6,std::fabs(0.5*(sl + sr)));

                    //check channel
                    if( i == 0 && CHANW != nullptr)
                    {
                        if(CHANW->data[r][c] > 1e-10 && CHANH->data[r][c] > 1e-10)
                        {
                            float QN = q->data[r][c];
                            float AN = 0.0;
                            float PN = 0.0;
                            SolveManningRectangular(CHANH->data[r][c],std::max(0.001f,SLOPE->data[r][c]), CHANW->data[r][c], (N->data[r][c]*MANNING_MULT_CSCPROF), QN,AN,PN,10);

                            q_estI->data[r][c] = std::min(q->data[r][c],QN);
                            A_estI->data[r][c] = AN;
                            P_estI->data[r][c] = PN;
                        }
                    }

                    //now check neighboring cells
                    float QN = q->data[r][c];
                    float AN = 0.0;
                    float PN = 0.0;
                    SolveManningTrapezoid(s,std::max(0.001f,SLOPE->data[r][c]), dxi, N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,10);
                    
                    if(q->data[r][c] > 100)
                    {
                        //FF::cout << "do_trapezoid " << i << " " << r << " " << c << "  : " <<q->data[r][c] << " " <<  QN << " " << SLOPE->data[r][c] << " " << N->data[r][c] << " " << s << " " << dx << " " << AN << " " << AN/dx << std::endl;
                    }

                    q_estI->data[r][c] = std::min(q->data[r][c],QN);
                    A_estI->data[r][c] = AN;
                    P_estI->data[r][c] = PN;

                    //save left-overs
                    float z_avg = (zl + zr + z)/3.0;
                    DEMAVI2->data[r][c] = z_avg;
                }
            }
        }

        cTMapf * temp = DEMAVI;
        DEMAVI = DEMAVI2;
        DEMAVI2 = temp;
    }

    /*#pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++)
    {
        for(int c = 0; c < DEM->Cols(); c++)
        {
            if(!pcr::isMV(DEM->data[r][c]) && DIRW->data[r][c] > 1e-8)
            {
                DIRX->data[r][c] = DIRX->data[r][c]/DIRW->data[r][c];
                DIRY->data[r][c] = DIRY->data[r][c]/DIRW->data[r][c];
            }
        }
    }*/


    //delete DIRW;

    delete DEMAVI;
    delete DEMAVI2;
    /*delete WHAVI;
    delete SlopeX1I;
    delete SlopeX1IP;
    delete SlopeX2I;
    delete SlopeX2IP;
    delete SlopeY1I;
    delete SlopeY1IP;
    delete SlopeY2I;
    delete SlopeY2IP;

    delete VMAXX1I;
    delete VMAXX2I;
    delete VMAXY1I;
    delete VMAXY2I;*/
    delete q_estI;
    delete P_estI;

    return A_estI;

}


inline float middleOfThree(float a, float b, float c)
{
    // Checking for b
    if ((a < b && b < c) || (c < b && b < a))
       return b;
 
    // Checking for a
    else if ((b < a && a < c) || (c < a && a < b))
       return a;
 
    else
       return c;
}


inline std::vector<cTMapf *> GetFlowCrossSectionalVolumeFromRiverDischargePreDistribute(cTMapf * DEM,cTMapf * q,  cTMapf * N, cTMapf * CHANW, cTMapf * CHANH,cTMapf * SLOPE, int iters)
{
    cTMapf * WidthI = DEM->GetCopy0();
    cTMapf * q_estI = q->GetCopy0();
    cTMapf * A_estI = q->GetCopy0();
    cTMapf * P_estI = q->GetCopy0();


    float dx = std::fabs(DEM->cellSizeX());

    std::vector<cTMapf *> DEM_MAXI;
    for(int i = 0; i < iters; i++)
    {
        DEM_MAXI.push_back(DEM->GetCopy0());
    }

    std::vector<cTMapf *> DEM_MAXL;
    for(int i = 0; i < iters; i++)
    {
        DEM_MAXL.push_back(DEM->GetCopy0());
    }

    std::vector<cTMapf *> DEM_MAXR;
    for(int i = 0; i < iters; i++)
    {
        DEM_MAXR.push_back(DEM->GetCopy0());
    }

    std::vector<cTMapf *> SR;
    for(int i = 0; i < iters; i++)
    {
        SR.push_back(DEM->GetCopy0());
    }


    std::vector<cTMapf *> SL;
    for(int i = 0; i < iters; i++)
    {
        SL.push_back(DEM->GetCopy0());
    }


    std::vector<cTMapf *> DEMAVI;
    for(int i = 0; i < iters; i++)
    {
        DEMAVI.push_back(DEM->GetCopy0());
    }

    std::vector<cTMapf *> VolI;
    for(int i = 0; i < iters; i++)
    {
        VolI.push_back(DEM->GetCopy0());
    }

    std::vector<cTMapf *> WX1;
    for(int i = 0; i < iters; i++)
    {
        WX1.push_back(DEM->GetCopy0());
    }
 std::vector<cTMapf *> WX2;
    for(int i = 0; i < iters; i++)
    {
        WX2.push_back(DEM->GetCopy0());
    }
     std::vector<cTMapf *> WY1;
    for(int i = 0; i < iters; i++)
    {
        WY1.push_back(DEM->GetCopy0());
    }
     std::vector<cTMapf *> WY2;
    for(int i = 0; i < iters; i++)
    {
        WY2.push_back(DEM->GetCopy0());
    }
   
    for(int i = 0; i < iters; i++)
    {
        //we look 2^i pixels away for slope
        float dxi = dx;

        //current with in logaritmic search
        dxi = dx * (1.0);
        if(i > 0)
        {
            dxi =dx * ((float)(std::pow(3,i)));
        }
        int idist = 1;
        if(i>0)
        {
            idist = std::pow(3,i);
        }


        #pragma omp parallel for collapse(2)
        for (int r = 0; r < DEM->Rows(); r++)
        {
            for(int c = 0; c < DEM->Cols(); c++)
            {
                if(!pcr::isMV(DEM->data[r][c]))
                {
                    float dem = DEM->data[r][c];
                    float z= dem;
                    float zx1;
                    float zx2;
                    float zy1;
                    float zy2;

                    float zx1y1;
                    float zx1y2;
                    float zx2y1;
                    float zx2y2;


                    float zmax_x1 = dem;
                    float zmax_y1 = dem;
                    float zmax_x2 = dem;
                    float zmax_y2 = dem;

                    if(i == 0)
                    {
                        if(c > 0)
                        {
                            zx1 = pcr::isMV(DEM->data[r][c-1])?dem:DEM->data[r][c-1];
                        }else {
                            zx1 = dem;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            zx2 = pcr::isMV(DEM->data[r][c+1])?dem:DEM->data[r][c+1];
                        }else {
                            zx2 = dem;
                        }
                        if(r > 0)
                        {
                            zy1 = pcr::isMV(DEM->data[r-1][c])?dem:DEM->data[r-1][c];
                        }else {
                            zy1 = dem;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            zy2 = pcr::isMV(DEM->data[r+1][c])?dem:DEM->data[r+1][c];
                        }else {
                            zy2 = dem;
                        }

                        if(c > 0 && c < DEM->nrCols() - 1 && r > 0 && r < DEM->nrRows() - 1)
                        {
                            zx1y1 = pcr::isMV(DEM->data[r-1][c-1])?dem:DEM->data[r-1][c-1];
                            zx2y1 = pcr::isMV(DEM->data[r-1][c+1])?dem:DEM->data[r-1][c+1];
                            zx1y2 = pcr::isMV(DEM->data[r+1][c-1])?dem:DEM->data[r+1][c-1];
                            zx2y2 = pcr::isMV(DEM->data[r+1][c+1])?dem:DEM->data[r+1][c+1];

                        }

                        zmax_x1 = std::max(std::max(zx1y1,zx1),zx1y2);
                        zmax_x2 = std::max(std::max(zx2y1,zx2),zx2y2);
                        float zmax_x =std::max(std::max(z,zx1),zx2);
                        zmax_y1 = std::max(std::max(zx2y1,zy1),zx1y1);
                        zmax_y2 = std::max(std::max(zx1y2,zy2),zx2y2);
                        float zmax_y =std::max(std::max(z,zy1),zy2);

                        float zm_y = std::min(zmax_x,std::min(zmax_x1,zmax_x2));
                        float zm_x = std::min(zmax_y,std::min(zmax_y1,zmax_y2));

                        DEM_MAXI[i]->data[r][c] = std::max(zm_y,zm_x);
                    
                    }else
                    {
                        dem = DEMAVI[i-1]->data[r][c];

                        if(c > idist)
                        {
                            zx1 = pcr::isMV(DEMAVI[i-1]->data[r][c-idist])?dem:DEMAVI[i-1]->data[r][c-idist];
                        }else {
                            zx1 = dem;
                        }
                        if(c < DEM->nrCols()-1-idist)
                        {
                            zx2 = pcr::isMV(DEMAVI[i-1]->data[r][c+idist])?dem:DEMAVI[i-1]->data[r][c+idist];
                        }else {
                            zx2 = dem;
                        }
                        if(r > idist)
                        {
                            zy1 = pcr::isMV(DEMAVI[i-1]->data[r-idist][c])?dem:DEMAVI[i-1]->data[r-idist][c];
                        }else {
                            zy1 = dem;
                        }
                        if(r < DEM->nrRows()-1-idist)
                        {
                            zy2 = pcr::isMV(DEMAVI[i-1]->data[r+idist][c])?dem:DEMAVI[i-1]->data[r+idist][c];
                        }else {
                            zy2 = dem;

                        }


                         if(c > idist)
                        {
                            zmax_x1 = pcr::isMV(DEM_MAXI[i-1]->data[r][c-idist])?dem:DEM_MAXI[i-1]->data[r][c-idist];
                        }else {
                            zmax_x1 = dem;
                        }
                        if(c < DEM->nrCols()-1-idist)
                        {
                            zmax_x2 = pcr::isMV(DEM_MAXI[i-1]->data[r][c+idist])?dem:DEM_MAXI[i-1]->data[r][c+idist];
                        }else {
                            zmax_x2 = dem;
                        }
                        if(r > idist)
                        {
                            zmax_y1 = pcr::isMV(DEM_MAXI[i-1]->data[r-idist][c])?dem:DEM_MAXI[i-1]->data[r-idist][c];
                        }else {
                            zmax_y1 = dem;
                        }
                        if(r < DEM->nrRows()-1-idist)
                        {
                            zmax_y2 = pcr::isMV(DEM_MAXI[i-1]->data[r+idist][c])?dem:DEM_MAXI[i-1]->data[r+idist][c];
                        }else {
                            zmax_y2 = dem;

                        }


                        if(c >  idist && c < DEM->nrCols() -1-idist && r >  idist && r < DEM->nrRows()-1-idist)
                        {
                            zx1y1 = pcr::isMV(DEM_MAXI[i-1]->data[r- idist][c- idist])?dem:DEM_MAXI[i-1]->data[r- idist][c- idist];
                            zx2y1 = pcr::isMV(DEM_MAXI[i-1]->data[r- idist][c+ idist])?dem:DEM_MAXI[i-1]->data[r- idist][c+ idist];
                            zx1y2 = pcr::isMV(DEM_MAXI[i-1]->data[r+ idist][c- idist])?dem:DEM_MAXI[i-1]->data[r+ idist][c- idist];
                            zx2y2 = pcr::isMV(DEM_MAXI[i-1]->data[r+ idist][c+ idist])?dem:DEM_MAXI[i-1]->data[r+ idist][c+ idist];

                        }

                        zmax_x1 = std::max(std::max(zx1y1,zx1),zx1y2);
                        zmax_x2 = std::max(std::max(zx2y1,zx2),zx2y2);
                        float zmax_x =std::max(std::max(z,zx1),zx2);
                        zmax_y1 = std::max(std::max(zx2y1,zy1),zx1y1);
                        zmax_y2 = std::max(std::max(zx1y2,zy2),zx2y2);
                        float zmax_y =std::max(std::max(z,zy1),zy2);

                        float zm_y = std::min(zmax_x,std::min(zmax_x1,zmax_x2));
                        float zm_x = std::min(zmax_y,std::min(zmax_y1,zmax_y2));

                        DEM_MAXI[i]->data[r][c] = std::max(zm_y,zm_x);

                        //DEM_MAXI[i]->data[r][c] = std::max(std::max(zmax_x1,zmax_x2),std::max(zmax_y1,zmax_y2));

                    }

                    float sx1 =(zx1 -z)/dxi;
                    float sy1 =(zy1 -z)/dxi;
                    float sx2 =(zx2 -z)/dxi;
                    float sy2 =(zy2 -z)/dxi;

                    //get direction of steepest descent
                    float ux = (sx1<sx2?-1.0 * sx1:1.0*sx2);
                    float uy = (sy1<sy2?-1.0 * sy1:1.0*sy2);
                    float norm = std::sqrt(ux*ux+uy*uy);
                    if(norm < 1e-12)
                    {
                        ux = 1.0;
                        uy = 0.0;
                        norm = 1.0;
                    }
                    ux= ux/std::max(1e-12f,norm);
                    uy= uy/std::max(1e-12f,norm);

                    //get left and right slope values (multiply vectors by 2d rotation matrix of 90 and -90 degrees)
                    float ulx = -uy;
                    float uly = ux;
                    float urx = uy;
                    float ury = -ux;

                    float sl= ((ulx > 0.0?std::fabs(ulx)*sx2: std::fabs(ulx)*sx1)) + std::fabs((uly > 0.0?std::fabs(uly)*sy2: std::fabs(uly)*sy1));
                    float sr= ((urx > 0.0?std::fabs(urx)*sx2: std::fabs(urx)*sx1)) + std::fabs((ury > 0.0?std::fabs(ury)*sy2: std::fabs(ury)*sy1));
                    float zl = z + dxi * sl;
                    float zr = z + dxi * sr;
                    float zml = ((ulx > 0.0? ulx * zmax_x2 : -ulx * zmax_x1) + (uly > 0.0? uly * zmax_y2 : -uly * zmax_y1))/std::max(1e-10f,std::fabs(uly) + std::fabs(ulx));
                    float zmr = ((urx > 0.0? urx * zmax_x2 : -urx * zmax_x1) + (ury > 0.0? ury * zmax_y2 : -ury * zmax_y1))/std::max(1e-10f,std::fabs(ury) + std::fabs(urx));

                    SL[i]->data[r][c] = std::max(0.001f,sl);
                    SR[i]->data[r][c] = std::max(0.001f,sr);

                    DEM_MAXR[i]->data[r][c] = zmr;
                    DEM_MAXL[i]->data[r][c] = zml;



                    WX1[i]->data[r][c] = ulx/(std::fabs(std::fabs(ulx) + std::fabs(uly)));
                    WX2[i]->data[r][c] = urx/(std::fabs(std::fabs(urx) + std::fabs(ury)));
                    WY1[i]->data[r][c] = uly/(std::fabs(std::fabs(ulx) + std::fabs(uly)));
                    WY2[i]->data[r][c] = ury/(std::fabs(std::fabs(urx) + std::fabs(ury)));


                    //save left-overs
                    float z_avg = (zx1 + zx2 + z + zy1 + zy2 )/5.0;
                    DEMAVI[i]->data[r][c] = z_avg;

                }
            }
        }
    }


        #pragma omp parallel for collapse(2)
        for (int r = 0; r < DEM->Rows(); r++)
        {
            for(int c = 0; c < DEM->Cols(); c++)
            {
                if(!pcr::isMV(DEM->data[r][c]) && (!(std::fabs(q_estI->data[r][c] - q->data[r][c]) < 1e-4)))
                {

                    float slope = std::max(0.001f,SLOPE->data[r][c]);
                    float dem = DEM->data[r][c];
                    int i_left = 0;
                    int i_right = 0;
                    bool left_blocked = false;
                    bool right_blocked = false;
                    float dx_middle = dx;

                    float h_left = 0.0 ;
                    float h_right = 0.0;

                    while(i_left < iters && i_right < iters)//(i_left < iters && i_right < iters)
                    {
                        if((std::fabs(q_estI->data[r][c] - q->data[r][c]) < 1e-4))
                        {
                            break;
                        }


                        //we look 3^i pixels away for slope
                        float dxi_left = dx;
                        if(i_left > 0)
                        {
                            dxi_left =dx * ((float)(std::pow(3,i_left)));
                        }
                        int idist_left = 1;
                        if(i_left>0)
                        {
                            idist_left = std::pow(3,i_left);
                        }
                        float dxi_right = dx;
                        if(i_right > 0)
                        {
                            dxi_right =dx * ((float)(std::pow(3,i_right)));
                        }
                        int idist_right = 1;
                        if(i_right>0)
                        {
                            idist_right = std::pow(3,i_right);
                        }

                        float dxi_left_p = dxi_left;
                        if(i_left > 0)
                        {
                            dxi_left_p =dx * ((float)(std::pow(3,i_left-1)));
                        }else
                        {
                            dxi_left_p = dx;
                        }
                        int idist_left_p = idist_left;
                        if(i_left>0)
                        {
                            idist_left_p = std::pow(3,i_left-1);
                        }else
                        {
                            idist_left_p = 0;
                        }
                        float dxi_right_p = dxi_right;
                        if(i_right > 0)
                        {
                            dxi_right_p =dx * ((float)(std::pow(3,i_right-1)));
                        }else{
                            dxi_right_p = dx;
                        }
                        int idist_right_p = idist_right_p;
                        if(i_right>0)
                        {
                            idist_right_p = std::pow(3,i_right-1);
                        }else{
                            idist_right_p = 0;
                        }

                        int i_left_p = std::max(0,i_left-1);
                        int i_right_p = std::max(0,i_right-1);

                        //get current left and right slope, and current left and right barrier height, if any
                        float sl= SL[i_left]->data[r][c];
                        float sr= SR[i_right]->data[r][c];

                        float zml = DEM_MAXL[i_left]->data[r][c];
                        float zmr = DEM_MAXR[i_right]->data[r][c];

                        //how much does the barrier extend above our average piece-wise linear approximation?
                        float demav_left = i_left > 0?DEMAVI[i_left-1]->data[r][c]:dem;
                        float demav_right = i_right > 0? DEMAVI[i_right-1]->data[r][c] : dem;
                        float barrier_left = std::max(0.0f,zml - (demav_left));
                        float barrier_right = std::max(0.0f,zmr - (demav_right));

                        if(q->data[r][c] > 1000)
                        {
                            //FF::cout << "bh "  << q->data[r][c] << " " <<   i_left << " " << i_right <<  " " << dem << " " << demav_left << " " << demav_right << " " <<  barrier_left << " " << barrier_right << std::endl;
                        }

                        //check channel
                        if( i_left == 0 && i_right == 0 && CHANW != nullptr)
                        {
                            if(CHANW->data[r][c] > 1e-10 && CHANH->data[r][c] > 1e-10)
                            {
                                float QN = q->data[r][c];
                                float AN = 0.0;
                                float PN = 0.0;
                                SolveManningRectangular(CHANH->data[r][c],slope, CHANW->data[r][c], (N->data[r][c]*MANNING_MULT_CSCPROF), QN,AN,PN,10);

                                q_estI->data[r][c] = std::min(q->data[r][c],QN);
                                A_estI->data[r][c] = AN;
                                P_estI->data[r][c] = PN;

                                //add volume here
                                VolI[0]->data[r][c] += AN/dx;
                                
                            }
                        }

                        //we do a symmetrical trapezoid
                        sl = 0.5 * sl * 4.0/3.0;
                        sr = 0.5 * sr * 4.0/3.0;
                        float s = std::max(1e-3, std::fabs(0.5*(sl + sr)));

                        bool left_is_blocking = false;
                        bool right_is_blocking = false;

                        //if not 0 iteration, check if barriers contain?
                        if(i_left > 0 && i_right > 0)
                        {
                            //rectangular channel Q

                            if(barrier_left > 0.0)
                            {
                                float AN = A_estI->data[r][c] + 2.0 * barrier_left * (dxi_left_p + dxi_right_p);
                                float PN = P_estI->data[r][c] + 2.0 * barrier_left;

                                float QN = AN * std::pow(AN/PN,2.0/3.0)*std::sqrt(slope)/(N->data[r][c]*MANNING_MULT_CSCPROF);
                                if(QN > q->data[r][c])
                                {
                                    left_is_blocking = true;
                                }
                            }

                            if(barrier_right > 0.0)
                            {
                                float AN = A_estI->data[r][c] + 2.0 * barrier_right * (dxi_left_p + dxi_right_p);
                                float PN = P_estI->data[r][c] + 2.0 * barrier_right;

                                float QN = AN * std::pow(AN/PN,2.0/3.0)*std::sqrt(slope)/(N->data[r][c]*MANNING_MULT_CSCPROF);
                                if(QN > q->data[r][c])
                                {
                                    right_is_blocking = true;
                                }

                            }
                        }

                        //if so, set area and break
                        if(left_is_blocking && right_is_blocking)
                        {

                            //solve rectangular channel
                            float QN = q->data[r][c];
                            float AN = A_estI->data[r][c];
                            
                            float AN_O = AN;
                            float PN = P_estI->data[r][c];
                            SolveManningRectangular(std::max(barrier_left, barrier_right),slope, dxi_left_p + dxi_right_p, N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,10);
                            q_estI->data[r][c] = std::min(q->data[r][c],QN);
                            A_estI->data[r][c] = AN;
                            P_estI->data[r][c] = PN;

                            //area to volume
                            float vol = (AN - AN_O)/(dxi_left_p + dxi_right_p);

                            if(q->data[r][c] > 1000.0)
                            {
                                FF::cout << "blockb " << barrier_left << " " << barrier_right << " " << vol << " " << dxi_left_p + dxi_right_p << " " << dx_middle << q_estI->data[r][c] << " " << q->data[r][c] << " " << QN << std::endl;
                            }

                            //divide over left and right
                            //float w_l = dxi_left_p /std::max(0.0001f,(dxi_left_p + dxi_right_p));
                            //float w_r = dxi_right_p /std::max(0.0001f,(dxi_left_p + dxi_right_p));

                            //we look 3^i pixels away for slope
                            
                            //left x weight
                            //left x offset
                            if((c + (WX1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left_p > -1) && (c + (WX1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left_p < DEM->nrCols()))
                            {
                                VolI[i_left_p]->data[r][c + (WX1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left_p] += (dx/dxi_left_p) *vol * std::fabs(WX1[i_left_p]->data[r][c]);
                                //VolI[i_left_p]->data[r][c] += (dx/dxi_left_p) *vol * std::fabs(WX1[i_left_p]->data[r][c]);
                            }
                            if((c + (WX2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right_p > -1) && (c +  (WX2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right_p < DEM->nrCols()))
                            {
                                VolI[i_right_p]->data[r][c + (WX2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right_p] += (dx/dxi_right_p) *vol * std::fabs(WX2[i_right_p]->data[r][c]);
                                //VolI[i_right_p]->data[r][c] += (dx/dxi_right_p) *vol * std::fabs(WX2[i_right_p]->data[r][c]);
                            }

                            if((r + (WY1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left_p > -1) && (r + (WY1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left_p < DEM->nrRows()))
                            {
                                VolI[i_left_p]->data[r + (WY1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left_p][c] += (dx/dxi_left_p) *vol * std::fabs(WY1[i_left_p]->data[r][c]);
                                //VolI[i_left_p]->data[r][c] += (dx/dxi_left_p) *vol * std::fabs(WY1[i_left_p]->data[r][c]);
                            }
                            if((r + (WY2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right_p > -1) && (r +  (WY2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right_p < DEM->nrRows()))
                            {
                                VolI[i_right_p]->data[r + (WY2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right_p][c] += (dx/dxi_right_p) *vol * std::fabs(WY2[i_right_p]->data[r][c]);
                                //VolI[i_right_p]->data[r][c] += (dx/dxi_right_p) *vol * std::fabs(WY2[i_right_p]->data[r][c]);
                            }

                            break;


                        }

                        //pick a side, do one-sided
                        //prefer non-blocking
                        //then prefer lower i
                        bool do_left = true;
                        if(left_is_blocking || (barrier_left > barrier_right))
                        {
                            do_left = false;
                        }

                        //run one-sided trapezoid 
                        if(i_left == 0 && i_right == 0)
                        {

                            float QN = q->data[r][c];
                            float AN = A_estI->data[r][c];
                            float AO = AN;
                            float PN = P_estI->data[r][c];
                            SolveManningTrapezoid(s,slope, dxi_left,N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,10);
                            q_estI->data[r][c] = std::min(q->data[r][c],QN);
                            A_estI->data[r][c] = AN;
                            P_estI->data[r][c] = PN;

                            i_left = i_left + 1;
                            i_right = i_right + 1;
                            h_left += s * dxi_left;
                            h_right += s * dxi_left;
                            dx_middle += 2.0 * dx;                           
                        
                            float vol = (AN-AO)/(3.0*dx);
                            //add volume to cascades
                            //get all weights 

                            //normalize to 1

                            //add them

                            if((c + (WX1[i_left]->data[r][c]>0.0?1:-1) * idist_left > -1) && (c + (WX1[i_left]->data[r][c]>0.0?1:-1) * idist_left < DEM->nrCols()))
                            {
                                VolI[0]->data[r][c + (WX1[i_left]->data[r][c]>0.0?1:-1) * 1] += (dx/dxi_left) *vol * std::fabs(WX1[i_left]->data[r][c]);
                            }
                            if((c + (WX2[i_right]->data[r][c]>0.0?1:-1) * idist_right > -1) && (c +  (WX2[i_right]->data[r][c]>0.0?1:-1) * idist_right < DEM->nrCols()))
                            {
                                VolI[0]->data[r][c + (WX2[i_right]->data[r][c]>0.0?1:-1) * 1] += (dx/dxi_right) *vol * std::fabs(WX2[i_right]->data[r][c]);
                            }

                            if((r + (WY1[i_left]->data[r][c]>0.0?1:-1) * idist_left > -1) && (r + (WY1[i_left]->data[r][c]>0.0?1:-1) * idist_left < DEM->nrRows()))
                            {
                                VolI[0]->data[r + (WY1[i_left]->data[r][c]>0.0?1:-1) * 1][c] += (dx/dxi_left) *vol * std::fabs(WY1[i_left]->data[r][c]);
                            }
                            if((r + (WY2[i_right]->data[r][c]>0.0?1:-1) * idist_right > -1) && (r +  (WY2[i_right]->data[r][c]>0.0?1:-1) * idist_right < DEM->nrRows()))
                            {
                                VolI[0]->data[r + (WY2[i_right]->data[r][c]>0.0?1:-1) * 1][c] += (dx/dxi_right) *vol * std::fabs(WY2[i_right]->data[r][c]);
                            }


                        }else
                        {
                           
                           if(do_left)
                           {
                                float QN = q->data[r][c];
                                float AN = A_estI->data[r][c];
                                float AO = AN;
                                float PN = P_estI->data[r][c];
                                SolveManningTrapezoidOS(sl,slope, dxi_left, dx_middle,  N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,10,1.0);
                                q_estI->data[r][c] = std::min(q->data[r][c],QN);
                                A_estI->data[r][c] = AN;
                                P_estI->data[r][c] = PN;

                                h_left += s * dxi_left;
                                h_right += s * dxi_left;
                        
                                float vol = (AN-AO)/(dxi_right_p + dxi_left);
                                
                                dx_middle += dxi_left;   

                                if((c + (WX1[i_left]->data[r][c]>0.0?1:-1) * idist_left > -1) && (c + (WX1[i_left]->data[r][c]>0.0?1:-1) * idist_left < DEM->nrCols()))
                                {
                                    VolI[i_left]->data[r][c + (WX1[i_left]->data[r][c]>0.0?1:-1) * idist_left] += (dx/dxi_left) *vol * std::fabs(WX1[i_left]->data[r][c]);
                                }
                                if((c + (WX2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right > -1) && (c +  (WX2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right < DEM->nrCols()))
                                {
                                    //VolI[i_right_p]->data[r][c + (WX2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right] += (dx/dxi_right_p) *vol * std::fabs(WX2[i_right_p]->data[r][c]);
                                    //VolI[i_right_p]->data[r][c] += (dx/dxi_right_p) *vol * std::fabs(WX2[i_right_p]->data[r][c]);
                                }

                                if((r + (WY1[i_left]->data[r][c]>0.0?1:-1) * idist_left > -1) && (r + (WY1[i_left]->data[r][c]>0.0?1:-1) * idist_left < DEM->nrRows()))
                                {
                                    VolI[i_left]->data[r + (WY1[i_left]->data[r][c]>0.0?1:-1) * idist_left][c] += (dx/dxi_left) *vol * std::fabs(WY1[i_left]->data[r][c]);
                                }
                                if((r + (WY2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right > -1) && (r +  (WY2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right < DEM->nrRows()))
                                {
                                    //VolI[i_right_p]->data[r + (WY2[i_right_p]->data[r][c]>0.0?1:-1) * idist_right][c] += (dx/dxi_right_p) *vol * std::fabs(WY2[i_right_p]->data[r][c]);
                                    //VolI[i_right_p]->data[r][c] += (dx/dxi_right_p) *vol * std::fabs(WY2[i_right_p]->data[r][c]);
                                }

                                
                                i_left = i_left + 1;    

                           }else{
                                float QN = q->data[r][c];
                                float AN = A_estI->data[r][c];
                                float AO = AN;
                                float PN = P_estI->data[r][c];
                                SolveManningTrapezoidOS(sr,slope, dxi_right, dx_middle,  N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,10,1.0);
                                q_estI->data[r][c] = std::min(q->data[r][c],QN);
                                A_estI->data[r][c] = AN;
                                P_estI->data[r][c] = PN;


                                h_left += s * dxi_right;    
                                h_right += s * dxi_right;                   
                        
                                float vol = (AN-AO)/(dxi_left_p + dxi_right);
                                
                                dx_middle += dxi_right;   

                                if((c + (WX1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left > -1) && (c + (WX1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left < DEM->nrCols()))
                                {
                                    //VolI[i_left_p]->data[r][c + (WX1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left] += (dx/dxi_left_p) * vol * std::fabs(WX1[i_left_p]->data[r][c]);

                                    //VolI[i_left_p]->data[r][c] += (dx/dxi_left_p) * vol * std::fabs(WX1[i_left_p]->data[r][c]);
                                    

                                }
                                if((c + (WX2[i_right]->data[r][c]>0.0?1:-1) * idist_right > -1) && (c +  (WX2[i_right]->data[r][c]>0.0?1:-1) * idist_right < DEM->nrCols()))
                                {
                                    VolI[i_right]->data[r][c + (WX2[i_right]->data[r][c]>0.0?1:-1) * idist_right] += (dx/dxi_right) * vol * std::fabs(WX2[i_right]->data[r][c]);
                                }

                                if((r + (WY1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left > -1) && (r + (WY1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left < DEM->nrRows()))
                                {
                                    //VolI[i_left_p]->data[r + (WY1[i_left_p]->data[r][c]>0.0?1:-1) * idist_left][c] += (dx/dxi_left_p) *vol * std::fabs(WY1[i_left_p]->data[r][c]);

                                    //VolI[i_left_p]->data[r][c] += (dx/dxi_left_p) *vol * std::fabs(WY1[i_left_p]->data[r][c]);
                                    
                                }
                                if((r + (WY2[i_right]->data[r][c]>0.0?1:-1) * idist_right > -1) && (r +  (WY2[i_right]->data[r][c]>0.0?1:-1) * idist_right < DEM->nrRows()))
                                {
                                    VolI[i_right]->data[r + (WY2[i_right]->data[r][c]>0.0?1:-1) * idist_right][c] += (dx/dxi_right) * vol * std::fabs(WY2[i_right]->data[r][c]);
                                }


                                i_right = i_right + 1;     

                           }

                        }
                    }
                }
            }
        }

    //now iterate backwards to re-distribute the water 
    for(int i = iters-1; i > 0; i--)
    {
        int idist = 1;
        if(i>0)
        {
            idist = std::pow(3,i-1);
        }

        #pragma omp parallel for collapse(2)
            for (int r = 0; r < DEM->Rows(); r++)
            {
                for(int c = 0; c < DEM->Cols(); c++)
                {
                    //if(!pcr::isMV(DEM->data[r][c])
                    {
                        //get water height
                        float vol = VolI[i]->data[r][c];
                        //rescale

                        //distribut over 9 sub-areas
                        if(r > -1 && r < DEM->Rows())
                        {
                            VolI[i-1]->data[r][c] += vol;
                        }
                        if(r-idist > -1 && r-idist < DEM->Rows())
                        {
                            VolI[i-1]->data[r-idist][c] += vol;
                        }
                        if(r+idist > -1 && r+idist < DEM->Rows())
                        {
                            VolI[i-1]->data[r+idist][c] += vol;
                        }
                        if(c+idist > -1 && c+idist < DEM->Cols())
                        {
                            VolI[i-1]->data[r][c+idist] += vol;
                        }
                        if(c-idist > -1 && c-idist < DEM->Cols())
                        {
                            VolI[i-1]->data[r][c-idist] += vol;
                        }
                        if(r-idist > -1 && r-idist < DEM->Rows() && c+idist > -1 && c+idist < DEM->Cols())
                        {
                            VolI[i-1]->data[r-idist][c + idist] += vol;
                        }
                        if(r+idist > -1 && r+idist < DEM->Rows() && c+idist > -1 && c+idist < DEM->Cols())
                        {
                            VolI[i-1]->data[r+idist][c + idist] += vol;
                        }
                        if(c-idist > -1 && c-idist < DEM->Cols() && r+idist > -1 && r+idist < DEM->Rows())
                        {
                            VolI[i-1]->data[r-idist][c-idist] += vol;
                        }
                        if(c-idist > -1 && c-idist < DEM->Cols() && r+idist > -1 && r+idist < DEM->Rows()) 
                        {
                            VolI[i-1]->data[r+idist][c-idist] += vol;
                        }
                    }
                }
            }
    }

    //delete all the maps we produced
    for(int i = 0; i < iters; i++)
    {
        delete DEM_MAXI[i];
    }
    for(int i = 0; i < iters; i++)
    {
        delete DEM_MAXL[i];
    }
    for(int i = 0; i < iters; i++)
    {
        delete DEM_MAXR[i];
    }
    for(int i = 0; i < iters; i++)
    {
       delete SR[i];
    }
    for(int i = 0; i < iters; i++)
    {
        delete SL[i];
    }
    for(int i = 0; i < iters; i++)
    {
        delete DEMAVI[i];
    }

    for(int i = 0; i < iters; i++)
    {
        delete WX1[i];
    }
    for(int i = 0; i < iters; i++)
    {
        delete WX2[i];
    }
    for(int i = 0; i < iters; i++)
    {
        delete WY1[i];
    }
    for(int i = 0; i < iters; i++)
    {
        delete WY2[i];
    }
    

    for(int i = 1; i < iters; i++)
    {
        delete VolI[i];
    }

    delete q_estI;
    delete P_estI;

    
    delete WidthI;


    return {A_estI,VolI[0]};

}



inline std::vector<cTMapf *> GetFlowCrossSectionalVolumeFromRiverDischargePreDistribute2(cTMapf * DEM,cTMapf * q,  cTMapf * N, cTMapf * CHANW, cTMapf * CHANH,cTMapf * SLOPE, int iters, float quality_p, cTMapf * FlowWL = nullptr, cTMapf * FlowWR = nullptr)
{
    cTMapf * test = q->GetCopy0();

    cTMapf * q_estI = q->GetCopy0();
    cTMapf * A_estI = q->GetCopy0();
    cTMapf * P_estI = q->GetCopy0();
    
    float dx = std::fabs(DEM->cellSizeX());

    std::vector<cTMapf *> DEM_MAXI;
    std::vector<cTMapf *> DEM_MAXL;
    std::vector<cTMapf *> DEM_MAXR;
    std::vector<cTMapf *> SR;
    std::vector<cTMapf *> SL;
    std::vector<cTMapf *> DEMAVI;
    std::vector<cTMapf *> VolI;
    std::vector<cTMapf *> WX1;
    std::vector<cTMapf *> WX2;
    std::vector<cTMapf *> WY1;
    std::vector<cTMapf *> WY2;

    int i_cascade_scale = 1;
    for(int i = 0; i < iters; i++)
    {
        if(i>0)
        {
            i_cascade_scale *= 3;
            DEM_MAXI.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            DEM_MAXL.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            DEM_MAXR.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            SR.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            SL.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            DEMAVI.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            VolI.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            WX1.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            WX2.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            WY1.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            WY2.push_back(DEM->GetCopy0Cascade(i_cascade_scale));
            FF::cout <<  i << " cascade : " << DEM_MAXI[i]->Cols() << " " << DEM_MAXI[i]->Rows()<<std::endl;
        }else{
            DEM_MAXI.push_back(DEM->GetCopy());
            DEM_MAXL.push_back(DEM->GetCopy());
            DEM_MAXR.push_back(DEM->GetCopy());
            SR.push_back(DEM->GetCopy0());
            SL.push_back(DEM->GetCopy0());
            DEMAVI.push_back(DEM->GetCopy());
            VolI.push_back(DEM->GetCopy0());
            WX1.push_back(DEM->GetCopy0());
            WX2.push_back(DEM->GetCopy0());
            WY1.push_back(DEM->GetCopy0());
            WY2.push_back(DEM->GetCopy0());
            FF::cout <<  i << " cascade : " << DEM_MAXI[i]->Cols() << " " << DEM_MAXI[i]->Rows()<<std::endl;
        }
    }

    FF::cout << "do dem prepare " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        //we look 2^i pixels away for slope
        float dxi = dx;

        //current with in logaritmic search
        dxi = dx * (1.0);
        if(i > 0)
        {
            dxi =dx * ((float)(std::pow(3,i)));
        }
        int idist = 1;
        if(i>0)
        {
            idist = std::pow(3,i);
        }

        cTMapf * mask = DEM_MAXI[i];
        cTMapf * mask_p = mask;
        if(i > 0)
        {
            mask_p = DEM_MAXI[i-1];
        }

        #pragma omp parallel for collapse(2)
        for (int r = 0; r < mask->Rows(); r++)
        {
            for(int c = 0; c < mask->Cols(); c++)
            {
                //first get the DEM values
                
                {
                    cTMapf * DEMAVIP;
                    cTMapf * DEMMAXIP;
                    if( i == 0)
                    {
                        DEMAVIP = DEM;
                    }else
                    {
                        DEMAVIP = DEMAVI[i-1];
                    }

                    if( i == 0)
                    {
                        DEMMAXIP = DEM;
                    }else
                    {
                        DEMMAXIP= DEM_MAXI[i-1];
                    }

                    int r_ip = r *3 + 1;
                    int c_ip = c *3 + 1;

                    if( i == 0)
                    {
                        r_ip = r;
                        c_ip = c;
                    }
                    float z  = 0.0;
                    pcr::setMV(z);
                    if(r_ip > -1 && c_ip > -1 && r_ip < mask_p->nrRows() && c_ip < mask_p->nrCols())
                    {
                        z = DEMAVIP->data[r_ip][c_ip];
                    }
                    float zm = z;
                    if(r_ip > -1 && c_ip > -1 && r_ip < mask_p->nrRows() && c_ip < mask_p->nrCols())
                    {
                        zm = DEMMAXIP->data[r_ip][c_ip];
                    }
                    if(!pcr::isMV(z))
                    {
                        

                        float zx1 = z;
                        float zx2 = z;
                        float zy1 = z;
                        float zy2 = z;

                        float zx1y1 = z;
                        float zx1y2 = z;
                        float zx2y1 = z;
                        float zx2y2 = z;

                        float zmx1 = z;
                        float zmx2 = z;
                        float zmy1 = z;
                        float zmy2 = z;

                        float zmx1y1 = z;
                        float zmx1y2 = z;
                        float zmx2y1 = z;
                        float zmx2y2 = z;

                        if(r_ip-1 > -1 && c_ip > -1 && r_ip-1 < mask_p->nrRows() && c_ip < mask_p->nrCols())
                        {
                            zy1 = DEMAVIP->data[r_ip-1][c_ip];
                            if(pcr::isMV(zy1))
                            {
                                zy1 = z;
                            }
                        }
                        if(r_ip+1 > -1 && c_ip > -1 && r_ip+1 < mask_p->nrRows() && c_ip < mask_p->nrCols())
                        {
                            zy2 = DEMAVIP->data[r_ip+1][c_ip];
                            if(pcr::isMV(zy2))
                            {
                                zy2 = z;
                            }
                        }
                        if(r_ip > -1 && c_ip-1 > -1 && r_ip < mask_p->nrRows() && c_ip-1 < mask_p->nrCols())
                        {
                            zx1 = DEMAVIP->data[r_ip][c_ip-1];
                            if(pcr::isMV(zx1))
                            {
                                zx1 = z;
                            }
                        }
                        if(r_ip > -1 && c_ip+1 > -1 && r_ip < mask_p->nrRows() && c_ip+1 < mask_p->nrCols())
                        {
                            zx2 = DEMAVIP->data[r_ip][c_ip+1];
                            if(pcr::isMV(zx2))
                            {
                                zx2 = z;
                            }
                        }

                        if(r_ip > 0 && c_ip > 0 && r_ip < mask_p->nrRows()-1 && c_ip < mask_p->nrCols()-1)
                        {
                            zx1y1 = DEMAVIP->data[r_ip-1][c_ip-1];
                            if(pcr::isMV(zx1y1))
                            {
                                zx1y1 = z;
                            }
                            zx1y2 = DEMAVIP->data[r_ip+1][c_ip-1];
                            if(pcr::isMV(zx1y2))
                            {
                                zx1y2 = z;
                            }
                            zx2y1 = DEMAVIP->data[r_ip-1][c_ip+1];
                            if(pcr::isMV(zx2y1))
                            {
                                zx2y1 = z;
                            }
                            zx2y2 = DEMAVIP->data[r_ip+1][c_ip+1];
                            if(pcr::isMV(zx2y2))
                            {
                                zx2y2 = z;
                            }
                        }

                        if(r_ip-1 > -1 && c_ip > -1 && r_ip-1 < mask_p->nrRows() && c_ip < mask_p->nrCols())
                        {
                            zmy1 = DEMMAXIP->data[r_ip-1][c_ip];
                            if(pcr::isMV(zmy1))
                            {
                                zmy1 = z;
                            }
                        }
                        if(r_ip+1 > -1 && c_ip > -1 && r_ip+1 < mask_p->nrRows() && c_ip < mask_p->nrCols())
                        {
                            zmy2 = DEMMAXIP->data[r_ip+1][c_ip];
                            if(pcr::isMV(zmy2))
                            {
                                zmy2 = z;
                            }
                        }
                        if(r_ip > -1 && c_ip-1 > -1 && r_ip < mask_p->nrRows() && c_ip-1 < mask_p->nrCols())
                        {
                            zmx1 = DEMMAXIP->data[r_ip][c_ip-1];
                            if(pcr::isMV(zmx1))
                            {
                                zmx1 = z;
                            }
                        }
                        if(r_ip > -1 && c_ip+1 > -1 && r_ip < mask_p->nrRows() && c_ip+1 < mask_p->nrCols())
                        {
                            zmx2 = DEMMAXIP->data[r_ip][c_ip+1];
                            if(pcr::isMV(zmx2))
                            {
                                zmx2 = z;
                            }
                        }

                        if(r_ip > 0 && c_ip > 0 && r_ip < mask_p->nrRows()-1 && c_ip < mask_p->nrCols()-1)
                        {
                            zmx1y1 = DEMMAXIP->data[r_ip-1][c_ip-1];
                            if(pcr::isMV(zmx1y1))
                            {
                                zmx1y1 = z;
                            }
                            zmx1y2 = DEMMAXIP->data[r_ip+1][c_ip-1];
                            if(pcr::isMV(zmx1y2))
                            {
                                zmx1y2 = z;
                            }
                            zmx2y1 = DEMMAXIP->data[r_ip-1][c_ip+1];
                            if(pcr::isMV(zmx2y1))
                            {
                                zmx2y1 = z;
                            }
                            zmx2y2 = DEMMAXIP->data[r_ip+1][c_ip+1];
                            if(pcr::isMV(zmx2y2))
                            {
                                zmx2y2 = z;
                            }
                        }


                        if(i !=0)
                        {
                            DEMAVI[i]->data[r][c] = (zx1 + zx2 + zy1 + zy2 + zx1y1 + zx1y2 + zx2y1 + zx2y2 + z)/9.0f;
                        }
                        if(i !=0)
                        {
                            DEM_MAXI[i]->data[r][c] = std::max(std::max(std::max(std::max(std::max(std::max(std::max(std::max(zm,zmx1),zmx2),zmy1),zmy2),zmx1y1),zmx1y2),zmx2y1),zmx2y2);
                        }

                        float zmax_x1 = std::max(zmx1,std::max(zmx1y1,zmx1y2));
                        float zmax_x2 = std::max(zmx2,std::max(zmx2y1,zmx2y2));
                        float zmax_y1 = std::max(zmy1,std::max(zmx1y1,zmx2y1));
                        float zmax_y2 = std::max(zmy2,std::max(zmx1y2,zmx2y2));

                        float sx1 =(zx1 -z)/dxi;
                        float sy1 =(zy1 -z)/dxi;
                        float sx2 =(zx2 -z)/dxi;
                        float sy2 =(zy2 -z)/dxi;

                        //get direction of steepest descent
                        float ux = (sx1<sx2?-1.0 * sx1:1.0*sx2);
                        float uy = (sy1<sy2?-1.0 * sy1:1.0*sy2);
                        float norm = std::sqrt(ux*ux+uy*uy);
                        if(norm < 1e-12)
                        {
                            ux = 1.0;
                            uy = 0.0;
                            norm = 1.0;
                        }
                        ux= ux/std::max(1e-12f,norm);
                        uy= uy/std::max(1e-12f,norm);

                        //get left and right slope values (multiply vectors by 2d rotation matrix of 90 and -90 degrees)
                        float ulx = -uy;
                        float uly = ux;
                        float urx = uy;
                        float ury = -ux;

                        float sl= ((ulx > 0.0?std::fabs(ulx)*sx2: std::fabs(ulx)*sx1)) + std::fabs((uly > 0.0?std::fabs(uly)*sy2: std::fabs(uly)*sy1));
                        float sr= ((urx > 0.0?std::fabs(urx)*sx2: std::fabs(urx)*sx1)) + std::fabs((ury > 0.0?std::fabs(ury)*sy2: std::fabs(ury)*sy1));
                        float zl = z + dxi * sl;
                        float zr = z + dxi * sr;
                        float zml = ((ulx > 0.0? ulx * zmax_x2 : -ulx * zmax_x1) + (uly > 0.0? uly * zmax_y2 : -uly * zmax_y1))/std::max(1e-10f,std::fabs(uly) + std::fabs(ulx));
                        float zmr = ((urx > 0.0? urx * zmax_x2 : -urx * zmax_x1) + (ury > 0.0? ury * zmax_y2 : -ury * zmax_y1))/std::max(1e-10f,std::fabs(ury) + std::fabs(urx));

                        SL[i]->data[r][c] = std::max(0.001f,sl);
                        SR[i]->data[r][c] = std::max(0.001f,sr);

                        DEM_MAXR[i]->data[r][c] = zmr;
                        DEM_MAXL[i]->data[r][c] = zml;

                        WX1[i]->data[r][c] = ulx/std::max(0.001f,std::fabs(std::fabs(ulx) + std::fabs(uly)));
                        WX2[i]->data[r][c] = urx/std::max(0.001f,std::fabs(std::fabs(urx) + std::fabs(ury)));
                        WY1[i]->data[r][c] = uly/std::max(0.001f,std::fabs(std::fabs(ulx) + std::fabs(uly)));
                        WY2[i]->data[r][c] = ury/std::max(0.001f,std::fabs(std::fabs(urx) + std::fabs(ury)));

                    }
                }
            }
        }
    }

    SLOPE = SLOPE->GetCopy0();
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++) {
        for(int c = 0; c < DEM->Cols(); c++) {

            if (!pcr::isMV(DEM->data[r][c])) {
                //use min(iters-1,3);
                int i_target =3;
                if (dx > 50) {
                    i_target = 2;
                }if (dx > 100) {
                    i_target = 1;
                }if (dx > 200) {
                    i_target = 0;
                }

                int i_left = std::min(iters-1,i_target);

                //get coordinates

                int r_il = std::max(0,std::min(DEMAVI[i_left]->Rows()-1,(int)(r/(i_left == 0? 1: std::pow(3,std::max(0,i_left))))));
                int c_il = std::max(0,std::min(DEMAVI[i_left]->Cols()-1,(int)(c/(i_left == 0? 1: std::pow(3,std::max(0,i_left))))));

                float z= DEMAVI[i_left]->data[r_il][c_il];
                if (pcr::isMV(z)) {
                    z = DEM->data[r][c];
                }

                float z_x1 = z;
                float z_x2 = z;
                float z_y1 = z;
                float z_y2 = z;

                if (r_il > 0) {
                    if (!pcr::isMV(DEMAVI[i_left]->data[r_il-1][c_il])) {
                        z_y1 = DEMAVI[i_left]->data[r_il-1][c_il];
                    }
                }
                if (c_il > 0) {
                    if (!pcr::isMV(DEMAVI[i_left]->data[r_il][c_il-1])) {
                        z_x1 = DEMAVI[i_left]->data[r_il][c_il-1];
                    }
                }
                if (r_il < DEMAVI[i_left]->Rows() -1) {
                    if (!pcr::isMV(DEMAVI[i_left]->data[r_il+1][c_il])) {
                        z_y2 = DEMAVI[i_left]->data[r_il+1][c_il];
                    }
                }
                if (c_il > DEMAVI[i_left]->Cols() - 1) {
                    if (!pcr::isMV(DEMAVI[i_left]->data[r_il][c_il+1])) {
                        z_x2 = DEMAVI[i_left]->data[r_il][c_il+1];
                    }
                }
                std::array<float, 4> values = {z_x1, z_x2, z_y1, z_y2};
                // 2. Sort the array in ascending order
                std::sort(values.begin(), values.end());
                float smallest = values[0];
                float second_smallest = values[1];
                float dxi = dx;
                if(i_left > 0)
                {
                    dxi =dx * ((float)(std::pow(3,i_left-1)));
                }else
                {
                    dxi = dx;
                }
                SLOPE->data[r][c] = std::max(1e-2f,std::fabs(smallest - second_smallest)/(2.0f*dxi));
            }

        }
    }
    FF::cout << "do actual distribution " << std::endl;
        
        #pragma omp parallel for collapse(2)
        for (int r = 0; r < DEM->Rows(); r++)
        {
            for(int c = 0; c < DEM->Cols(); c++)
            {
                if(!pcr::isMV(DEM->data[r][c]) && (!(std::fabs(q_estI->data[r][c] - q->data[r][c]) < 1e-4)))
                {
                    float slope = std::max(0.001f,SLOPE->data[r][c]);
                    float dem = DEM->data[r][c];
                    int i_left = 0;
                    int i_right = 0;
                    bool left_blocked = false;
                    bool right_blocked = false;
                    float dx_middle = dx;

                    float h_left = 0.0 ;
                    float h_right = 0.0;
                    if(r == 261 && c== 105)
                    {
                        //std::cout << "sus pixel " << q->data[r][c] << std::endl;
                    }

                    while(i_left < iters && i_right < iters)//(i_left < iters && i_right < iters)
                    {
                        if((std::fabs(q_estI->data[r][c] - q->data[r][c]) < 1e-4))
                        {
                            break;
                        }

                        //we look 3^i pixels away for slope
                        float dxi_left = dx;
                        if(i_left > 0)
                        {
                            dxi_left =dx * ((float)(std::pow(3,i_left)));
                        }
                        int idist_left = 1;
                        if(i_left>0)
                        {
                            idist_left = std::pow(3,i_left);
                        }
                        float dxi_right = dx;
                        if(i_right > 0)
                        {
                            dxi_right =dx * ((float)(std::pow(3,i_right)));
                        }
                        int idist_right = 1;
                        if(i_right>0)
                        {
                            idist_right = std::pow(3,i_right);
                        }

                        float dxi_left_p = dxi_left;
                        if(i_left > 0)
                        {
                            dxi_left_p =dx * ((float)(std::pow(3,i_left-1)));
                        }else
                        {
                            dxi_left_p = dx;
                        }
                        int idist_left_p = idist_left;
                        if(i_left>0)
                        {
                            idist_left_p = std::pow(3,i_left-1);
                        }else
                        {
                            idist_left_p = 0;
                        }
                        float dxi_right_p = dxi_right;
                        if(i_right > 0)
                        {
                            dxi_right_p =dx * ((float)(std::pow(3,i_right-1)));
                        }else{
                            dxi_right_p = dx;
                        }
                        int idist_right_p = idist_right;
                        if(i_right>0)
                        {
                            idist_right_p = std::pow(3,i_right-1);
                        }else{
                            idist_right_p = 0;
                        }

                        if(FlowWL != nullptr)
                        {
                            FlowWL->data[r][c] = std::max(dx,dxi_left);
                        }
                           if(FlowWR != nullptr)
                        {
                            FlowWR->data[r][c] = std::max(dx,dxi_right);
                        }                 


                        int i_left_p = std::max(0,i_left-1);
                        int i_right_p = std::max(0,i_right-1);
                        
                        int r_il = std::max(0,std::min(DEMAVI[i_left]->Rows()-1,(int)(r/(i_left == 0? 1: std::pow(3,std::max(0,i_left))))));
                        int c_il = std::max(0,std::min(DEMAVI[i_left]->Cols()-1,(int)(c/(i_left == 0? 1: std::pow(3,std::max(0,i_left))))));
                        int r_ir = std::max(0,std::min(DEMAVI[i_right]->Rows()-1,(int)(r/(i_right == 0? 1: std::pow(3,std::max(0,i_right))))));
                        int c_ir = std::max(0,std::min(DEMAVI[i_right]->Cols()-1,(int)(c/(i_right == 0? 1: std::pow(3,std::max(0,i_right))))));
                        int r_ilp = std::max(0,std::min(DEMAVI[i_left_p]->Rows()-1,(int)(r/(i_left < 2? 1: std::pow(3,std::max(0,i_left-1))))));
                        int c_ilp = std::max(0,std::min(DEMAVI[i_left_p]->Cols()-1,(int)(c/(i_left < 2? 1: std::pow(3,std::max(0,i_left-1))))));
                        int r_irp = std::max(0,std::min(DEMAVI[i_right_p]->Rows()-1,(int)(r/(i_right < 2? 1: std::pow(3,std::max(0,i_right-1))))));
                        int c_irp = std::max(0,std::min(DEMAVI[i_right_p]->Cols()-1,(int)(c/(i_right < 2? 1: std::pow(3,std::max(0,i_right-1))))));
                        
                        //get current left and right slope, and current left and right barrier height, if any
                        float sl= SL[i_left]->data[r_il][c_il];
                        float sr= SR[i_right]->data[r_ir][c_ir];
                        if(pcr::isMV(sl))
                        {
                            sl = 0.001f;
                        }
                        if(pcr::isMV(sr))
                        {
                            sr = 0.001f;
                        }

                        float zml = std::max(
                            DEM_MAXL[i_left]->data[std::max(0,std::min(WY1[i_left]->Rows()-1,r_il + (WY1[i_left]->data[r_il][c_il]>0.0?1:-1)))][c_il]
                            ,std::max(
                                DEM_MAXL[i_left]->data[r_il][std::max(0,std::min(WY1[i_left]->Cols()-1,c_il + (WX1[i_left]->data[r_il][c_il]>0.0?1:-1)))]
                                ,DEM_MAXL[i_left]->data[r_il][c_il])
                                );
                        float zmr = std::max(
                                    DEM_MAXR[i_right]->data[std::max(0,std::min(WY1[i_right]->Rows()-1,r_ir + (WY2[i_right]->data[r_ir][c_ir]>0.0?1:-1)))][c_ir]
                                ,std::max(
                                    DEM_MAXR[i_right]->data[r_ir][std::max(0,std::min(WX2[i_right]->Cols()-1,c_ir + (WX2[i_right]->data[r_ir][c_ir]>0.0?1:-1)))]
                                ,DEM_MAXR[i_right]->data[r_ir][c_ir])
                                );
                        if(pcr::isMV(zml))
                        {
                            zml = dem+ 0.001f;
                        }
                        if(pcr::isMV(zmr))
                        {
                            zmr = dem + 0.001f;
                        }

                        //how much does the barrier extend above our average piece-wise linear approximation?
                        float demav_left = i_left > 0?DEMAVI[i_left]->data[r_il][c_il]:dem;
                        float demav_right = i_right > 0? DEMAVI[i_right]->data[r_ir][c_ir] : dem;
                        float barrier_left = std::max(0.0f,zml - (dem + h_left));
                        float barrier_right = std::max(0.0f,zmr - (dem + h_right));

                        if(q->data[r][c] > 50)
                        {
                            //FF::cout << i_left << " " << i_right  << "  " <<  i_left_p << " " << i_right_p << " " <<  r << " " << c << " " << r_il << " " << c_il << " " << r_ir << " " << c_ir << " prev " << " " << r_ilp << " " << c_ilp << " " << r_irp << " " << c_irp << std::endl;
                            //FF::cout << "bh "  << q->data[r][c] << " " <<   i_left << " " << i_right <<  " dem: " << dem << " dev:" << demav_left << " " << demav_right << " b:" <<  barrier_left << " " << barrier_right <<  " z:" << zml << " " << zmr <<  "  wh:" << h_left << std::endl;
                        }

                        //check channel
                        if( i_left == 0 && i_right == 0 && CHANW != nullptr)
                        {
                            if(CHANW->data[r][c] > 1e-10 && CHANH->data[r][c] > 1e-10)
                            {
                                float QN = q->data[r][c];
                                float AN = 0.0;
                                float PN = 0.0;
                                SolveManningRectangular(CHANH->data[r][c],slope, CHANW->data[r][c], N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,FF_QH_ITERS);

                                q_estI->data[r][c] = std::min(q->data[r][c],QN);
                                A_estI->data[r][c] = AN;
                                P_estI->data[r][c] = PN;

                                //add volume here
                                VolI[0]->data[r][c] += AN/dx;
                                
                            }
                            
                            if((std::fabs(q_estI->data[r][c] - q->data[r][c]) < 1e-4))
                            {
                                break;
                            }
                        }


                        //we do a symmetrical trapezoid
                        sl = 0.5 * sl * 4.0/3.0;
                        sr = 0.5 * sr * 4.0/3.0;
                        float s = std::max(1e-3, std::fabs(0.5*(sl + sr)));

                        bool left_is_blocking = false;
                        bool right_is_blocking = false;

                        //if not 0 iteration, check if barriers contain?
                        if(i_left > 0 || i_right > 0)
                        {
                            //rectangular channel Q

                            if(barrier_left > 0.0)
                            {
                                float AN = A_estI->data[r][c] + 2.0 * barrier_left * (dxi_left_p + dxi_right_p);
                                float PN = P_estI->data[r][c] + 2.0 * barrier_left;

                                float QN = AN * std::pow(AN/PN,2.0/3.0)*std::sqrt(slope)/(N->data[r][c]*MANNING_MULT_CSCPROF);
                                if(QN > q->data[r][c])
                                {
                                    left_is_blocking = true;
                                }
                            }

                            if(barrier_right > 0.0)
                            {
                                float AN = A_estI->data[r][c] + 2.0 * barrier_right * (dxi_left_p + dxi_right_p);
                                float PN = P_estI->data[r][c] + 2.0 * barrier_right;

                                float QN = AN * std::pow(AN/PN,2.0/3.0)*std::sqrt(slope)/(N->data[r][c]*MANNING_MULT_CSCPROF);
                                if(QN > q->data[r][c])
                                {
                                    right_is_blocking = true;
                                }

                            }
                        }

                        //FF::cout << "do distribution" <<std::endl;

                        //if so, set area and break
                        if(left_is_blocking && right_is_blocking && (i_left > 0 || i_right > 0))
                        {

                            //solve rectangular channel
                            float QN = q->data[r][c];
                            float AN = A_estI->data[r][c];
                            
                            float AN_O = AN;
                            float PN = P_estI->data[r][c];
                            SolveManningRectangular(std::max(barrier_left, barrier_right),slope, dxi_left_p + dxi_right_p, N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,FF_QH_ITERS);
                            q_estI->data[r][c] = std::min(q->data[r][c],QN);
                            A_estI->data[r][c] = AN;
                            P_estI->data[r][c] = PN;

                            test->data[r][c] += AN;

                            //area to volume
                            float vol = (AN - AN_O)/(dxi_left_p + dxi_right_p);

                            if(q->data[r][c] > 50.0)
                            {
                                //FF::cout << "blockb " << barrier_left << " " << barrier_right << " " << vol << " " << dxi_left_p + dxi_right_p << " " << dx_middle << q_estI->data[r][c] << " " << q->data[r][c] << " " << QN << std::endl;
                            }

                            //FF::cout << "blocked " << std::endl;

                            if(r == 261 && c== 105)
                            {
                                //std::cout << "sus pixelbb " << AN << " " << test->data[r][c] << " " << q_estI->data[r][c] << " " << vol << std::endl;

                            }
                            //divide over left and right
                            //float w_l = dxi_left_p /std::max(0.0001f,(dxi_left_p + dxi_right_p));
                            //float w_r = dxi_right_p /std::max(0.0001f,(dxi_left_p + dxi_right_p));

                            //we look 3^i pixels away for slope
                            
                            //left x weight
                            //left x offset
                            
                            if((c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) > -1) && (c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)  < DEMAVI[i_left_p]->nrCols()))
                            {
                                VolI[i_left_p]->data[r_ilp][c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)] += (dx/dxi_left_p) *vol * std::fabs(WX1[i_left_p]->data[r_ilp][c_ilp]);
                            }
                            if((c_irp + (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) > -1) && (c_irp +  (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) < DEMAVI[i_right_p]->nrCols()))
                            {
                                VolI[i_right_p]->data[r_irp][c_irp + (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1)] += (dx/dxi_right_p) *vol * std::fabs(WX2[i_right_p]->data[r_irp][c_irp]);
                            }
                            if((r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) > -1) && (r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) < DEMAVI[i_left_p]->nrRows()))
                            {
                                VolI[i_left_p]->data[r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)][c_ilp] += (dx/dxi_left_p) *vol * std::fabs(WY1[i_left_p]->data[r_ilp][c_ilp]);
                            }
                            if((r_irp + (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) > -1) && (r_irp +  (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) < DEMAVI[i_right_p]->nrRows()))
                            {
                                VolI[i_right_p]->data[r_irp + (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1)][c_irp] += (dx/dxi_right_p) *vol * std::fabs(WY2[i_right_p]->data[r_irp][c_irp]);
                            }

                            break;
                        }else if(std::min(barrier_left, barrier_right) > 1e-12f  && (i_left > 0 || i_right > 0))
                        {
                            //solve rectangular channel
                            float QN = q->data[r][c];
                            float AN = A_estI->data[r][c];

                            float AN_O = AN;
                            float PN = P_estI->data[r][c];
                            SolveManningRectangular(std::min(barrier_left, barrier_right),slope, dxi_left_p + dxi_right_p, N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,FF_QH_ITERS);
                            q_estI->data[r][c] = std::min(q->data[r][c],QN);
                            A_estI->data[r][c] = AN;
                            P_estI->data[r][c] = PN;

                            test->data[r][c] += AN;

                            //area to volume
                            float vol = (AN - AN_O)/(dxi_left_p + dxi_right_p);
                            if(r == 261 && c== 105)
                            {
                                //std::cout << "sus pixelbnew " << AN << " " << test->data[r][c] << " " << q_estI->data[r][c] << " " << vol << std::endl;

                            }
                            if(q->data[r][c] > 50.0)
                            {
                                //FF::cout << "blockb " << barrier_left << " " << barrier_right << " " << vol << " " << dxi_left_p + dxi_right_p << " " << dx_middle << q_estI->data[r][c] << " " << q->data[r][c] << " " << QN << std::endl;
                            }

                            h_left += std::min(barrier_left, barrier_right);
                            h_right += std::min(barrier_left, barrier_right);
                            //FF::cout << "blocked " << std::endl;

                            //divide over left and right
                            //float w_l = dxi_left_p /std::max(0.0001f,(dxi_left_p + dxi_right_p));
                            //float w_r = dxi_right_p /std::max(0.0001f,(dxi_left_p + dxi_right_p));

                            //we look 3^i pixels away for slope

                            //left x weight
                            //left x offset

                            if((c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) > -1) && (c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)  < DEMAVI[i_left_p]->nrCols()))
                            {
                                VolI[i_left_p]->data[r_ilp][c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)] += (dx/dxi_left_p) *vol * std::fabs(WX1[i_left_p]->data[r_ilp][c_ilp]);
                            }
                            if((c_irp + (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) > -1) && (c_irp +  (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) < DEMAVI[i_right_p]->nrCols()))
                            {
                                VolI[i_right_p]->data[r_irp][c_irp + (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1)] += (dx/dxi_right_p) *vol * std::fabs(WX2[i_right_p]->data[r_irp][c_irp]);
                            }
                            if((r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) > -1) && (r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) < DEMAVI[i_left_p]->nrRows()))
                            {
                                VolI[i_left_p]->data[r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)][c_ilp] += (dx/dxi_left_p) *vol * std::fabs(WY1[i_left_p]->data[r_ilp][c_ilp]);
                            }
                            if((r_irp + (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) > -1) && (r_irp +  (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) < DEMAVI[i_right_p]->nrRows()))
                            {
                                VolI[i_right_p]->data[r_irp + (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1)][c_irp] += (dx/dxi_right_p) *vol * std::fabs(WY2[i_right_p]->data[r_irp][c_irp]);
                            }
                        }

                        //pick a side, do one-sided
                        //prefer non-blocking
                        //then prefer lower i
                        bool do_left = true;
                        if(left_is_blocking || (barrier_left > barrier_right))
                        {
                            do_left = false;
                        }

                        //run one-sided trapezoid 
                        if(i_left == 0 && i_right == 0)
                        {

                            float QN = q->data[r][c];
                            float AN = A_estI->data[r][c];
                            float AO = AN;
                            float PN = P_estI->data[r][c];
                            SolveManningTrapezoid(s,slope, dxi_left,N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,FF_QH_ITERS);
                            q_estI->data[r][c] = std::min(q->data[r][c],QN);
                            A_estI->data[r][c] = AN;
                            P_estI->data[r][c] = PN;

                            test->data[r][c] += AN;
                            i_left = i_left + 1;
                            i_right = i_right + 1;
                            h_left += s * dxi_left;
                            h_right += s * dxi_right;
                            dx_middle += 2.0 * dx;                           
                        
                            float vol = (AN-AO)/(dxi_left + dxi_right);

                            if(r == 261 && c== 105)
                            {
                                //std::cout << "sus pixelfirst " << AN << " " <<  test->data[r][c]<< "  " << vol<< std::endl;

                            }


                            if(q->data[r][c] >25.0)
                            {
                                //FF::cout << "middle " << r << " " << c << " " << q->data[r][c] << " " << slope << " " << s << " "<< dxi_left<< " " << N->data[r][c]<< " " << vol << " " << AN << " " << AO << " " << dx <<  std::endl;
                                //FF::cout << r_il << " " << c_il << " " << r_ir << " " << c_ir << std::endl;
                            }
                
                            if((c + (WX1[0]->data[r][c]>0.0?1:-1) > -1) && (c + (WX1[0]->data[r][c]>0.0?1:-1)  < DEMAVI[0]->nrCols()))
                            {
                                VolI[0]->data[r][c] += (dx/dxi_left) *vol * 0.5*  std::fabs(WX1[0]->data[r][c]);
                            }
                            if((c + (WX2[0]->data[r][c]>0.0?1:-1) > -1) && (c +  (WX2[0]->data[r][c]>0.0?1:-1) < DEMAVI[0]->nrCols()))
                            {
                                VolI[0]->data[r][c ] += (dx/dxi_right) *vol *0.5*  std::fabs(WX2[0]->data[r][c]);
                            }
                            if((r + (WY1[0]->data[r][c]>0.0?1:-1) > -1) && (r + (WY1[0]->data[r][c]>0.0?1:-1) < DEMAVI[0]->nrRows()))
                            {
                                VolI[0]->data[r ][c] += (dx/dxi_left) *vol *0.5*  std::fabs(WY1[0]->data[r][c]);
                            }
                            if((r + (WY2[0]->data[r][c]>0.0?1:-1) > -1) && (r +  (WY2[0]->data[r][c]>0.0?1:-1) < DEMAVI[0]->nrRows()))
                            {
                                VolI[0]->data[r ][c] += (dx/dxi_right) *vol *0.5*  std::fabs(WY2[0]->data[r][c]);
                            }


                        }else
                        {
                           
                           if(do_left)
                           {
                                float QN = q->data[r][c];
                                float AN = A_estI->data[r][c];
                                float AO = AN;
                                float PN = P_estI->data[r][c];
                                SolveManningTrapezoidOS(sl,slope, dxi_left, dx_middle,  N->data[r][c]*MANNING_MULT_CSCPROF, QN,AN,PN,FF_QH_ITERS,1.0);
                                q_estI->data[r][c] = std::min(q->data[r][c],QN);
                                A_estI->data[r][c] = AN;
                                P_estI->data[r][c] = PN;

                                test->data[r][c] += AN;

                                h_left += s * dxi_left;
                                h_right += s * dxi_left;
                        
                                float vol = (AN-AO)/(dxi_right_p + dxi_left);
                                if(r == 261 && c== 105)
                                {
                                    //std::cout << "sus pixell " << AN <<  " " << test->data[r][c]<<" " << q_estI->data[r][c] <<  " " << vol << std::endl;

                                }
                                dx_middle += dxi_left;   

                                if(q->data[r][c] >25.0)
                                {
                                    //FF::cout << "middle " << r << " " << c << " " << q->data[r][c] << " " << slope << " " << s << " "<< dxi_left<< " " << N->data[r][c]<< " " << vol << " " << AN << " " << AO << " " << dx <<  std::endl;
                                    //FF::cout << "left "<< i_left << " " << i_right << " " << i_left_p << " " << i_right_p <<  std::endl;
                                    //FF::cout << r_il << " " << c_il << " " << r_ir << " " << c_ir << " " << r_ilp << " " << c_ilp << std::endl;
                                }

                                if((c_il + (WX1[i_left]->data[r_il][c_il]>0.0?1:-1) > -1) && (c_il + (WX1[i_left]->data[r_il][c_il]>0.0?1:-1)  < DEMAVI[i_left]->nrCols()))
                                {
                                    VolI[i_left]->data[r_il][c_il + (WX1[i_left]->data[r_il][c_il]>0.0?1:-1)] += (dx/dxi_left) *vol * std::fabs(WX1[i_left]->data[r_il][c_il]);
                                }
                                if((c_irp + (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) > -1) && (c_irp +  (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) < DEMAVI[i_right_p]->nrCols()))
                                {
                                    VolI[i_right_p]->data[r_irp][c_irp + (WX2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1)] += (dx/dxi_right_p) *vol * std::fabs(WX2[i_right_p]->data[r_irp][c_irp]);
                                }
                                if((r_il + (WY1[i_left]->data[r_il][c_il]>0.0?1:-1) > -1) && (r_il + (WY1[i_left]->data[r_il][c_il]>0.0?1:-1) < DEMAVI[i_left]->nrRows()))
                                {
                                    VolI[i_left]->data[r_il + (WY1[i_left]->data[r_il][c_il]>0.0?1:-1)][c_il] += (dx/dxi_left) *vol * std::fabs(WY1[i_left]->data[r_il][c_il]);
                                }
                                if((r_irp + (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) > -1) && (r_irp +  (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1) < DEMAVI[i_right_p]->nrRows()))
                                {
                                    VolI[i_right_p]->data[r_irp + (WY2[i_right_p]->data[r_irp][c_irp]>0.0?1:-1)][c_irp] += (dx/dxi_right_p) *vol * std::fabs(WY2[i_right_p]->data[r_irp][c_irp]);
                                }

                                i_left = i_left + 1;    

                           }else{
                                float QN = q->data[r][c];
                                float AN = A_estI->data[r][c];
                                float AO = AN;
                                float PN = P_estI->data[r][c];
                                SolveManningTrapezoidOS(sr,slope, dxi_right, dx_middle,  (N->data[r][c]*MANNING_MULT_CSCPROF), QN,AN,PN,FF_QH_ITERS,1.0);
                                q_estI->data[r][c] = std::min(q->data[r][c],QN);
                                A_estI->data[r][c] = AN;
                                P_estI->data[r][c] = PN;

                                test->data[r][c] += AN;
                                h_left += s * dxi_right;    
                                h_right += s * dxi_right;                   
                        
                                float vol = (AN-AO)/(dxi_left_p+dxi_right);

                                if(r == 261 && c== 105)
                                {
                                    //std::cout << "sus pixelr " << AN <<  " " << test->data[r][c]<<" " << q_estI->data[r][c] <<  " " << vol << std::endl;

                                }
                                dx_middle += dxi_right; 

                                if(q->data[r][c] >25.0)
                                {
                                    //FF::cout << "middle " << r << " " << c << " " << q->data[r][c] << " " << slope << " " << s << " "<< dxi_left<< " " << N->data[r][c]<< " " << vol << " " << AN << " " << AO << " " << dx <<  std::endl;
                                    //FF::cout << "right "<< i_left << " " << i_right << " " << i_left_p << " " << i_right_p <<  std::endl;
                                    //FF::cout << r_il << " " << c_il << " " << r_ir << " " << c_ir << " " << r_ilp << " " << c_ilp << std::endl;
                                }  

                                if((c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) > -1) && (c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)  < DEMAVI[i_left_p]->nrCols()))
                                {
                                    VolI[i_left_p]->data[r_ilp][c_ilp + (WX1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)] += (dx/dxi_left_p) *vol * std::fabs(WX1[i_left_p]->data[r_ilp][c_ilp]);
                                }
                                if((c_ir + (WX2[i_right]->data[r_ir][c_ir]>0.0?1:-1) > -1) && (c_ir +  (WX2[i_right]->data[r_ir][c_ir]>0.0?1:-1) < DEMAVI[i_right]->nrCols()))
                                {
                                    VolI[i_right]->data[r_ir][c_ir + (WX2[i_right]->data[r_ir][c_ir]>0.0?1:-1)] += (dx/dxi_right) *vol * std::fabs(WX2[i_right]->data[r_ir][c_ir]);
                                }
                                 if((r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) > -1) && (r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1) < DEMAVI[i_left_p]->nrRows()))
                                {
                                    VolI[i_left_p]->data[r_ilp + (WY1[i_left_p]->data[r_ilp][c_ilp]>0.0?1:-1)][c_ilp] += (dx/dxi_left_p) *vol * std::fabs(WY1[i_left_p]->data[r_ilp][c_ilp]);
                                }
                                if((r_ir + (WY2[i_right]->data[r_ir][c_ir]>0.0?1:-1) > -1) && (r_ir +  (WY2[i_right]->data[r_ir][c_ir]>0.0?1:-1) < DEMAVI[i_right]->nrRows()))
                                {
                                    VolI[i_right]->data[r_ir + (WY2[i_right]->data[r_ir][c_ir]>0.0?1:-1)][c_ir] += (dx/dxi_right) *vol * std::fabs(WY2[i_right]->data[r_ir][c_ir]);
                                }

                                i_right = i_right + 1;     
                           }
                        }
                    }
                }
            }
        }

    //writeRaster(*test,"test.tif");

    FF::cout << "get average flow heights " << std::endl;
    /*for (int r = 0; r < DEMAVI[0]->Rows(); r++)
    {
        for(int c = 0; c < DEMAVI[0]->Cols(); c++)
        {
            if(( r/36 < DEMAVI[4]->Rows()-1) && (c/36 < DEMAVI[4]->Cols() -1))
            {
                VolI[0]->data[r][c] = DEM_MAXI[4]->data[r/36][c/36];
            }//else{
             //   VolI[0]->data[r][c] = 0.0;
            //}
        }
    }*/


    FF::cout << "do dem prepare " << std::endl;
    for(int i = 1; i < iters; i++)
    {
        
        FF::cout << "add wh av to dem av " << i << std::endl;

        cTMapf * mask = DEM_MAXI[i];
        cTMapf * mask_p = mask;
        if(i > 0)
        {
            mask_p = DEM_MAXI[i-1];
        }

        #pragma omp parallel for collapse(2)
        for (int r = 0; r < mask->Rows(); r++)
        {
            for(int c = 0; c < mask->Cols(); c++)
            {
                //first get the DEM values
                {
                    cTMapf * DEMAVIP;
                    if( i == 1 || i == 0)
                    {
                        DEMAVIP = VolI[0];
                    }else
                    {
                        DEMAVIP = WX1[i-1];
                    }

                    int r_ip = r *3 + 1;
                    int c_ip = c *3 + 1;

                    if( i == 0)
                    {
                        r_ip = r;
                        c_ip = c;
                    }
                    float z  = 0.0;
                    pcr::setMV(z);
                    if(r_ip > -1 && c_ip > -1 && r_ip < mask_p->nrRows() && c_ip < mask_p->nrCols())
                    {
                        z = DEMAVIP->data[r_ip][c_ip];
                    }
                    
                    if(!pcr::isMV(z))
                    {
                        float zx1 = z;
                        float zx2 = z;
                        float zy1 = z;
                        float zy2 = z;

                        float zx1y1 = z;
                        float zx1y2 = z;
                        float zx2y1 = z;
                        float zx2y2 = z;

                        if(r_ip-1 > -1 && c_ip > -1 && r_ip-1 < mask_p->nrRows() && c_ip < mask_p->nrCols())
                        {
                            zy1 = DEMAVIP->data[r_ip-1][c_ip];
                            if(pcr::isMV(zy1))
                            {
                                zy1 = z;
                            }
                        }
                        if(r_ip+1 > -1 && c_ip > -1 && r_ip+1 < mask_p->nrRows() && c_ip < mask_p->nrCols())
                        {
                            zy2 = DEMAVIP->data[r_ip+1][c_ip];
                            if(pcr::isMV(zy2))
                            {
                                zy2 = z;
                            }
                        }
                        if(r_ip > -1 && c_ip-1 > -1 && r_ip < mask_p->nrRows() && c_ip-1 < mask_p->nrCols())
                        {
                            zx1 = DEMAVIP->data[r_ip][c_ip-1];
                            if(pcr::isMV(zx1))
                            {
                                zx1 = z;
                            }
                        }
                        if(r_ip > -1 && c_ip+1 > -1 && r_ip < mask_p->nrRows() && c_ip+1 < mask_p->nrCols())
                        {
                            zx2 = DEMAVIP->data[r_ip][c_ip+1];
                            if(pcr::isMV(zx2))
                            {
                                zx2 = z;
                            }
                        }

                        if(r_ip > 0 && c_ip > 0 && r_ip < mask_p->nrRows()-1 && c_ip < mask_p->nrCols()-1)
                        {
                            zx1y1 = DEMAVIP->data[r_ip-1][c_ip-1];
                            if(pcr::isMV(zx1y1))
                            {
                                zx1y1 = z;
                            }
                            zx1y2 = DEMAVIP->data[r_ip+1][c_ip-1];
                            if(pcr::isMV(zx1y2))
                            {
                                zx1y2 = z;
                            }
                            zx2y1 = DEMAVIP->data[r_ip-1][c_ip+1];
                            if(pcr::isMV(zx2y1))
                            {
                                zx2y1 = z;
                            }
                            zx2y2 = DEMAVIP->data[r_ip+1][c_ip+1];
                            if(pcr::isMV(zx2y2))
                            {
                                zx2y2 = z;
                            }
                        }

                        if(i !=0)
                        {
                            WX1[i]->data[r][c] = ((zx1 + zx2 + zy1 + zy2 + zx1y1 + zx1y2 + zx2y1 + zx2y2 + z)/9.0f);
                        }
                    }
                }
            }
        }

        if(i > 0)
        {
            #pragma omp parallel for collapse(2)
            for (int r = 0; r < mask->Rows(); r++)
            {
                for(int c = 0; c < mask->Cols(); c++)
                {
                    if(!pcr::isMV(DEMAVI[i]->data[r][c]))
                    {
                        DEMAVI[i]->data[r][c] = DEMAVI[i]->data[r][c] + WX1[i]->data[r][c];
            
                    }

                }
            }
        }

        if(i > 0)
        {
            #pragma omp parallel for collapse(2)
            for (int r = 0; r < mask->Rows(); r++)
            {
                for(int c = 0; c < mask->Cols(); c++)
                {
                    WX1[i]->data[r][c] = WX1[i]->data[r][c] + VolI[i]->data[r][c];

                }
            }
        }
                
    }

    /*for (int r0 = 0; r0 < DEMAVI[0]->Rows(); r0++)
    {
        for(int c0 = 0; c0 < DEMAVI[0]->Cols(); c0++)
        {
            if(r0/9 < DEM_MAXI[2]->Rows() && c0/9 < DEM_MAXI[2]->Cols())
            {

                VolI[0]->data[r0][c0] = DEM_MAXI[2]->data[r0/9][c0/9];
            }
        }
    }*/


    //number of iterations to use for each layer of the cascading structure

    
    std::vector<int> iters_p = {(int)(quality_p * 15.0),(int)(quality_p *30.0),(int)(quality_p *30),(int)(quality_p*30),(int)(quality_p * 15),(int)(quality_p*15),(int)(quality_p * 5),6,6,6,6};
    //now iterate backwards to re-distribute the water 
    for(int i = iters-1; i > 0; i--)
    {
        if(i < iters -1)
        {
            //now do a small diffusive max wave
            FF::cout << "diffusive max wave cascade: " << i <<  " " << iters_p[i] << std::endl;
            FF::cout << DEMAVI[i]->Rows() << " " << DEMAVI[i]->Cols() << " " << DEM_MAXI[i]->Rows() << " " << DEM_MAXI[i]->Cols() << " " << VolI[i]->Rows() << " " << VolI[i]->Cols() << std::endl;
            cTMapf * res = DiffusiveMaxWaveCascade(DEMAVI[i],DEM_MAXI[i],VolI[i],iters_p[i],0.25);
            FF::cout << "replace data " << i << std::endl;

            //replace the local values
            for (int r0 = 0; r0 < DEMAVI[i]->Rows(); r0++)
            {
                for(int c0 = 0; c0 < DEMAVI[i]->Cols(); c0++)
                {
                    VolI[i]->data[r0][c0] = res->data[r0][c0];
                }
            }
            FF::cout << "done: " << i << std::endl;

            delete res;
        }

        int idist = 1;
        //#pragma omp parallel for collapse(2)
            for (int r0 = 0; r0 < DEMAVI[i]->Rows(); r0++)
            {
                for(int c0 = 0; c0 < DEMAVI[i]->Cols(); c0++)
                {
                    //if(!pcr::isMV(DEM->data[r][c])
                    {
                        //get water height
                        float vol = VolI[i]->data[r0][c0];
                        //rescale

                        int r = (r0 * 3) + 1;
                        int c = (c0 * 3) + 1;

                        //distribut over 9 sub-areas
                        if((r > -1) && (r < DEMAVI[i-1]->Rows()) && (c > -1) && (c < DEMAVI[i-1]->Cols()))
                        {
                            VolI[i-1]->data[r][c] += vol;
                        }
                        if((r-1 > -1) && (r-1 < DEMAVI[i-1]->Rows()) && (c > -1) && (c < DEMAVI[i-1]->Cols()))
                        {
                            VolI[i-1]->data[r-1][c] += vol;
                        }
                        if((r+1 > -1) && (r+1 < DEMAVI[i-1]->Rows()) && (c > -1) && (c < DEMAVI[i-1]->Cols()))
                        {
                            VolI[i-1]->data[r+1][c] += vol;
                        }
                        if((c+1 > -1) && (c+1 < DEMAVI[i-1]->Cols()) && (r > -1) && (r < DEMAVI[i-1]->Rows())) 
                        {
                            VolI[i-1]->data[r][c+1] += vol;
                        }
                        if((c-1 > -1) && (c-1 < DEMAVI[i-1]->Cols()) && (r > -1) && (r < DEMAVI[i-1]->Rows()))
                        {
                            VolI[i-1]->data[r][c-1] += vol;
                        }
                        if((r-1 > -1) && (r-1 < DEMAVI[i-1]->Rows()) && (c+1 > -1) && (c+1 < DEMAVI[i-1]->Cols()))
                        {
                            VolI[i-1]->data[r-1][c + 1] += vol;
                        }
                        if((r+1 > -1) && (r+1 < DEMAVI[i-1]->Rows()) && (c+1 > -1) && (c+1 < DEMAVI[i-1]->Cols()))
                        {
                            VolI[i-1]->data[r+1][c + 1] += vol;
                        }
                        if((c-1 > -1) && (c-1 < DEMAVI[i-1]->Cols()) && (r+1 > -1) && (r+1 < DEMAVI[i-1]->Rows()))
                        {
                            VolI[i-1]->data[r-1][c-1] += vol;
                        }
                        if((c-1 > -1) && (c-1 < DEMAVI[i-1]->Cols()) && (r+1 > -1) && (r+1 < DEMAVI[i-1]->Rows())) 
                        {
                            VolI[i-1]->data[r+1][c-1] += vol;
                        }
                    }
                }
            }
    }

    FF::cout << "delete layers " << std::endl;

    //delete all the maps we produced
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete DEM_MAXI[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete DEM_MAXL[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete DEM_MAXR[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
       delete SR[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete SL[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete DEMAVI[i];
    }
    FF::cout << "del1 " << std::endl;
    
    for(int i = 0; i < iters; i++)
    {
        FF::cout << " i " << i << WX1[i] << std::endl;
        delete WX1[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete WX2[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete WY1[i];
    }
    FF::cout << "del1 " << std::endl;
    for(int i = 0; i < iters; i++)
    {
        delete WY2[i];
    }
    FF::cout << "del1 " << std::endl;


    delete q_estI;
    delete P_estI;
    delete A_estI;
    delete SLOPE;
    //we replaced slope with a temporary copy, so we own that memory and should delete it.

    return VolI;

}
