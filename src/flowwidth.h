#ifndef OSS_FLOWWIDTH_H
#define OSS_FLOWWIDTH_H

#include "map.h"

inline std::vector<cTMapf *> GetFlowWidthDirectionFromRiverVolume(cTMapf * DEM,cTMapf * h, int iters)
{

    cTMapf * DEMAVI = DEM->GetCopy();
    cTMapf * DEMAVI2 = DEM->GetCopy();
    cTMapf * WidthI = DEM->GetCopy0();
    cTMapf * vol_leftI = h->GetCopy();

    cTMapf * DIRX = DEM->GetCopy0();
    cTMapf * DIRY = DEM->GetCopy0();
    cTMapf * DIRW = DEM->GetCopy0();

    float dx = std::fabs(DEM->cellSizeX());

    float dxr = 1.0/dx;
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++)
    {
        for(int c = 0; c < DEM->Cols(); c++)
        {
            if(!pcr::isMV(DEM->data[r][c]) && vol_leftI->data[r][c] > 1e-6)
            {
                vol_leftI->data[r][c] = h->data[r][c]*dx;
            }

        }
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
                if(!pcr::isMV(DEM->data[r][c]) && vol_leftI->data[r][c] > 1e-4 * dxi)
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

                    if( i == 0)
                    {

                        sl = sl * 4.0/3.0;
                        sr = sr * 4.0/3.0;
                        float s = std::max(1e-6,std::fabs(0.5*(sl + sr)));

                        //get width
                        float width = 2.0*std::sqrt( vol_leftI->data[r][c])/std::sqrt(s);
                        WidthI->data[r][c] = std::min(dxi * 3.0f,width);

                        float volh = WidthI->data[r][c] *  0.5*WidthI->data[r][c] * 0.5 * s;

                        DIRX->data[r][c] += volh * ux;
                        DIRY->data[r][c] += volh * uy;
                        DIRW->data[r][c] += volh;

                        vol_leftI->data[r][c] = std::max(0.0f,vol_leftI->data[r][c] - volh);
                    }else
                    {
                        float s = std::max(1e-6,std::fabs(0.5*(sl + sr)));

                        //get width
                        float width = std::sqrt(dxi*dxi * s + 4.0 * vol_leftI->data[r][c])/std::sqrt(s);
                        float widtho = WidthI->data[r][c];

                        WidthI->data[r][c] = std::min(dxi * 3.0f,width);
                        float volh = (WidthI->data[r][c]-dxi) *s *0.5* (WidthI->data[r][c]-dxi) *0.5 + dxi * (s*0.5*((WidthI->data[r][c]-dxi) ));

                        DIRX->data[r][c] += volh * ux;
                        DIRY->data[r][c] += volh * uy;
                        DIRW->data[r][c] += volh;

                        vol_leftI->data[r][c] = std::max(0.0f,vol_leftI->data[r][c] - volh);

                    }

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

    #pragma omp parallel for collapse(2)
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
    }


    delete DIRW;

    delete DEMAVI;
    delete DEMAVI2;
    delete vol_leftI;


    return {WidthI,DIRX,DIRY};

}

inline cTMapf * GetFlowWidthFromRiverVolume(cTMapf * DEM,cTMapf * h, int iters)
{

    cTMapf * DEMAVI = DEM->GetCopy();
    cTMapf * DEMAVI2 = DEM->GetCopy();
    cTMapf * WidthI = DEM->GetCopy0();
    cTMapf * vol_leftI = h->GetCopy();

    float dx = std::fabs(DEM->cellSizeX());

    float dxr = 1.0/dx;
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++)
    {
        for(int c = 0; c < DEM->Cols(); c++)
        {
            if(!pcr::isMV(DEM->data[r][c]) && vol_leftI->data[r][c] > 1e-6)
            {
                vol_leftI->data[r][c] = h->data[r][c]*dx;
            }

        }
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
                if(!pcr::isMV(DEM->data[r][c]) && vol_leftI->data[r][c] > 1e-4 * dxi)
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

                    if( i == 0)
                    {

                        sl = sl * 4.0/3.0;
                        sr = sr * 4.0/3.0;
                        float s = std::max(1e-6,std::fabs(0.5*(sl + sr)));

                        //get width
                        float width = 2.0*std::sqrt( vol_leftI->data[r][c])/std::sqrt(s);
                        WidthI->data[r][c] = std::min(dxi * 3.0f,width);

                        float volh = WidthI->data[r][c] *  0.5*WidthI->data[r][c] * 0.5 * s;

                        vol_leftI->data[r][c] = std::max(0.0f,vol_leftI->data[r][c] - volh);
                    }else
                    {
                        float s = std::max(1e-6,std::fabs(0.5*(sl + sr)));

                        //get width
                        float width = std::sqrt(dxi*dxi * s + 4.0 * vol_leftI->data[r][c])/std::sqrt(s);
                        float widtho = WidthI->data[r][c];

                        WidthI->data[r][c] = std::min(dxi * 3.0f,width);
                        float volh = (WidthI->data[r][c]-dxi) *s *0.5* (WidthI->data[r][c]-dxi) *0.5 + dxi * (s*0.5*((WidthI->data[r][c]-dxi) ));
                        vol_leftI->data[r][c] = std::max(0.0f,vol_leftI->data[r][c] - volh);

                    }

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


    delete DEMAVI;
    delete DEMAVI2;
    delete vol_leftI;


    return WidthI;


}


inline void GetAvSoilPropsFromRiverVolume(cTMapf * DEM,cTMapf * soilp, std::vector<cTMapf * > &output_soilprop, std::vector<cTMapf *> &output_area, int iters, bool do_area =true)
{

    cTMapf * DEMAVI = DEM->GetCopy();
    cTMapf * DEMAVI2 = DEM->GetCopy();

    cTMapf * SOILPAVI = soilp->GetCopy();
    cTMapf * SOILPAVI2 = soilp->GetCopy();

    cTMapf * WidthI = DEM->GetCopy0();

    float dx = std::fabs(DEM->cellSizeX());

    float dxr = 1.0/dx;

    //put initial maps

    for(int i = 0; i < iters; i++)
    {
        output_soilprop.push_back(DEM->GetCopy0());
        if(do_area)
        {
            output_area.push_back(DEM->GetCopy0());

        }
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


                    float sp= soilp->data[r][c];
                    float spx1;
                    float spx2;
                    float spy1;
                    float spy2;

                    if(i == 0)
                    {
                        if(c > 0)
                        {
                            zx1 = pcr::isMV(DEM->data[r][c-1])?dem:DEM->data[r][c-1];
                            spx1 = pcr::isMV(soilp->data[r][c-1])?sp:soilp->data[r][c-1];
                        }else {
                            zx1 = dem;
                            spx1 = sp;
                        }
                        if(c < DEM->nrCols()-1)
                        {
                            zx2 = pcr::isMV(DEM->data[r][c+1])?dem:DEM->data[r][c+1];
                            spx2 = pcr::isMV(soilp->data[r][c+1])?sp:soilp->data[r][c+1];
                        }else {
                            zx2 = dem;
                            spx2 = sp;
                        }
                        if(r > 0)
                        {
                            zy1 = pcr::isMV(DEM->data[r-1][c])?dem:DEM->data[r-1][c];
                            spy1 = pcr::isMV(soilp->data[r-1][c])?sp:soilp->data[r-1][c];
                        }else {
                            zy1 = dem;
                            spy1 = sp;
                        }
                        if(r < DEM->nrRows()-1)
                        {
                            zy2 = pcr::isMV(DEM->data[r+1][c])?dem:DEM->data[r+1][c];
                            spy2 = pcr::isMV(soilp->data[r+1][c])?sp:soilp->data[r+1][c];
                        }else {
                            zy2 = dem;
                            spy2= sp;

                        }
                    }else
                    {
                        dem = DEMAVI->data[r][c];

                        if(c > idist)
                        {
                            zx1 = pcr::isMV(DEMAVI->data[r][c-idist])?dem:DEMAVI->data[r][c-idist];
                            spx1 = pcr::isMV(SOILPAVI->data[r][c-idist])?sp:SOILPAVI->data[r][c-idist];
                        }else {
                            zx1 = dem;
                            spx1 = sp;
                        }
                        if(c < DEM->nrCols()-1-idist)
                        {
                            zx2 = pcr::isMV(DEMAVI->data[r][c+idist])?dem:DEMAVI->data[r][c+idist];
                            spx2 = pcr::isMV(SOILPAVI->data[r][c+idist])?sp:SOILPAVI->data[r][c+idist];
                        }else {
                            zx2 = dem;
                            spx2 = sp;
                        }
                        if(r > idist)
                        {
                            zy1 = pcr::isMV(DEMAVI->data[r-idist][c])?dem:DEMAVI->data[r-idist][c];
                            spy1 = pcr::isMV(SOILPAVI->data[r-idist][c])?sp:SOILPAVI->data[r-idist][c];

                        }else {
                            zy1 = dem;
                            spy1 =sp;
                        }
                        if(r < DEM->nrRows()-1-idist)
                        {
                            zy2 = pcr::isMV(DEMAVI->data[r+idist][c])?dem:DEMAVI->data[r+idist][c];
                            spy2 = pcr::isMV(SOILPAVI->data[r+idist][c])?sp:SOILPAVI->data[r+idist][c];
                        }else {
                            zy2 = dem;
                            spy2 = sp;

                        }
                    }



                    if(do_area)
                    {

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

                        sl = sl * 4.0/3.0;
                        sr = sr * 4.0/3.0;
                        float s = std::max(1e-6,std::fabs(0.5*(sl + sr)));

                        //save left-overs
                        float z_avg = (zl + zr + z)/3.0;
                        DEMAVI2->data[r][c] = z_avg;
                        SOILPAVI2->data[r][c] = ((spx1 + spx2 + spy1 + spy2)/2.0 + sp)/3.0;

                        if(i >0)
                        {
                            output_area[i]->data[r][c]  = output_area[i-1]->data[r][c] + s * dxi * dxi;
                        }else
                        {
                            output_area[i]->data[r][c] = s * dxi * dxi;
                        }
                    }
                    output_soilprop[i]->data[r][c] = SOILPAVI2->data[r][c];
                }
            }
        }

        cTMapf * temp = DEMAVI;
        DEMAVI = DEMAVI2;
        DEMAVI2 = temp;

        temp = SOILPAVI;
        SOILPAVI = SOILPAVI2;
        SOILPAVI2 = temp;
    }


    delete DEMAVI;
    delete DEMAVI2;
    delete SOILPAVI;
    delete SOILPAVI2;

    delete WidthI;


}

#endif // OSS_FLOWWIDTH_H
