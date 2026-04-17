#ifndef OSS_DEMCHANNEL_H
#define OSS_DEMCHANNEL_H

#include "map.h"
#include "spreadflow.h"
#include <cmath>
#include <algorithm>
#include "spread.h"

static const int UF_DIRECTION_X = 1;
static const int UF_DIRECTION_Y = 2;
static const int UF_DIRECTION_XY = 3;
static const int UF_DERIVATIVE_LR = 0;
static const int UF_DERIVATIVE_L = 1;
static const int UF_DERIVATIVE_R = 2;

static bool UF_OUTORMV(cTMapf * map, int r, int c)
{
    if(r > -1 && c > -1 && r < map->nrRows() && c < map->nrCols())
    {
        return pcr::isMV(map->data[r][c]);
    }
    return true;
}

static double UF2D_Derivative_scaled(cTMapf * _dem, int r, int c, int direction,double scale, int calculationside)
{
    cTMapf * _in = _dem;
    float _dx = _dem->cellSizeX();
    if(UF_OUTORMV(_dem,r,c))
    {
        return 0;
    }
    if(direction == UF_DIRECTION_X)
    {
        if(calculationside == UF_DERIVATIVE_LR)
        {

            double dx1 = !(UF_OUTORMV(_dem,r,c+1))? scale *_in->data[r][c+1] -scale *_in->Drc :0.0;
            double dx2 = !(UF_OUTORMV(_dem,r,c-1))? scale *_in->Drc - scale *_in->data[r][c-1] :0.0;
            return (dx1 + dx2)/(2.0*_dx);

        }else if(calculationside == UF_DERIVATIVE_L)
        {
            return (!UF_OUTORMV(_dem,r,c-1))? (scale *_in->Drc -scale *_in->data[r][c-1])/_dx :0.0;

        }else if(calculationside == UF_DERIVATIVE_R)
        {
            return (!UF_OUTORMV(_dem,r,c+1))? (scale *_in->data[r][c+1] -scale *_in->Drc)/_dx :0.0;

        }
    }
    if(direction == UF_DIRECTION_Y)
    {
        if(calculationside == UF_DERIVATIVE_LR)
        {
            double dy1 = (!UF_OUTORMV(_dem,r+1,c))? scale *_in->data[r+1][c] -scale *_in->Drc :0.0;
            double dy2 = (!UF_OUTORMV(_dem,r-1,c))? scale *_in->Drc - scale *_in->data[r-1][c] :0.0;
            return (dy1 + dy2)/(2.0*_dx);

        }else if(calculationside == UF_DERIVATIVE_L)
        {
            return (!UF_OUTORMV(_dem,r-1,c))? (scale *_in->Drc - scale *_in->data[r-1][c])/_dx :0.0;

        }else if(calculationside == UF_DERIVATIVE_R)
        {
            return (!UF_OUTORMV(_dem,r+1,c))? (scale *_in->data[r+1][c] -scale *_in->Drc)/_dx :0.0;
        }
    }
    return 0;

}


inline std::vector<cTMapf*> GetDInfiniteNetwork(cTMapf * DEM, float delta, int iters)
{
    cTMapf * DEM2 = SpreadMonotonicReconstruct(DEM,delta,iters);//SpreadDirectionalAbsMaxMD(Seeds,DEM,GradX1,GradX2,GradY1,GradY2);
    cTMapf * GradX1 = DEM->GetCopy0();
    cTMapf * GradX2 = DEM->GetCopy0();
    cTMapf * GradY1 = DEM->GetCopy0();
    cTMapf * GradY2 = DEM->GetCopy0();

    //get new normalized flux network
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++) {
        for (int c = 0; c < DEM->Cols(); c++) {
            if (!pcr::isMV(DEM->data[r][c])) {

                float gx1 =std::min(0.0,-UF2D_Derivative_scaled(DEM2,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_L));
                float gx2 =std::max(0.0,-UF2D_Derivative_scaled(DEM2,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_R));
                float gy1 =std::min(0.0,-UF2D_Derivative_scaled(DEM2,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_L));
                float gy2 =std::max(0.0,-UF2D_Derivative_scaled(DEM2,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_R));

                float gtot = std::max(1e-12f,std::fabs(gx1) + std::fabs(gx2) + std::fabs(gy1) + std::fabs(gy2));

                GradX1->data[r][c] = gx1/gtot;
                GradX2->data[r][c] = gx2/gtot;
                GradY1->data[r][c] = gy1/gtot;
                GradY2->data[r][c] = gy2/gtot;

            } else
            {
                pcr::setMV(GradX1->data[r][c]);
                pcr::setMV(GradX2->data[r][c]);
                pcr::setMV(GradY1->data[r][c]);
                pcr::setMV(GradY2->data[r][c]);
            }
        }
    }

    return {GradX1,GradX2,GradY1,GradY2};
}

/**
 * @brief Applies Conditional LDD4 Correction to direct flow correctly across DEM.
 *
 * Flow accumulation is performed to estimate local drainage, then cells that pass
 * the drainage channel threshold are corrected to enforce an LDD (Local Drain Direction)
 * 4-connectivity consistency, forcing flow into the dominant gradient direction.
 *
 * @param DEM The digital elevation model map.
 * @param GradX1 Left gradient component.
 * @param GradX2 Right gradient component.
 * @param GradY1 Top gradient component.
 * @param GradY2 Bottom gradient component.
 * @param iters Number of iterations for the flow spread.
 * @param dx Cell size resolution.
 * @param min_drainchancor Threshold area above which correction is applied.
 * @param RAININ External rain/water source input map.
 */
inline void ApplyConditionalLDD4Correction(cTMapf* DEM, cTMapf* GradX1, cTMapf* GradX2, cTMapf* GradY1, cTMapf* GradY2, int iters, float dx, float min_drainchancor = 0.0, cTMapf* RAININ = 0)
{
    cTMapf * rtemp = RAININ ? RAININ : GradX1->GetCopy1();

    // Accumulate the flow to determine drainage strength
    cTMapf * qdrain = SpreadFlowMD2(rtemp, GradX1, GradX2, GradY1, GradY2, iters);

    double cellareakm = (dx * dx) / 1000000.0;
    for (int r = 0; r < GradX1->Rows(); r++) {
        for (int c = 0; c < GradX1->Cols(); c++) {
            if (!pcr::isMV(DEM->data[r][c])) {

                // If drainage area is above threshold, correct the direction to LDD-4
                if(qdrain->data[r][c] * cellareakm > min_drainchancor)
                {
                    float gx1 = std::fabs(GradX1->data[r][c]);
                    float gx2 = std::fabs(GradX2->data[r][c]);
                    float gy1 = std::fabs(GradY1->data[r][c]);
                    float gy2 = std::fabs(GradY2->data[r][c]);

                    GradX1->data[r][c] = 0.0;
                    GradX2->data[r][c] = 0.0;
                    GradY1->data[r][c] = 0.0;
                    GradY2->data[r][c] = 0.0;

                    int retldd = 0;
                    if(gx1 > 0.0) {
                        if(gy1 > 0.0) {
                            if(gy1 > gx1) retldd = 2;
                            else retldd = 4;
                        } else if(gy2 > 0.0) {
                            if(gy2 > gx1) retldd = 8;
                            else retldd = 4;
                        } else {
                            retldd = 4;
                        }
                    } else if(gx2 > 0.0) {
                        if(gy1 > 0.0) {
                            if(gy1 > gx2) retldd = 2;
                            else retldd = 6;
                        } else if(gy2 > 0.0) {
                            if(gy2 > gx2) retldd = 8;
                            else retldd = 6;
                        } else {
                            retldd = 6;
                        }
                    } else if(gy1 > 0.0) {
                        retldd = 2;
                    } else if(gy2 > 0.0) {
                        retldd = 8;
                    }

                    if(retldd == 8) GradY2->data[r][c] = 1.0;
                    else if(retldd == 6) GradX2->data[r][c] = 1.0;
                    else if(retldd == 2) GradY1->data[r][c] = -1.0;
                    else if(retldd == 4) GradX1->data[r][c] = -1.0;
                }
            }
        }
    }

    if(!RAININ) delete rtemp;
    delete qdrain;
}

#endif
