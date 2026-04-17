#ifndef FULLFLOOD_KINEMATIC_H
#define FULLFLOOD_KINEMATIC_H

#include "../algorithms/rasteralgorithms/ldd.h"

static int ch_dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
static int ch_dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
#define GRAV 9.81


#define LISEM_KIN_MIN_FLUX 1e-7
#define LISEM_KIN_MAX_ITERS 100

inline static double IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX)
{
    /* Using Newton-Raphson Method */
    double  ab_pQ, deltaTX, C;  //auxillary vars
    int   count;
    double Qkx; //iterated discharge, becomes Qnew
    double fQkx = 1.0; //function
    double dfQkx;  //derivative
    const double _epsilon = 1e-12;
    const double beta = 0.6;

    /* common terms */
     // ab_pQ = alpha*beta*pow(((Qold+Qin)/2),beta-1);
    // derivative of diagonal average (space-time)

    deltaTX = deltaT/deltaX;
    C = deltaTX*Qin + alpha*pow(Qold,beta) + deltaT*q;
    //dt/dx*Q = m3/s*s/m=m2; a*Q^b = A = m2; q*dt = s*m2/s = m2
    //C is unit volume of water
    // can be negative because of q

    // if C < 0 than all infiltrates, return 0, if all fluxes 0 then return
    if (C < 0)
    {
        return(0);
    }

    // pow function sum flux must be > 0
    if (Qold+Qin > 0)
    {
        ab_pQ = alpha*beta*pow(((Qold+Qin)/2),beta-1);
        // derivative of diagonal average (space-time), must be > 0 because of pow function
        Qkx = (deltaTX * Qin + Qold * ab_pQ + deltaT * q) / (deltaTX + ab_pQ);
        // explicit first guess Qkx, VERY important
        Qkx   = std::max(Qkx, 0.0);
    }
    else
        Qkx =  0;

    Qkx   = std::isnan(Qkx) ? 0.0 : std::max(Qkx, 0.0);
    if (Qkx < LISEM_KIN_MIN_FLUX)
        return(0);

    // avoid spurious iteration
    count = 0;
    do {
        fQkx  = deltaTX * Qkx + alpha * pow(Qkx, beta) - C;   /* Current k */ //m2
        dfQkx = deltaTX + alpha * beta * pow(Qkx, beta - 1);  /* Current k */
        Qkx   -= fQkx / dfQkx;                                /* Next k */

        Qkx   = std::isnan(Qkx) ? 0.0 : std::max(Qkx, 0.0);
        count++;
    } while(fabs(fQkx) > _epsilon && count < LISEM_KIN_MAX_ITERS);

    return Qkx;
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::Kinematic(int pitRowNr, int pitColNr, cTMapf *_LDD,cTMapf *_Q, cTMapf *_Qn, cTMapf *_q, cTMapf *_Alpha, cTMapf *_DX,cTMapf *_Vol,cTMapf *_StorVol)
 * @brief Spatial implementation of the kinematic wave
 *
 * Kinematic wave spatial part, used for slope, channel and tiledrain system:
 * Kinematic is called for each pit (i.e. a cell with value 5 in the LDD):
 * A linked list of all cells connected to the pit is made, after that it 'walks' through the list
 * calculating the water fluxes from upstream to downstream.
 *
 * @param pitRowNr : row nr of the current pit (can be more than one, LISEM loops through all the pits)
 * @param pitColNr : col nr of the current pit (can be more than one, LISEM loops through all the pits)
 * @param _LDD : LDD map, can be slope, channel or tiledrain system
 * @param _Q : incoming discharge at the start of the timestep (m3/s)
 * @param _Qn : Outgoing new discharge at the end of the timestep (m3/s)
 * @param _q : infiltration surplus from infiltration functions, 0 or negative value
 * @param _Alpha : a in A=aQ^b
 * @param _DX : dx corrected for slope
 * @param _Vol : water volume in cell (m3)
 * @param _StorVol : water volume in buffers (m3)
 * @return new water discharge
 *
 * @see TWorld::IterateToQnew
 * @see TWorld::LDD
 */
inline static void Kinematic(int pitRowNr, int pitColNr, cTMapf *_LDD,
                       cTMapf *_Q, cTMapf *_Qn,
                       cTMapf *_q, cTMapf *_Alpha, float _DX,
                        double _dt)//,cTMapf *_StorVol)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    /// Linked list of cells in order of LDD flow network, ordered from pit upwards
    LDD_LINKEDLIST *list = nullptr, *temp = nullptr;
    list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

    list->prev = nullptr;
    /// start gridcell: outflow point of area
    list->rowNr = pitRowNr;
    list->colNr = pitColNr;

    //    if (SwitchErosion)
    //        _Qsn->fill(0);
    // set output sed flux to 0
    //   _Qn->setMV();
    // flag all Qn gridcell with MV
    // does not work with multtple pits because earlier iterations are set to zero

    while (list != nullptr)
    {
        int i = 0;
        bool  subCachDone = true; // are sub-catchment cells done ?
        int rowNr = list->rowNr;
        int colNr = list->colNr;

        /** put all points that have to be calculated to calculate the current point in the list,
         before the current point */
        for (i=1; i<=9; i++)
        {
            int r, c;
            int ldd = 0;

            // this is the current cell
            if (i==5)
                continue;

            r = rowNr+dy[i];
            c = colNr+dx[i];

            if (INSIDE(_LDD,r, c) && !pcr::isMV(_LDD->Drc))
                ldd = (int) _LDD->Drc;
            else
                continue;

            // check if there are more cells upstream, if not subCatchDone remains true
            if (pcr::isMV(_Qn->Drc) &&
                    FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                    INSIDE(_LDD,r, c))
            {
                temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                temp->prev = list;
                list = temp;
                list->rowNr = r;
                list->colNr = c;
                subCachDone = false;
            }
        }

        // all cells above a cell are linked in a "sub-catchment or branch
        // continue with water and sed calculations
        // rowNr and colNr are the last upstream cell linked
        if (subCachDone)
        {
            double Qin=0.0;//, Sin=0.0;

            // for all incoming cells of a cell, sum Q and Sed and put in Qin and Sin
            // note these are values of Qn and Qsn so the new flux
            for (i=1;i<=9;i++)
            {
                int r, c, ldd = 0;

                if (i==5)  // Skip current cell
                    continue;

                r = rowNr+dy[i];
                c = colNr+dx[i];

                if (INSIDE(_LDD,r, c) && !pcr::isMV(_LDD->Drc))
                    ldd = (int) _LDD->Drc;
                else
                    continue;

                if (INSIDE(_LDD,r, c) &&
                        FLOWS_TO(ldd, r,c,rowNr, colNr) &&
                        !pcr::isMV(_LDD->Drc) )
                {
                    Qin += _Qn->Drc;
                    //QinKW->data[rowNr][colNr] = Qin;

                    /*if (SwitchErosion)
                        Sin += _Qsn->Drc;*/
                }
            }



            _Qn->data[rowNr][colNr] = IterateToQnew(Qin, _Q->data[rowNr][colNr], _q->data[rowNr][colNr],
                                                    _Alpha->data[rowNr][colNr], _dt, _DX);

            if(std::isnan(_Qn->data[rowNr][colNr]))
            {
                //FF::cout << " iterate to new q: " << Qin << " " << _Q->data[rowNr][colNr] << " " << _q->data[rowNr][colNr] << " " << _Alpha->data[rowNr][colNr] << " " <<_dt << " " << _DX<<std::endl;

            }

                                                    //FF::cout << _Qn->data[rowNr][colNr] << std::endl;
            // Newton Rapson iteration for water of current cell

              /* cell rowN, colNr is now done */

            temp=list;
            list=list->prev;
            free(temp);
            // go to the previous cell in the list

        }/* eof subcatchment done */
    } /* eowhile list != nullptr */
}


inline static void FlowKinematic(cTMapf * LDD, cTMapf * Slope, cTMapf * Width, float man, cTMapf * H, cTMapf * Q, float dt, float dx,cTMapf * Alpha,cTMapf * Qn,cTMapf * infil)
{
    //FF::cout << "Flow Kinematic " << std::endl;
    //convert height to discharge

    for(int r = 0; r < Alpha->nrRows(); r++)
    {
        for(int c = 0; c < Alpha->nrCols(); c++)
        {
            Q->data[r][c] = 0.0;
            Alpha->data[r][c] = 0.0;
            Qn->data[r][c] = 0.0;
            infil->data[r][c] = 0.0;
        }
    }

    for(int r = 0; r < Alpha->nrRows(); r++)
    {
        for(int c = 0; c < Alpha->nrCols(); c++)
        {
            if(!pcr::isMV(LDD->data[r][c]) && Width->data[r][c] > 1e-6)
            {

                double Perim, Radius, Area;
                const double beta = 0.6;
                const double _23 = 2.0/3.0;
                double beta1 = 1.0/beta;
                double wh = H->data[r][c];
                double FW = Width->data[r][c];
                double grad = sqrt(sin(atan(std::max(Slope->data[r][c],1e-12f))));

                {
                    Perim = FW + 2.0*wh;
                    Area = FW*wh;
                }

                Radius = (Perim > 0 ? Area/Perim : 0);

                Alpha->Drc = std::pow(man/grad * std::pow(Perim, _23),beta);

                if (Alpha->Drc > 0)
                    Q->Drc = std::pow(Area/Alpha->Drc, beta1);
                else
                    Q->Drc = 0;

                Q->Drc = std::max(1e-6f,Q->Drc);

                pcr::setMV(Qn->Drc);
            }

        }

    }
    //FF::cout << "Flow Kinematic2" << std::endl;

    //now run the kinematic algorithm to get new q

    for(int r = 0; r < LDD->nrRows(); r++)
    {
        for(int c = 0; c < LDD->nrCols(); c++)
        {
            if(!pcr::isMV(LDD->data[r][c]))
            {

                if((int) LDD->data[r][c] == 5)
                {

                    Kinematic(r,c,LDD,Q,Qn,infil,Alpha,dx,dt);

                }
            }
        }
    }
 //FF::cout << "Flow Kinematic3" << std::endl;
    //now convert back to flow height

    for(int r = 0; r < LDD->nrRows(); r++)
    {
        for(int c = 0; c < LDD->nrCols(); c++)
        {
            if(!pcr::isMV(LDD->data[r][c]) && Width->data[r][c] > 1e-6)
            {

                Q->Drc = (Alpha->Drc*pow(Qn->Drc, 0.6))/Width->data[r][c];

                if(std::isnan(Q->Drc))
                {
                    //FF::cout << "set h " << Alpha->Drc << " " << Qn->Drc << " " << Width->data[r][c] << " " << Slope->data[r][c] << " " << H->data[r][c] << " " << dt << " " << dx << std::endl;

                    Q->Drc = 0;
                }

            }
        }
    }

}
#endif // FULLFLOOD_KINEMATIC_H
