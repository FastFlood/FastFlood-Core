#ifndef FULLFLOOD_COURANT_H
#define FULLFLOOD_COURANT_H

inline void SetMinDT3Casc(double tglob, double dttot, double _dt,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, double dtmin, double dtleft, float * ref_sp = nullptr)
{
    //FF::cout << "MinDT  " << r << " " << c << " " << l <<std::endl;

    if(!(r > -1 && c > -1 && r < flowstates[l][0]->nrRows() && c < flowstates[l][0]->nrCols()))
    {
        return;
    }

    bool is_cache = false;

    if(l == 0 && flowstates[l][17] != nullptr)
    {
        if(flowstates[l][17]->data[r][c] > 0.5f)
        {
            is_cache = true;
        }
    }

   if(!is_cache && flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
   {
    double dtmin = _dt;
      for(int r2 = 0; r2 < 3; r2++)
      {
          for (int c2 = 0; c2 < 3; c2++)
          {
              SetMinDT3Casc(tglob,dttot,_dt,flowstates, r * 3 + r2,c*3 + c2, l+1,dtmin,dtleft,ref_sp);
          }
      }

      return ;

   }else{

    //do calculations here

    bool has_neighbor_dt = false;
    bool has_self_dt = tglob + 1.0 * _dt > (flowstates[l][12]->data[r][c]);

    auto dtlambda = [&](int r2,int c2,int l2,float dy, int ca, int ca2)
    {
        if(!(r2 > -1 && c2 > -1 && r2 < flowstates[l2][0]->nrRows() && c2 < flowstates[l2][0]->nrCols()))
        {
            return;
        }

        //is this cell active itself?
        bool has_self_dt_nb = tglob + 1.0 * _dt > (flowstates[l2][12]->data[r2][c2]);
        if(has_self_dt_nb)
        {
            has_neighbor_dt = true;
        }
    };

    ForEach3CascNeighborCachedSuperCells(_dt,1,flowstates, r, c, l, r,c+1,l,1,0,false, dtlambda,5,0);
    ForEach3CascNeighborCachedSuperCells(_dt,1,flowstates, r, c, l, r,c-1,l,-1,0,false, dtlambda,5,0);
    ForEach3CascNeighborCachedSuperCells(_dt,1,flowstates, r, c, l, r+1,c,l,0,1,false, dtlambda,6,0);
    ForEach3CascNeighborCachedSuperCells(_dt,1,flowstates, r, c, l, r-1,c,l,0,-1,false, dtlambda,6,0);

    float t_initialguess = flowstates[l][13]->data[r][c];
    if(has_self_dt)
    {
        flowstates[l][13]->data[r][c] = std::min(dttot,std::max((double)dtmin,std::min(dttot - tglob,(double)flowstates[l][13]->data[r][c])));
        *ref_sp = 1.0;
    }else if(has_neighbor_dt)
    {
        flowstates[l][13]->data[r][c] = std::max(dtmin,std::min(dttot-tglob,(double)flowstates[l][12]->data[r][c]) - std::min(dttot,tglob));
        *ref_sp = 1.0;
    }else
    {
        flowstates[l][13]->data[r][c] = 0.0;
    }

    if((r == 15 && c== 15 && l ==0) || (r == 15 && c== 15 && l ==0))
    {
        //FF::cout << "TIMESTEP CALC FIN SET "<< r << " "<< c << " "<< flowstates[l][13]->data[r][c] << " tglob: " << tglob  << " dttot " << dttot << " tg " << t_initialguess << std::endl;

    }


    //flowstates[l][13]->data[r][c] = std::min(dtmin,dtleft);

    //FF::cout << "timestep final: " << tglob << "  " << _dt << "  " <<  flowstates[l][12]->data[r][c] << " " <<  flowstates[l][13]->data[r][c]  << "  " << has_self_dt << " " << has_neighbor_dt << std::endl;
    //flowstates[l][47]->data[r][c] = flowstates[l][13]->data[r][c];

    return ;
   }

  return ;
}


inline double GetMinDT3Casc(double tglob,double courant, double _dt,float dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, double dtmin, double dtleft, float * ref_tnextmin)
{
  //FF::cout << "MinDT  " << r << " " << c << " " << l <<std::endl;
    if(!(r > -1 && c > -1 && r < flowstates[l][0]->nrRows() && c < flowstates[l][0]->nrCols()))
    {
        return _dt;
    }

    //
    bool is_cache = false;

    if(l == 0 && flowstates[l][17] != nullptr)
    {
        if(flowstates[l][17]->data[r][c] > 0.5f)
        {
            is_cache = true;
        }
    }

   if(!is_cache && flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
   {

       if(l == 0)
       {

       }

      double dtmin = _dt;
      for(int r2 = 0; r2 < 3; r2++)
      {
          for (int c2 = 0; c2 < 3; c2++)
          {
              double dtminh = GetMinDT3Casc(tglob,courant, _dt,dx/ 3.0,flowstates, r * 3 + r2,c*3 + c2, l+1,dtmin,dtleft, ref_tnextmin);
              dtmin = std::min(dtminh,dtmin);
          }
      }

      return dtmin;

   }else{

    double dt_req_final = _dt;
    //do calculations here
    if((r == 688 && c== 591 && l ==3)  || (r == 687 && c== 591 && l ==3))
    {

        //FF::cout << "start dt is  " << dt_req_final <<std::endl;
    }
    auto dtlambda = [&](int r2,int c2,int l2,float dy, int ca, int ca2)
    {
        if(!(r2 > -1 && c2 > -1 && r2 < flowstates[l2][0]->nrRows() && c2 < flowstates[l2][0]->nrCols()))
        {
            return;
        }
        float z = flowstates[l][0]->data[r][c];
        float z_x1 = z;

        float h = flowstates[l][4]->data[r][c];

        float h_x1 = h;

        float vx = flowstates[l][5]->data[r][c];

        float vy = flowstates[l][6]->data[r][c];


        if(pcr::isMV(z))
        {
          h = 0;
          vx= 0;
          vy = 0;
        }

        float vh = 0.0;
        if(r2 > -1 && c2 > -1 && r2 < flowstates[l2][0]->Rows() && c2 < flowstates[l2][0]->Cols() )
        {
          z_x1 = flowstates[l2][0]->data[r2][c2];
          h_x1 = flowstates[l2][4]->data[r2][c2];
          vh = flowstates[l2][ca]->data[r2][c2];
        }


        float zc_x1 = z_x1 < -99995? z : z_x1;


        float hw_x1 = h > h_x1? std::max(0.0f,std::min(h,h+z - (h_x1 +z_x1))) : 0.0;
        float vmax = 50.0;

        vx = std::max(std::sqrt(9.81f * std::fabs(hw_x1)),std::min(vmax,std::max(-vmax,vx)));
        vy = std::max(std::sqrt(9.81f * std::fabs(hw_x1)),std::min(vmax,std::max(-vmax,vy)));

        double dt_req = courant *(0.5*(dx+dy))/( std::min(1000.0f,std::max(std::max(std::fabs(vh),std::max(0.01f,9.81f * std::fabs(h))),(sqrt(vx*vx + vy * vy)))));
        dt_req_final = std::min(dt_req_final,dt_req);



        if((r == 688 && c== 591 && l ==3)  || (r == 687 && c== 591 && l ==3) || (r2 == 688 && c2== 591 && l2 ==3))
        {
            //FF::cout << "TIMESTEP CALC "<< std::endl;
            //FF::cout << "fsolve " << r << " " << c << " " << l << " loc2: " << r2 << " " << c2 << " " << l2 << " h:" <<  h << " hx1:" << h_x1 << " dt: " << dt_req <<  "  dtnow:L " << dt_req_final << std::endl;

        }


    };

    ForEach3CascNeighborCached(_dt,dx,flowstates, r, c, l, r,c+1,l,1,0,false, dtlambda,5,0);
    ForEach3CascNeighborCached(_dt,dx,flowstates, r, c, l, r,c-1,l,-1,0,false, dtlambda,5,0);
    ForEach3CascNeighborCached(_dt,dx,flowstates, r, c, l, r+1,c,l,0,1,false, dtlambda,6,0);
    ForEach3CascNeighborCached(_dt,dx,flowstates, r, c, l, r-1,c,l,0,-1,false, dtlambda,6,0);


    //potentially set this value
    if(flowstates[l][21] != nullptr)
    {
        double dtleftn = std::max(0.0,(_dt)- flowstates[l][12]->data[r][c]);
        flowstates[l][13]->data[r][c] = std::min(dtleftn,std::max(dtmin,dt_req_final));

        if((r == 15 && c== 15 && l ==0) || (r == 15 && c== 15 && l ==0))
        {
            //FF::cout << "TIMESTEP CALC FIN  "<< r << " "<< c << " "<< flowstates[l][13]->data[r][c] << " dtmin " << dtmin << " there "<< flowstates[l][12]->data[r][c] << " _dt " << _dt <<  " dtleft " << dtleftn << std::endl;

        }



        if(ref_tnextmin != nullptr){
            *ref_tnextmin = std::min(flowstates[l][13]->data[r][c],*ref_tnextmin);
        }

        //FF::cout << "timestep " << tglob << "  " << dtleft << " " << dtleftn  << "  "  <<  flowstates[l][12]->data[r][c] << " " <<  flowstates[l][13]->data[r][c] << "  " << dt_req_final << " " << dx << "  "  << dtleft << std::endl;
        //flowstates[l][47]->data[r][c] = dt_req_final;

        //flowstates[l][47]->data[r][c] = flowstates[l][12]->data[r][c];
    }

    return dt_req_final;
   }

  return _dt;
}
#endif
