#ifndef FULLFLOOD_GRIDADAPT_H
#define FULLFLOOD_GRIDADAPT_H

inline float GetMaxWaterHeightNB(std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l)
{
  if(r < 0 || c < 0 || r > flowstates[l][0]->Rows()-1 || c > flowstates[l][0]->Cols()-1)
  {
    return 0.0;
  }
 if(flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
    {
        float maxh = 0.0;
        for(int r2 = 0; r2 < 3; r2++)
        {
            for (int c2 = 0; c2 < 3; c2++)
            {
                maxh = std::max(maxh,GetMaxWaterHeightNB(flowstates,r * 3 + r2,c*3 + c2,l+1));
            }
        }
        return maxh;
    }else{

        bool do_super = false;
        int r2_ll = r/3;
        int c2_ll = c/3;
        if(l > 0 )
        {
            do_super = flowstates[l-1][2]->data[r2_ll][c2_ll] < 0.5;
        }

        //should we check the supergrid cell?
        if(do_super)
        {
          return GetMaxWaterHeightNB(flowstates,r2_ll,c2_ll,l-1);

        }else{

          if(r > -1 && c> 0 && r < flowstates[l][4]->Rows() && c < flowstates[l][4]->Cols())
          {
            return flowstates[l][4]->data[r][c];

          }else{
            return 0.0;
          }

        }
    }
    //backup, never reached
    return 0.0;

}

inline int GetMaxSubGridLevel(std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l)
{
    if(r < 0 || c < 0 || r > flowstates[l][0]->Rows()-1 || c > flowstates[l][0]->Cols()-1)
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
                maxl = std::max(maxl,GetMaxSubGridLevel(flowstates,r * 3 + r2,c*3 + c2,l+1));
            }
        }
        return maxl;
    }else{
      return l;
    }
}
inline void CheckSubGridDo(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates,float courant, int r, int c, int l, float scale_req =1.0)
{

    //only if this is not a cache-active gridcell



  if(r < 0 || c < 0 || r > flowstates[l][0]->Rows()-1 || c > flowstates[l][0]->Cols()-1)
  {
    return;
  }
  if(l >= flowstates.size()-1)
  {
    return;
  }
   if(pcr::isMV(flowstates[l][2]->data[r][c]))
   {
    return;
   }

   if(flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
   {
      for(int r2 = 0; r2 < 3; r2++)
      {
          for (int c2 = 0; c2 < 3; c2++)
          {
              CheckSubGridDo(_dt,dx / 3.0,flowstates,courant,r * 3 + r2,c*3 + c2,l+1, scale_req);
          }
      }
   }else if(l+1 < flowstates.size())
   {

        int max_c1 = GetMaxSubGridLevel(flowstates,r,c+1,l);
        max_c1 = std::max(max_c1,GetMaxSubGridLevel(flowstates,r,c-1,l));
        max_c1 = std::max(max_c1,GetMaxSubGridLevel(flowstates,r+1,c,l));
        max_c1 = std::max(max_c1,GetMaxSubGridLevel(flowstates,r-1,c,l));
      //ratio of flow height vs pixel size
      //at least max level of neighbors - 1
       //FF::cout << "sga  " << l << " " << flowstates[l][4]->data[r][c] << " " << 27.0/dx << " " << dx << std::endl;
      //check max neighbor flow height

      float h_c1 = GetMaxWaterHeightNB(flowstates,r,c+1,l);
      h_c1 = std::max(h_c1,GetMaxWaterHeightNB(flowstates,r,c-1,l));
      h_c1 = std::max(h_c1,GetMaxWaterHeightNB(flowstates,r+1,c,l));
      h_c1 = std::max(h_c1,GetMaxWaterHeightNB(flowstates,r-1,c,l));


      //FF::cout << "DO NOT activate sub-grid, wh "<< std::max(h_c1,flowstates[l][4]->data[r][c]) << " crit "  << (flowstates[l][11]->data[r][c] * scale_req/std::sqrt(dx)) << " " << scale_req/std::sqrt(dx) << " " << scale_req << " " << flowstates[l][11]->data[r][c] << "  lreq "<< (max_c1 > l +1) <<std::endl;


      if(std::max(h_c1,flowstates[l][4]->data[r][c]) > (std::min(std::min(1.0f,(0.05f * scale_req + 0.35f) * 10.0f / dx),flowstates[l][11]->data[r][c] * scale_req/std::sqrt(dx))) || max_c1 > l +1) //(3.0/std::sqrt(dx))
      {

        //FF::cout << "activate sub-grid, wh "<< std::max(h_c1,flowstates[l][4]->data[r][c]) << " crit "  << (flowstates[l][11]->data[r][c] * scale_req/std::sqrt(dx)) << " " << scale_req/std::sqrt(dx) << " " << scale_req << " " << flowstates[l][11]->data[r][c] << "  lreq "<< (max_c1 > l +1) <<std::endl;
        //activate sub-grid


        float h = flowstates[l][4]->data[r][c];
        float vx = flowstates[l][5]->data[r][c];
        float vy = flowstates[l][6]->data[r][c];
        flowstates[l][4]->data[r][c] = 0;
        flowstates[l][5]->data[r][c] = 0;
        flowstates[l][6]->data[r][c] = 0;
        flowstates[l][7]->data[r][c] = 0;
        flowstates[l][8]->data[r][c] = 0;
        flowstates[l][9]->data[r][c] = 0;

        flowstates[l][2]->data[r][c] = 1.0;

        for(int r2 = 0; r2 < 3; r2++)
        {
            for (int c2 = 0; c2 < 3; c2++)
            {
               int r3 = r*3 + r2;
               int c3 = c*3 + c2;

               if(r3 > -1 && c3 > -1 && r3 < flowstates[l+1][4]->Rows() && c3 < flowstates[l+1][4]->Cols())
               {
                flowstates[l+1][4]->data[r3][c3] = h;
                flowstates[l+1][5]->data[r3][c3] = vx;
                flowstates[l+1][6]->data[r3][c3] = vy;
                flowstates[l+1][7]->data[r3][c3] = h;
                flowstates[l+1][8]->data[r3][c3] = vx;
                flowstates[l+1][9]->data[r3][c3] = vy;
               }
            }
        }

      }
   }

}
static inline void FlowSubgriddify(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates,float courant = 0.1, float scale_req =1.0)
{
    //detect when we want to do subgridding
    //strong gradient in flow surface?
    //strong discharge?


    cTMapf * DEM = flowstates[0][0];
    for(int r = 0; r < DEM->nrRows(); r++)
    {
        for (int c = 0; c < DEM->nrCols(); c++)
        {
            if(flowstates[0][17] != nullptr)
            {
                if(flowstates[0][17]->data[r][c] > 0.5)
                {
                    continue;
                }
            }

            if(!pcr::isMV(DEM->data[r][c]))
            {
                //for all pixels in this supergrid pixel
                CheckSubGridDo(_dt,dx,flowstates,courant,r,c,0, scale_req);
            }
        }
    }

    //detect when we dont want subgridding any longer
}

inline bool HasSubGrid(std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l)
{
 if(r < 0 || c < 0 || r > flowstates[l][0]->Rows()-1 || c > flowstates[l][0]->Cols()-1)
  {
    return false;
  }

  if(flowstates[l][2]->data[r][c] > 0.5)
   {
    return true;
   }else{
    return false;
   }
}

inline void CheckSubGridRemove(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates,float courant, int r, int c, int l,float scale_req)
{

  if(r < 0 || c < 0 || r > flowstates[l][0]->Rows()-1 || c > flowstates[l][0]->Cols()-1)
  {
    return;
  }
  if(l > flowstates.size()-1)
  {
    return;
  }
   if(pcr::isMV(flowstates[l][2]->data[r][c]))
   {
    return;
   }
   if(flowstates[l][2]->data[r][c] < 0.5)
   {
    return;
   }

   if(flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
   {
      bool has_any_subgrid = false;
      for(int r2 = 0; r2 < 3; r2++)
      {
          for (int c2 = 0; c2 < 3; c2++)
          {

              CheckSubGridRemove(_dt,dx / 3.0,flowstates,courant,r * 3 + r2,c*3 + c2,l+1,scale_req);
              has_any_subgrid = has_any_subgrid || HasSubGrid(flowstates,r * 3 + r2,c*3 + c2,l+1);
          }
      }

      if(!has_any_subgrid)
      {
        //try to see if we need to remove this subgrid?

        int max_c1 = GetMaxSubGridLevel(flowstates,r,c+1,l);
        max_c1 = std::max(max_c1,GetMaxSubGridLevel(flowstates,r,c-1,l));
        max_c1 = std::max(max_c1,GetMaxSubGridLevel(flowstates,r+1,c,l));
        max_c1 = std::max(max_c1,GetMaxSubGridLevel(flowstates,r-1,c,l));

        float h_c1 = GetMaxWaterHeightNB(flowstates,r,c+1,l);
        h_c1 = std::max(h_c1,GetMaxWaterHeightNB(flowstates,r,c-1,l));
        h_c1 = std::max(h_c1,GetMaxWaterHeightNB(flowstates,r+1,c,l));
        h_c1 = std::max(h_c1,GetMaxWaterHeightNB(flowstates,r-1,c,l));


        if(!(std::max(h_c1,flowstates[l][4]->data[r][c]) > (std::min(std::min(1.0f,(0.05f * scale_req + 0.35f) * 10.0f / dx),flowstates[l][11]->data[r][c] * scale_req/std::sqrt(dx))) || max_c1 > l +1)) //(3.0/std::sqrt(dx))
        {

            //remove it
            flowstates[l][2]->data[r][c] = 0;

            //get average flow properties
            float h_av = 0.0; //flowstates[l][4]->data[r][c];
            float vx_av = 0.0; //flowstates[l][5]->data[r][c];
            float vy_av = 0.0; //flowstates[l][6]->data[r][c];

            float h_n = 0;

            for(int r2 = 0; r2 < 3; r2++)
            {
                for (int c2 = 0; c2 < 3; c2++)
                {
                  int r3 =  r * 3 + r2;
                  int c3 = c*3 + c2;

                  if(r3 > -1 && c3 > -1 && r3 < flowstates[l+1][4]->Rows() && c3 < flowstates[l+1][4]->Cols())
                  {

                    h_n += 1.0;
                    h_av +=flowstates[l+1][4]->data[r3][c3];
                    vx_av +=flowstates[l+1][5]->data[r3][c3];
                    vy_av +=flowstates[l+1][6]->data[r3][c3];

                    flowstates[l+1][4]->data[r3][c3] = 0.0;
                    flowstates[l+1][5]->data[r3][c3] = 0.0;
                    flowstates[l+1][6]->data[r3][c3] = 0.0;
                    flowstates[l+1][7]->data[r3][c3] = 0.0;
                    flowstates[l+1][8]->data[r3][c3] = 0.0;
                    flowstates[l+1][9]->data[r3][c3] = 0.0;

                  }
                }
            }


            h_av = h_av/std::max(1.0f,h_n);
            vx_av = vx_av/std::max(1.0f,h_n);
            vy_av = vy_av/std::max(1.0f,h_n);

            flowstates[l][4]->data[r][c] = h_av;
            flowstates[l][5]->data[r][c] = vx_av;
            flowstates[l][6]->data[r][c] = vy_av;
            flowstates[l][7]->data[r][c] = h_av;
            flowstates[l][8]->data[r][c] = vx_av;
            flowstates[l][9]->data[r][c] = vy_av;
        }
      }
   }

}

static inline void FlowUnSubgriddify(double _dt,float dx, std::vector<std::vector<cTMapf *>> &flowstates,float courant = 0.1, float scale_req = 1.0)
{
    //detect when we want to do subgridding
    //strong gradient in flow surface?
    //strong discharge?


    cTMapf * DEM = flowstates[0][0];
    for(int r = 0; r < DEM->nrRows(); r++)
    {
        for (int c = 0; c < DEM->nrCols(); c++)
        {
            if(flowstates[0][17] != nullptr)
            {
                if(flowstates[0][17]->data[r][c] > 0.5)
                {
                    continue;
                }
            }
            if(!pcr::isMV(DEM->data[r][c]))
            {
                //for all pixels in this supergrid pixel
                CheckSubGridRemove(_dt,dx,flowstates,courant,r,c,0,scale_req);
            }
        }
    }



}

#endif
