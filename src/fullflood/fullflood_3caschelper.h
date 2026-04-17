#ifndef FULLFLOOD_3CASCHELPER_H
#define FULLFLOOD_3CASCHELPER_H


template<typename F>
void ForEach3CascEdge(double _dt,float _dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, int r2, int c2, int l2, int dx, int dy, bool include_edge, const F & lambda,int ca, int ca2, int r0, int c0)
{
   //check
    //FF::cout << "for reach edge " << r2 << "  " << c2 << "  "<<  l2 << " " << std::endl;
    if(r2 < 0 || c2 <0 || r2 >= flowstates[l2][0]->nrRows() || c2 >= flowstates[l2][0]->nrCols())
    {
        return;
    }
    bool do_this = true;

    if(r2+dy > -1 && c2+dx > -1 && r2+dy < flowstates[l2][0]->Rows() && c2+dx < flowstates[l2][0]->Cols())
    {
        if(pcr::isMV(flowstates[l2][0]->data[r2+dy][c2+dx]) && !include_edge)
        {
            do_this=false;
        }
    }
     if(do_this)
     {
        if(flowstates[l2][2]->data[r2][c2] > 0.5 && l < flowstates.size() -1)
        {
            //do the edges in the specific direction
            ForEach3CascEdge(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1 + dy-dx, c2*3+1 +dx-dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2, r0*3+1, c0*3+1);
            ForEach3CascEdge(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1 + dy   , c2*3+1 +dx   , l2+1,dx,dy,include_edge, (lambda),ca, ca2, r0*3+1, c0*3+1);
            ForEach3CascEdge(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1 + dy+dx, c2*3+1 +dx+dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2, r0*3+1, c0*3+1);

        }else
        {
            if(r2+dy > -1 && c2+dx > -1 && r2+dy < flowstates[l2][0]->Rows() && c2+dx < flowstates[l2][0]->Cols())
            {
                //only do the corners once
                if(std::abs(r2 - r0) == std::abs(c2-c0))
                {
                    //if(dx == 0)
                    {
                        lambda(r2,c2,l2,_dx,ca, ca2);
                    }
                }else
                {
                    lambda(r2,c2,l2,_dx,ca, ca2);
                }
            }else if(include_edge)
            {
               //get z's, h's, v's
                if(std::abs(r2 - r0) == std::abs(c2-c0))
                {
                    //if(dx == 0)
                    {
                        lambda(r2,c2,l2,_dx,ca, ca2);
                    }
                }else
                {
                    lambda(r2,c2,l2,_dx,ca, ca2);
                }
            }
        }
    }

}

template<typename F>
void ForEach3CascEdgeP(double _dt,float _dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, int r2, int c2, int l2, bool include_edge, const F & lambda)
{
    std::vector<int> dxi = {1,-1,0,0};
    std::vector<int> dyi = {0,0,1,-1};
    std::vector<int> cai = {1,1,2,2};
    std::vector<int> ca2i = {1,-1,1,-1};

    if(!(flowstates[l2][2]->data[r2][c2] > 0.5 && l < flowstates.size() -1))
    {
        lambda(r2,c2,l2,_dx,0,0);
        return;
    }

    for(int i = 0; i < 4; i++)
    {
        int dx = dxi[i];
        int dy = dyi[i];
        int ca = cai[i];
        int ca2 = ca2i[i];


        bool do_this = true;

        if(r2+dy > -1 && c2+dx > -1 && r2+dy < flowstates[l2][0]->Rows() && c2+dx < flowstates[l2][0]->Cols())
        {
            if(pcr::isMV(flowstates[l2][0]->data[r2+dy][c2+dx]) && !include_edge)
            {
                do_this=false;
            }
        }
         if(do_this)
         {
             if(flowstates[l2][2]->data[r2][c2] > 0.5 && l < flowstates.size() -1)
             {
                 //do the edges in the specific direction
                 ForEach3CascEdge(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1 + dy-dx, c2*3+1 +dx-dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2, r*3+1, c*3+1);
                 ForEach3CascEdge(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1 + dy   , c2*3+1 +dx   , l2+1,dx,dy,include_edge, (lambda),ca, ca2, r*3+1, c*3+1);
                 ForEach3CascEdge(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1 + dy+dx, c2*3+1 +dx+dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2, r*3+1, c*3+1);

             }else
             {
                 if(r2+dy > -1 && c2+dx > -1 && r2+dy < flowstates[l2][0]->Rows() && c2+dx < flowstates[l2][0]->Cols())
                 {
                     lambda(r2,c2,l2,_dx,ca, ca2);
                 }else if(include_edge)
                 {
                    //get z's, h's, v's
                    lambda(r2,c2,l2,_dx,ca, ca2);
                 }

             }
         }

    }
}


template<typename F>
void ForEach3CascNeighborCachedSuperCells(double _dt,float _dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, int r2, int c2, int l2, int dx, int dy, bool include_edge, const F & lambda,int ca, int ca2)
{
   //check
   //FF::cout << "MinDTnb  " << r2 << " " << c2 << " " << l2 <<std::endl;
   if(r2 > -1 && c2 > -1 && r2 < flowstates[l2][0]->Rows() && c2 < flowstates[l2][0]->Cols())
   {
    if(!pcr::isMV(flowstates[l2][0]->data[r2][c2]) || include_edge)
    {
        //check if this neighbor is part of a cache pixel
        if(flowstates[0][17] != nullptr)
        {
            int div = std::max(1,(int)std::pow(3,l2-1));
            int rb = r2/div;
            int cb = c2/div;

            if(rb < flowstates[0][17]->nrRows() && cb < flowstates[0][17]->nrCols())
            {
                if(flowstates[0][17]->data[rb][cb] > 0.5f)
                {
                    lambda(rb,cb,0,_dx*div,ca, ca2);
                    return;
                }
            }
        }
        //check if we have the a neighbor which is in another superpixel, that is not subgridded
        bool do_super = false;
        int r2_ll = r2/3;
        int c2_ll = c2/3;
        if(l2 > 0 )
        {
            do_super = flowstates[l2-1][2]->data[r2_ll][c2_ll] < 0.5;
        }
        if(do_super) //should we use the supergrid for this neighbor
        {
          ForEach3CascNeighbor(_dt,_dx*3.0,flowstates,r,c,l, (r2/3), (c2/3), l2-1,dx,dy, include_edge, (lambda),ca, ca2);
        }
        //is the neighbor not subgridded, or are we in the last grid
        else if(!(flowstates[l2][2]->data[r2][c2] > 0.5)  || !(l2 < flowstates.size()-1))
        {
          //call lambda
          lambda(r2,c2,l2,_dx,ca, ca2);
        }else{ //if the neighbor is further subgridded
          //call this same function three times, for each neighboring sub-pixel at the boundary
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy+dx, c2*3+1-dx+dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2);
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy   , c2*3+1-dx,    l2+1,dx,dy,include_edge, (lambda),ca, ca2);
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy-dx, c2*3+1-dx-dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2);
        }
    }
   }else if(include_edge)
   {
      //get z's, h's, v's
      lambda(r2,c2,l2,_dx,ca, ca2);
   }

}

template<typename F>
void ForEach3CascNeighborCached(double _dt,float _dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, int r2, int c2, int l2, int dx, int dy, bool include_edge, const F & lambda,int ca, int ca2)
{
   //check
   //FF::cout << "MinDTnb  " << r2 << " " << c2 << " " << l2 <<std::endl;
   if(r2 > -1 && c2 > -1 && r2 < flowstates[l2][0]->Rows() && c2 < flowstates[l2][0]->Cols())
   {
    if(!pcr::isMV(flowstates[l2][0]->data[r2][c2]) || include_edge)
    {
        //check if this neighbor is part of a cache pixel
        if(flowstates[0][17] != nullptr)
        {
            int div = std::max(1,(int)std::pow(3,l2-1));
            int rb = r2/div;
            int cb = c2/div;

            if(rb < flowstates[0][17]->nrRows() && cb < flowstates[0][17]->nrCols())
            {
                if(flowstates[0][17]->data[rb][cb] > 0.5f)
                {
                    ForEach3CascNeighbor(_dt,_dx*3.0,flowstates,r,c,l, rb, cb, 0 /*l0*/,dx,dy, include_edge, (lambda),ca, ca2);
                    return;
                }
            }
        }
        //check if we have the a neighbor which is in another superpixel, that is not subgridded
        bool do_super = false;
        int r2_ll = r2/3;
        int c2_ll = c2/3;
        if(l2 > 0 )
        {
            do_super = flowstates[l2-1][2]->data[r2_ll][c2_ll] < 0.5;
        }
        if(do_super) //should we use the supergrid for this neighbor
        {
          ForEach3CascNeighbor(_dt,_dx*3.0,flowstates,r,c,l, (r2/3), (c2/3), l2-1,dx,dy, include_edge, (lambda),ca, ca2);
        }
        //is the neighbor not subgridded, or are we in the last grid
        else if(!(flowstates[l2][2]->data[r2][c2] > 0.5)  || !(l2 < flowstates.size()-1))
        {
          //call lambda
          lambda(r2,c2,l2,_dx,ca, ca2);
        }else{ //if the neighbor is further subgridded
          //call this same function three times, for each neighboring sub-pixel at the boundary
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy+dx, c2*3+1-dx+dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2);
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy   , c2*3+1-dx,    l2+1,dx,dy,include_edge, (lambda),ca, ca2);
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy-dx, c2*3+1-dx-dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2);
        }
    }
   }else if(include_edge)
   {
      //get z's, h's, v's
      lambda(r2,c2,l2,_dx,ca, ca2);
   }

}

template<typename F>
void ForEach3CascNeighbor(double _dt,float _dx,std::vector<std::vector<cTMapf *>> &flowstates, int r, int c, int l, int r2, int c2, int l2, int dx, int dy, bool include_edge, const F & lambda,int ca, int ca2)
{
   //check
   //FF::cout << "MinDTnb  " << r2 << " " << c2 << " " << l2 <<std::endl;
   if(r2 > -1 && c2 > -1 && r2 < flowstates[l2][0]->Rows() && c2 < flowstates[l2][0]->Cols())
   {
    if(!pcr::isMV(flowstates[l2][0]->data[r2][c2]) || include_edge)
    {
        //check if we have the a neighbor which is in another superpixel, that is not subgridded
        bool do_super = false;
        int r2_ll = r2/3;
        int c2_ll = c2/3;
        if(l2 > 0 )
        {
            do_super = flowstates[l2-1][2]->data[r2_ll][c2_ll] < 0.5;
        }
        if(do_super) //should we use the supergrid for this neighbor
        {
          ForEach3CascNeighbor(_dt,_dx*3.0,flowstates,r,c,l, (r2/3), (c2/3), l2-1,dx,dy, include_edge, (lambda),ca, ca2);
        }
        //is the neighbor not subgridded, or are we in the last grid
        else if(!(flowstates[l2][2]->data[r2][c2] > 0.5)  || !(l2 < flowstates.size()-1))
        {
          //call lambda
          lambda(r2,c2,l2,_dx,ca, ca2);
        }else{ //if the neighbor is further subgridded
          //call this same function three times, for each neighboring sub-pixel at the boundary
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy+dx, c2*3+1-dx+dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2);
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy   , c2*3+1-dx,    l2+1,dx,dy,include_edge, (lambda),ca, ca2);
          ForEach3CascNeighbor(_dt,_dx/3.0,flowstates,r,c,l, r2*3+1-dy-dx, c2*3+1-dx-dy, l2+1,dx,dy,include_edge, (lambda),ca, ca2);
        }
    }
   }else if(include_edge)
   {
      //get z's, h's, v's
      lambda(r2,c2,l2,_dx,ca, ca2);
   }

}

template<typename F>
void ForEachSubCell(int r, int c, int l, std::vector<std::vector<cTMapf *>> &flowstates,float dx, const F & lambda)
{
    if(!(r > -1 && c > -1 && r < flowstates[l][0]->nrRows() && c < flowstates[l][0]->nrCols()))
    {
        return;
    }
   if(flowstates[l][2]->data[r][c] > 0.5f && l < flowstates.size()-1)
   {
      for(int r2 = 0; r2 < 3; r2++)
      {
          for (int c2 = 0; c2 < 3; c2++)
          {
              ForEachSubCell(r * 3 + r2,c*3 + c2, l+1,flowstates,dx/3.0,lambda);
          }
      }

      return;

   }else{

       lambda(r,c,l,flowstates,dx);
   }



}

#endif // FULLFLOOD_3CASCHELPER_H
