#ifndef FULLFLOOD_GRIDRESAMPLE_H
#define FULLFLOOD_GRIDRESAMPLE_H

#include "../algorithms/rasteralgorithms/gridresample.h"

inline void CreateMap3CascadeAvMax(cTMapf * DEM, cTMapf * MAN, cTMapf * RAIN, cTMapf * KSAT, cTMapf * R3DEM, cTMapf * R3DEMM, cTMapf * R3MAN, cTMapf * R3RAIN, cTMapf * R3KSAT)
{
  //FF::cout <<"3-av cascade calc: " << DEM->Rows() << " " << DEM->Cols() << "  " << R3DEM->Rows() << " " << R3DEM->Cols() << std::endl;
  for (int r = 0; r < R3DEM->Rows(); r++) {
        for (int c = 0; c < R3DEM->Cols(); c++) {

          //FF::cout << "rc " << r << "  "  << c<<std::endl;
        int r_ip = r *3 + 1;
        int c_ip = c *3 + 1;

        float z = 0.0;
        float man = 0.06;
        float rain = 0.00;
        float ksat = 0.00;
        if(r_ip > -1 && c_ip > -1 && r_ip < DEM->nrRows() && c_ip < DEM->nrCols())
        {
            z = DEM->data[r_ip][c_ip];
            man = MAN->data[r_ip][c_ip];
            rain = RAIN->data[r_ip][c_ip];
            ksat = KSAT->data[r_ip][c_ip];
        }

        if(pcr::isMV(z))
        {
            pcr::setMV(R3DEM->data[r][c]);

            R3MAN->data[r][c] = 0.06;

            int n_znb = 0;
            float z_total = 0;
            //get backup z
            if(r_ip-1 > -1 && c_ip > -1 && r_ip-1 < DEM->nrRows() && c_ip < DEM->nrCols())
            {
                float zy1 = DEM->data[r_ip-1][c_ip];
                if(!pcr::isMV(zy1))
                {
                    z_total += zy1;
                    n_znb += 1;
                }
            }
            if(r_ip+1 > -1 && c_ip > -1 && r_ip+1 < DEM->nrRows() && c_ip < DEM->nrCols())
            {
                float zy2 = DEM->data[r_ip+1][c_ip];
                if(!pcr::isMV(zy2))
                {
                    z_total += zy2;
                    n_znb += 1;
                }
            }
            if(r_ip > -1 && c_ip-1 > -1 && r_ip < DEM->nrRows() && c_ip-1 < DEM->nrCols())
            {
                float zx1 = DEM->data[r_ip][c_ip-1];
                if(!pcr::isMV(zx1))
                {
                    z_total += zx1;
                    n_znb += 1;
                }
            }
            if(r_ip > -1 && c_ip+1 > -1 && r_ip < DEM->nrRows() && c_ip+1 < DEM->nrCols())
            {
                float zx2 = DEM->data[r_ip][c_ip+1];
                if(!pcr::isMV(zx2))
                {
                    z_total += zx2;
                    n_znb += 1;
                }
            }

            if(n_znb > 0)
            {
                z = z_total/std::max(1.0f,((float)(n_znb)));
            }else{
                continue;
            }

        }

        //FF::cout << 1<<std::endl;

          //get 3x3 subgrid values
          float zx1 = z;
          float zx2 = z;
          float zy1 = z;
          float zy2 = z;

          float zx1y1 = z;
          float zx1y2 = z;
          float zx2y1 = z;
          float zx2y2 = z;

          float nx1 = man;
          float nx2 = man;
          float ny1 = man;
          float ny2 = man;

          float nx1y1 = man;
          float nx1y2 = man;
          float nx2y1 = man;
          float nx2y2 = man;

          float rx1 = rain;
          float rx2 = rain;
          float ry1 = rain;
          float ry2 = rain;

          float rx1y1 = rain;
          float rx1y2 = rain;
          float rx2y1 = rain;
          float rx2y2 = rain;

          float kx1 = ksat;
          float kx2 = ksat;
          float ky1 = ksat;
          float ky2 = ksat;

          float kx1y1 = ksat;
          float kx1y2 = ksat;
          float kx2y1 = ksat;
          float kx2y2 = ksat;

          if(r_ip-1 > -1 && c_ip > -1 && r_ip-1 < DEM->nrRows() && c_ip < DEM->nrCols())
            {
                zy1 = DEM->data[r_ip-1][c_ip];
                ny1 = MAN->data[r_ip-1][c_ip];
                ky1 = KSAT->data[r_ip-1][c_ip];
                ry1 = RAIN->data[r_ip-1][c_ip];
                if(pcr::isMV(zy1))
                {
                    zy1 = z;
                    ny1 = man;
                    ry1 = rain;
                }
            }
            if(r_ip+1 > -1 && c_ip > -1 && r_ip+1 < DEM->nrRows() && c_ip < DEM->nrCols())
            {
                zy2 = DEM->data[r_ip+1][c_ip];
                ny2 = MAN->data[r_ip+1][c_ip];
                ry2 = RAIN->data[r_ip+1][c_ip];
                ky2 = KSAT->data[r_ip+1][c_ip];
                if(pcr::isMV(zy2))
                {
                    zy2 = z;
                    ny2 = man;
                    ry2 = rain;
                }
            }
            if(r_ip > -1 && c_ip-1 > -1 && r_ip < DEM->nrRows() && c_ip-1 < DEM->nrCols())
            {
                zx1 = DEM->data[r_ip][c_ip-1];
                nx1 = MAN->data[r_ip][c_ip-1];
                rx1 = RAIN->data[r_ip][c_ip-1];
                kx1 = KSAT->data[r_ip][c_ip-1];
                if(pcr::isMV(zx1))
                {
                    zx1 = z;
                    nx1 = man;
                    rx1 = rain;
                }
            }
            if(r_ip > -1 && c_ip+1 > -1 && r_ip < DEM->nrRows() && c_ip+1 < DEM->nrCols())
            {
                zx2 = DEM->data[r_ip][c_ip+1];
                nx2 = MAN->data[r_ip][c_ip+1];
                rx2 = RAIN->data[r_ip][c_ip+1];
                kx2 = KSAT->data[r_ip][c_ip+1];
                if(pcr::isMV(zx2))
                {
                    zx2 = z;
                    nx2 = man;
                    rx2 = rain;
                }
            }

        //FF::cout << 2<<std::endl;
            if(r_ip-1 > 0 && c_ip-1 > 0 && r_ip-1 < DEM->nrRows() && c_ip-1 < DEM->nrCols())
            {
                zx1y1 = DEM->data[r_ip-1][c_ip-1];
                nx1y1 = MAN->data[r_ip-1][c_ip-1];
                rx1y1 = RAIN->data[r_ip-1][c_ip-1];
                kx1y1 = KSAT->data[r_ip-1][c_ip-1];
                if(pcr::isMV(zx1y1))
                {
                    zx1y1 = z;
                    nx1y1 = man;
                    rx1y1 = rain;
                }
            }
            if(r_ip+1 > 0 && c_ip-1 > 0 && r_ip+1 < DEM->nrRows() && c_ip-1 < DEM->nrCols())
            {

                zx1y2 = DEM->data[r_ip+1][c_ip-1];
                nx1y2 = MAN->data[r_ip+1][c_ip-1];
                rx1y2 = RAIN->data[r_ip+1][c_ip-1];
                kx1y2 = KSAT->data[r_ip+1][c_ip-1];
                if(pcr::isMV(zx1y2))
                {
                    zx1y2 = z;
                    nx1y2 = man;
                    rx1y2 = rain;
                }
            }
            if(r_ip-1 > 0 && c_ip+1 > 0 && r_ip-1 < DEM->nrRows() && c_ip+1 < DEM->nrCols())
            {
                zx2y1 = DEM->data[r_ip-1][c_ip+1];
                nx2y1 = MAN->data[r_ip-1][c_ip+1];
                rx2y1 = RAIN->data[r_ip-1][c_ip+1];
                kx2y1 = KSAT->data[r_ip-1][c_ip+1];
                if(pcr::isMV(zx2y1))
                {
                    zx2y1 = z;
                    nx2y1 = man;
                    rx2y1 = rain;
                }
            }
            if(r_ip+1 > 0 && c_ip+1 > 0 && r_ip+1 < DEM->nrRows() && c_ip+1 < DEM->nrCols())
            {
                zx2y2 = DEM->data[r_ip+1][c_ip+1];
                nx2y2 = MAN->data[r_ip+1][c_ip+1];
                rx2y2 = RAIN->data[r_ip+1][c_ip+1];
                kx2y2 = KSAT->data[r_ip+1][c_ip+1];
                if(pcr::isMV(zx2y2))
                {
                    zx2y2 = z;
                    nx2y2 = man;
                    rx2y2 = rain;
                }
            }

        //FF::cout << 3<<std::endl;
            R3DEM->data[r][c] = (zx1 + zx2 + zy1 + zy2 + zx1y1 + zx1y2 + zx2y1 + zx2y2 + z)/9.0f;
            R3MAN->data[r][c] = (nx1 + nx2 + ny1 + ny2 + nx1y1 + nx1y2 + nx2y1 + nx2y2 + man)/9.0f;
            R3RAIN->data[r][c] = (rx1 + rx2 + ry1 + ry2 + rx1y1 + rx1y2 + rx2y1 + rx2y2 + rain)/9.0f;
            R3KSAT->data[r][c] = (kx1 + kx2 + ky1 + ky2 + kx1y1 + kx1y2 + kx2y1 + kx2y2 + ksat)/9.0f;
            R3DEMM->data[r][c] = std::max(std::max(std::max(std::max(std::max(std::max(std::max(std::max(z,zx1),zx2),zy1),zy2),zx1y1),zx1y2),zx2y1),zx2y2);
        }
   }

}

#endif // FULLFLOOD_GRIDRESAMPLE_H
