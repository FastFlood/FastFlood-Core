
inline static std::vector<cTMapf*> AS_SlopeStabilityISS(cTMapf * DEM, cTMapf * soildepth, cTMapf * coh, cTMapf * ifa, cTMapf * dens, cTMapf * WH, cTMapf * PGA, float sfmax, bool do_failure, float coh_mult, float ifa_mult)
{

    cTMapf * SF = DEM->GetCopy0();
    cTMapf * FD = DEM->GetCopy0();

    float dx = DEM->cellSize();
    #pragma omp parallel for collapse(2)
    for(int r = 0; r < DEM->data.nr_rows();r++)
    {
        for(int c = 0; c < DEM->data.nr_cols();c++)
        {
            if(!pcr::isMV(DEM->data[r][c]))
            {
                cTMapf * Other = DEM;

                float x_11 = GetMapValue_OUTORMV3x3Av(Other,r-1,c-1);
                float x_12 = GetMapValue_OUTORMV3x3Av(Other,r-1,c);
                float x_13 = GetMapValue_OUTORMV3x3Av(Other,r-1,c+1);
                float x_21 = GetMapValue_OUTORMV3x3Av(Other,r,c-1);
                float x_22 = GetMapValue_OUTORMV3x3Av(Other,r,c);
                float x_23 = GetMapValue_OUTORMV3x3Av(Other,r,c+1);
                float x_31 = GetMapValue_OUTORMV3x3Av(Other,r+1,c-1);
                float x_32 = GetMapValue_OUTORMV3x3Av(Other,r+1,c);
                float x_33 = GetMapValue_OUTORMV3x3Av(Other,r+1,c+1);

                float dzdx = (x_13 + 2.0f * x_23 + x_33 - x_11 - 2.0f* x_21 - x_31)/(8.0f*dx);
                float dzdy = (x_31 + 2.0f * x_32 + x_33 - x_11 - 2.0f * x_12 - x_13)/(8.0f*dx);

                float pga = PGA->data[r][c];
                float slope = std::sqrt(dzdx*dzdx + dzdy*dzdy);
                float sins = sin(atan(slope));
                float coss = cos(atan(slope));
                float tanphi = ifa_mult* tan(ifa->data[r][c]);
                float den = coh_mult*coh->data[r][c] -coss * sins*(pga/9.81)*
                        (
                            dens->data[r][c] *(soildepth->data[r][c] - WH->data[r][c]) + 1000.0f *(WH->data[r][c])
                        ) + tanphi*(soildepth->data[r][c] * (dens->data[r][c]) - WH->data[r][c] * 1000.0f) * coss*coss;
                float div =  coss*coss*(pga/9.81)*(dens->data[r][c] * soildepth->data[r][c] + 1000.0 *(WH->data[r][c])) + (soildepth->data[r][c] * (dens->data[r][c]) + WH->data[r][c] * 1000.0f) * sins*coss;

                FD->data[r][c] = 0.0;
                if(div > 1e-12)
                {
                    SF->data[r][c] = den/div;

                    if(SF->data[r][c] < 1.0)
                    {

                        float fheight =  0.0f;

                        float swh = WH->data[r][c];
                        float sd = soildepth->data[r][c];

                        //get all relevant vvariables for calculation of stable depth
                        float wf = swh/sd;
                        float swf = swh/sd;
                        float wd = 1000.0f;
                        float hf = std::max(0.5f * sd,sd - slope * dx);
                        float cohf = coh_mult*coh->data[r][c];

                        float h0 = hf;
                        float hdx = 0.5f * dx;
                        float sde = dens->data[r][c];

                        float sc = coh_mult* coh->data[r][c] + (pga/9.81) *(wd * wf + (1.0 - wf) *dens->data[r][c])*coss*sins *tanphi -  (pga/9.81) *(wd * wf + (1.0 - wf) *dens->data[r][c])*coss*coss;

                        //inversion of safety factor equation (solved using trigonometric identies of squares in a sin/squares in a cosine)
                        //reduces to quadratic equation, of which two solutions can be found

                            double h1 = soildepth->data[r][c];
                            double h2 = soildepth->data[r][c];

                            //if(!(SwitchSeismic || SwitchUpslopeForcing || SwitchDownslopeForcing))
                            {
                                //get all relevant vvariables for calculation of stable depth
                                double cif = cos(ifa->data[r][c]);
                                double sif = sin(ifa->data[r][c]);
                                double d = dx;
                                double dx2 = dx*dx;
                                double pp = 0.01f;
                                double pc = 0.01;;
                                double sd = dens->data[r][c];
                                double ws = 0.01f;

                                //threshold safety factor value
                                //include sf calibration (compensation factor to make sure that the initial state is stable)
                                double sf = 1.0;// / _SFCalibration->Drc;

                                //solution for stable depth at which safety factor equals a threshold value
                                double t1 = (2.0 *hf *pc *cif + d* pp* sf*cif + 2.0 *hf *sc *cif -
                                             d *hf *sd *sf * cif + 2.0 *hf *ws *cif -
                                             dx2 *sd *sif + dx2 *wd *wf *sif);
                                double t21 = (-2.0 *hf *pc *cif - d*pp*sf*cif - 2.0 *hf *sc *cif +
                                        d *hf *sd *sf*cif - 2.0 *hf *ws *cif + dx2 *sd *sif -
                                        dx2 *wd *wf *sif);
                                double t22 = (pc*cif + sc*cif - d*sf*sd*cif +
                                        ws*cif)*(dx2*pc*cif + hf*hf*pc*cif + d*hf*pp*sf*cif +
                                        dx2*sc*cif + hf*hf*sc*cif + dx2*ws*cif +
                                        hf*hf*ws*cif + dx2*pp*sif);

                                //in the end, it is nothin more than a really complex implementation of the abc-formula
                                // x1,2 = (-b +- Sqrt(b^2-4ac))/2a
                                //Two solutions for this formula, in our case: 1 positive, 1 negative, so we pick the positive value

                                double t2 = sqrt(t21*t21 - 4.0*t22);
                                double t3 = (2.0*(pc*cif + sc*cif - d*sd*sf*cif + ws*cif));

                                //positive and negative solution
                                h1 = (t1 + t2)/t3;
                                h2 = (t1 - t2)/t3;

                                h1 = std::isnan(h1)? 0.0:h1;
                                h2 = std::isnan(h2)? 0.0:h2;

                                //qDebug() << h1 << h2 << _DFForcing->Drc << _DFForcingUp->Drc;
                            }


                        float hnew_select1 = std::max(h1,h2);
                        float hnew_select2 = std::min(h1,h2);

                        float hnew_select = hnew_select2;

                        float hnew = std::min(sd,std::max(sd * 0.5f,hnew_select));


                        fheight = std::min(sd,std::max(0.0f,sd - hnew));

                        FD->data[r][c] = fheight;
                    }

                }else {
                    pcr::setMV(SF->data[r][c]);
                }
            }


        }
    }

    return {SF,FD};
}




inline static std::vector<cTMapf*> AS_SlopeStabilityIFM(cTMapf * DEMORG, cTMapf * soildepthorg, cTMapf * coh, cTMapf * ifa, cTMapf * dens, cTMapf * WH, float sfmax,int iter, bool do_failure, float coh_mult, float ifa_mult)
{

    cTMapf * soildepth = soildepthorg->GetCopy();
    cTMapf * DEM = DEMORG->GetCopy();
    cTMapf * SF = DEM->GetCopy0();
    cTMapf * FD = DEM->GetCopy0();
    cTMapf * FDT = DEM->GetCopy0();

    for(int i = 0; i < iter; i++)
    {
        float dx = DEM->cellSize();
        #pragma omp parallel for collapse(2)
        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                if(!pcr::isMV(DEM->data[r][c]))
                {
                    cTMapf * Other = DEM;

                    float x_11 = GetMapValue_OUTORMV3x3Av(Other,r-1,c-1);
                    float x_12 = GetMapValue_OUTORMV3x3Av(Other,r-1,c);
                    float x_13 = GetMapValue_OUTORMV3x3Av(Other,r-1,c+1);
                    float x_21 = GetMapValue_OUTORMV3x3Av(Other,r,c-1);
                    float x_22 = GetMapValue_OUTORMV3x3Av(Other,r,c);
                    float x_23 = GetMapValue_OUTORMV3x3Av(Other,r,c+1);
                    float x_31 = GetMapValue_OUTORMV3x3Av(Other,r+1,c-1);
                    float x_32 = GetMapValue_OUTORMV3x3Av(Other,r+1,c);
                    float x_33 = GetMapValue_OUTORMV3x3Av(Other,r+1,c+1);

                    float dzdx = (x_13 + 2.0f * x_23 + x_33 - x_11 - 2.0f* x_21 - x_31)/(8.0f*dx);
                    float dzdy = (x_31 + 2.0f * x_32 + x_33 - x_11 - 2.0f * x_12 - x_13)/(8.0f*dx);

                    float slope = std::sqrt(dzdx*dzdx + dzdy*dzdy);
                    float sins = sin(atan(slope));
                    float coss = cos(atan(slope));
                    float tanphi = ifa_mult* tan(ifa->data[r][c]);
                    float den = coh_mult*coh->data[r][c] + tanphi*(soildepth->data[r][c] * (dens->data[r][c]) - WH->data[r][c] * 1000.0f) * coss*coss;
                    float div = (soildepth->data[r][c] * (dens->data[r][c]) + WH->data[r][c] * 1000.0f) * sins*coss;

                    FD->data[r][c] = 0.0;
                    if(div > 1e-12)
                    {
                        if(i == 0)
                        {
                            SF->data[r][c] = den/div;
                        }else
                        {

                            SF->data[r][c] = std::min(SF->data[r][c],den/div);
                        }

                        if(den/div < 1.0)
                        {

                            float fheight =  0.0f;

                            float swh = WH->data[r][c];
                            float sd = soildepth->data[r][c];

                            //get all relevant vvariables for calculation of stable depth
                            float wf = swh/sd;
                            float swf = swh/sd;
                            float wd = 1000.0f;
                            float hf = std::max(0.5f * sd,sd - slope * dx);
                            float cohf = coh_mult *coh->data[r][c];

                            float h0 = hf;
                            float hdx = 0.5f * dx;
                            float sde = dens->data[r][c];
                            float h1 = sd;
                            float h2 = sd;
                            float sc = 0.0f;

                            //inversion of safety factor equation (solved using trigonometric identies of squares in a sin/squares in a cosine)
                            //reduces to quadratic equation, of which two solutions can be found
                            {


                                    h1 = ((-cohf)*h0*sde + h0*sde*sfmax + dx*sde*tanphi -
                                          cohf*h0*swf*wd + h0*sfmax*swf*wd + dx*swf*tanphi*wd -

                                          sqrt(-4.0f*dx*sde*swh*
                                             tanphi*((-cohf)*sde + sde*sfmax - cohf*swf*wd + sfmax*swf*wd) +
                                                    pow(cohf*h0*sde - h0*sde*sfmax - dx*sde*tanphi +
                                               cohf*h0*swf*wd - h0*sfmax*swf*wd - dx*swf*tanphi*wd,2.0f)))/
                                             (2.0f*((-cohf)*sde + sde*sfmax - cohf*swf*wd + sfmax*swf*wd));
                                    h2 = ((-cohf)*h0*sde + h0*sde*sfmax + dx*sde*tanphi -
                                          cohf*h0*swf*wd + h0*sfmax*swf*wd + dx*swf*tanphi*wd +

                                          sqrt(-4.0f*dx*sde*swh*
                                             tanphi*((-cohf)*sde + sde*sfmax - cohf*swf*wd + sfmax*swf*wd) +
                                                    pow(cohf*h0*sde - h0*sde*sfmax - dx*sde*tanphi +
                                               cohf*h0*swf*wd - h0*sfmax*swf*wd - dx*swf*tanphi*wd,2.0f)))/
                                             (2.0f*((-cohf)*sde + sde*sfmax - cohf*swf*wd + sfmax*swf*wd));

                                    h1 = std::isnan(h1)? sd:h1;
                                    h2 = std::isnan(h2)? sd:h2;

                            }

                            float hnew_select1 = std::max(h1,h2);
                            float hnew_select2 = std::min(h1,h2);

                            float hnew_select = hnew_select1;

                            float hnew = std::min(sd,std::max(sd * 0.5f,hnew_select));


                            fheight = std::min(sd,std::max(0.0f,sd - hnew));

                            FD->data[r][c] = fheight;
                        }else
                        {

                            FD->data[r][c] = 0.0;
                        }

                    }else {
                        FD->data[r][c] = 0.0;
                    }
                }


            }
        }

        bool has_failure= false;

        //#pragma omp parallel for collapse(2)
        for(int r = 0; r < DEM->data.nr_rows();r++)
        {
            for(int c = 0; c < DEM->data.nr_cols();c++)
            {
                if(!pcr::isMV(DEM->data[r][c]))
                {
                    if(FD->data[r][c] > 0.05)
                    {
                        has_failure = true;

                        DEM->data[r][c] -= FD->data[r][c];
                        soildepth->data[r][c] -= FD->data[r][c];
                        FDT->data[r][c] += FD->data[r][c];
                        
                    }
                }
            }
        }
        if(!has_failure)
        {
            break;
        }
    }

    delete soildepth;
    delete DEM;
    delete FD;

    return {SF,FDT};

}