

inline std::vector<cTMapf *> DinfiniteLandslideMomentumDepositionSolve(std::vector<cTMapf*> dinf,cTMapf * DEM, cTMapf * manning, cTMapf * source_solid, cTMapf * source_fluid, cTMapf * density, cTMapf * acc_org, cTMapf * width_org, cTMapf * ifamap, cTMapf * dist_variance,cTMapf * DEPRESSIONID, std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, int iter_max = 0)
{

    //for each pixel
    //check if material
    //get flow width/height combo
    //first, check if accaleration is positive
    //get local peak volume
    //calculate duration of travel and momentum and mass going to next cell
    //second, if negative accaleration
    //reduce momentum going to next cell
    //do we have no momentum left?
    //calculate the amount of material needed as deposition to have the depositions reach the mid-point of the next cell
    //solve for slope to get a zero accaleration? or solve for slope to get accelaration so that dx *(1+0.5) is still travelled
    //deposit this material and propagate reduced material and momentum

    cTMapf * map = DEM->GetCopy0();
    cTMapf * materialf = source_solid->GetCopy0();
    cTMapf * materials = source_fluid->GetCopy0();
    cTMapf * deposition = source_solid->GetCopy0();
    cTMapf * momentum = DEM->GetCopy0();
    cTMapf * velavg = DEM->GetCopy0();
    cTMapf * velacc = DEM->GetCopy0();

    cTMapf * fracx1 = dinf[0];
    cTMapf * fracx2 = dinf[1];
    cTMapf * fracy1 = dinf[2];
    cTMapf * fracy2 = dinf[3];

    float in_dx = DEM->cellSize();

 
    float dx = map->cellSize();
    cTMapf * slope = DEM->GetCopy0();
    //get average distance, N and velocity
    #pragma omp parallel for collapse(2)
    for (int r = 0; r < DEM->Rows(); r++) {
        for (int c = 0; c < DEM->Cols(); c++) {
            if (!pcr::isMV(DEM->data[r][c])) {
                float SlopeX = UF2D_Derivative(DEM,r,c,UF_DIRECTION_X,1.0,UF_DERIVATIVE_LR);
                float SlopeY = UF2D_Derivative(DEM,r,c,UF_DIRECTION_Y,1.0,UF_DERIVATIVE_LR);
                float Slope = std::fabs(SlopeX) + std::fabs(SlopeY);
                slope->data[r][c] = std::max(1e-4f,Slope);
            }
        }
    }



#define NUM_ITERS_LANDSLIDE_MOMENTUM 10

    //we keep iterating through this algorithm untill there is no change left to make
    bool change = true;
    bool first = true;

    float iters_inv = 1.0/((float)(NUM_ITERS_LANDSLIDE_MOMENTUM));

    int iter = 0;

    while(change && ((iter_max <= 0) || (iter_max > 0 && iter < iter_max)))
    {
        iter ++;
        //change = false;

        for(int k = 0; k < DEPRESSIONCVOLS.size(); k++)
        {
            DEPRESSIONCVOLS[k] = 0.0;
        }

        if(iter%2 == 0)
        {

            //first we move in right-lower direction
            for(int r = 0; r < map->data.nr_rows();r++)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = DEM->data[r][c];
                    if(!pcr::isMV(v_points))
                    {


                        double w_x1 = 0.0;
                        double w_x2 = 0.0;
                        double w_y1 = 0.0;
                        double w_y2 = 0.0;

                        double wf_x1 = 0.0;
                        double wf_x2 = 0.0;
                        double wf_y1 = 0.0;
                        double wf_y2 = 0.0;

                        double ws_x1 = 0.0;
                        double ws_x2 = 0.0;
                        double ws_y1 = 0.0;
                        double ws_y2 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                float h = momentum->data[r][c-1];
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * h;
                                float hf = materialf->data[r][c-1];
                                wf_x1 = std::fabs(fracx2->data[r][c-1]) * hf;
                                float hs = materials->data[r][c-1];
                                ws_x1 = std::fabs(fracx2->data[r][c-1]) * hs;
                            }
                        }
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10) && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                float h = momentum->data[r-1][c];
                                w_y1 = std::fabs(fracy2->data[r-1][c]) *  h;
                                float hf = materialf->data[r-1][c];
                                wf_y1 = std::fabs(fracy2->data[r-1][c]) * hf;
                                float hs = materials->data[r-1][c];
                                ws_y1 = std::fabs(fracy2->data[r-1][c]) * hs;
                            }
                        }
                        if(c < map->nrCols()-1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                float h = momentum->data[r][c+1];
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * h;
                                float hf = materialf->data[r][c+1];
                                wf_x2 = std::fabs(fracx1->data[r][c+1]) * hf;
                                float hs = materials->data[r][c+1];
                                ws_x2 = std::fabs(fracx1->data[r][c+1]) * hs;
                            }
                        }
                        if(r < map->nrRows()-1)
                        {
                            if((std::fabs(fracy1->data[r+1][c]) > 1e-10) && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                float h = momentum->data[r+1][c];
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * h;
                                float hf = materialf->data[r+1][c];
                                wf_y2 = std::fabs(fracy1->data[r+1][c]) * hf;
                                float hs = materials->data[r+1][c];
                                ws_y2 = std::fabs(fracy1->data[r+1][c]) * hs;
                            }
                        }



                        //get local properties
                        float momentum_in = w_x1 + w_x2 + w_y1 + w_y2;
                        float ifa = ifamap->data[r][c];
                        float s = slope->data[r][c];
                        float water = wf_x1 + wf_x2 + wf_y1 + wf_y2 + source_fluid->data[r][c] * in_dx* in_dx;
                        float solid = ws_x1 + ws_x2 + ws_y1 + ws_y2 + source_solid->data[r][c] * in_dx* in_dx;
                        float n = manning->data[r][c];
                        float dens = density->data[r][c];

                        //check potential inflow into depression 
                        if(DEPRESSIONID != nullptr)
                        {
                            //chcek if there is a local depression here
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    
                                    float totalflow = water + solid;
                                    float frac_depr = std::max(0.0f,std::min(1.0f,vol/std::max(0.00001f,totalflow)));
                                    water = water * (1.0-frac_depr);
                                    solid = solid * (1.0-frac_depr);

                                    //FF::cout << "depr: " << DEPRESSIONID->data[r][c] << " " << r << " " << c << " " << frac_depr << " " << vol << " " << totalflow << std::endl;
                                    
                                    //otherwise do deposition
                                    DEPRESSIONCVOLS.at(id) += frac_depr * totalflow;

                                    //fractionally remove momentum
                                    momentum_in = momentum_in * (1.0-frac_depr);
                                }
                            }
                        }

                        float material_total = water + solid;

                        //estimate height
                        float blocklength = dist_variance->data[r][c];

                        float flowheight = material_total/std::max(blocklength*dx,0.001f);
                        //get flow width and correct flow height
                        float ratio_vol = material_total/(std::max(0.001f,acc_org->data[r][c])*dx*dx);
                        float flowwidth = ratio_vol * width_org->data[r][c];
                        flowheight = flowheight * dx/std::max(dx,flowwidth);

                        //get solid volumetric content ratio
                        float sratio = std::min(1.0f,solid / std::max(0.001f,solid+water));

                        float weight = (dens * sratio + 1000.0 * (1.0-sratio)) * flowheight;
                        float basal_porepressure = flowheight * 1000.0 * (1.0-sratio);

                        //get coefficient for the differential equation
                        float acc = 9.81 * (s - std::tan(ifa)*std::max(0.0f,(1.0f - basal_porepressure/weight)));

                        float vel_in = momentum_in/(std::max(0.00001f,water*1000.0f+solid*dens));
                        float a = acc;
                        float b = std::max(9.81f*n*n/std::max(0.01f,flowheight),0.01f);

                        //are we increasing or decreasing momentum?

                        float deposit = 0.0;
                        float vel_final = vel_in;

                        //if accalerating (at least, relative to rest)
                        if(acc > 0.0)
                        {
                            //get balance velocity
                            float vel_bal = std::max(0.01f,std::sqrt(a/b));

                            bool slow = vel_bal < vel_in;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = std::min(vel_in,vel_bal);
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }else
                                {
                                    v1 = std::max(std::min(vel_bal,v1),vel_in);
                                }

                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;

                        //if decellaratingto a level that it might stop moving
                        }else
                        {
                            float vel_bal = 0.0;

                            bool slow = true;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = 0.0;
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            float v00 = vel_in;

                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                                        float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }
                                //we have stalled, check if we can continue by depositing partially, otherwise just continue zero momentum and mass
                                if(v1 < 1e-6)
                                {
                                    //we want to deposit as much as we need to maintain momentum, meaning our slope angle (due to deposition) should change to be equal to the basal friction angle
                                    //(std::tan(ifa) - s)*dx * dx is what we deposit for this
                                    deposit  = std::max(0.0f,(std::tan(ifa) - s)*dx*dx*(dx+width_org->data[r][c]));
                                    v0 = v00;
                                    break;
                                }


                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;
                        }


                        //deal with deposits
                        if(deposit > 0.0)
                        {
                            float vfrac_dep = std::min(1.0f,std::max(0.0f,deposit/(std::max(0.001f,water+solid))));
                            water = water * (1.0-vfrac_dep);
                            solid = solid * (1.0-vfrac_dep);
                        }

                         velavg->data[r][c] = flowwidth;;
                        velacc->data[r][c] = a;
                        //velavg->data[r][c] = vel_in;
                        //velacc->data[r][c] = a;
                        //deal with calculating outgoing mass and momentum
                        deposition->data[r][c] = deposit;
                        momentum->data[r][c] = vel_final * (1000.0 * water  + dens * solid);
                        materialf->data[r][c] = water;
                        materials->data[r][c] = solid;


                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = DEM->data[r][c];
                    if(!pcr::isMV(v_points))
                    {

                        double w_x1 = 0.0;
                        double w_x2 = 0.0;
                        double w_y1 = 0.0;
                        double w_y2 = 0.0;

                        double wf_x1 = 0.0;
                        double wf_x2 = 0.0;
                        double wf_y1 = 0.0;
                        double wf_y2 = 0.0;

                        double ws_x1 = 0.0;
                        double ws_x2 = 0.0;
                        double ws_y1 = 0.0;
                        double ws_y2 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                float h = momentum->data[r][c-1];
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * h;
                                float hf = materialf->data[r][c-1];
                                wf_x1 = std::fabs(fracx2->data[r][c-1]) * hf;
                                float hs = materials->data[r][c-1];
                                ws_x1 = std::fabs(fracx2->data[r][c-1]) * hs;
                            }
                        }
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10) && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                float h = momentum->data[r-1][c];
                                w_y1 = std::fabs(fracy2->data[r-1][c]) *  h;
                                float hf = materialf->data[r-1][c];
                                wf_y1 = std::fabs(fracy2->data[r-1][c]) * hf;
                                float hs = materials->data[r-1][c];
                                ws_y1 = std::fabs(fracy2->data[r-1][c]) * hs;
                            }
                        }
                        if(c < map->nrCols()-1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                float h = momentum->data[r][c+1];
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * h;
                                float hf = materialf->data[r][c+1];
                                wf_x2 = std::fabs(fracx1->data[r][c+1]) * hf;
                                float hs = materials->data[r][c+1];
                                ws_x2 = std::fabs(fracx1->data[r][c+1]) * hs;
                            }
                        }
                        if(r < map->nrRows()-1)
                        {
                            if((std::fabs(fracy1->data[r+1][c]) > 1e-10) && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                float h = momentum->data[r+1][c];
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * h;
                                float hf = materialf->data[r+1][c];
                                wf_y2 = std::fabs(fracy1->data[r+1][c]) * hf;
                                float hs = materials->data[r+1][c];
                                ws_y2 = std::fabs(fracy1->data[r+1][c]) * hs;
                            }
                        }



                        //get local properties
                        float momentum_in = w_x1 + w_x2 + w_y1 + w_y2;
                        float ifa = ifamap->data[r][c];
                        float s = slope->data[r][c];
                        float water = wf_x1 + wf_x2 + wf_y1 + wf_y2 + source_fluid->data[r][c] * in_dx* in_dx;
                        float solid = ws_x1 + ws_x2 + ws_y1 + ws_y2 + source_solid->data[r][c] * in_dx* in_dx;
                        float n = manning->data[r][c];
                        float dens = density->data[r][c];


//check potential inflow into depression 
                        if(DEPRESSIONID != nullptr)
                        {
                            //chcek if there is a local depression here
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    
                                    float totalflow = water + solid;
                                    float frac_depr = std::max(0.0f,std::min(1.0f,vol/std::max(0.00001f,totalflow)));
                                    water = water * (1.0-frac_depr);
                                    solid = solid * (1.0-frac_depr);

                                    //FF::cout << "depr: " << r << " " << c << " " << frac_depr << " " << vol << " " << totalflow << std::endl;
                                    //otherwise do deposition
                                    DEPRESSIONCVOLS.at(id) += frac_depr * totalflow;

                                    //fractionally remove momentum
                                    momentum_in = momentum_in * (1.0-frac_depr);
                                }
                            }
                        }

                        float material_total = water + solid;

                        //estimate height
                        float blocklength = dist_variance->data[r][c];

                        float flowheight = material_total/std::max(blocklength*dx,0.001f);
                        //get flow width and correct flow height
                        float ratio_vol = material_total/(std::max(0.001f,acc_org->data[r][c])*dx*dx);
                        float flowwidth = ratio_vol * width_org->data[r][c];
                        flowheight = flowheight * dx/std::max(dx,flowwidth);

                        //get solid volumetric content ratio
                        float sratio = std::min(1.0f,solid / std::max(0.001f,solid+water));

                        float weight = (dens * sratio + 1000.0 * (1.0-sratio)) * flowheight;
                        float basal_porepressure = flowheight * 1000.0 * (1.0-sratio);

                        //get coefficient for the differential equation
                        float acc = 9.81 * (s - std::tan(ifa)*std::max(0.0f,(1.0f - basal_porepressure/weight)));

                        float vel_in = momentum_in/(std::max(0.00001f,water*1000.0f+solid*dens));
                        float a = acc;
                        float b = std::max(9.81f*n*n/std::max(0.01f,flowheight),0.01f);

                        //are we increasing or decreasing momentum?

                        float deposit = 0.0;
                        float vel_final = vel_in;

                        //if accalerating (at least, relative to rest)
                        if(acc > 0.0)
                        {
                            //get balance velocity
                            float vel_bal = std::max(0.01f,std::sqrt(a/b));

                            bool slow = vel_bal < vel_in;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = std::min(vel_in,vel_bal);
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }else
                                {
                                    v1 = std::max(std::min(vel_bal,v1),vel_in);
                                }

                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;

                        //if decellaratingto a level that it might stop moving
                        }else
                        {
                            float vel_bal = 0.0;

                            bool slow = true;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = 0.0;
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            float v00 = vel_in;

                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                                        float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }
                                //we have stalled, check if we can continue by depositing partially, otherwise just continue zero momentum and mass
                                if(v1 < 1e-6)
                                {
                                    //we want to deposit as much as we need to maintain momentum, meaning our slope angle (due to deposition) should change to be equal to the basal friction angle
                                    //(std::tan(ifa) - s)*dx * dx is what we deposit for this
                                    deposit  = std::max(0.0f,(std::tan(ifa) - s)*dx*dx*(dx+width_org->data[r][c]));
                                    v0 = v00;
                                    break;
                                }


                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;
                        }


                        //deal with deposits
                        if(deposit > 0.0)
                        {
                            float vfrac_dep = std::min(1.0f,std::max(0.0f,deposit/(std::max(0.001f,water+solid))));
                            water = water * (1.0-vfrac_dep);
                            solid = solid * (1.0-vfrac_dep);
                        }

                        velavg->data[r][c] = flowwidth;;
                        velacc->data[r][c] = a;
                        //velavg->data[r][c] = vel_in;
                        //velacc->data[r][c] = a;
                        //deal with calculating outgoing mass and momentum
                        deposition->data[r][c] = deposit;
                        momentum->data[r][c] = vel_final * (1000.0 * water  + dens * solid);
                        materialf->data[r][c] = water;
                        materials->data[r][c] = solid;


                    }


                }
            }
        }else
        {

            //first we move in right-lower direction
            for(int r = 0; r < map->data.nr_rows();r++)
            {
                for(int c = map->data.nr_cols()-1; c > -1 ;c--)
                {
                    float v_points = DEM->data[r][c];
                    if(!pcr::isMV(v_points))
                    {


                        double w_x1 = 0.0;
                        double w_x2 = 0.0;
                        double w_y1 = 0.0;
                        double w_y2 = 0.0;

                        double wf_x1 = 0.0;
                        double wf_x2 = 0.0;
                        double wf_y1 = 0.0;
                        double wf_y2 = 0.0;

                        double ws_x1 = 0.0;
                        double ws_x2 = 0.0;
                        double ws_y1 = 0.0;
                        double ws_y2 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                float h = momentum->data[r][c-1];
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * h;
                                float hf = materialf->data[r][c-1];
                                wf_x1 = std::fabs(fracx2->data[r][c-1]) * hf;
                                float hs = materials->data[r][c-1];
                                ws_x1 = std::fabs(fracx2->data[r][c-1]) * hs;
                            }
                        }
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10) && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                float h = momentum->data[r-1][c];
                                w_y1 = std::fabs(fracy2->data[r-1][c]) *  h;
                                float hf = materialf->data[r-1][c];
                                wf_y1 = std::fabs(fracy2->data[r-1][c]) * hf;
                                float hs = materials->data[r-1][c];
                                ws_y1 = std::fabs(fracy2->data[r-1][c]) * hs;
                            }
                        }
                        if(c < map->nrCols()-1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                float h = momentum->data[r][c+1];
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * h;
                                float hf = materialf->data[r][c+1];
                                wf_x2 = std::fabs(fracx1->data[r][c+1]) * hf;
                                float hs = materials->data[r][c+1];
                                ws_x2 = std::fabs(fracx1->data[r][c+1]) * hs;
                            }
                        }
                        if(r < map->nrRows()-1)
                        {
                            if((std::fabs(fracy1->data[r+1][c]) > 1e-10) && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                float h = momentum->data[r+1][c];
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * h;
                                float hf = materialf->data[r+1][c];
                                wf_y2 = std::fabs(fracy1->data[r+1][c]) * hf;
                                float hs = materials->data[r+1][c];
                                ws_y2 = std::fabs(fracy1->data[r+1][c]) * hs;
                            }
                        }



                        //get local properties
                        float momentum_in = w_x1 + w_x2 + w_y1 + w_y2;
                        float ifa = ifamap->data[r][c];
                        float s = slope->data[r][c];
                        float water = wf_x1 + wf_x2 + wf_y1 + wf_y2 + source_fluid->data[r][c] * in_dx* in_dx;
                        float solid = ws_x1 + ws_x2 + ws_y1 + ws_y2 + source_solid->data[r][c] * in_dx* in_dx;
                        float n = manning->data[r][c];
                        float dens = density->data[r][c];

//check potential inflow into depression 
                        if(DEPRESSIONID != nullptr)
                        {
                            //chcek if there is a local depression here
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    
                                    float totalflow = water + solid;
                                     float frac_depr = std::max(0.0f,std::min(1.0f,vol/std::max(0.00001f,totalflow)));
                                    water = water * (1.0-frac_depr);
                                    solid = solid * (1.0-frac_depr);

                                    //FF::cout << "depr: " << r << " " << c << " " << frac_depr << " " << vol << " " << totalflow << std::endl;
                                    
                                    //otherwise do deposition
                                    DEPRESSIONCVOLS.at(id) += frac_depr * totalflow;

                                    //fractionally remove momentum
                                    momentum_in = momentum_in * (1.0-frac_depr);
                                }
                            }
                        }

                        float material_total = water + solid;

                        //estimate height
                        float blocklength = dist_variance->data[r][c];

                        float flowheight = material_total/std::max(blocklength*dx,0.001f);
                        //get flow width and correct flow height
                        float ratio_vol = material_total/(std::max(0.001f,acc_org->data[r][c])*dx*dx);
                        float flowwidth = ratio_vol * width_org->data[r][c];
                        flowheight = flowheight * dx/std::max(dx,flowwidth);

                        //get solid volumetric content ratio
                        float sratio = std::min(1.0f,solid / std::max(0.001f,solid+water));

                        float weight = (dens * sratio + 1000.0 * (1.0-sratio)) * flowheight;
                        float basal_porepressure = flowheight * 1000.0 * (1.0-sratio);

                        //get coefficient for the differential equation
                        float acc = 9.81 * (s - std::tan(ifa)*std::max(0.0f,(1.0f - basal_porepressure/weight)));

                        float vel_in = momentum_in/(std::max(0.00001f,water*1000.0f+solid*dens));
                        float a = acc;
                        float b = std::max(9.81f*n*n/std::max(0.01f,flowheight),0.01f);

                        //are we increasing or decreasing momentum?

                        float deposit = 0.0;
                        float vel_final = vel_in;

                        //if accalerating (at least, relative to rest)
                        if(acc > 0.0)
                        {
                            //get balance velocity
                            float vel_bal = std::max(0.01f,std::sqrt(a/b));

                            bool slow = vel_bal < vel_in;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = std::min(vel_in,vel_bal);
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }else
                                {
                                    v1 = std::max(std::min(vel_bal,v1),vel_in);
                                }

                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;

                        //if decellaratingto a level that it might stop moving
                        }else
                        {
                            float vel_bal = 0.0;

                            bool slow = true;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = 0.0;
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            float v00 = vel_in;

                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                                        float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }
                                //we have stalled, check if we can continue by depositing partially, otherwise just continue zero momentum and mass
                                if(v1 < 1e-6)
                                {
                                    //we want to deposit as much as we need to maintain momentum, meaning our slope angle (due to deposition) should change to be equal to the basal friction angle
                                    //(std::tan(ifa) - s)*dx * dx is what we deposit for this
                                    deposit  = std::max(0.0f,(std::tan(ifa) - s)*dx*dx*(dx+width_org->data[r][c]));
                                    v0 = v00;
                                    break;
                                }


                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;
                        }


                        //deal with deposits
                        if(deposit > 0.0)
                        {
                            float vfrac_dep = std::min(1.0f,std::max(0.0f,deposit/(std::max(0.001f,water+solid))));
                            water = water * (1.0-vfrac_dep);
                            solid = solid * (1.0-vfrac_dep);
                        }
                         velavg->data[r][c] = flowwidth;;
                        velacc->data[r][c] = a;
                        //velavg->data[r][c] = vel_in;
                        //velacc->data[r][c] = a;
                        //deal with calculating outgoing mass and momentum
                        deposition->data[r][c] = deposit;
                        momentum->data[r][c] = vel_final * (1000.0 * water  + dens * solid);
                        materialf->data[r][c] = water;
                        materials->data[r][c] = solid;


                    }
                }
            }

            //then we move in left-upper direction
            for(int r = map->data.nr_rows()-1; r > -1 ;r--)
            {
                for(int c = 0; c < map->data.nr_cols();c++)
                {
                    float v_points = DEM->data[r][c];
                    if(!pcr::isMV(v_points))
                    {


                        double w_x1 = 0.0;
                        double w_x2 = 0.0;
                        double w_y1 = 0.0;
                        double w_y2 = 0.0;

                        double wf_x1 = 0.0;
                        double wf_x2 = 0.0;
                        double wf_y1 = 0.0;
                        double wf_y2 = 0.0;

                        double ws_x1 = 0.0;
                        double ws_x2 = 0.0;
                        double ws_y1 = 0.0;
                        double ws_y2 = 0.0;
                        if(c > 0)
                        {
                            if((std::fabs(fracx2->data[r][c-1]) > 1e-10) && !pcr::isMV(DEM->data[r][c-1]))
                            {
                                float h = momentum->data[r][c-1];
                                w_x1 = std::fabs(fracx2->data[r][c-1]) * h;
                                float hf = materialf->data[r][c-1];
                                wf_x1 = std::fabs(fracx2->data[r][c-1]) * hf;
                                float hs = materials->data[r][c-1];
                                ws_x1 = std::fabs(fracx2->data[r][c-1]) * hs;
                            }
                        }
                        if(r > 0)
                        {
                            if((std::fabs(fracy2->data[r-1][c]) > 1e-10) && !pcr::isMV(DEM->data[r-1][c]))
                            {
                                float h = momentum->data[r-1][c];
                                w_y1 = std::fabs(fracy2->data[r-1][c]) *  h;
                                float hf = materialf->data[r-1][c];
                                wf_y1 = std::fabs(fracy2->data[r-1][c]) * hf;
                                float hs = materials->data[r-1][c];
                                ws_y1 = std::fabs(fracy2->data[r-1][c]) * hs;
                            }
                        }
                        if(c < map->nrCols()-1)
                        {
                            if((std::fabs(fracx1->data[r][c+1]) > 1e-10) && !pcr::isMV(DEM->data[r][c+1]))
                            {
                                float h = momentum->data[r][c+1];
                                w_x2 = std::fabs(fracx1->data[r][c+1]) * h;
                                float hf = materialf->data[r][c+1];
                                wf_x2 = std::fabs(fracx1->data[r][c+1]) * hf;
                                float hs = materials->data[r][c+1];
                                ws_x2 = std::fabs(fracx1->data[r][c+1]) * hs;
                            }
                        }
                        if(r < map->nrRows()-1)
                        {
                            if((std::fabs(fracy1->data[r+1][c]) > 1e-10) && !pcr::isMV(DEM->data[r+1][c]))
                            {
                                float h = momentum->data[r+1][c];
                                w_y2 = std::fabs(fracy1->data[r+1][c]) * h;
                                float hf = materialf->data[r+1][c];
                                wf_y2 = std::fabs(fracy1->data[r+1][c]) * hf;
                                float hs = materials->data[r+1][c];
                                ws_y2 = std::fabs(fracy1->data[r+1][c]) * hs;
                            }
                        }



                        //get local properties
                        float momentum_in = w_x1 + w_x2 + w_y1 + w_y2;
                        float ifa = ifamap->data[r][c];
                        float s = slope->data[r][c];
                        float water = wf_x1 + wf_x2 + wf_y1 + wf_y2 + source_fluid->data[r][c] * in_dx* in_dx;
                        float solid = ws_x1 + ws_x2 + ws_y1 + ws_y2 + source_solid->data[r][c] * in_dx* in_dx;
                        float n = manning->data[r][c];
                        float dens = density->data[r][c];

//check potential inflow into depression 
                        if(DEPRESSIONID != nullptr)
                        {
                            //chcek if there is a local depression here
                            if(!pcr::isMV(DEPRESSIONID->data[r][c]))
                            {
                                float frac_depr = 0.0;
                                int id = DEPRESSIONID->data[r][c];
                                if(!((id >= DEPRESSIONVOLS.size()) || id < 0))
                                {
                                    float vol = DEPRESSIONVOLS.at(id);
                                    
                                    float totalflow = water + solid;
                                     float frac_depr = std::max(0.0f,std::min(1.0f,vol/std::max(0.00001f,totalflow)));
                                    water = water * (1.0-frac_depr);
                                    solid = solid * (1.0-frac_depr);

                                    //FF::cout << "depr: " << r << " " << c << " " << frac_depr << " " << vol << " " << totalflow << std::endl;
                                    
                                    //otherwise do deposition
                                    DEPRESSIONCVOLS.at(id) += frac_depr * totalflow;

                                    //fractionally remove momentum
                                    momentum_in = momentum_in * (1.0-frac_depr);
                                }
                            }
                        }

                        float material_total = water + solid;

                        //estimate height
                        float blocklength = dist_variance->data[r][c];

                        float flowheight = material_total/std::max(blocklength*dx,0.001f);
                        //get flow width and correct flow height
                        float ratio_vol = material_total/(std::max(0.001f,acc_org->data[r][c])*dx*dx);
                        float flowwidth = ratio_vol * width_org->data[r][c];
                        flowheight = flowheight * dx/std::max(dx,flowwidth);

                        //get solid volumetric content ratio
                        float sratio = std::min(1.0f,solid / std::max(0.001f,solid+water));

                        float weight = (dens * sratio + 1000.0 * (1.0-sratio)) * flowheight;
                        float basal_porepressure = flowheight * 1000.0 * (1.0-sratio);

                        //get coefficient for the differential equation
                        float acc = 9.81 * (s - std::tan(ifa)*std::max(0.0f,(1.0f - basal_porepressure/weight)));

                        float vel_in = momentum_in/(std::max(0.00001f,water*1000.0f+solid*dens));
                        float a = acc;
                        float b = std::max(9.81f*n*n/std::max(0.01f,flowheight),0.01f);

                        //are we increasing or decreasing momentum?

                        float deposit = 0.0;
                        float vel_final = vel_in;

                        //if accalerating (at least, relative to rest)
                        if(acc > 0.0)
                        {
                            //get balance velocity
                            float vel_bal = std::max(0.01f,std::sqrt(a/b));

                            bool slow = vel_bal < vel_in;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = std::min(vel_in,vel_bal);
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }else
                                {
                                    v1 = std::max(std::min(vel_bal,v1),vel_in);
                                }

                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;

                        //if decellaratingto a level that it might stop moving
                        }else
                        {
                            float vel_bal = 0.0;

                            bool slow = true;

                            //now since we always stay between vel_in and vel_bal, we can adjust the timestep to this.
                            //make sure it takes at most 20 steps to cover the entire pixel
                            float vel_min = 0.0;
                            float dt_start = iters_inv*dx/std::max(0.01f,vel_in);
                            float dt_min = iters_inv*dx/std::max(0.01f,vel_min);

                            //numerically integrate position and velocity, untill reaching dx distance
                            //prevent oscillating by limiting to balance velocity
                            float x0 = 0.0;
                            float v0 = vel_in;
                            float v00 = vel_in;

                            for(int i = 0; i < NUM_ITERS_LANDSLIDE_MOMENTUM; i++)
                            {
                                float dt = std::max(dt_min,iters_inv*dx/std::max(0.01f,v0));
                                                        float v1 = v0 + (a - b * v0 * v0)*dt;

                                if(slow)
                                {
                                    v1 = std::min(std::max(vel_bal,v1),vel_in);
                                }
                                //we have stalled, check if we can continue by depositing partially, otherwise just continue zero momentum and mass
                                if(v1 < 1e-6)
                                {
                                    //we want to deposit as much as we need to maintain momentum, meaning our slope angle (due to deposition) should change to be equal to the basal friction angle
                                    //(std::tan(ifa) - s)*dx * dx is what we deposit for this
                                    deposit  = std::max(0.0f,(std::tan(ifa) - s)*dx*dx*(dx+width_org->data[r][c]));
                                    v0 = v00;
                                    break;
                                }


                                float x1 = x0 + v0 * dt;
                                if(x1 > dx)
                                {
                                    //we skipped the cell end, so lets interpolate linearly to get an approximation for the final velocity

                                    float lfrac = (dx - x0)/std::max(0.01f,(x1-x0));
                                    v0 = lfrac * v1 + (1.0-lfrac) * v0;
                                    break;
                                }else
                                {
                                    v0 = v1;
                                    x0 = x1;
                                }
                            }
                            vel_final = v0;
                        }


                        //deal with deposits
                        if(deposit > 0.0)
                        {
                            float vfrac_dep = std::min(1.0f,std::max(0.0f,deposit/(std::max(0.001f,water+solid))));
                            water = water * (1.0-vfrac_dep);
                            solid = solid * (1.0-vfrac_dep);
                        }
                        velavg->data[r][c] = flowwidth;;
                        velacc->data[r][c] = a;
                        //velavg->data[r][c] = vel_in;
                        //velacc->data[r][c] = a;
                        //deal with calculating outgoing mass and momentum
                        deposition->data[r][c] = deposit;
                        momentum->data[r][c] = vel_final * (1000.0 * water  + dens * solid);
                        materialf->data[r][c] = water;
                        materials->data[r][c] = solid;


                    }


                }
            }

        }
        first = false;

    }

    delete slope;
    delete map;

    return {materialf, materials, velavg,velacc,deposition,momentum};

}
