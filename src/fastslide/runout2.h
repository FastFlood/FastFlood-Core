#ifndef RUNOUT2_H
#define RUNOUT2_H

#define NUM_ITERS_LANDSLIDE_MOMENTUM 10

#include "../model/fastflood/flowwidth.h"

inline void ls_momentumsolve(int r, int c, float dx, cTMapf * DEM,cTMapf * map, cTMapf * manning,  cTMapf * momentum, cTMapf * fracx1,cTMapf * fracx2,cTMapf * fracy1, cTMapf * fracy2,cTMapf * materialf, cTMapf * materials, cTMapf * DEPRESSIONID,
cTMapf * source_solid,cTMapf * source_fluid,cTMapf * ifamap,cTMapf * rocksizemap, cTMapf * density, cTMapf * dist_variance,cTMapf * width_org,cTMapf * acc_org,cTMapf * velavg,cTMapf * velacc,cTMapf * slope,
std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, cTMapf *deposition, cTMapf * HFLOOD)
{


    float in_dx = dx;

    float iters_inv = 1.0/((float)(NUM_ITERS_LANDSLIDE_MOMENTUM));

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

    float wh_flood = 0.0;
    if(HFLOOD != nullptr)
    {
        wh_flood = HFLOOD->data[r][c];
    }

    
    if(solid + water + wh_flood < 1e-8)
    {
        velavg->data[r][c] = 0.0;
        velacc->data[r][c] = 0.0;
        //velavg->data[r][c] = vel_in;
        //velacc->data[r][c] = a;
        //deal with calculating outgoing mass and momentum
        deposition->data[r][c] = 0.0;
        momentum->data[r][c] = 0.0;
        materialf->data[r][c] = water;
        materials->data[r][c] = solid;
        return;
    }


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
                float vol = std::max(0.0f,DEPRESSIONVOLS.at(id)-DEPRESSIONCVOLS.at(id));
                
                float totalflow = water + solid;
                float frac_depr = std::max(0.0f,std::min(1.0f,vol/std::max(0.00001f,totalflow)));
                water = water * (1.0-frac_depr);
                solid = solid * (1.0-frac_depr);

                //FF::cout << "depr: " << DEPRESSIONID->data[r][c] << " " << r << " " << c << " " << frac_depr << " " << vol << " " << totalflow << std::endl;
                
                //otherwise do deposition
                DEPRESSIONCVOLS.at(id) =  std::min(DEPRESSIONVOLS.at(id),DEPRESSIONCVOLS.at(id) + frac_depr * totalflow);

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
    float flowwidth = ratio_vol * width_org->data[r][c] + wh_flood;
    flowheight = flowheight * dx/std::max(dx,flowwidth);

    //get solid volumetric content ratio
    float sratio = std::min(1.0f,solid / std::max(0.001f,solid+water));
    float sratio_fl = std::min(1.0f,flowheight * sratio  / std::max(0.001f,flowheight + wh_flood));
    float fratio = std::min(1.0f,water / std::max(0.001f,solid+water));

    //estimate viscosity
    float visc = std::min(100.0,0.032 * std::exp(10.0 * sratio_fl));

    //based on jamming ratio distance? https://www.pnas.org/doi/10.1073/pnas.2209109119
    //better relation with grain sizes in there, upper bound 100? 
    //visc = std::min(100.0,0.032 * std::exp(10.0 * sratio));

    //based on simple volumetric concentration? https://www.engr.colostate.edu/~pierre/ce_old/resume/Paperspdf/O'Brien-Julien-ASCE88.pdf
    //Traditional and well-tested, upper bound 100
    //estimate jamming ratio (0.499 + polydispersity (m3/m22) * 0.0019) Mk = the kth moment of the size distribution ( https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.040302)
    //estimate maximimal packing density of the grains (jamming ratio * 1.05)
    //visc = 0.01* delta_jamming ^2  (simple fit from graph in paper)

    //further simplificatition based on two-parameter model (DebrisIntermixing: https://gmd.copernicus.org/articles/9/2909/2016/gmd-9-2909-2016.pdf)

    //overview HEC: https://www.hec.usace.army.mil/confluence/rasdocs/rasmuddebris/non-newtonian-user-s-manual/introduction-taxonomy-and-rheology-of-debris-flows
    //paper behind HEC: Comparing single-phase, non-Newtonian approaches with experimental results: Validating flume-scale mud and debris flow in HEC-RAS

    float weight = (dens * sratio + 1000.0 * (1.0-sratio)) * flowheight;
    float basal_porepressure = flowheight * 1000.0 * (1.0-sratio) + wh_flood * 1000.0;

    float gamma = dens/1000.0;
    float invgamma = 1000.0/std::max(0.001f,dens);

    //get coefficient for the differential equation
    float acc = 9.81 * gamma* (s - sratio *std::tan(ifa)*std::max(0.0f,(1.0f - basal_porepressure/weight)));

    float vel_in = momentum_in/(std::max(0.00001f,water*1000.0f+solid*dens));
    float a = acc;
    float b = invgamma * std::max(9.81f*n*n/std::max(0.01f,flowheight + wh_flood),0.01f);
    float cc = invgamma * visc * 3.0 * fratio*(1.0/std::max(0.01f,flowheight + wh_flood));


    //drag force?? dependent on velocity differences between flood waters and debris mixture
    //, so not balance any more ... that sucks

    //are we increasing or decreasing momentum?

    float deposit = 0.0;
    float vel_final = vel_in;

    //if accalerating (at least, relative to rest)
    if(water + solid > 1e-6)
    {
        if(acc > 0.0)
        {
            //get balance velocity
            float vel_bal = std::max(0.01f,(-cc + std::sqrt(cc*cc + 4.0f * a*b))/(2.0f*b));

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
                float v1 = v0 + (a -cc * v0 - b * v0 * v0)*dt;

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
    }else{
        vel_final = 0.0;
    }


    //deal with deposits
    if(deposit > 0.0)
    {
        float vfrac_dep = std::min(1.0f,std::max(0.0f,deposit/(std::max(0.001f,water+solid))));
        water = water * (1.0-vfrac_dep);
        solid = solid * (1.0-vfrac_dep);
    }

        velavg->data[r][c] = vel_final;
    velacc->data[r][c] = a;
    //velavg->data[r][c] = vel_in;
    //velacc->data[r][c] = a;
    //deal with calculating outgoing mass and momentum
    deposition->data[r][c] = deposit;
    momentum->data[r][c] = vel_final * (1000.0 * water  + dens * solid);
    materialf->data[r][c] = water;
    materials->data[r][c] = solid;

}


inline void ls_momentumerosionsolve(int r, int c, float dx, cTMapf * DEM,cTMapf * map, cTMapf * manning,  cTMapf * momentum, cTMapf * fracx1,cTMapf * fracx2,cTMapf * fracy1, cTMapf * fracy2,cTMapf * materialf, cTMapf * materials, cTMapf * DEPRESSIONID,
cTMapf * source_solid,cTMapf * source_fluid,cTMapf * ifamap,cTMapf * rocksizemap, cTMapf * density, cTMapf * dist_variance,cTMapf * width_org,cTMapf * acc_org,cTMapf * velavg,cTMapf * velacc,cTMapf * slope,
std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, cTMapf *deposition, cTMapf * HFLOOD, cTMapf * COH, cTMapf *IFAB, cTMapf * SD, std::vector<cTMapf *> &CMAPS_AREA, std::vector<cTMapf *> &CMAPS_COH,std::vector<cTMapf *> &CMAPS_IFAB,std::vector<cTMapf *> &CMAPS_SD, cTMapf * EROSION, float erosion_mult)
{


    float in_dx = dx;

    float iters_inv = 1.0/((float)(NUM_ITERS_LANDSLIDE_MOMENTUM));

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

    float wh_flood = 0.0;
    if(HFLOOD != nullptr)
    {
        wh_flood = HFLOOD->data[r][c];
    }

     if(solid + water + wh_flood < 1e-8)
    {
        velavg->data[r][c] = 0.0;
        velacc->data[r][c] = 0.0;
        //velavg->data[r][c] = vel_in;
        //velacc->data[r][c] = a;
        //deal with calculating outgoing mass and momentum
        deposition->data[r][c] = 0.0;
        momentum->data[r][c] = 0.0;
        materialf->data[r][c] = water;
        materials->data[r][c] = solid;
        return;
    }


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
                float vol = std::max(0.0f,DEPRESSIONVOLS.at(id)-DEPRESSIONCVOLS.at(id));
                
                float totalflow = water + solid;
                float frac_depr = std::max(0.0f,std::min(1.0f,vol/std::max(0.00001f,totalflow)));
                water = water * (1.0-frac_depr);
                solid = solid * (1.0-frac_depr);

                //FF::cout << "depr: " << DEPRESSIONID->data[r][c] << " " << r << " " << c << " " << frac_depr << " " << vol << " " << totalflow << std::endl;
                
                //otherwise do deposition
                DEPRESSIONCVOLS.at(id) =  std::min(DEPRESSIONVOLS.at(id),DEPRESSIONCVOLS.at(id) + frac_depr * totalflow);

                //fractionally remove momentum
                momentum_in = momentum_in * (1.0-frac_depr);
            }
        }
    }

    float material_total = water + solid;

    //estimate height
    float blocklength = dist_variance->data[r][c];

    float flowheight = material_total/std::max(blocklength*dx,dx*dx);
    float cross_area = flowheight * dx;
    //get flow width and correct flow height
    float ratio_vol = material_total/(std::max(0.001f,acc_org->data[r][c])*dx*dx);
    float flowwidth = ratio_vol * width_org->data[r][c] + wh_flood;
    flowheight = flowheight * dx/std::max(dx,flowwidth);

    //get solid volumetric content ratio
    float sratio = std::min(1.0f,solid / std::max(0.001f,solid+water));
    float sratio_fl = std::min(1.0f,flowheight * sratio  / std::max(0.001f,flowheight + wh_flood));
    float fratio = std::min(1.0f,water / std::max(0.001f,solid+water));

    //estimate viscosity
    float visc = std::min(100.0,0.032 * std::exp(10.0 * sratio_fl));

    //based on jamming ratio distance? https://www.pnas.org/doi/10.1073/pnas.2209109119
    //better relation with grain sizes in there, upper bound 100? 
    //visc = std::min(100.0,0.032 * std::exp(10.0 * sratio));

    //based on simple volumetric concentration? https://www.engr.colostate.edu/~pierre/ce_old/resume/Paperspdf/O'Brien-Julien-ASCE88.pdf
    //Traditional and well-tested, upper bound 100
    //estimate jamming ratio (0.499 + polydispersity (m3/m22) * 0.0019) Mk = the kth moment of the size distribution ( https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.040302)
    //estimate maximimal packing density of the grains (jamming ratio * 1.05)
    //visc = 0.01* delta_jamming ^2  (simple fit from graph in paper)

    //further simplificatition based on two-parameter model (DebrisIntermixing: https://gmd.copernicus.org/articles/9/2909/2016/gmd-9-2909-2016.pdf)

    //overview HEC: https://www.hec.usace.army.mil/confluence/rasdocs/rasmuddebris/non-newtonian-user-s-manual/introduction-taxonomy-and-rheology-of-debris-flows
    //paper behind HEC: Comparing single-phase, non-Newtonian approaches with experimental results: Validating flume-scale mud and debris flow in HEC-RAS

    float weight = (dens * sratio + 1000.0 * (1.0-sratio)) * flowheight;
    float basal_porepressure = flowheight * 1000.0 * (1.0-sratio) + wh_flood * 1000.0;

    float gamma = dens/1000.0;
    float invgamma = 1000.0/std::max(0.001f,dens);

    //get coefficient for the differential equation
    float acc = 9.81 * gamma* (s - sratio *std::tan(ifa)*std::max(0.0f,(1.0f - basal_porepressure/weight)));

    float vel_in = momentum_in/(std::max(0.00001f,water*1000.0f+solid*dens));
    float a = acc;
    float b = invgamma * std::max(9.81f*n*n/std::max(0.01f,flowheight + wh_flood),0.01f);
    float cc = invgamma * visc * 3.0 * fratio*(1.0/std::max(0.01f,flowheight + wh_flood));


    //drag force?? dependent on velocity differences between flood waters and debris mixture
    //, so not balance any more ... that sucks

    //are we increasing or decreasing momentum?

    float deposit = 0.0;
    float vel_final = vel_in;

float t_erosion = 0.0;
    //if accalerating (at least, relative to rest)
    if(water + solid > 1e-6)
    {
        if(acc > 0.0)
        {
            //get balance velocity
            float vel_bal = std::max(0.01f,(-cc + std::sqrt(cc*cc + 4.0f * a*b))/(2.0f*b));


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
                float v1 = v0 + (a -cc * v0 - b * v0 * v0)*dt;

    t_erosion += dt;
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

                t_erosion += dt;

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
                    t_erosion = 0.0;
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
    }else{
        vel_final = 0.0;
    }

    //now check erosion

    //get average soil properties based on flow width

    //get current flow cross-sectional area

    //linearly interpolate between options in cascading map data

    float coh_here = 0.0;
    float ifab_here = 1.2;
    float sd_here = 0.0;

    float area_p = 0.0;
    for(int i = 0; i < CMAPS_AREA.size(); i++)
    {
        if(cross_area < CMAPS_AREA[i]->data[r][c])
        {
            float area_frac = std::max(0.0f,std::min(1.0f,(cross_area - area_p)/std::max(0.001f,(CMAPS_AREA[i]->data[r][c] -area_p))));
            if(i == 0)
            {
                coh_here = (1.0-area_frac) * COH->data[r][c] + area_frac * CMAPS_COH[i]->data[r][c];
                ifab_here = (1.0-area_frac) * IFAB->data[r][c] +area_frac * CMAPS_IFAB[i]->data[r][c];
                sd_here = (1.0-area_frac) * SD->data[r][c] + area_frac * CMAPS_SD[i]->data[r][c];
            }else
            {
                coh_here = (1.0-area_frac) * CMAPS_COH[i]->data[r][c] + area_frac * CMAPS_COH[i]->data[r][c];
                ifab_here = (1.0-area_frac) * CMAPS_IFAB[i]->data[r][c] +area_frac * CMAPS_IFAB[i]->data[r][c];
                sd_here = (1.0-area_frac) * CMAPS_SD[i]->data[r][c] + area_frac * CMAPS_SD[i]->data[r][c];


            }
        }else if(i == CMAPS_AREA.size()-1)
        {
            coh_here = CMAPS_COH[i]->data[r][c];
            ifab_here = CMAPS_IFAB[i]->data[r][c];
            sd_here = CMAPS_SD[i]->data[r][c];
        }else
        {
            area_p = CMAPS_AREA[i]->data[r][c];
        }
    }

    //calculate erosion of block flowing over soil

    //calculate chezy coefficient
    float chezy = (1.0/n) * std::pow(flowheight + wh_flood,0.166);
    //lower-flow boundary mobility factor for fluid and solid

    //gamma
    float vhere = vel_final;
    float zeta = std::sin(std::atan(s));//slope angle
    float sqrte = 0.025;
    float e_eff = 0.79 * std::exp(-0.85*coh_here);
    float v_waterman = (1.0/n) * std::pow(wh_flood,2.0/3.0) *sqrt(s);


    if(wh_flood > 1e-6)
    {
        t_erosion = 0.5 * t_erosion + std::min(2.5f * dx,dx/std::max(0.1f,v_waterman));
    }

    float ef = erosion_mult * e_eff * chezy * 0.5 *(0.1 + 0.1 * sratio_fl)* (std::max(0.0,1.0 - 2.5 * sratio_fl)) * (wh_flood + flowheight) * (vhere*flowheight + wh_flood * v_waterman)/std::max(0.001f,wh_flood + flowheight);
    
    //limit fluid flow erosion to 25 percent of flow volume
    float conc_base = 0.6;
    float conc_crit = 0.6;
    
    float es = std::max(0.0f,(float)(erosion_mult *e_eff * (1.0 - gamma)*(ifa - ifab_here *conc_base/std::max(0.1f,sratio))*0.035*cos(zeta)*sqrt(9.81 * flowheight * sratio)*(1.0 + conc_base/std::max(0.1f,sratio))));
    
    float erosion_total_f = std::min(dx*wh_flood*dx * 0.35f,std::min(t_erosion*(ef+es),sd_here) * dx * std::min(flowwidth,100.0f));
    float erosion_total_s = std::min(std::max(0.0f,conc_base * (ifab_here/std::max(0.01f,ifa)) - sratio)* dx*flowheight*std::max(dx,std::min(flowwidth,100.0f)), std::min(t_erosion*(es),sd_here) * dx * std::min(flowwidth,100.0f));

    float erosion_total = erosion_total_f + erosion_total_s;

    //add materials
    water = water + conc_base * erosion_total;
    solid = solid + (1.0-conc_base) * erosion_total;

    //alter momentum balance
    //deal with deposits
    EROSION->data[r][c] = erosion_total;// / (std::max(0.00001f,e_eff*flowheight));
    
    if(deposit > 0.0)
    {
        float vfrac_dep = std::min(1.0f,std::max(0.0f,deposit/(std::max(0.001f,water+solid))));
        water = water * (1.0-vfrac_dep);
        solid = solid * (1.0-vfrac_dep);
    }

        velavg->data[r][c] = vel_final;
    velacc->data[r][c] = a;
    //velavg->data[r][c] = vel_in;
    //velacc->data[r][c] = a;
    //deal with calculating outgoing mass and momentum
    deposition->data[r][c] = deposit;
    momentum->data[r][c] = vel_final * (1000.0 * water  + dens * solid);
    materialf->data[r][c] = water;
    materials->data[r][c] = solid;

}


inline std::vector<cTMapf *> DinfiniteLandslideMomentumDepositionSolve2(std::vector<cTMapf*> dinf,cTMapf * DEM, cTMapf * manning, cTMapf * source_solid, cTMapf * source_fluid, cTMapf * density, cTMapf * acc_org, cTMapf * width_org, cTMapf * ifamap,cTMapf * rocksizemap, cTMapf * HFLOOD, cTMapf * dist_variance,cTMapf * DEPRESSIONID, std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, int iter_max = 0)
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



    //we keep iterating through this algorithm untill there is no change left to make
    bool change = true;
    bool first = true;

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

                        ls_momentumsolve(r, c, dx, DEM, map, manning, momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition, HFLOOD);





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


                        ls_momentumsolve(r, c, dx, DEM,map,manning,momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
                        source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition,HFLOOD);



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


                           ls_momentumsolve(r, c,dx, DEM, map, manning,momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition, HFLOOD);

     

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

                        ls_momentumsolve(r, c, dx,DEM, map, manning,momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
                        source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition,HFLOOD);



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






inline std::vector<cTMapf *> DinfiniteLandslideMomentumDepositionErosionSolve2(std::vector<cTMapf*> dinf,cTMapf * DEM, cTMapf * manning, cTMapf * source_solid, cTMapf * source_fluid, cTMapf * density, cTMapf * acc_org, cTMapf * width_org, cTMapf * ifamap,cTMapf * rocksizemap, cTMapf * HFLOOD, cTMapf * dist_variance,cTMapf * DEPRESSIONID, std::vector<float> & DEPRESSIONVOLS, std::vector<float> & DEPRESSIONCVOLS, cTMapf * COH, cTMapf * IFAB, cTMapf * SD, cTMapf * ERODE, float erosion_mult, int iter_max = 0)
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


    //for erosion, we want the exponentially-cascaded structure for soil properties
    std::vector<cTMapf*> maps_soilp_area;
    std::vector<cTMapf*> maps_soilp_coh;
    std::vector<cTMapf*> maps_soilp_ifab;
    std::vector<cTMapf*> maps_soilp_sd;

    GetAvSoilPropsFromRiverVolume(DEM,COH,maps_soilp_coh,maps_soilp_area,5,true);
    GetAvSoilPropsFromRiverVolume(DEM,IFAB,maps_soilp_ifab,maps_soilp_area,5,false);
    GetAvSoilPropsFromRiverVolume(DEM,SD,maps_soilp_sd,maps_soilp_area,5,false);


    //we keep iterating through this algorithm untill there is no change left to make
    bool change = true;
    bool first = true;

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

                        ls_momentumerosionsolve(r, c, dx, DEM, map, manning, momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition, HFLOOD,COH, IFAB, SD, maps_soilp_area,maps_soilp_coh,maps_soilp_ifab,maps_soilp_sd,ERODE, erosion_mult);





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


                        ls_momentumerosionsolve(r, c, dx, DEM,map,manning,momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
                        source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition,HFLOOD,COH, IFAB, SD,maps_soilp_area, maps_soilp_coh,maps_soilp_ifab,maps_soilp_sd,ERODE,erosion_mult);



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


                           ls_momentumerosionsolve(r, c,dx, DEM, map, manning,momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition, HFLOOD,COH, IFAB, SD,maps_soilp_area, maps_soilp_coh,maps_soilp_ifab,maps_soilp_sd,ERODE, erosion_mult);

     

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

                        ls_momentumerosionsolve(r, c, dx,DEM, map, manning,momentum, fracx1,fracx2,fracy1, fracy2,materialf, materials, DEPRESSIONID,
                        source_solid,source_fluid,ifamap,rocksizemap, density, dist_variance,width_org,acc_org,velavg,velacc,slope,DEPRESSIONVOLS, DEPRESSIONCVOLS, deposition,HFLOOD,COH, IFAB, SD,maps_soilp_area, maps_soilp_coh,maps_soilp_ifab,maps_soilp_sd,ERODE, erosion_mult);



                    }


                }
            }

        }
        first = false;

    }

    for(int i = 0; i < maps_soilp_area.size(); i++)
    {
        delete maps_soilp_area[i];
    }
    for(int i = 0; i < maps_soilp_coh.size(); i++)
    {
        delete maps_soilp_coh[i];
    }
    for(int i = 0; i < maps_soilp_ifab.size(); i++)
    {
        delete maps_soilp_ifab[i];
    }
    for(int i = 0; i < maps_soilp_sd.size(); i++)
    {
        delete maps_soilp_sd[i];
    }

    delete slope;
    delete map;

    return {materialf, materials, velavg,velacc,deposition,momentum};

}

#endif