#ifndef FASTWIND_H
#define FASTWIND_H

#include <fastmath.h>

#include "../geospatial/map/map.h"
#include "../common/threading/parallel.h"


/*
 -key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -sim -dem bolund_dem.tif -vxbcond 10 -vybcond 0 -dout dem.tif -bhout bh.tif -vtotout vtot.tif -vxyzout v.tif -gout g.tif -pout p.tif


-key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -sim -dem dem_block.tif -man rl_block.tif -bh bh_block.tif -vxbcond 10 -vybcond 0 -dout dem.tif -bhout bh.tif -vtotout vtot.tif -vxyzout v.tif -gout g.tif -pout p.tif
-key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -sim -autoscenario -61.559959018 15.688957979 -61.165043425 15.169320632 50 high 270 -dout dem.tif -bhout bh.tif -vtotout vtot.tif -vxyzout v.tif -gout g.tif -pout p.tif
-key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -sim -autoscenario -61.399107658355845 15.321686105862176 -61.36651054564369 15.29024521473749 50 high 270 -dout dem.tif -bhout bh.tif -vtotout vtot.tif -vxyzout v.tif -gout g.tif -pout p.tif

-autohazard flashfloodfluvialflood 100 medium
-key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -sim -autoscenario -61.391806431 15.302905255 -61.377766048 15.29216288 50 medium 270 -dout dem.tif -bhout bh.tif -vtotout vtot.tif -vxyzout v.tif -gout g.tif -pout p.tif
-key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -sim -autohazard -61.39863489324715 15.320432145233246 -61.366037975773295 15.288991254108561 50 medium 270 -dout dem.tif -bhout bh.tif -vtotout vtot.tif -vxyzout v.tif -gout g.tif -pout p.tif
-key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -key Xx1BfVGbNP76RwYlmQtR8os4halVCpxH -sim -autohazard -61.39863489324715 15.320432145233246 -61.366037975773295 15.288991254108561 50 wind medium 270 -dout dem.tif -bhout bh.tif -vtotout vtot.tif -vxyzout v.tif -gout g.tif -pout p.tif
*/

#define FW_MAX_SLOPE 1.0f

static const bool BLOCK_BUILDINGS = true;

/**
 * @brief Helper to dispatch row iterations across the thread pool via efficient block chunking.
 * Each chunk runs contiguous rows directly on a thread without work-stealing or per-row atomics.
 * Uses a local completion counter so concurrent directional simulations never block each other
 * on global thread pool barriers.
 */
template <typename F>
inline void ParallelOrSerialFor(int start, int end, F&& func, int min_rows_for_parallel = 16) {
    int total = end - start;
    if (total <= 0) return;
    int n_threads = static_cast<int>(GetTPool().get_active_threads());
    if (total < min_rows_for_parallel || n_threads <= 1) {
        for (int i = start; i < end; ++i) {
            func(i);
        }
        return;
    }

    int n_chunks = std::min(total, n_threads);
    std::atomic<int> remaining(n_chunks - 1);

    for (int k = 1; k < n_chunks; ++k) {
        int r_start = start + (total * k) / n_chunks;
        int r_end = start + (total * (k + 1)) / n_chunks;

        GetTPool().push([&func, r_start, r_end, &remaining]() {
            for (int r = r_start; r < r_end; ++r) {
                func(r);
            }
            remaining.fetch_sub(1, std::memory_order_release);
        });
    }

    // Calling thread computes chunk 0 directly
    int r0_start = start;
    int r0_end = start + total / n_chunks;
    for (int r = r0_start; r < r0_end; ++r) {
        func(r);
    }

    while (remaining.load(std::memory_order_acquire) > 0) {
        std::this_thread::yield();
    }
}

/**
 * @brief Resamples a source grid to a destination grid.
 *
 * @param src The source map to resample from.
 * @param dst The destination map to resample to.
 * @param method The resampling method (1 for bilinear interpolation).
 */
inline void ResampleGrid(cTMapf * src, cTMapf * dst, int method) {
    if (method == 1) { // Bilinear interpolation
        ParallelOrSerialFor(0, dst->nrRows(), [&](int r) {
            for (int c = 0; c < dst->nrCols(); ++c) {
                // Calculate world coordinates of the center of the destination cell
                double x = dst->west() + (c + 0.5) * dst->cellSizeX();
                // Note: y-coordinate calculation depends on whether cellSizeY is positive or negative.
                // The SampleLinear function should handle this correctly.
                double y = dst->north() + (r + 0.5) * dst->cellSizeY();

                // Sample from the source map using bilinear interpolation
                dst->data[r][c] = src->SampleLinear(x, y);
            }
        });
    }
}

//FastWind function,
//does a mass balancing, physics corrected
//rapid diagnostic wind model
//data needed, extreme wind speeds ERA5 /ERA5-LAND + some statistical/AI refinement based on station wind data
//then local modelling based on approach below
//the gust output is the total velocity including gustiness, not the added gust velocity
//so it can be directly used
//returns 1 for success
//and negative for error
inline int FastWind(
        //raw high resolution elevation data in meter
        //must not be 0
        //must be somewhat smooth(continuous) (though we can handle a lot)
        //all other maps must in principle be in this same rows/cols/dx/dy/ulx/uly/brx/bry
    cTMapf * in_DEM,
        //building height in meters
        //we need to seperate as the topography-following coordinates
        //need a contiuous surface
    cTMapf * in_BH,
        //frictional length for the terrain, for energy loss
        //must not be 0
    const cTMapf * in_N,
        //the layer heights, would suggest something like
        //layers_z = {2,4,8,20,100,500}
    const std::vector<float> &layers_z,
        //the layer lateral pixel scales, would suggest something like
        //layers_z = {1,1,3,9,27,81}
        //but it needs to be in cascades of scale 1 or 3, and first one has to be 1
        //((layers_pxlatscale[i]/layers_pxlatscale[i-1]) == 3 || ((layers_pxlatscale[i]/layers_pxlatscale[i-1]) == 1
        //though the function checks this for you
    const std::vector<int> &layers_pxlatscale,
        //the reference velocity without topographical effects
        //forced at boundary conditions
    const cTMapf* VEL_HIGH_X,
    const cTMapf* VEL_HIGH_Y,
        //output, velocity maps, are going to be in coarsers
        //resolution, lower layer at high resolution first,
        //then coarser lower resolution layers
        //same numnber of maps as count of layers_z provided to function
        //these maps are allocated within this function, ownership returned to caller
        //these references to vectors of maps should be empty at the start, though we
        //will call clear before adding the actual output
    std::vector<std::shared_ptr<cTMapf>> &out_velx,
        //for y component
    std::vector<std::shared_ptr<cTMapf>> &out_vely,
        //for z component
    std::vector<std::shared_ptr<cTMapf>> &out_velz,
        //output pressure layers,
        std::vector<std::shared_ptr<cTMapf>> &out_P,
        //output gust speed layers
         std::vector<std::shared_ptr<cTMapf>> &out_gust,

         bool output_test,

        //iterations for the shadow sweeping
        int n_iters_sweepshadow = 6,
        //reference height for the input wind field
         float ref_wind_height = 2.0,
        //iterations for the potential flow
        int n_iters_flow = 250,
        //iterations for the turbulence propagation
        int n_iters_turb = 3,
        //Velocity shadowing coefficient
        float velocity_deficit_coefficient = 0.4,
        //critical shadowing angle
        float velocity_deficit_angle_strength = 0.5,
        //Wake expension rate
        float wake_diffusion_rate = 0.25,
        //vertical weight ratio multiplier
        float alpha_sq = 25.0f,
        // Reference roughness length for the input wind field (VEL_HIGH_X/Y)
        float ref_roughness_length = 0.25,
        // Von Karman constant for turbulence calculations
        float von_karman_constant = 0.41f
    ) {
    //Mass consistent Diagnostic Wind Model


    //allocate variables

    std::vector<std::shared_ptr<cTMapf>> MAXBH;//max building height underneath this cascade voxel
    std::vector<std::shared_ptr<cTMapf>> BBF; //building block fraction

    std::vector<std::shared_ptr<cTMapf>> WAKE_DIST;
    std::vector<std::shared_ptr<cTMapf>> WAKE_DIST_MAX;
    std::vector<std::shared_ptr<cTMapf>> WAKE_LAT_DIST;
    std::vector<std::shared_ptr<cTMapf>> WAKE_LAT_DIST_MAX;

    std::vector<std::shared_ptr<cTMapf>> PROF;
    std::vector<std::shared_ptr<cTMapf>> SHADE;
    std::vector<std::shared_ptr<cTMapf>> TKE;
    std::vector<std::shared_ptr<cTMapf>> RL;
    std::vector<std::shared_ptr<cTMapf>> N;
    std::vector<std::shared_ptr<cTMapf>> DEM;
    std::vector<std::shared_ptr<cTMapf>> PSI;
    std::vector<std::shared_ptr<cTMapf>> SINTF;
    std::vector<std::shared_ptr<cTMapf>> COSTF;
    std::vector<std::shared_ptr<cTMapf>> ux;
    std::vector<std::shared_ptr<cTMapf>> uy;
    std::vector<std::shared_ptr<cTMapf>> uz;
    std::vector<std::shared_ptr<cTMapf>> ux_initial;
    std::vector<std::shared_ptr<cTMapf>> uy_initial;
    std::vector<std::shared_ptr<cTMapf>> uz_initial;

    if (layers_z.empty() || layers_pxlatscale.empty()) {
        return -1;
    }

    //define grid sizes
    std::vector<float> layers_dz = {};
    std::vector<float> layers_dx = {};
    std::vector<int> layers_scaledown = {};
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        if (i > 0) {
            if (layers_z[i]-layers_z[i-1] < 1e-6) {
                return -1;
            }

            layers_dz.push_back(layers_z[i]-layers_z[i-1]);
        }else {
            layers_dz.push_back(layers_z[i]-0);
        }
    }



    //then setup 3d grid stack for variables
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {

        if (layers_pxlatscale[i] < 1) {
            return -1;
        }
        int cscale = 1;
        cTMapf* ref = in_DEM;
        cTMapf* refn = const_cast<cTMapf*>(in_N);
        if (i == 0) {
            cscale = 1;
            ref = in_DEM;
            refn = const_cast<cTMapf*>(in_N);
        }else
        {
            ref = DEM[i-1].get();
            refn = N[i-1].get();
            if (layers_pxlatscale[i] > layers_pxlatscale[i-1]) {
                cscale = layers_pxlatscale[i]/layers_pxlatscale[i-1];
            }else if(layers_pxlatscale[i] == layers_pxlatscale[i-1]) {
                cscale = 1;
            }else {
                return -1;
            }
        }

        std::cout << "cscale :  " << cscale << std::endl;
        layers_scaledown.push_back(cscale);

        cTMapf* map_casc = ref->GetCopy0Cascade(cscale);
        cTMapf* map_cascn = map_casc->GetCopy0();
        CreateMapCascadeAvg(ref,map_casc,cscale);
        CreateMapCascadeAvg(refn,map_cascn,cscale);
        DEM.push_back(std::shared_ptr<cTMapf>(map_casc));
        RL.push_back(map_casc->GetCopySP0());
        PSI.push_back(map_casc->GetCopySP0());
        N.push_back(std::shared_ptr<cTMapf>(map_cascn));
        ux.push_back(map_casc->GetCopySP0());
        uy.push_back(map_casc->GetCopySP0());
        uz.push_back(map_casc->GetCopySP0());
        ux_initial.push_back(map_casc->GetCopySP0());
        uy_initial.push_back(map_casc->GetCopySP0());
        uz_initial.push_back(map_casc->GetCopySP0());
        SINTF.push_back(map_casc->GetCopySP0());
        COSTF.push_back(map_casc->GetCopySP0());
        SHADE.push_back(map_casc->GetCopySP1());
        TKE.push_back(map_casc->GetCopySP0());

        PROF.push_back(map_casc->GetCopySP1());

        WAKE_DIST.push_back(map_casc->GetCopySP0());
        WAKE_DIST_MAX.push_back(map_casc->GetCopySP0());
        WAKE_LAT_DIST.push_back(map_casc->GetCopySP0());
        WAKE_LAT_DIST_MAX.push_back(map_casc->GetCopySP0());

        //track max building height per level
        //and fractional blocking per level
        BBF.push_back(map_casc->GetCopySP0());
        MAXBH.push_back(map_casc->GetCopySP0());

        if (i == 0)
        {
            for (int r= 0; r < BBF[i]->nrRows(); r++) {
                for (int c = 0; c < BBF[i]->nrCols(); c++) {
                    if (pcr::isMV(in_BH->data[r][c])) {
                        BBF[i]->data[r][c] = 0.0f;
                        MAXBH[i]->data[r][c] = 0.0f;
                    }else {
                        BBF[i]->data[r][c] = std::min(1.0f,in_BH->data[r][c]/layers_dz[i]);
                        MAXBH[i]->data[r][c] = in_BH->data[r][c];
                    }

                }
            }
        }else {
            for (int r = 0; r < BBF[i]->nrRows(); r++) {
                for (int c = 0; c < BBF[i]->nrCols(); c++) {
                    float max_bh = 0.0f;
                    float frac_cover = 0.0f;
                    //do we need to heck 3x3 subgrid
                    if (layers_scaledown[i] == 3) {
                        for ( int dr = 0; dr < 3; dr++) {
                            for ( int dc = 0; dc < 3; dc++) {
                                int r2 = r*3+dr;
                                int c2 = c*3+dc;
                                if (r2 < MAXBH[i-1]->nrRows() && c2 < MAXBH[i-1]->nrCols()) {
                                    max_bh = std::max(max_bh,MAXBH[i-1]->data[r2][c2]) ;
                                    frac_cover += (1.0f/9.0f) * std::min(1.0f,std::max(0.0f,(MAXBH[i-1]->data[r2][c2] - layers_z[i-1]))/layers_dz[i]);
                                }
                            }
                        }

                        MAXBH[i]->data[r][c] = max_bh;
                        BBF[i]->data[r][c] = frac_cover;

                    }else {
                        //otherwise just continue with a 1x1 voxel underneath
                        if (r < MAXBH[i-1]->nrRows() && c < MAXBH[i-1]->nrCols()){
                            MAXBH[i]->data[r][c] = MAXBH[i-1]->data[r][c];
                            //get fractional presence for remaining building height withing this voxel
                            BBF[i]->data[r][c] = std::min(1.0f,std::max(0.0f,(MAXBH[i]->data[r][c] - layers_z[i-1]))/layers_dz[i]);
                        }
                    }
                }
            }
        }
    }

    //set lateral grid spacing for each layer
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        layers_dx.push_back(DEM[i]->cellSizeX());
        std::cout << "layer " << i << "  dx : " <<  layers_dx[i] <<  "  dz : " <<  layers_dz[i] << std::endl;
    }

    //shielding analysis
    //sweep through in direction of wind input
    //and track for each layer shadowing effect in this direction

    //this helper function requires that there is onyl 1 distance for either i, r or c, with i2, r2 and c2
    //we assume here that the user ignores the scales between layers. So just uses r, c for r2 and c2,
    //even if the next layer is 3 times finer or courser in scale
    std::function<float(std::vector<std::shared_ptr<cTMapf>>&, int, int, int, int, int, int)> getvoxelmasked;
    getvoxelmasked = [&](std::vector<std::shared_ptr<cTMapf>> & VOX, int i, int r, int c, int i2, int r2, int c2) -> float {
        //try to get the voxel value from the struct
        //we might need to do 3x3 pixel average from layer below
        //or look at larger gridcells above
        if (i2 < 0 || i2 >= static_cast<int>(VOX.size()) ||( (i2 == i) && (r2 < 0 || r2 >= VOX[i2]->nrRows() || c2 < 0 || c2 >= VOX[i2]->nrCols()))) {
            return VOX[i]->data[r][c];
        }
        if (i2 == i-1) {

            if (layers_scaledown[i] == 3) {
                float av = 0.0f;
                float w = 0.0f;
                for (int dr = 0; dr < 3; dr++) {
                    for (int dc = 0; dc < 3; dc++) {
                        if (r2*3+dr < VOX[i2]->nrRows() && c2*3+dc < VOX[i2]->nrCols()) {
                            if (!pcr::isMV(VOX[i2]->data[r2*3+dr ][c2*3+dc])) {
                                av += (1.0f/9.0f)* VOX[i2]->data[r2*3+dr ][c2*3+dc];
                                w+= 1.0f/9.0f;
                            }
                        }
                    }
                }
                if ( w < 1e-6f) {
                    return VOX[i]->data[r][c];
                }

                return av/std::max(w,1e-6f);
            }else {
                return VOX[i2]->data[r2][c2];
            }
        }else if (i2 == i+1)
        {
            if (layers_scaledown[i2] == 3) {
                // 1. Map fine cell center to continuous coarse coordinates
                float r_c_exact = (r2 + 0.5f) / 3.0f - 0.5f;
                float c_c_exact = (c2 + 0.5f) / 3.0f - 0.5f;

                // 2. Base top-left cell for the central 2x2 of the 4x4 stencil
                int r0 = (int)std::floor(r_c_exact);
                int c0 = (int)std::floor(c_c_exact);

                // 3. Fractional weights [0.0, 1.0)
                float wr = r_c_exact - (float)r0;
                float wc = c_c_exact - (float)c0;

                // 4. Helper to fetch coarse values safely (clamped to boundaries)
                auto get_c = [&](int row, int col) {
                    int cr = std::max(0, std::min(VOX[i2]->nrRows() - 1, row));
                    int cc = std::max(0, std::min(VOX[i2]->nrCols() - 1, col));
                    return VOX[i2]->data[cr][cc];
                };

                // 5. 1D Catmull-Rom cubic interpolation helper
                auto cubic_interp = [](float p0, float p1, float p2, float p3, float t) {
                    return p1 + 0.5f * t * (p2 - p0 + t * (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3 + t * (3.0f * (p1 - p2) + p3 - p0)));
                };

                // 6. Evaluate the 4x4 grid: interpolate horizontally across the 4 rows
                float col_vals[4];
                for (int dr = -1; dr <= 2; ++dr) {
                    float p0 = get_c(r0 + dr, c0 - 1);
                    float p1 = get_c(r0 + dr, c0);
                    float p2 = get_c(r0 + dr, c0 + 1);
                    float p3 = get_c(r0 + dr, c0 + 2);
                    col_vals[dr + 1] = cubic_interp(p0, p1, p2, p3, wc);
                }

                // 7. Interpolate vertically across the 4 results
                float val = cubic_interp(col_vals[0], col_vals[1], col_vals[2], col_vals[3], wr);

                return val;

            } else {
                // Fallback for non-scaled cascades
                if (r2 < 0 || r2 >= VOX[i2]->nrRows() || c2 < 0 || c2 >= VOX[i2]->nrCols()) {
                    return VOX[i]->data[r][c];
                }
                return VOX[i2]->data[r2][c2];
            }
        }else if (i == i2){
            return VOX[i2]->data[r2][c2];
        }

        return VOX[i]->data[r][c];
    };

    //this helper function requires that there is onyl 1 distance for either i, r or c, with i2, r2 and c2
    //we assume here that the user ignores the scales between layers. So just uses r, c for r2 and c2,
    //even if the next layer is 3 times finer or courser in scale
    std::function<float(std::vector<std::shared_ptr<cTMapf>>&, int, int, int, int, int, int)> getvoxelmaskedzerob;
    getvoxelmaskedzerob = [&](std::vector<std::shared_ptr<cTMapf>> & VOX, int i, int r, int c, int i2, int r2, int c2) -> float {
        //try to get the voxel value from the struct
        //we might need to do 3x3 pixel average from layer below
        //or look at larger gridcells above
        if (i2 < 0) {
            return 0.0f;
        }

        if (i2 >= static_cast<int>(VOX.size())) {
            return getvoxelmaskedzerob(VOX, i, r, c, i - 1, r, c);
        }

        // Only check lateral out-of-bounds using direct r2, c2 if we are on the same layer
        if (i2 == i && (r2 < 0 || r2 >= VOX[i2]->nrRows() || c2 < 0 || c2 >= VOX[i2]->nrCols())) {
            return VOX[i]->data[r][c];
        }
        if (i2 == i-1) {

            if (layers_scaledown[i] == 3) {
                float av = 0.0f;
                float w = 0.0f;
                for (int dr = 0; dr < 3; dr++) {
                    for (int dc = 0; dc < 3; dc++) {
                        if (r*3+dr < VOX[i2]->nrRows() && c*3+dc < VOX[i2]->nrCols())
                            {
                            if (!pcr::isMV(VOX[i2]->data[r*3+dr ][c*3+dc])) {
                                av += (1.0f/9.0f)* VOX[i2]->data[r*3+dr ][c*3+dc];
                                w+= 1.0f/9.0f;
                            }
                        }
                    }
                }
                if ( w < 1e-6f) {
                    return VOX[i]->data[r][c];
                }
                return av/std::max(w,1e-6f);
            }else {
                return VOX[i2]->data[r2][c2];
            }
        }else if (i2 == i+1)
        {
            if (layers_scaledown[i2] == 3) {
                int r_coarse = r/3;
                int c_coarse = c/3;
                if (r_coarse < 0 || r_coarse > VOX[i2]->nrRows()-1 || c_coarse < 0 || c_coarse > VOX[i2]->nrCols()-1) {
                    return VOX[i]->data[r][c];
                }
                return VOX[i2]->data[r_coarse][c_coarse];
            }else {

                return VOX[i2]->data[r2][c2];
            }
        }else if (i == i2){
            return VOX[i2]->data[r2][c2];
        }

        return VOX[i]->data[r][c];
    };

    auto fw_slopex =[&](int i, int r, int c) {

        r = std::max(0,std::min(r,DEM[i]->nrRows()-1));
        c = std::max(0,std::min(c,DEM[i]->nrCols()-1));
        float z = DEM[i]->data[r][c];
        float z1 = z; float z2 = z;
        if (c > 0){ z1 = DEM[i]->data[r][c-1];};
        if (c < DEM[i]->nrCols()-1){ z2 = DEM[i]->data[r][c+1];};
        return std::min(FW_MAX_SLOPE,std::max(-FW_MAX_SLOPE,(z2-z1)/(2.0f*layers_dx[i])));
    };

    auto fw_slopey = [&](int i, int r, int c) {

        r = std::max(0,std::min(r,DEM[i]->nrRows()-1));
        c = std::max(0,std::min(c,DEM[i]->nrCols()-1));
        float z = DEM[i]->data[r][c];
        float z1 = z; float z2 = z;
        if (r > 0){ z1 = DEM[i]->data[r-1][c];};
        if (r < DEM[i]->nrRows()-1){ z2 = DEM[i]->data[r+1][c];};
        return std::min(FW_MAX_SLOPE,std::max(-FW_MAX_SLOPE,(z2-z1)/(2.0f*layers_dx[i])));
    };

    //pre-compute sin and cos of terrain angle
    for (int i = static_cast<int>(layers_z.size()) - 1; i > -1; i--) {
        float dx = layers_dx[i];
        cTMapf* map = DEM[i].get(); //our reference domain for this layer
        for (int r = 0; r < map->nrRows(); r++) {
            for (int c = 0; c < map->nrCols(); c++) {

                float dz_x = fw_slopex(i,r,c) * dx;
                float dz_y = fw_slopey(i,r,c) * dx;
                float total_dz = std::fabs(dz_x) + std::fabs(dz_y);
                float dist = M_SQRT2*dx;

                float inflow_lat = std::fabs(fast_cos_atanf(total_dz,dist));
                float inflow_z = std::fabs(fast_sin_atanf(total_dz,dist));

                float inflow_t = std::max(1e-6f, inflow_z + inflow_lat);
                inflow_lat /= inflow_t;
                inflow_z /= inflow_t;

                COSTF[i]->data[r][c] = inflow_lat;
                SINTF[i]->data[r][c] = inflow_z;

            }
        }
        writeRaster(*SINTF[i], "inflowz_" +to_string<int>(i) + ".tif");
    }

    std::shared_ptr<cTMapf> seperation = DEM[0]->GetCopySP0();
    std::shared_ptr<cTMapf> sx = DEM[0]->GetCopySP0();
    std::shared_ptr<cTMapf> sy = DEM[0]->GetCopySP0();
    std::shared_ptr<cTMapf> sdir = DEM[0]->GetCopySP0();

    std::cout << "do shadow sweep " << n_iters_sweepshadow << std::endl;
    // 2. Shadow sweep
    // We use terrain and buildings to look at obstacles blocking wind flow,
    // Using intiial estimated of direction through the boundary conditions
    for (int iterz = 0; iterz < n_iters_sweepshadow; iterz++) {

        for (int i = static_cast<int>(layers_z.size()) - 1; i > -1; i--) {
        //for each direaction (combine 2 adjecent directions as well in the sweep to optimize for diagonal movement at the same time)

            for (int iter = 0; iter < 1; iter++) {
                auto do_shadow_progress = [&](int r, int c, float vx, float vy, float dx) {
                    // 1. Calculate the local reference velocity magnitude
                    float U_ref = std::max(1e-3f, std::sqrt(vx * vx + vy * vy));

                    // 2. Fetch slopes and upstream shadows
                    float slope_x = fw_slopex(i, r, c);
                    float slope_y = fw_slopey(i, r, c);
                    float dz_x = slope_x * dx;
                    float dz_y = slope_y * dx;

                    float shadowfax_upstream_x = vx > 0.0f? getvoxelmasked(SHADE,i,r,c,i,r,c-1) : getvoxelmasked(SHADE,i,r,c,i,r,c+1);
                    float shadowfax_upstream_y = vy > 0.0f? getvoxelmasked(SHADE,i,r,c,i,r-1,c) : getvoxelmasked(SHADE,i,r,c,i,r+1,c);
                    float shadowfax_upstream_z_x = dz_x * vx < 0.0f? getvoxelmaskedzerob(SHADE,i,r,c,i-1,r,c) : getvoxelmasked(SHADE,i,r,c,i+1,r,c);
                    float shadowfax_upstream_z_y = dz_y * vy < 0.0f? getvoxelmaskedzerob(SHADE,i,r,c,i-1,r,c) : getvoxelmasked(SHADE,i,r,c,i+1,r,c);

                    float vx_l = -vy; float vy_l = vx;
                    float vx_r = vy;  float vy_r = -vx;
                    float shadowfax_upstream_l = (std::fabs(vx_l) *(vx_l > 0.0f? getvoxelmasked(SHADE,i,r,c,i,r,c-1) : getvoxelmasked(SHADE,i,r,c,i,r,c+1)) + std::fabs(vy_l)* (vy_l > 0.0f? getvoxelmasked(SHADE,i,r,c,i,r-1,c) : getvoxelmasked(SHADE,i,r,c,i,r+1,c)))/(std::fabs(vx_l) + std::fabs(vy_l));
                    float shadowfax_upstream_r = (std::fabs(vx_r) *(vx_r > 0.0f? getvoxelmasked(SHADE,i,r,c,i,r,c-1) : getvoxelmasked(SHADE,i,r,c,i,r,c+1)) + std::fabs(vy_r)* (vy_r >  0.0f? getvoxelmasked(SHADE,i,r,c,i,r-1,c) : getvoxelmasked(SHADE,i,r,c,i,r+1,c)))/(std::fabs(vx_r) + std::fabs(vy_r));

                    float bbfax_upstream_x = 0.5 * (vx > 0.0f? getvoxelmasked(BBF,i,r,c,i,r,c-1) : getvoxelmasked(BBF,i,r,c,i,r,c+1));
                    float bbfax_upstream_y = 0.5 * (vy > 0.0f? getvoxelmasked(BBF,i,r,c,i,r-1,c) : getvoxelmasked(BBF,i,r,c,i,r+1,c));
                    float bbfax_upstream_z_x = 0.5 * (dz_x * vx < 0.0f? getvoxelmasked(BBF,i,r,c,i-1,r,c) : getvoxelmasked(BBF,i,r,c,i+1,r,c));
                    float bbfax_upstream_z_y = 0.5 * (dz_y * vy < 0.0f? getvoxelmasked(BBF,i,r,c,i-1,r,c) : getvoxelmasked(BBF,i,r,c,i+1,r,c));

                    float inflow_lat = COSTF[i]->data[r][c];
                    float inflow_z = SINTF[i]->data[r][c];
                    float w_lat = std::max(1e-6f, std::fabs(vx) + std::fabs(vy));
                    float w_vert = std::max(1e-6f, std::fabs(dz_x) + std::fabs(dz_y));


                    // --- NEW: PHYSICALLY BASED TERRAIN SEPARATION ---
                    // Calculate the directional gradient of the terrain strictly along the wind vector
                    float directional_slope = (vx * slope_x + vy * slope_y) / U_ref;

                    // Dynamic critical slope: Lowers at high wind speeds (flow separates easier due to inertia)
                    // Base is 0.32 (~18 degrees). Decreases by 0.01 per m/s of wind speed.
                    float critical_slope = std::max(0.0f,0.25f - (0.02f * U_ref));
                    critical_slope = std::max(0.05f, std::min(0.32f, critical_slope)); // Clamp to realistic aerodynamic bounds

                    float separation_severity = (-directional_slope - critical_slope) / -directional_slope;
                    float terrain_shadow = 1.0f; // 1.0 means no shadow
                    // If slope is negative (downhill) and steeper than the critical threshold, flow separates
                    if (directional_slope < 0 && i == 0) { //-critical_slope
                        // "the fraction out of total thats over the angle should be the relative applied shading"
                        // Interpret "total" as the directional slope. This normalizes severity.
                        // Generate shadow deficit based on severity and your input coefficient
                        terrain_shadow = 1.0f - velocity_deficit_coefficient  *(1.0f- std::exp(- (separation_severity * velocity_deficit_angle_strength)));
                        terrain_shadow = std::min(1.0f,std::max(0.25f, terrain_shadow)); // Cap max wake deficit at 75%
                    }else {
                        terrain_shadow = 0.0f;
                    }

                    if (i == 0) {
                        sx->data[r][c] = slope_x;
                    }
                    if (i == 0) {
                        sy->data[r][c] = slope_y;
                    }
                    if (i == 0) {
                        sdir->data[r][c] = separation_severity;
                    }
                    if (i == 0) {
                        seperation->data[r][c] = terrain_shadow;
                    }
                    shadowfax_upstream_z_x = dz_x * vx < 0.0f? (terrain_shadow * shadowfax_upstream_x + (1.0f-terrain_shadow) * shadowfax_upstream_z_x) : shadowfax_upstream_z_x;
                    shadowfax_upstream_z_y = dz_y * vy < 0.0f? (terrain_shadow * shadowfax_upstream_y + (1.0f-terrain_shadow) * shadowfax_upstream_z_y) : shadowfax_upstream_z_y;

                    // 3. Calculate Base Advected Shadow (from buildings and upstream obstacles)
                    float shadowfax_advected = inflow_lat * (std::fabs(vx) * shadowfax_upstream_x * (1.0f-bbfax_upstream_x) + std::fabs(vy) * shadowfax_upstream_y* (1.0f-bbfax_upstream_y))/w_lat
                                             + inflow_z * (std::fabs(dz_x) * shadowfax_upstream_z_x * (1.0f-bbfax_upstream_z_x) + std::fabs(dz_y) * shadowfax_upstream_z_y* (1.0f-bbfax_upstream_z_y))/w_vert;


                    // Combine the advected shadow with any new leeward terrain separation
                    // (Takes the minimum value, as 0.0 is total shadow and 1.0 is free stream)
                    float combined_shadow = shadowfax_advected;

                    // --- NEW: PHYSICALLY BASED WAKE DIFFUSION ---
                    // Mechanical turbulence scales with velocity. Faster wind = faster wake recovery.
                    float dynamic_diffusion = wake_diffusion_rate * (1.0f + 0.1f * U_ref);
                    float fac_wakediffusive = std::exp(-(dynamic_diffusion * dx)/dx);

                    // Final Shadow Output
                    float shadowfax_new = fac_wakediffusive * combined_shadow + (1.0f - fac_wakediffusive) * 0.5f * (shadowfax_upstream_l + shadowfax_upstream_r);

                    SHADE[i]->data[r][c] = std::max(0.0f, std::min(1.0f, std::min(SHADE[i]->data[r][c], shadowfax_new)));
                };

                float dx = layers_dx[i];
                cTMapf* map = DEM[i].get(); //our reference domain for this layer
                for (int r = 0; r < map->nrRows(); r++) {
                    for (int c = 0; c < map->nrCols(); c++) {

                        int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                        int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                        //sample reference velocity here
                        float vx = VEL_HIGH_X->data[r_vel][c_vel];
                        float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                        //look in direction opposite to sweep direction
                        //if velocity doesnt have a negative dot product to sweep direction
                        if (vx > 0.0f && vy > 0.0f) {
                            do_shadow_progress(r,c,vx,vy,dx);
                        }
                    }
                }
                for (int r = map->nrRows()-1; r > -1; r--) {
                    for (int c = 0; c < map->nrCols(); c++) {
                        int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                        int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                        //sample reference velocity here
                        float vx = VEL_HIGH_X->data[r_vel][c_vel];
                        float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                        //look in direction opposite to sweep direction
                        //if velocity doesnt have a negative dot product to sweep direction
                        if (vx > 0.0f && vy <= 0.0f) {
                            do_shadow_progress(r,c,vx,vy,dx);
                        }
                    }
                }
                for (int r = 0; r < map->nrRows(); r++) {
                    for (int c = map->nrCols()-1; c > -1 ; c--) {
                        int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                        int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                        //sample reference velocity here
                        float vx = VEL_HIGH_X->data[r_vel][c_vel];
                        float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                        //look in direction opposite to sweep direction
                        //if velocity doesnt have a negative dot product to sweep direction
                        if (vx <= 0.0f && vy > 0.0f) {
                            do_shadow_progress(r,c,vx,vy,dx);
                        }
                    }
                }
                for (int r = map->nrRows()-1; r > -1; r--) {
                    for (int c = map->nrCols()-1; c > -1 ; c--) {
                        int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                        int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                        //sample reference velocity here
                        float vx = VEL_HIGH_X->data[r_vel][c_vel];
                        float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                        //look in direction opposite to sweep direction
                        //if velocity doesnt have a negative dot product to sweep direction
                        if (vx <= 0.0f && vy <= 0.0f) {
                            do_shadow_progress(r,c,vx,vy,dx);
                        }
                    }
                }
            }

        }
    }

    if (output_test)
    {
        writeRaster(*seperation,"seperation.tif");
        writeRaster(*sx,"sx.tif");
        writeRaster(*sy,"sy.tif");
        writeRaster(*sdir,"sdir.tif");
        writeRaster(*DEM[0],"sdem.tif");

        for (int i = 0; i < static_cast<int>(SHADE.size()); i++) {
            writeRaster(*SHADE[i], "shadow_" + std::to_string(i) + ".tif");
        }

        for (int i = 0; i < static_cast<int>(BBF.size()); i++) {
            writeRaster(*BBF[i], "bbf_" + std::to_string(i) + ".tif");
        }

    }

    //and estimate roughness attached to these grid scales.
    //roughness length is about one tenth of the high variations in the landscape
    //for each layer
    for (int i = static_cast<int>(layers_z.size()) - 1; i > -1; i--) {
        float dx = layers_dx[i];
        int scale = layers_pxlatscale[i];
        cTMapf* map = DEM[i].get(); //our reference domain for this layer
        for (int r = 0; r < map->nrRows(); r++) {
            for (int c = 0; c < map->nrCols(); c++) {

                // Roughness from Manning's n
                float rl_n = N[i]->data[r][c] * 25.0f; // Correlation: 0.2 Manning -> 15m roughness

                // Roughness from building height
                float rl_bh = MAXBH[i]->data[r][c] * 0.1f;

                // Roughness from sub-grid topography
                float rl_topo = 0.0f;
                float mean_dem = 0.0f;
                float mean_sq_dem = 0.0f;
                int count = 0;
                for (int dr = 0; dr < scale; ++dr) {
                    for (int dc = 0; dc < scale; ++dc) {
                        int r_fine = r * scale + dr;
                        int c_fine = c * scale + dc;
                        if (r_fine < in_DEM->nrRows() && c_fine < in_DEM->nrCols()) {
                            float dem_val = in_DEM->data[r_fine][c_fine];
                            mean_dem += dem_val;
                            mean_sq_dem += dem_val * dem_val;
                            count++;
                        }
                    }
                }
                if (count > 1) {
                    mean_dem /= count;
                    // Prevent NaN from floating point inaccuracies when variance is near zero
                    float variance = mean_sq_dem / count - mean_dem * mean_dem;
                    float std_dev = std::sqrt(std::max(0.0f, variance));
                    rl_topo = std_dev * 0.1f;
                }

                // Combine roughness lengths, capping at a reasonable maximum
                RL[i]->data[r][c] = std::min(15.0f, rl_n + rl_bh + rl_topo);
            }
        }

        // Apply 9x9 window average smoothing to roughness length to avoid blocky transitions
        cTMapf* rl_map = RL[i].get();
        cTMapf* rl_temp = rl_map->GetCopy0();
        int rows = rl_map->nrRows();
        int cols = rl_map->nrCols();
        const int rad = 4; // 9x9 window: [-4, +4]

        ParallelOrSerialFor(0, rows, [&](int r) {
            for (int c = 0; c < cols; c++) {
                float sum = 0.0f;
                int count = 0;
                for (int dc = -rad; dc <= rad; dc++) {
                    int cc = c + dc;
                    if (cc >= 0 && cc < cols) {
                        sum += rl_map->data[r][cc];
                        count++;
                    }
                }
                rl_temp->data[r][c] = sum / (float)count;
            }
        });

        ParallelOrSerialFor(0, rows, [&](int r) {
            for (int c = 0; c < cols; c++) {
                float sum = 0.0f;
                int count = 0;
                for (int dr = -rad; dr <= rad; dr++) {
                    int rr = r + dr;
                    if (rr >= 0 && rr < rows) {
                        sum += rl_temp->data[rr][c];
                        count++;
                    }
                }
                rl_map->data[r][c] = sum / (float)count;
            }
        });
        delete rl_temp;
    }
    if (output_test) {
        writeRaster(*RL[0],"rl_0.tif");
    }

    // 6. Internal Boundary Layer (IBL) Friction Propagation
    // We advect a velocity multiplier (PROF) downstream. This allows wind at 2m to
    // speed up over water and slow down over forests gradually, propagating upward.

    // Define the reference roughness length for the input wind field.
    float z0_ref_input = ref_roughness_length;
    float z0_rl_disappear = 2.5f;
    // Calculate the log term for the reference input velocity (denominator for PROF)
    float log_ref_input_denom = std::log(std::max(ref_wind_height, z0_ref_input * 1.5f) / z0_ref_input);


    ///TODO:
    //Monin-Obukhov stability corrections
    //relevant for low wind speed events,
    //where therman instability is influential for the profile shape.
    //need thermal parameters to estimate this well

    // Initialize the profile multiplier for each cell based on its local roughness length.
    // This creates a spatially varying initial state for the advection scheme.
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = PROF[i].get();
        float z = layers_z[i];

        ParallelOrSerialFor(0, map->nrRows(), [&](int r) {
            for (int c = 0; c < map->nrCols(); c++) {
                // Use the local roughness length (z0) for this cell.
                float rl = RL[i]->data[r][c];
                if ( z > z0_rl_disappear*RL[i]->data[r][c] ) {
                    rl = ref_roughness_length;
                }else {
                    float fac = std::max(0.0f,std::min(1.0f,(z)/(z0_rl_disappear*RL[i]->data[r][c])));
                    rl = fac * ref_roughness_length + (1.0f-fac)*RL[i]->data[r][c];
                }

                float z0_local = std::max(0.001f, rl);

                // Log-law numerator: based on current layer height (z) and local roughness.
                float z_eff = std::max(z, z0_local * 1.5f);
                float log_at_z_local_num = std::log(z_eff / z0_local);
                float prof_log = log_at_z_local_num / std::max(1e-6f, log_ref_input_denom);

                // Sample reference velocity magnitude to evaluate wind speed regime
                int r_vel = std::max(0, std::min(VEL_HIGH_X->nrRows() - 1, static_cast<int>(((map->north() + r * map->cellSizeY()) - VEL_HIGH_X->north()) / VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0, std::min(VEL_HIGH_X->nrCols() - 1, static_cast<int>(((map->west() + c * map->cellSizeX()) - VEL_HIGH_X->west()) / VEL_HIGH_X->cellSizeX())));
                float vx_ref = VEL_HIGH_X->data[r_vel][c_vel];
                float vy_ref = VEL_HIGH_Y->data[r_vel][c_vel];
                float u_ref = std::sqrt(vx_ref * vx_ref + vy_ref * vy_ref);

                // For hurricane wind speeds, strong atmospheric mixing flattens vertical profile
                // Linear/flatter profile: starts at 1.0 at reference height and grows mildly with height
                float prof_linear = 1.0f + 0.20f * std::max(0.0f, z - ref_wind_height) / std::max(50.0f, z + 50.0f);

                // Smoothly blend between log profile (normal winds <= 17 m/s) and flat profile (hurricane >= 33 m/s)
                float w_hurr = std::max(0.0f, std::min(1.0f, (u_ref - 17.0f) / 16.0f));
                w_hurr = w_hurr * w_hurr * (3.0f - 2.0f * w_hurr); // Smoothstep

                float prof_final = (1.0f - w_hurr) * prof_log + w_hurr * prof_linear;
                PROF[i]->data[r][c] = std::max(0.2f, std::min(2.5f, prof_final));
            }
        });

        writeRaster(*RL[i],"rl_"+to_string<int>(i)+".tif");
        writeRaster(*PROF[i],"profmult_"+to_string<int>(i)+".tif");
    }

    for (int iter = 0; iter < 2; iter++) {
        for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
            cTMapf* map = PROF[i].get();
            float dx = layers_dx[i];
            float z = layers_z[i];

            auto do_prof_progress = [&](int r, int c, float vx, float vy) {
                float U_mag = std::fabs(vx) + std::fabs(vy);

                if (U_mag < 1e-6f) return; // Prevent division by zero

                // 1. Get advected upstream multiplier
                float prof_up_x = vx > 0.0f ? getvoxelmasked(PROF, i, r, c, i, r, c - 1) : getvoxelmasked(PROF, i, r, c, i, r, c + 1);
                float prof_up_y = vy > 0.0f ? getvoxelmasked(PROF, i, r, c, i, r - 1, c) : getvoxelmasked(PROF, i, r, c, i, r + 1, c);
                float prof_advected = (std::fabs(vx) * prof_up_x + std::fabs(vy) * prof_up_y) / U_mag;

                // 2. Calculate local TARGET equilibrium multiplier based on local roughness.
                float z0_local = std::max(0.001f, RL[i]->data[r][c]);
                float z_eff = std::max(z, z0_local * 1.5f);
                float log_at_z_local_num = std::log(z_eff / z0_local);

                // Target profile multiplier: log(z / z0_local) / log(ref_wind_height / z0_ref_input)
                float target_prof = log_at_z_local_num / std::max(1e-6f, log_ref_input_denom);

                float u_ref = std::sqrt(vx * vx + vy * vy);
                float prof_linear = 1.0f + 0.20f * std::max(0.0f, z - ref_wind_height) / std::max(50.0f, z + 50.0f);
                float w_hurr = std::max(0.0f, std::min(1.0f, (u_ref - 17.0f) / 16.0f));
                w_hurr = w_hurr * w_hurr * (3.0f - 2.0f * w_hurr);
                target_prof = (1.0f - w_hurr) * target_prof + w_hurr * prof_linear;
                target_prof = std::max(0.2f, std::min(2.5f, target_prof)); // Clamp physical extremes

                // 3. Calculate relaxation factor (how fast the wind adapts to the new surface).
                // IBL grows at a ~1:10 slope. Higher layers require vastly more fetch to adjust.
                float adjustment_length = std::max(10.0f, 50.0f * z);
                float relaxation_factor = 1.0f - std::exp(-dx / adjustment_length);

                // 4. Blend the advected wind profile with the local terrain's target profile
                PROF[i]->data[r][c] = prof_advected * (1.0f - relaxation_factor) + target_prof * relaxation_factor;
            };

            // Perform sweeps in all 4 quadrants to follow the wind vector
            for (int r_s = 0; r_s < map->nrRows(); r_s++) {
                for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                    int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                    int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                    float vx = VEL_HIGH_X->data[r_vel][c_vel];
                    float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                    if (vx > 0.0f && vy > 0.0f) do_prof_progress(r_s, c_s, vx, vy);
                }
            }
            for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
                for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                    int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                    int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                    float vx = VEL_HIGH_X->data[r_vel][c_vel];
                    float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                    if (vx > 0.0f && vy < 0.0f) do_prof_progress(r_s, c_s, vx, vy);
                }
            }
            for (int r_s = 0; r_s < map->nrRows(); r_s++) {
                for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                    int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                    int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                    float vx = VEL_HIGH_X->data[r_vel][c_vel];
                    float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                    if (vx < 0.0f && vy > 0.0f) do_prof_progress(r_s, c_s, vx, vy);
                }
            }
            for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
                for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                    int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                    int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                    float vx = VEL_HIGH_X->data[r_vel][c_vel];
                    float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                    if (vx < 0.0f && vy < 0.0f) do_prof_progress(r_s, c_s, vx, vy);
                }
            }
        }
    }

    // Construct the Initial Velocity Field (V0)
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = ux_initial[i].get();
        float dx = layers_dx[i];
        ParallelOrSerialFor(0, map->nrRows(), [&](int r) {
            for (int c = 0; c < map->nrCols(); c++) {
                if (BBF[i]->data[r][c] > 0.9f) {
                    ux_initial[i]->data[r][c] = 0.0f;
                    uy_initial[i]->data[r][c] = 0.0f;
                    uz_initial[i]->data[r][c] = 0.0f;
                } else {
                    int r_vel = std::max(0, std::min(VEL_HIGH_X->nrRows() - 1, static_cast<int>(((map->north() + r * map->cellSizeY()) - VEL_HIGH_X->north()) / VEL_HIGH_X->cellSizeY())));
                    int c_vel = std::max(0, std::min(VEL_HIGH_X->nrCols() - 1, static_cast<int>(((map->west() + c * map->cellSizeX()) - VEL_HIGH_X->west()) / VEL_HIGH_X->cellSizeX())));
                    float vx_ref = VEL_HIGH_X->data[r_vel][c_vel];
                    float vy_ref = VEL_HIGH_Y->data[r_vel][c_vel];

                    float shade_factor = SHADE[i]->data[r][c];
                    float prof_factor = PROF[i]->data[r][c];

                    float u0_x = vx_ref * shade_factor * prof_factor;
                    float u0_y = vy_ref * shade_factor * prof_factor;
                    float u0_mag = std::sqrt(u0_x * u0_x + u0_y * u0_y);

                    // Smoothly limit extreme/supersonic velocities entering the solver
                    // Preserves wind direction while softly compressing unphysical extremes
                    const float V_SOFT_LIMIT = 65.0f; // m/s (~234 km/h)
                    const float V_MAX_CAP = 90.0f;    // Asymptotic upper ceiling (~324 km/h)
                    if (u0_mag > V_SOFT_LIMIT) {
                        float excess = u0_mag - V_SOFT_LIMIT;
                        float max_excess = V_MAX_CAP - V_SOFT_LIMIT;
                        float limited_mag = V_SOFT_LIMIT + max_excess * std::tanh(excess / max_excess);
                        float scale = limited_mag / u0_mag;
                        u0_x *= scale;
                        u0_y *= scale;
                    }

                    ux_initial[i]->data[r][c] = u0_x;
                    uy_initial[i]->data[r][c] = u0_y;

                    // The initial vertical velocity MUST be zero for a diagnostic model.
                    // This creates the mass imbalance that the potential solver is designed to fix.
                    uz_initial[i]->data[r][c] = 0.0f;
                }
            }
        });
    }

   //now, we want to add rotors
   //first do a sweep in the direction of the Uref (1 sweep for each cascade is enough)
   //make a map with the distance along the Uref direction of the wake
   //assume wake starts if shade < 0.5 or so (define a threshold constant)
   //in the sweep, distance is 0 if no wake, otherwise, increase by dx for every cell
   //as long as the shadow is there.
   //if building fraction is greater than 0.5, have 0 as distance (no wake)
   //and afterward it starts again
    const float WAKE_SHADOW_THRESHOLD = 0.5f;
    const float WAKE_BUILDING_THRESHOLD = 0.5f;

    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = WAKE_DIST[i].get();
        float dx = layers_dx[i];

        auto do_wake_dist_progress = [&](int r, int c, float vx, float vy) {
            if (SHADE[i]->data[r][c] < WAKE_SHADOW_THRESHOLD && BBF[i]->data[r][c] < WAKE_BUILDING_THRESHOLD) {
                float u_mag = std::fabs(vx) + std::fabs(vy);
                if (u_mag > 1e-6f) {
                    float dist_up_x = vx > 0.0f ? getvoxelmasked(WAKE_DIST, i, r, c, i, r, c - 1) : getvoxelmasked(WAKE_DIST, i, r, c, i, r, c + 1);
                    float dist_up_y = vy > 0.0f ? getvoxelmasked(WAKE_DIST, i, r, c, i, r - 1, c) : getvoxelmasked(WAKE_DIST, i, r, c, i, r + 1, c);
                    float advected_dist = (std::fabs(vx) * dist_up_x + std::fabs(vy) * dist_up_y) / u_mag;
                    WAKE_DIST[i]->data[r][c] = advected_dist + dx;
                } else {
                    WAKE_DIST[i]->data[r][c] = 0.0f;
                }
            } else {
                WAKE_DIST[i]->data[r][c] = 0.0f;
            }
        };

        // Sweeps for WAKE_DIST
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx > 0.0f && vy > 0.0f) do_wake_dist_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx > 0.0f && vy < 0.0f) do_wake_dist_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx < 0.0f && vy > 0.0f) do_wake_dist_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx < 0.0f && vy < 0.0f) do_wake_dist_progress(r_s, c_s, vx, vy);
            }
        }
    }


    //also do a sweep for lateral distance within the wakes
    //so in a sweep direction 90 degrees rotated (perpendicular)
    //do the same distance tracking with the same conditions
    //and consider building presence
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = WAKE_LAT_DIST[i].get();
        float dx = layers_dx[i];

        auto do_wake_lat_dist_progress = [&](int r, int c, float vx, float vy) {
            if (SHADE[i]->data[r][c] < WAKE_SHADOW_THRESHOLD && BBF[i]->data[r][c] < WAKE_BUILDING_THRESHOLD) {
                float pvx = -vy;
                float pvy = vx;
                float u_mag = std::fabs(pvx) + std::fabs(pvy);
                if (u_mag > 1e-6f) {
                    float dist_up_x = pvx > 0.0f ? getvoxelmasked(WAKE_LAT_DIST, i, r, c, i, r, c - 1) : getvoxelmasked(WAKE_LAT_DIST, i, r, c, i, r, c + 1);
                    float dist_up_y = pvy > 0.0f ? getvoxelmasked(WAKE_LAT_DIST, i, r, c, i, r - 1, c) : getvoxelmasked(WAKE_LAT_DIST, i, r, c, i, r + 1, c);
                    float advected_dist = (std::fabs(pvx) * dist_up_x + std::fabs(pvy) * dist_up_y) / u_mag;
                    WAKE_LAT_DIST[i]->data[r][c] = advected_dist + dx;
                } else {
                    WAKE_LAT_DIST[i]->data[r][c] = 0.0f;
                }
            } else {
                WAKE_LAT_DIST[i]->data[r][c] = 0.0f;
            }
        };

        // Sweeps for WAKE_LAT_DIST
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx > 0.0f && pvy > 0.0f) do_wake_lat_dist_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx > 0.0f && pvy < 0.0f) do_wake_lat_dist_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx < 0.0f && pvy > 0.0f) do_wake_lat_dist_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx < 0.0f && pvy < 0.0f) do_wake_lat_dist_progress(r_s, c_s, vx, vy);
            }
        }
    }


    //now sweep in the opposite direction, and create a new map
    //this way, we can
    //sweep into the wake regions the maximum distance
    //and finally we can know where we are relatively
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = WAKE_DIST_MAX[i].get();
        // Initialize WAKE_DIST_MAX with WAKE_DIST
        for(int r = 0; r < map->nrRows(); r++) {
            for(int c = 0; c < map->nrCols(); c++) {
                WAKE_DIST_MAX[i]->data[r][c] = WAKE_DIST[i]->data[r][c];
            }
        }

        auto do_wake_dist_max_progress = [&](int r, int c, float vx, float vy) {
            if (WAKE_DIST[i]->data[r][c] > 1e-6f) { // only propagate within a wake
                float u_mag = std::fabs(vx) + std::fabs(vy);
                if (u_mag > 1e-6f) {
                    float dist_down_x = vx > 0.0f ? getvoxelmasked(WAKE_DIST_MAX, i, r, c, i, r, c + 1) : getvoxelmasked(WAKE_DIST_MAX, i, r, c, i, r, c - 1);
                    float dist_down_y = vy > 0.0f ? getvoxelmasked(WAKE_DIST_MAX, i, r, c, i, r + 1, c) : getvoxelmasked(WAKE_DIST_MAX, i, r, c, i, r - 1, c);
                    float advected_max_dist = (std::fabs(vx) * dist_down_x + std::fabs(vy) * dist_down_y) / u_mag;
                    WAKE_DIST_MAX[i]->data[r][c] = std::max(WAKE_DIST_MAX[i]->data[r][c], advected_max_dist);
                }
            } else {
                WAKE_DIST_MAX[i]->data[r][c] = 0.0f;
            }
        };

        // Reverse sweeps for WAKE_DIST_MAX
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx > 0.0f && vy > 0.0f) do_wake_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx > 0.0f && vy < 0.0f) do_wake_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx < 0.0f && vy > 0.0f) do_wake_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                if (vx < 0.0f && vy < 0.0f) do_wake_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
    }


    //also do this reverse sweep for the lateral distance
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = WAKE_LAT_DIST_MAX[i].get();
        // Initialize WAKE_LAT_DIST_MAX with WAKE_LAT_DIST
        for(int r = 0; r < map->nrRows(); r++) {
            for(int c = 0; c < map->nrCols(); c++) {
                WAKE_LAT_DIST_MAX[i]->data[r][c] = WAKE_LAT_DIST[i]->data[r][c];
            }
        }

        auto do_wake_lat_dist_max_progress = [&](int r, int c, float vx, float vy) {
            if (WAKE_LAT_DIST[i]->data[r][c] > 1e-6f) { // only propagate within a wake
                float pvx = -vy;
                float pvy = vx;
                float u_mag = std::fabs(pvx) + std::fabs(pvy);
                if (u_mag > 1e-6f) {
                    float dist_down_x = pvx > 0.0f ? getvoxelmasked(WAKE_LAT_DIST_MAX, i, r, c, i, r, c + 1) : getvoxelmasked(WAKE_LAT_DIST_MAX, i, r, c, i, r, c - 1);
                    float dist_down_y = pvy > 0.0f ? getvoxelmasked(WAKE_LAT_DIST_MAX, i, r, c, i, r + 1, c) : getvoxelmasked(WAKE_LAT_DIST_MAX, i, r, c, i, r - 1, c);
                    float advected_max_dist = (std::fabs(pvx) * dist_down_x + std::fabs(pvy) * dist_down_y) / u_mag;
                    WAKE_LAT_DIST_MAX[i]->data[r][c] = std::max(WAKE_LAT_DIST_MAX[i]->data[r][c], advected_max_dist);
                }
            } else {
                WAKE_LAT_DIST_MAX[i]->data[r][c] = 0.0f;
            }
        };

        // Reverse sweeps for WAKE_LAT_DIST_MAX
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx > 0.0f && pvy > 0.0f) do_wake_lat_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx > 0.0f && pvy < 0.0f) do_wake_lat_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx < 0.0f && pvy > 0.0f) do_wake_lat_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
        for (int r_s = 0; r_s < map->nrRows(); r_s++) {
            for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                int r_vel = std::max(0,std::min(VEL_HIGH_X->nrRows()-1,static_cast<int>(((map->north() + r_s * map->cellSizeY())-VEL_HIGH_X->north())/VEL_HIGH_X->cellSizeY())));
                int c_vel = std::max(0,std::min(VEL_HIGH_X->nrCols()-1,static_cast<int>(((map->west() + c_s * map->cellSizeX())-VEL_HIGH_X->west())/VEL_HIGH_X->cellSizeX())));
                float vx = VEL_HIGH_X->data[r_vel][c_vel];
                float vy = VEL_HIGH_Y->data[r_vel][c_vel];
                float pvx = -vy; float pvy = vx;
                if (pvx < 0.0f && pvy < 0.0f) do_wake_lat_dist_max_progress(r_s, c_s, vx, vy);
            }
        }
    }



    //finally we have maps for each cascade for
    //wake space distance
    //wake space maximum distance
    //wake space lateral distance
    //wake space maximum lateral distance
    if (output_test)
    {
        for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
            writeRaster(*WAKE_DIST[i], "wakedist_" + to_string<int>(i) + ".tif");
            writeRaster(*WAKE_LAT_DIST[i], "wakelatdist_" + to_string<int>(i) + ".tif");
            writeRaster(*WAKE_DIST_MAX[i], "wakedistmax_" + to_string<int>(i) + ".tif");
            writeRaster(*WAKE_LAT_DIST_MAX[i], "wakelatdistmax_" + to_string<int>(i) + ".tif");
        }
    }


    //Run our flow equations in 3d
    //in topography following coordinates
    //with varying lateral and vertical grid cell size (both coarser as you go up
    //potential flow with turbulence tracking and energy/loss advection
    //first assume perfect flow lateral to terrain,
    //which means potential function only varies laterally
    //solve potential flow on coursest top layer, then pre-assume solution for other layers
    //but interpolated to somewhat higher resolution
    //then final solve in 3d which will require lots fewer iterations


    // 3. Initialization: Coarse-to-fine initialization
     // Initialize the potential field (PSI) for each layer directly.
     // This provides a good initial guess for the iterative solver.
     for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
         cTMapf* map = PSI[i].get();
         for (int r = 0; r < map->nrRows(); r++) {
             for (int c = 0; c < map->nrCols(); c++) {
                 PSI[i]->data[r][c] = 0.0f;
             }
         }
     }
    if (output_test) {
        writeRaster(*PSI[PSI.size()-1],"psitemp.tif");

        writeRaster(*PSI[0],"psitempb.tif");
    }

    int top_layer = static_cast<int>(PSI.size())-1;
     // 4. Main Iterative Solver: Red-Black Gauss-Seidel with Over-relaxation
     float omega = 1.7f; // Over-relaxation factor (must be > 1.0 for acceleration)

     for (int iter = 0; iter < n_iters_flow; iter++) {
         // Two passes: one for "white" cells (r+c is even), one for "black" cells (r+c is odd).
         // This allows for safe parallelization as a cell's neighbors in one pass are from the other color.
         for (int color = 0; color < 2; color++) {
             for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
                 cTMapf* map = PSI[i].get();
                 float dx = layers_dx[i];
                 float dy = layers_dx[i]; // Assuming square cells
                 float dz = layers_dz[i];

                 ParallelOrSerialFor(0, map->nrRows(), [&](int r) {
                        for (int c = 0; c < map->nrCols(); c++) {
                            if ((r + c) % 2 != color) continue; // Process only cells of the current color

                            // Pin potential at one corner to make the solution unique
                            if (i == top_layer && r == 0 && c == 0) {
                                PSI[i]->data[r][c] = 0.0f;
                                continue;
                            }

                            // Get neighbor potential values
                            float phi_e = getvoxelmasked(PSI, i, r, c, i, r, c + 1); // East
                            float phi_w = getvoxelmasked(PSI, i, r, c, i, r, c - 1); // West
                            float phi_n = getvoxelmasked(PSI, i, r, c, i, r - 1, c); // North
                            float phi_s = getvoxelmasked(PSI, i, r, c, i, r + 1, c); // South
                            float phi_u = getvoxelmasked(PSI, i, r, c, i + 1, r, c); // Up
                            float phi_d = getvoxelmasked(PSI, i, r, c, i - 1, r, c); // Down

                            // Handle lateral boundary conditions (Neumann: d(phi)/dn = 0 for correction potential)
                            if (c == 0)                   phi_w = phi_e;
                            if (c == map->nrCols() - 1)   phi_e = phi_w;
                            if (r == 0)                   phi_n = phi_s;
                            if (r == map->nrRows() - 1)   phi_s = phi_n;
                            if (i == 0)                   phi_d = PSI[i]->data[r][c];
                            if (i == top_layer)           phi_u = PSI[i]->data[r][c];

                            if (BLOCK_BUILDINGS) {
                                float bbf_e = getvoxelmasked(BBF, i, r, c, i, r, c + 1);
                                float bbf_w = getvoxelmasked(BBF, i, r, c, i, r, c - 1);
                                float bbf_n = getvoxelmasked(BBF, i, r, c, i, r - 1, c);
                                float bbf_s = getvoxelmasked(BBF, i, r, c, i, r + 1, c);
                                float bbf_u = getvoxelmasked(BBF, i, r, c, i + 1, r, c);
                                float bbf_d = getvoxelmaskedzerob(BBF, i, r, c, i - 1, r, c);

                                if (bbf_e > 0.9f) phi_e = PSI[i]->data[r][c];
                                if (bbf_w > 0.9f) phi_w = PSI[i]->data[r][c];
                                if (bbf_n > 0.9f) phi_n = PSI[i]->data[r][c];
                                if (bbf_s > 0.9f) phi_s = PSI[i]->data[r][c];
                                if (bbf_u > 0.9f) phi_u = PSI[i]->data[r][c];
                                if (bbf_d > 0.9f) phi_d = PSI[i]->data[r][c];
                            }

                            // Calculate Divergence of Initial Velocity Field (V0)
                            float u0_e = getvoxelmasked(ux_initial, i, r, c, i, r, c + 1);
                            float u0_w = getvoxelmasked(ux_initial, i, r, c, i, r, c - 1);
                            float v0_n = getvoxelmasked(uy_initial, i, r, c, i, r - 1, c);
                            float v0_s = getvoxelmasked(uy_initial, i, r, c, i, r + 1, c);
                            float w0_u = getvoxelmasked(uz_initial, i, r, c, i + 1, r, c);
                            float w0_d = getvoxelmasked(uz_initial, i, r, c, i - 1, r, c);

                            // Standard Cartesian divergence terms
                            float div_cartesian = (u0_e - u0_w) / (2.0f * dx) +      // d(u)/dx
                                                  (v0_n - v0_s) / (2.0f * dy) +      // d(v)/dy (y decreases as r increases)
                                                  (w0_u - w0_d) / (2.0f * dz);      // d(w)/dz

                            // --- Metric terms for topography-following coordinates ---
                            // These terms account for the vertical flux induced by horizontal flow over sloped terrain.
                            // This is crucial for capturing pressure effects from topography.
                            
                            // Get values from layers above and below, applying boundary conditions
                            int i_up = (i == top_layer) ? i : i + 1;
                            int i_down = (i == 0) ? i : i - 1;
                            
                            int r_up = r;
                            int c_up = c;
                            int r_down = r;
                            int c_down = c;

                            if (i > 0 && layers_scaledown[i] == 3) {
                                r_down = r * 3;
                                c_down = c * 3;
                            }
                            if (i < top_layer && layers_scaledown[i+1] == 3) {
                                r_up = r / 3;
                                c_up = c / 3;
                            }


                            // Term 1: -d(u * slope_x)/dz, with safe boundary handling
                            float slope_x_u = fw_slopex(i_up, r_up, c_up);
                            float u0_up = getvoxelmasked(ux_initial, i, r, c, i_up, r, c);
                            float slope_x_d = fw_slopex(i_down, r_down, c_down);
                            float u0_down = getvoxelmasked(ux_initial, i, r, c, i_down, r, c);
                            float metric_term_x = (u0_up * slope_x_u - u0_down * slope_x_d) / (2.0f * dz);

                            // Term 2: -d(v * slope_y)/dz, with safe boundary handling
                            float slope_y_u = fw_slopey(i_up, r_up, c_up);
                            float v0_up = getvoxelmasked(uy_initial, i, r, c, i_up, r, c);
                            float slope_y_d = fw_slopey(i_down, r_down, c_down);
                            float v0_down = getvoxelmasked(uy_initial, i, r, c, i_down, r, c);
                            float metric_term_y = (v0_up * slope_y_u - v0_down * slope_y_d) / (2.0f * dz);

                            float divergence = div_cartesian - metric_term_x - metric_term_y;

                            // --- BUILDING BLOCKAGE EFFECTS ---
                            float bbf_e = getvoxelmasked(BBF, i, r, c, i, r, c + 1);
                            float bbf_w = getvoxelmasked(BBF, i, r, c, i, r, c - 1);
                            float bbf_n = getvoxelmasked(BBF, i, r, c, i, r - 1, c);
                            float bbf_s = getvoxelmasked(BBF, i, r, c, i, r + 1, c);
                            float bbf_u = getvoxelmasked(BBF, i, r, c, i + 1, r, c);
                            float bbf_d = getvoxelmaskedzerob(BBF, i, r, c, i - 1, r, c);

                            float p_e = 1.0f, p_w = 1.0f, p_n = 1.0f, p_s = 1.0f, p_u = 1.0f, p_d = 1.0f;

                            if (BLOCK_BUILDINGS) {
                                // Enforce Neumann boundary condition: d(phi)/dn = 0 (since U_initial . n is 0 at solid walls)
                                if (bbf_e > 0.9f) {
                                    phi_e = PSI[i]->data[r][c];
                                }
                                if (bbf_w > 0.9f) {
                                    phi_w = PSI[i]->data[r][c];
                                }
                                if (bbf_n > 0.9f) {
                                    phi_n = PSI[i]->data[r][c];
                                }
                                if (bbf_s > 0.9f) {
                                    phi_s = PSI[i]->data[r][c];
                                }
                                if (bbf_u > 0.9f) phi_u = PSI[i]->data[r][c]; // No vertical flow through ceiling/roof
                                if (bbf_d > 0.9f) phi_d = PSI[i]->data[r][c]; // No vertical flow through floor
                            } else {
                                // Variable Porosity
                                float bbf_c = BBF[i]->data[r][c];
                                auto to_theta = [](float bbf) { return std::max(0.001f, 1.0f - 0.1f* bbf); };
                                float theta_c = to_theta(bbf_c);
                                float theta_e = to_theta(bbf_e);
                                float theta_w = to_theta(bbf_w);
                                float theta_n = to_theta(bbf_n);
                                float theta_s = to_theta(bbf_s);
                                float theta_u = to_theta(bbf_u);
                                float theta_d = to_theta(bbf_d);

                                auto h_mean = [](float a, float b) {
                                    if (a + b < 1e-9f) return 0.0f;
                                    return 2.0f * a * b / (a + b);
                                };
                                p_e = h_mean(theta_c, theta_e);
                                p_w = h_mean(theta_c, theta_w);
                                p_n = h_mean(theta_c, theta_n);
                                p_s = h_mean(theta_c, theta_s);
                                p_u = h_mean(theta_c, theta_u);
                                p_d = h_mean(theta_c, theta_d);
                            }

                            // --- Transformed Laplacian for Terrain-Following Coordinates ---
                            // The standard Laplacian weights (w_x, w_y, w_z) must be modified
                            // to account for the non-orthogonal grid.

                            float w_x = 1.0f / (dx * dx);
                            float w_y = 1.0f / (dy * dy);
                            // The vertical weight, scaled by alpha_sq, balances horizontal vs. vertical flow.
                            float w_z = alpha_sq * 1.0f / (dz * dz);

                            // Get local slopes for the metric terms in the Laplacian
                            float slope_x = fw_slopex(i, r, c);
                            float slope_y = fw_slopey(i, r, c);

                            // --- Cross-derivative terms using only axial neighbors ---
                            // Approximate diagonal phi values by averaging axial neighbors.
                            // This respects the limitations of getvoxelmasked.
                            float phi_ue = 0.5f * (phi_u + phi_e); float phi_uw = 0.5f * (phi_u + phi_w);
                            float phi_de = 0.5f * (phi_d + phi_e); float phi_dw = 0.5f * (phi_d + phi_w);
                            float phi_un = 0.5f * (phi_u + phi_n); float phi_us = 0.5f * (phi_u + phi_s);
                            float phi_dn = 0.5f * (phi_d + phi_n); float phi_ds = 0.5f * (phi_d + phi_s);

                            // Add metric terms to the numerator and denominator
                            // The full Laplacian term for the cross-derivative is -2 * slope_x * d2(phi)/dx'dz'
                            // Our finite difference approximation for d2(phi)/dx'dz' is (phi_ue - phi_uw - phi_de + phi_dw) / (4*dx*dz).
                            // The sign is inverted here so it can be cleanly subtracted in the numerator.
                            float cross_term_x = 2.0f * slope_x * (phi_ue - phi_uw - phi_de + phi_dw) / (4.0f * dx * dz);
                            float cross_term_y = 2.0f * slope_y * (phi_un - phi_us - phi_dn + phi_ds) / (4.0f * dy * dz);

                            float numer = w_x * (p_e * phi_e + p_w * phi_w) +
                                          w_y * (p_n * phi_n + p_s * phi_s) +
                                          w_z * (p_u * phi_u + p_d * phi_d) * (1.0f + slope_x * slope_x + slope_y * slope_y) -
                                          cross_term_x - cross_term_y - divergence;

                            float denom = (p_e + p_w) * w_x + (p_n + p_s) * w_y + (p_u + p_d) * w_z * (1.0f + slope_x * slope_x + slope_y * slope_y);

                            float phi_new = numer / denom;

                            //forces stability but limits output variability
                            //float max_neighbor = std::max({phi_e, phi_w, phi_n, phi_s, phi_u, phi_d});
                            //float min_neighbor = std::min({phi_e, phi_w, phi_n, phi_s, phi_u, phi_d});
                            //phi_new = std::max(min_neighbor, std::min(max_neighbor, phi_new));

                            if (std::isnormal(phi_new)) {
                                // This clamping step was preventing the solver from creating necessary pressure peaks (local maxima/minima).
                                // Removing it is essential for the solver to function correctly.
                               PSI[i]->data[r][c] = (1.0f - omega) * PSI[i]->data[r][c] + omega * phi_new;
                            }
                        }
                    });
                 if ( iter < 5 && color == 1) {
                     if (output_test) {
                         writeRaster(*PSI[i], "psi_bottom_" +to_string<int>(i) + "_"+ to_string<int>(iter) + ".tif");
                     }
                 }
             }
         }

     }
    if (output_test) {
        for (int i = 0; i < static_cast<int>(PSI.size());i++)
        {
            writeRaster(*PSI[i], "psi_bottom_" + to_string<int>(i) +"_final.tif");
        }
    }
    // 5. Final Velocity Recovery
     for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
         cTMapf* map = PSI[i].get();
         float dx = layers_dx[i];
         float dy = layers_dx[i]; // Assuming square cells
         float dz = layers_dz[i];

         ParallelOrSerialFor(0, map->nrRows(), [&](int r) {
            for (int c = 0; c < map->nrCols(); c++) {
                float phi_e = getvoxelmasked(PSI, i, r, c, i, r, c + 1);
                float phi_w = getvoxelmasked(PSI, i, r, c, i, r, c - 1);
                float phi_n = getvoxelmasked(PSI, i, r, c, i, r - 1, c);
                float phi_s = getvoxelmasked(PSI, i, r, c, i, r + 1, c);
                float phi_u = getvoxelmasked(PSI, i, r, c, i + 1, r, c);
                float phi_d = getvoxelmasked(PSI, i, r, c, i - 1, r, c);

                if (c == 0)                   phi_w = phi_e;
                if (c == map->nrCols() - 1)   phi_e = phi_w;
                if (r == 0)                   phi_n = phi_s;
                if (r == map->nrRows() - 1)   phi_s = phi_n;
                if (i == 0)                   phi_d = phi_u;
                if (i == top_layer)           phi_u = phi_d;

                // --- Transformed Gradient for Velocity Recovery ---
                // The gradient operator must also be transformed to match the solver's coordinate system.
                float slope_x = fw_slopex(i, r, c);
                float slope_y = fw_slopey(i, r, c);

                // Partial derivatives of phi in computational space (x', y', z')
                float dphi_dx_prime = (phi_e - phi_w) / (2.0f * dx);
                float dphi_dy_prime = (phi_n - phi_s) / (2.0f * dy); // Corrected sign for y-axis (y decreases as r increases)
                float dphi_dz_prime = (phi_u - phi_d) / (2.0f * dz);

                // Physical velocity corrections include metric terms
                float grad_phi_x = dphi_dx_prime - slope_x * dphi_dz_prime;
                float grad_phi_y = dphi_dy_prime - slope_y * dphi_dz_prime;
                // The physical vertical velocity (w) is the sum of the initial guess (0) and the correction.
                // The correction for w in terrain-following coordinates is more complex than just dphi_dz_prime,
                // but for the final velocity, we add the full terrain-following effect.

                ux[i]->data[r][c] = ux_initial[i]->data[r][c] - grad_phi_x;
                uy[i]->data[r][c] = uy_initial[i]->data[r][c] - grad_phi_y;
                // The final vertical velocity is the terrain-following component plus the pressure-driven vertical correction.
                uz[i]->data[r][c] = (ux[i]->data[r][c] * slope_x + uy[i]->data[r][c] * slope_y) - dphi_dz_prime;

                if (BBF[i]->data[r][c] > 0.9f) {
                    ux[i]->data[r][c] = 0.0f;
                    uy[i]->data[r][c] = 0.0f;
                    uz[i]->data[r][c] = 0.0f;
                }
            }
        });
     }

    // Turbulence (TKE) Production, Propagation, and Dissipation
    // This is done after the main velocity field is calculated.

    // 1. TKE Production
    // Turbulence is primarily produced by shear and interaction with surface roughness.
    // Here, we'll use a simplified production model where TKE is generated based on
    // the local roughness length and the square of the wind speed.
    const float C_mu = 0.09f; // Empirical constant for TKE production

    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = TKE[i].get();
        ParallelOrSerialFor(0, map->nrRows(), [&](int r) {
            for (int c = 0; c < map->nrCols(); c++) {
                // --- TKE Production from Shear ---
                // TKE_prod = l^2 * (dU/dz)^2, where l is mixing length and dU/dz is shear.

                // 1. Calculate vertical wind shear (dU/dz)
                float vel_mag_up = 0.0f;
                if (i < static_cast<int>(layers_z.size()) - 1) {
                    float ux_up = getvoxelmasked(ux, i, r, c, i + 1, r, c);
                    float uy_up = getvoxelmasked(uy, i, r, c, i + 1, r, c);
                    vel_mag_up = std::sqrt(ux_up * ux_up + uy_up * uy_up);
                }

                float vel_mag_current = std::sqrt(ux[i]->data[r][c] * ux[i]->data[r][c] + uy[i]->data[r][c] * uy[i]->data[r][c]);
                float delta_U = std::fabs(vel_mag_up - vel_mag_current);
                float delta_Z = (i < static_cast<int>(layers_z.size()) - 1) ? (layers_z[i+1] - layers_z[i]) : layers_dz[i];
                float shear = delta_U / std::max(1e-3f, delta_Z);

                // 2. Calculate mixing length (l)
                // l = k * z in the surface layer, where k is Von Karman constant.
                // We cap it to avoid excessive values high up.
                float z_eff = layers_z[i]; // Effective height for this layer
                float mixing_length = von_karman_constant * z_eff;

                // 3. Calculate TKE production
                // The formula is C_mu * l^2 * (dU/dz)^2, which is dimensionally correct (m^2/s^2).
                float tke_prod = C_mu * mixing_length * mixing_length * shear * shear;
                TKE[i]->data[r][c] = tke_prod;
            }
        });
    }

    // 2. TKE Propagation and Dissipation (Sweep)
    // Similar to the shadow sweep, we advect and decay TKE downstream, driven by the calculated wind field.
    for (int iter = 0; iter < n_iters_turb; iter++) {
        for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
            cTMapf* map = TKE[i].get();
            float dx = layers_dx[i];

            auto do_tke_progress = [&](int r, int c) {
                float vx = ux[i]->data[r][c];
                float vy = uy[i]->data[r][c];

                // Get upstream TKE values from lateral and vertical neighbors
                float tke_up_x = vx > 0.0f ? getvoxelmasked(TKE, i, r, c, i, r, c - 1) : getvoxelmasked(TKE, i, r, c, i, r, c + 1);
                float tke_up_y = vy > 0.0f ? getvoxelmasked(TKE, i, r, c, i, r - 1, c) : getvoxelmasked(TKE, i, r, c, i, r + 1, c);
                float tke_up_z = getvoxelmaskedzerob(TKE, i, r, c, i - 1, r, c); // Turbulence generally comes from below

                // Weighted average of upstream TKE based on wind direction
                float total_vel = std::fabs(vx) + std::fabs(vy);
                if (total_vel < 1e-6f) {
                    return; // Skip advection if velocity is zero to prevent division by zero
                }
                float tke_advected = (std::fabs(vx) * tke_up_x + std::fabs(vy) * tke_up_y) / std::max(1e-6f, total_vel);

                // Add vertical contribution (mixing from below)
                float vert_mix_factor = 0.1f;
                tke_advected = (1.0f - vert_mix_factor) * tke_advected + vert_mix_factor * tke_up_z;

                // Dissipation: TKE decays over distance. Use a similar exponential decay as the wake.
                float decay_factor = std::exp(-wake_diffusion_rate * dx * 0.5f); // Use a slightly different rate for TKE
                float tke_new = tke_advected * decay_factor;

                // The current cell's TKE is the max of its own production and the advected/dissipated TKE from upstream.
                TKE[i]->data[r][c] = std::max(TKE[i]->data[r][c], tke_new);
            };

            // Perform sweeps in all 4 quadrants
            for (int r_s = 0; r_s < map->nrRows(); r_s++) { // Top-left to bottom-right
                for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                    if (ux[i]->data[r_s][c_s] > 0.0f && uy[i]->data[r_s][c_s] > 0.0f) do_tke_progress(r_s, c_s);
                }
            }
            for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) { // Bottom-left to top-right
                for (int c_s = 0; c_s < map->nrCols(); c_s++) {
                    if (ux[i]->data[r_s][c_s] > 0.0f && uy[i]->data[r_s][c_s] < 0.0f) do_tke_progress(r_s, c_s);
                }
            }
            for (int r_s = 0; r_s < map->nrRows(); r_s++) { // Top-right to bottom-left
                for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                    if (ux[i]->data[r_s][c_s] < 0.0f && uy[i]->data[r_s][c_s] > 0.0f) do_tke_progress(r_s, c_s);
                }
            }
            for (int r_s = map->nrRows() - 1; r_s > -1; r_s--) { // Bottom-right to top-left
                for (int c_s = map->nrCols() - 1; c_s > -1; c_s--) {
                    if (ux[i]->data[r_s][c_s] < 0.0f && uy[i]->data[r_s][c_s] < 0.0f) do_tke_progress(r_s, c_s);
                }
            }
        }
    }

     // 7. Estimate Dynamic Pressure
     // From Bernoulli's principle, dynamic pressure P = 0.5 * rho * |U|^2
     // We assume rho (air density) is constant, so we calculate P/rho.
     float rho_air = 1.225f; // Standard air density at sea level, kg/m^3
     out_P.clear();
     for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
         cTMapf* map = ux[i].get();
         out_P.push_back(map->GetCopySP0());
         ParallelOrSerialFor(0, map->nrRows(), [&](int r) {
             for (int c = 0; c < map->nrCols(); c++) {
                 float vel_x = ux[i]->data[r][c];
                 float vel_y = uy[i]->data[r][c];
                 float vel_z = uz[i]->data[r][c];
                 float vel_sq = vel_x * vel_x + vel_y * vel_y + vel_z * vel_z;
                 out_P[i]->data[r][c] = 0.5f * rho_air * vel_sq;
             }
         });
     }


    // 8. Estimate Wind Gusts from TKE
    // Gust speed = Mean speed + n * sigma, where sigma is the std dev of velocity fluctuations.
    // TKE = 0.5 * (sigma_u^2 + sigma_v^2 + sigma_w^2). Assuming isotropic turbulence, TKE ~= 1.5 * sigma_u^2.
    // So, sigma_u ~= sqrt(TKE / 1.5). We'll use a gust factor of ~2-3.
    out_gust.clear();
    for (int i = 0; i < static_cast<int>(layers_z.size()); i++) {
        cTMapf* map = ux[i].get();
        out_gust.push_back(map->GetCopySP0());
        ParallelOrSerialFor(0, map->nrRows(), [&](int r) {
            for (int c = 0; c < map->nrCols(); c++) {
                float mean_vel_mag = std::sqrt(ux[i]->data[r][c] * ux[i]->data[r][c] + uy[i]->data[r][c] * uy[i]->data[r][c] + uz[i]->data[r][c] * uz[i]->data[r][c]);
                float tke_val = TKE[i]->data[r][c];
                float sigma = std::sqrt(std::max(0.0f, tke_val / 1.5f));
                float gust_addition = 2.5f * sigma; // Gust factor of 2.5
                out_gust[i]->data[r][c] = mean_vel_mag + gust_addition;
            }
        });
    }




     //and output the final grids
     out_velx.clear();
     out_vely.clear();
     out_velz.clear();
     for(int i = 0; i < static_cast<int>(layers_z.size()); i++)
     {
         out_velx.push_back(ux[i]);
         out_vely.push_back(uy[i]);
         out_velz.push_back(uz[i]);
     }

     return 0;

}



#endif