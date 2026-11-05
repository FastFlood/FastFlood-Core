//
// Created by basti on 31/01/2026.
//

#ifndef BACKWATER_H
#define BACKWATER_H


inline cTMapf * BackWaterEffect(cTMapf * DEM, cTMapf * Q, cTMapf* GradX1,  cTMapf* GradX2,  cTMapf* GradY1,  cTMapf* GradY2)
{
  //detect all of the places with backwater effect,
  //partial obstruction of flow cross sectional area
  //either due to bridge, ocean boundary condition or terrain itself.
  //then calculate required water height at that location,
  //move upslope and apply same partial blocking reduced by the
  //slope, as well as the volume loss within that step

  //bridges are input, easy to place herer
  //coastal boundary is input, we need map of water depth from that one.
  //terrain based is harder. this woudl require us to do the cascaded grid steps
  //but we need this output to do that one.
  //some other trick for terrain needed




}


#endif //BACKWATER_H
