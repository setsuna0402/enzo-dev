/***********************************************************************
/
/  Calculate the effect of ionizing radiation on the species
/  HI, HeI and HeII
/
/  written by: John Wise
/  date:     
/  modified1: John Regan
/             Moved to its own function and file from
/             WalkPhotonPackage
/  modified2: Ka Hou Leong
              
/
/  PURPOSE: Calculate the ionisations and photon absorbed by HI, HeI and HeII
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"

#define DEVCODE 0
int grid::RadiativeTransferIonization(PhotonPackageEntry **PP, FLOAT *dPi, int cellindex, 
				      int species, float tau, FLOAT photonrate, 
				      FLOAT *excessrate, float geo_correction,
				      const int *kphNum, int gammaNum)
{
  FLOAT dP1 = 0.0;
#if DEVCODE
  dPi[species] = (*PP)->Photons*(1-expf(-tau));
#else
  // at most use all photons for photo-ionizations
  if (tau > 29.0) //Completely Optically Thick  original: 20
    dPi[species] = (1.0+BFLOAT_EPSILON) * (*PP)->Photons;
  else if (tau > 1.e-4) //Exponential decline in photons
    dPi[species] = min((*PP)->Photons*(-expm1(-tau)), (*PP)->Photons);
    // dPi[species] = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
  else //Optically thin case
    dPi[species] = min((*PP)->Photons*tau, (*PP)->Photons);
#endif
  
  //dP1 is the number of absorptions
  dP1 = dPi[species] * geo_correction;
  // KH 2021/7/30 : try to use the whole photonpackage to ionising  
  // dP1 = (*PP)->Photons * geo_correction;

  // contributions to the photoionization rate is over whole timestep
  // Units = 1/(LengthUnits^3)*1/CodeTime
  // BaryonField[kphNum[species]] needs to be normalised
  // see Grid_FinalizeRadiationField.C
  BaryonField[kphNum[species]][cellindex] += dP1*photonrate;


  // the heating rate is just the number of photo ionizations (1/(LengthUnits^3))
  // times the excess energy units here are eV/CodeTime.
  // Units = Ev per time [Ev/TimeUnits/(LengthUnits^3)]
  // BaryonField[gammaNum] needs to be normalised
  // see Grid_FinalizeRadiationField.C
  BaryonField[gammaNum][cellindex] += dP1*excessrate[species];
/*
   * This scheme creates numerical errors
   * When Ionisation rate and density of neutral atoms are both very small
   * The result of normalisation (Ion rate = Ion rate / number density) 
   * becomes unpredictable
   * Also, setting ionisation rate to tiny_number, in fact, magnify the ionisation rate.
   * Notice that, if the radiation source is very dim and the system is neutral,
   * This scheme may be useful. 
*/
/*
#if !DEVCODE
    
   // Check to make sure we are not just dealing with very small numbers 
   // that could cause problems later on
   
  if(BaryonField[kphNum[species]][cellindex] < tiny_number) 
    BaryonField[kphNum[species]][cellindex] = tiny_number;
  if(BaryonField[gammaNum][cellindex] < tiny_number) 
    BaryonField[gammaNum][cellindex] = tiny_number;

#endif
*/
  return SUCCESS;
}
