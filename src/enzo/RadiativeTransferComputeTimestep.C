/***********************************************************************
/
/  PREPARE THE RADIATIVE TRANSFER MODULE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: October 2021, Edward Leong (add HeIII Ionisation restriction)
/
/ PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "Star.h"
#include "CommunicationUtilities.h"
#include "phys_constants.h"

#define MAX_NUMBER_CROSS_CELL 0.5

extern int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
//int LastTimestepUseHII = FALSE;
float LastPhotonDT[2] = {-1,-1};

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level)
{

  if (RadiativeTransferAdaptiveTimestep == FALSE && 
      dtPhoton != FLOAT_UNDEFINED)
    return SUCCESS;

  const int MaxStepsPerHydroStep = 8;
  float PhotonCourantFactor;

  switch (RadiativeTransferAdaptiveTimestep) {
  case 1:  PhotonCourantFactor = 1.0;  break;
  case 2:  PhotonCourantFactor = 0.2;  break;
  default: PhotonCourantFactor = 1.0;
  }

  // Restrict the increase in dtPhoton to this factor
  const float MaxDTChange = 30.0;
  // Restrict the decrease in dtPhoton to this factor
  const float MaxDTDecrease = 1e-10;

  LevelHierarchyEntry *Temp;
  bool InitialTimestep;
  int l, maxLevel, ncells_rad;
  FLOAT HydroTime;
  float ThisPhotonDT;
  const float unchangedLimit = 0.1*PhotonCourantFactor*huge_number;
  const float lowerLimit = MaxDTDecrease * LastPhotonDT[0];

  // Search for the maximum level with radiation
  maxLevel = -1;
  for (l = MAX_DEPTH_OF_HIERARCHY-1; l >= 0; l--) {
    ncells_rad = 0;
    for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) {
      ncells_rad = Temp->GridData->CountRadiationCells();
      if (ncells_rad > 0) break;
    }
    if (ncells_rad > 0) {
      maxLevel = l;
      break;
    }
  }
  maxLevel = CommunicationMaxValue(maxLevel);
  if (debug)
    fprintf(stdout, "EvolvePhotons: Maximum level with radiation = %d\n", maxLevel);
  // If no radiation, find the maximum level
  if (maxLevel < 0) {
    for (l = MAX_DEPTH_OF_HIERARCHY-1; l >= 0; l--) {
      if (LevelArray[l] != NULL) {
	maxLevel = l;
	break;
      }
    }
  }

  // Determine if this is the first timestep (not in restart)
  InitialTimestep = true;
  for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++)
    if (LevelCycleCount[l] > 0) {
      InitialTimestep = false;
      break;
    }

  dtPhoton = 10*huge_number;

  float AvgLastTimestep;
  FLOAT TimeNow = LevelArray[level]->GridData->ReturnTime();
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
	   &VelocityUnits, TimeNow);

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(TimeNow, &a, &dadt);
  float afloat = float(a);

// KH 2021/8/3: Original code sets a lowerLimit for dtPhoton based on last step. 
// It is generally fine. However, dtPhoton could be overestimated at the first-call of this function,
// causing the lowerLimit may be so high. (When the luminosity of the radiation source is QSO-like.)
// Here, I set the first two timesteps according to factor * cellwidth / lightspeed.
// Copy the relalavt code from grid::ComputePhotonTimestep()
// RadiativeTransferTimestepVelocityLimit in km/s
  float dx_level, dx_ratio;
  float dtPhotonSafety = 10*huge_number;
  float LightSpeed;
  LightSpeed = RadiativeTransferPropagationSpeedFraction * (clight/VelocityUnits);

/*
  if (RadiativeTransferTimestepVelocityLevel >= 0) {
    dx_level = TopGridDx[0] * POW(RefineBy, -RadiativeTransferTimestepVelocityLevel);
  }
*/

  if (RadiativeTransferHIIRestrictedTimestep || RadiativeTransferAdaptiveTimestep == 2 ||
      RadiativeTransferHeIIIRestrictedTimestep) {
    if (LastPhotonDT[0] < 0 || LastPhotonDT[1] < 0) {
        for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++) {
            for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) {
                ThisPhotonDT = afloat * MAX_NUMBER_CROSS_CELL * 
                Temp->GridData->GetCellWidth(0, 0) / LightSpeed;
/*
                 if (RadiativeTransferTimestepVelocityLevel >= 0) {
                    dx_ratio = dx_level / Temp->GridData->GetCellWidth(0, 0);
                    if (dx_ratio > 1)
                    ThisPhotonDT *= dx_ratio;
                 }
*/
                dtPhotonSafety = min(dtPhotonSafety, ThisPhotonDT);
            }
        }
        if (LastPhotonDT[0] > 0) dtPhotonSafety *= 2.0;
        dtPhoton = min(dtPhoton, dtPhotonSafety);
        fprintf(stderr, "dtPhotonSafety = %g ", dtPhotonSafety);
        fprintf(stderr, "MetaData->GlobalMaximumkphIfront = %g \n", MetaData->GlobalMaximumkphIfront);
        fprintf(stderr, "MetaData->GlobalMaximumkpHeIIfront = %g \n", MetaData->GlobalMaximumkpHeIIfront);
    }

    // Calculate timestep by limiting to a max change in HII
    if (RadiativeTransferHIIRestrictedTimestep) {
        for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++) {
	        for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) {
	            ThisPhotonDT = Temp->GridData->
	            ComputePhotonTimestepHII(DensityUnits, LengthUnits, VelocityUnits, 
				                        afloat, MetaData->GlobalMaximumkphIfront);
	            if (ThisPhotonDT > lowerLimit)
	            dtPhoton = min(dtPhoton, ThisPhotonDT);
	        }
        }
    }

    // KH 2021/10/17 :
    // modify RadiativeTransferHIIRestrictedTimestep -> RadiativeTransferHeIIIRestrictedTimestep
    // Later
    // Calculate timestep by limiting to a max change in HeIII
    if (RadiativeTransferHeIIIRestrictedTimestep) {
        for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++) {
	        for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) {
	            ThisPhotonDT = Temp->GridData->
	            ComputePhotonTimestepHeIII(DensityUnits, LengthUnits, VelocityUnits, 
				                        afloat, MetaData->GlobalMaximumkpHeIIfront);
	            if (ThisPhotonDT > lowerLimit)
	            dtPhoton = min(dtPhoton, ThisPhotonDT);
	        }
        }
    }

  // Calculate timestep by limiting to a max change in intensity
  // (proportional to the I-front speed)
    if (RadiativeTransferAdaptiveTimestep == 2)
      for (l = 0; l < MAX_DEPTH_OF_HIERARCHY-1; l++)
	for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) {
	  ThisPhotonDT = Temp->GridData->
	    ComputePhotonTimestepTau(DensityUnits, LengthUnits, VelocityUnits, 
				     afloat);
	  if (ThisPhotonDT > lowerLimit)
	    dtPhoton = min(dtPhoton, ThisPhotonDT);
	}

    dtPhoton = PhotonCourantFactor * CommunicationMinValue(dtPhoton);

    /* Use the average because the minimum ionization timescale can
       fluctuate significantly.  It gets even worse if the dtPhoton is
       allowed to vary a lot (>factor of a few). */

    if (debug)
      printf("dtPhoton=%g, LastPhotonDT=%g %g, dtHydro/dtPhoton=%g\n",
	     dtPhoton, LastPhotonDT[0], LastPhotonDT[1],
	     dtLevelAbove/dtPhoton); 

    if (RadiativeTransferHIIRestrictedTimestep || RadiativeTransferHeIIIRestrictedTimestep)
      if (LastPhotonDT[0] > 0 && LastPhotonDT[1] > 0 && 
	  dtPhoton < unchangedLimit) {
	AvgLastTimestep = sqrt(LastPhotonDT[0] * LastPhotonDT[1]);
	if (dtPhoton > MaxDTChange * AvgLastTimestep)
	  dtPhoton = MaxDTChange * AvgLastTimestep;
      }

  // KH 03/08/2021: testing shorter dtPhoton
  //  if (dtPhotonSafety >= unchangedLimit) dtPhoton = 0.1 * dtPhoton;
// KH 2021/8/3: Only update LastPhotonDT if GlobalMaximumkphIfront > tiny_number
    // if (dtPhoton < unchangedLimit) 
    if ((dtPhoton < unchangedLimit) && 
        (MetaData->GlobalMaximumkphIfront > tiny_number || MetaData->GlobalMaximumkpHeIIfront > tiny_number  || dtPhotonSafety >= unchangedLimit )) {
      // Store dtPhoton before modifying it based on the next topgrid timestep
      LastPhotonDT[1] = LastPhotonDT[0];
      LastPhotonDT[0] = dtPhoton;  
    }

  } // ENDIF
  
  // if we didn't find any cells that restrict timestep or the option
  // isn't requested, use hydro timestep on finest level
  if (dtPhoton >= unchangedLimit) {
    for (Temp = LevelArray[maxLevel]; Temp; Temp = Temp->NextGridThisLevel) {
      ThisPhotonDT = Temp->GridData->ComputePhotonTimestep();
      dtPhoton = min(dtPhoton, ThisPhotonDT);
    } // ENDFOR grids
    dtPhoton = CommunicationMinValue(dtPhoton);

    // Ensure that not too many photon timesteps are taken per hydro step
    HydroTime = LevelArray[maxLevel]->GridData->ReturnTime();
    dtPhoton = max(dtPhoton, 
		   (HydroTime - PhotonTime) / MaxStepsPerHydroStep);
    //LastTimestepUseHII = FALSE;
  } // ENDIF
  // KH 04/07/2021: testing shorter dtPhoton
  // dtPhoton = 0.1 * dtPhoton;
  /* Do not go past the level-0 time (+timestep) */

  float Saved_dtPhoton = dtPhoton;
  HydroTime = LevelArray[0]->GridData->ReturnTime();
  if (level == 0)
    HydroTime += LevelArray[0]->GridData->ReturnTimeStep();
  float dtTol = PFLOAT_EPSILON * HydroTime;
  if ((HydroTime+dtTol - PhotonTime) < dtPhoton) {
    if (debug) 
      printf("HydroTime = %"PSYM", PhotonTime = %"PSYM
	     ", dtPhoton = %g, dtPhoton0 = %g\n",
	     HydroTime, PhotonTime, dtPhoton, Saved_dtPhoton);
    dtPhoton = min(1.01 * (HydroTime - PhotonTime), Saved_dtPhoton);
    dtPhoton = max(dtPhoton, 1e-4*Saved_dtPhoton);
    //LastTimestepUseHII = FALSE;
  }

  //LastTimestepUseHII = CommunicationMaxValue(LastTimestepUseHII);

  //if (InitialTimestep && !MetaData->FirstTimestepAfterRestart)
  //  dtPhoton = min(dtPhoton, dtLevelAbove);

  if (RadiativeTransferAdaptiveTimestep == FALSE && debug)
    printf("RadiativeTransfer: Setting dtPhoton = %g = %g years\n",
	   dtPhoton, dtPhoton*TimeUnits/yr_s);
  
  return SUCCESS;

}
