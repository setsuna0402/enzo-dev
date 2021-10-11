#define DEBUG 1
#define GHOST_DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber

/***********************************************************************
/
/  GRID CLASS (FINALIZE THE PHOTO-RATES BY DIVIDING BY NUMBER OF PARTICLES)
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: 
/
************************************************************************/


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
#include "fortran.def"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::FinalizeRadiationFields(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

#ifdef TRANSFER

  int i, j, k, index, dim;
  float CellVolume = 1;

  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0];

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  /* Find radiative transfer fields. */

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits; 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  float DensityConversion = DensityUnits / mh;
  float factor = DensityConversion * CellVolume;

  // 13/09/2021: KH modify the ionisation field in ghost and see what will change.

#if GHOST_DEBUG
/*
  for (k = GridStartIndex[2] - 3; k < GridStartIndex[2]; k++)
    for (j = GridStartIndex[1]; j < GridStartIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	BaryonField[kphHINum][index] /= factor * BaryonField[HINum][index];
	BaryonField[gammaNum][index] /= factor * BaryonField[HINum][index]; //divide by N_HI = n_HI*(dx)^3
      } // ENDFOR i
    } // ENDFOR j
*/



  int index_test, temp_k, temp_j, temp_i; 
  int temp_mark_i, temp_mark_j, temp_mark_k; 
  // for (k = GridStartIndex[2] - 3; k <= GridEndIndex[2] + 3; k++) 
  for (k = 0; k < GridDimension[2]; k++) {
    temp_mark_k = 0;
    temp_k = k;
    if (k < GridStartIndex[2]) {
      temp_k = GridStartIndex[2];
      temp_mark_k = 1;
    }
    if (k > GridEndIndex[2]  ) {
      temp_k = GridEndIndex[2]  ;
      temp_mark_k = 1;
    }
    // for (j = GridStartIndex[1] - 3; j <= GridEndIndex[1] + 3; j++) 
    for (j = 0; j < GridDimension[1]; j++) {
      temp_mark_j = 0;
      temp_j = j;
      if (j < GridStartIndex[1]) {
        temp_j = GridStartIndex[1];
        temp_mark_j = 1;
      }
      if (j > GridEndIndex[1]  ) {
        temp_j = GridEndIndex[1]  ;
        temp_mark_j = 1;
      }
      //if (temp_mark == 1) {
      //  index_test = GRIDINDEX_NOGHOST(GridStartIndex[0] - 3, temp_j, temp_k);
      //} 

      // index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      // for (i = GridStartIndex[0] - 3; i <= GridEndIndex[0] + 3; i++) 
      for (i = 0; i < GridDimension[0]; i++) {
        index = GRIDINDEX_NOGHOST(i, j, k);
        temp_mark_i = 0;
        temp_i = i;
        if (i < GridStartIndex[0]) {
          temp_i = GridStartIndex[0];
          temp_mark_i = 1;
        }
        if (i > GridEndIndex[0]  ) {
          temp_i = GridEndIndex[0]  ;
          temp_mark_i = 1;
        }
/*
        if (temp_mark_i == 1 || temp_mark_j == 1 || temp_mark_k == 1) {
          index_test = GRIDINDEX_NOGHOST(temp_i, temp_j, temp_k); 
          BaryonField[kphHINum  ][index] = BaryonField[kphHINum  ][index_test];
	      BaryonField[gammaNum  ][index] = BaryonField[gammaNum  ][index_test]; 
          BaryonField[kphHeINum ][index] = BaryonField[kphHeINum ][index_test]; 
          BaryonField[kphHeIINum][index] = BaryonField[kphHeIINum][index_test];
        }
*/       
	    BaryonField[kphHINum  ][index] /= factor * BaryonField[HINum][index];
	    BaryonField[gammaNum  ][index] /= factor * BaryonField[HINum][index];
    	BaryonField[kphHeINum ][index] /= 0.25 * factor * BaryonField[HeINum ][index];
	    BaryonField[kphHeIINum][index] /= 0.25 * factor * BaryonField[HeIINum][index];
        if (temp_mark_i == 1 || temp_mark_j == 1 || temp_mark_k == 1) {
          index_test = GRIDINDEX_NOGHOST(temp_i, temp_j, temp_k); 
          BaryonField[kphHINum  ][index] = BaryonField[kphHINum  ][index_test];
	      BaryonField[gammaNum  ][index] = BaryonField[gammaNum  ][index_test]; 
          BaryonField[kphHeINum ][index] = BaryonField[kphHeINum ][index_test]; 
          BaryonField[kphHeIINum][index] = BaryonField[kphHeIINum][index_test];
          BaryonField[HINum     ][index] = BaryonField[HINum     ][index_test];
          BaryonField[HeINum    ][index] = BaryonField[HeINum    ][index_test];
          BaryonField[HeIINum   ][index] = BaryonField[HeIINum   ][index_test];
        }
     } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k
#endif

  // 13/09/2021: original codes.
#if DEBUG

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	BaryonField[kphHINum][index] /= factor * BaryonField[HINum][index];
	BaryonField[gammaNum][index] /= factor * BaryonField[HINum][index]; //divide by N_HI = n_HI*(dx)^3
      } // ENDFOR i
    } // ENDFOR j

  if (RadiativeTransferHydrogenOnly == FALSE)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  BaryonField[kphHeINum][index] /= 
	    0.25 * factor * BaryonField[HeINum][index];
	  BaryonField[kphHeIINum][index] /= 
	    0.25 * factor * BaryonField[HeIINum][index];
	} // ENDFOR i
      } // ENDFOR j
#endif  
   if (MultiSpecies > 1 && !RadiativeTransferFLD)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  if(!RadiativeTransferUseH2Shielding) {
	    BaryonField[kdissH2INum][index] /= 
	     1.0 * factor * BaryonField[H2INum][index];
	  }
	  BaryonField[kphHMNum][index] /= 
	    1.0 * factor * BaryonField[HMNum][index];
	  BaryonField[kdissH2IINum][index] /= 
	    1.0 * factor * BaryonField[H2IINum][index];
	} // ENDFOR i
      } // ENDFOR j

   if(this->IndexOfMaximumkph >= 0)
     this->MaximumkphIfront /= (factor * BaryonField[HINum][IndexOfMaximumkph]);

#endif /* TRANSFER */  
  
  return SUCCESS;
}
