/***********************************************************************
/
/  DETERMINES WHETHER THE STAR PARTICLE IS RADIATIVE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
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

bool Star::IsARadiationSource(FLOAT Time)
{

  int i;
  bool *rules, result = true;

  /* To add rules, you must also modify NumberOfRules. */

  const int NumberOfRules = 4;
  rules = new bool[NumberOfRules];

  for (i = 0; i < NumberOfRules; i++) 
    rules[i] = false; 

  /*******************************************************************
     Below are the multiple definitions for a radiation source.  If
     all of the rules are met, the star particle is a radiation
     source. 
  ********************************************************************/

  // Particles only marked for nothing or continuous supernova
  rules[0] = (FeedbackFlag == NO_FEEDBACK || 
	      FeedbackFlag == CONT_SUPERNOVA ||
	      FeedbackFlag == MBH_THERMAL ||
	      FeedbackFlag == MBH_JETS);
  
  // Living
  rules[1] = (Time >= BirthTime && Time <= BirthTime+LifeTime && type > 0);

  // Non-zero BH accretion (usually accretion_rate[] here is NULL - Ji-hoon Kim Sep.2009)
  // KH 2021, naccretions here is always equal to one. This may be a bug or a "developing" function. 
  if ((type == BlackHole || type == MBH) && naccretions > 0)
    rules[2] = (accretion_rate[0] > tiny_number); 
  else
    rules[2] = true;

  // Non-zero mass
  rules[3] = (Mass > tiny_number);

  /******************** END RULES ********************/
  // KH Leong : Check BirthTime and LifeTime
  // fprintf(stderr, "KH Function = %s, ", __FUNCTION__);
  // fprintf(stderr, "ID = %lld, ", Identifier);
  // fprintf(stderr, "star_type = %ld, ", type);
  // fprintf(stderr, "MBH = %ld, ", MBH);
  // fprintf(stderr, "rule_0 = %d, rule_1 = %d, rule_2 = %d, rule_3 = %d, ", rules[0], rules[1], rules[2], rules[3]);
  // fprintf(stderr, "naccretions = %d, ", naccretions);
  // if (naccretions > 0) fprintf(stderr, "accretion_rate[0] = %g, ", accretion_rate[0]);
  // fprintf(stderr, "last_accretion_rate = %g, ", last_accretion_rate);
  // fprintf(stderr, "Time = %"FSYM", ", Time);
  // fprintf(stderr, "BirthTime = %"FSYM", LifeTime = %"FSYM".\n", BirthTime, LifeTime);

  for (i = 0; i < NumberOfRules; i++)
    result &= rules[i];

  delete [] rules;

  return result;

}
