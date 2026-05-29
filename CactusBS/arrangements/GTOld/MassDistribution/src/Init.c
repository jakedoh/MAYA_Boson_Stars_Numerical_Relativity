#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "MassDistribution.h"

/*******************************
 ****** External Routines ******
 *******************************/
void MassDistribution_Init(CCTK_ARGUMENTS);

void MassDistribution_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* incl. total mass bin in initialization */
  for (int i=0; i<intvolumes*TOTAL_NUMBER_OF_BINS; i++) {
	MassPerSpecificEnergy[i] = -1e234;
  } 
  for (int i=0; i<intvolumes; i++) {
    HistogramMinimum[i] = +1e235;
    HistogramMaximum[i] = -1e235;
    file_created[i] = 0;
  }

  *have_work_to_do = 0;
}
