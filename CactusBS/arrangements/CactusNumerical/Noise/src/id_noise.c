/*
  IDRandom

  Denis Pollney <pollney@aei.mpg.de>
  14 February 2002

  Put random data into the evolution variables for robust stability
  testing a-la-winicour.

  $Header$
*/

#include <stdlib.h>

#include "noise.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void
add_noise_to_var (int idx, const char* optstring, void* cctkGH)
{
  DECLARE_CCTK_PARAMETERS;
  int i, j, k, ijk;
  CCTK_REAL* data;
  cGH* GH = cctkGH;

  CCTK_VInfo(CCTK_THORNSTRING, "Adding initial data noise to %s",
		  CCTK_VarName(idx));

  data = (CCTK_REAL*) CCTK_VarDataPtrI(GH, 0, idx);

  for (k=1; k< GH->cctk_lsh[2]-1; ++k)
    {
      for (j=1; j< GH->cctk_lsh[1]-1; ++j)
	{
	  for (i=1; i< GH->cctk_lsh[0]-1; ++i)
	    {
	      ijk = CCTK_GFINDEX3D(GH, i, j, k);

	      data[ijk] += RAND_VAL;
      	    }
	}
    }
}


void
id_noise(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (CCTK_TraverseString(id_vars, add_noise_to_var, cctkGH,
			  CCTK_GROUP_OR_VAR) < 0)
    {
      CCTK_WARN (1, "Failed to parse 'Noise::id_vars' parameter");
    }
}
