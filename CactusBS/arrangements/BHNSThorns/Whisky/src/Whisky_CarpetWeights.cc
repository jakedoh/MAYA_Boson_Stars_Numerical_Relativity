 /*@@
   @file     Whisky_CarpetWeights.cc
   @date     Wed Aug 20 23:08:28 2008
   @author   Luca Baiotti
   @desc 
   Copy the variable CarpetReduce::weight in a Whisky variable.
   Written with the help of Thomas.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "SpaceMask.h"

#include <stdio.h>
#include <stdlib.h>


#ifdef __cplusplus
  extern "C" {
#endif
    
    /* Scheduled functions */
    void Whisky_CarpetWeights(CCTK_ARGUMENTS);

#ifdef __cplusplus
  } /* extern "C" */
#endif


void Whisky_CarpetWeights(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int i;

  if (CCTK_IsThornActive("CarpetReduce")!=0){
    // Copy CarpetReduce::weight to Whisky::Whisky_Carpet
    const CCTK_REAL *CarpetWeights = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,"CarpetReduce::weight"));
    for (i=0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
      {
	Whisky_CarpetWeights[i] = CarpetWeights[i];
      }
  }
  else{
    // If CarpetReduce isn't active, set all weights to unity
    CCTK_INFO("We are not using Carpet, so I set Whisky_CarpetWeights equal to 1");
    for (i=0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
      {
	Whisky_CarpetWeights[i] = 1; //if we do not use Carpet then we set this variable to be equal to 1
      }
  }

  return;
}
