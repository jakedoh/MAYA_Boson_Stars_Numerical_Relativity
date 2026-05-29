 /*@@
   @file    InitSymBound.c
   @date    
   @author  Gabrielle Allen
   @desc
            Sets the symmetries for Wave Toy
   @enddesc
   @version $Header$
 @@*/


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "Symmetry.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusWave_WaveToyC_InitSymBound_c)

void WaveToyC_InitSymBound(CCTK_ARGUMENTS);

 /*@@
   @routine    WaveToyC_InitSymBound
   @date      
   @author     Gabrielle Allen
   @desc 
               Sets the symmetries for Wave Toy
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void WaveToyC_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH, sym,"wavetoy::phi");

  /* automatic boundary condition applied by driver */
  if (CCTK_IsFunctionAliased("Driver_SelectVarForBC")) {
    const char *bctype;

    bctype = NULL;
    if (CCTK_EQUALS(bound,"flat") || CCTK_EQUALS(bound,"static") ||
        CCTK_EQUALS(bound,"radiation") || CCTK_EQUALS(bound,"robin") ||
        CCTK_EQUALS(bound,"none"))
    {
      bctype = bound;
    }
    else if (CCTK_EQUALS(bound,"zero"))
    {
      bctype = "scalar";
    }

    /* Uses all default arguments, so invalid table handle -1 can be passed */
    if (bctype && Driver_SelectVarForBC (cctkGH, CCTK_ALL_FACES, 1, -1,
                                         "wavetoy::phi", bctype) < 0)
    {
      CCTK_ERROR("Failed to register bound BC for wavemol::scalarevolvemol_scalar!");
    }
  }

  return;
}
