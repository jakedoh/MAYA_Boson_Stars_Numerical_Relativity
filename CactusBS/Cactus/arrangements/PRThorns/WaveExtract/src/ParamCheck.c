 /*@@
   @file      ParamCheck.c
   @date      April 26 2002
   @author    Gabrielle Allen
   @desc 
              Check parameters for WaveExtract
   @enddesc 
   @version $Id: ParamCheck.c 57 2009-04-23 18:02:39Z schnetter $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(HerrmannCVS_WaveExtract_ParamCheck_c)

void WavExtr_ParamCheck(CCTK_ARGUMENTS);

 /*@@
   @routine    WavExtr_ParamCheck
   @date       April 26 2002
   @author     Gabrielle Allen
   @desc 
               Check parameters for WaveExtract
   @enddesc 
   @calls     
@@*/
void WavExtr_ParamCheck(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose>3)
    CCTK_INFO("Checking parameters");

  if(! (CCTK_EQUALS(metric_type, "physical") || 
        CCTK_EQUALS(metric_type, "static conformal")))
  {
    CCTK_PARAMWARN("WaveExtract only works currently with metric_type \"static conformal\" or \"physical\"");
  }
}
