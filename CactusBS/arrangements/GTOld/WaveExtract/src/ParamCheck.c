 /*@@
   @file      ParamCheck.c
   @date      April 26 2002
   @author    Gabrielle Allen
   @desc 
              Check parameters for WaveExtract
   @enddesc 
   @version $Id: ParamCheck.c,v 1.1.1.1 2004/06/10 15:08:38 herrmann Exp $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header: /numrelcvs/HerrmannCVS/WaveExtract/src/ParamCheck.c,v 1.1.1.1 2004/06/10 15:08:38 herrmann Exp $";

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

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>3)
    CCTK_INFO("Checking parameters");

  if(! (CCTK_EQUALS(metric_type, "physical") || 
        CCTK_EQUALS(metric_type, "static conformal")))
  {
    CCTK_PARAMWARN("WaveExtract only works currently with metric_type \"static conformal\" or \"physical\"");
  }

  if (out_offset != 0) {
    CCTK_WARN(1,"out_offset is deprecated - assuming 0");
  }

}
