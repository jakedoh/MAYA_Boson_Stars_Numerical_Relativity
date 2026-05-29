 /*@@
   @file      ParamCheck.c
   @date      April 26 2002
   @author    Gabrielle Allen
   @desc 
      Check parameters for Extract
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_Extract_ParamCheck_c)

void Extract_ParamCheck(CCTK_ARGUMENTS);

 /*@@
   @routine    Extract_ParamCheck
   @date       April 26 2002
   @author     Gabrielle Allen
   @desc 
      Check parameters for Extract
   @enddesc 
   @calls     
   @history  
 
   @endhistory 

@@*/

void Extract_ParamCheck(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(! (CCTK_EQUALS(metric_type, "physical") || 
        CCTK_EQUALS(metric_type, "static conformal")))
  {
    CCTK_PARAMWARN("Extract only works currently with metric_type \"static conformal\" or \"physical\"");
  }
}
