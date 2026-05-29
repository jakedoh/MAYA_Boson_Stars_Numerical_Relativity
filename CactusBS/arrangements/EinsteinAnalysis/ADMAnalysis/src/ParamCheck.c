 /*@@
   @file      ParamCheck.c
   @date      Thu Apr 25 19:02:51 2002
   @author    Tom Goodale
   @desc 
   Parameter checking stuff for ADMAnalysis
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMAnalysis_ParamCheck_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/
void ADMAnalysis_ParamCheck(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    ADMAnalysis_ParamCheck
   @date       Thu Apr 25 19:04:06 2002
   @author     Tom Goodale
   @desc 
   Scheduled routine to detect invalid parameter settings.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
 @@*/
void ADMAnalysis_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_ADMAnalysis_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  if(! CCTK_EQUALS(metric_type, "physical") &&
     ! CCTK_EQUALS(metric_type, "static conformal"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\" and \"static conformal\"");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
