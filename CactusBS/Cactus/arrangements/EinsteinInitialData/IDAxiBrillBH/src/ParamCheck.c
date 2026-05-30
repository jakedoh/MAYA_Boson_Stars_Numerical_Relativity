 /*@@
   @file      ParamCheck.c
   @date      Thu May  2 19:23:47 CEST 2002
   @author    David Rideout
   @desc 
   Check the parameters for IDAxiBrillBH
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_IDAxiBrillBH_ParamCheck_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void IDAxiBrillBH_ParamChecker(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

void IDAxiBrillBH_ParamChecker(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_IDAxiBrillBH_ParamChecker
  DECLARE_CCTK_PARAMETERS

  /* Do we know how to deal with this type of metric ? */
  if( ! CCTK_EQUALS(metric_type, "static conformal") &&
      ! CCTK_EQUALS(metric_type, "physical"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\" and \"static conformal\"");
  }

  /* Report on parameters */
  
  CCTK_INFO("Initial data will be axisymmetric BH+Brill Wave");
  CCTK_VInfo(CCTK_THORNSTRING,"  ... wave amplitude: %f",amp);
  CCTK_VInfo(CCTK_THORNSTRING,"  ... wave center (in eta coords): %f",eta0);
  CCTK_VInfo(CCTK_THORNSTRING,"  ... wave sigma: %f",sigma);
  CCTK_VInfo(CCTK_THORNSTRING,"  ... wave power of sin theta: %d",(int)n);
  CCTK_VInfo(CCTK_THORNSTRING,"  ... outer edge of eta grid: %f",etamax);

}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
