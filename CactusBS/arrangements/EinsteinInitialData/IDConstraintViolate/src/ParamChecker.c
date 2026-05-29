 /*@@
   @file      ParamChecker.F
   @date      May 2002
   @author    Gabrielle Allen
   @desc 
      Check parameters for constraint violating initial data
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(AEIDevelopment_IDConstraintViolate_ParamChecker_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void IDConstraintViolate_ParamChecker(CCTK_ARGUMENTS);

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
   @routine    IDConstraintViolate_ParamChecker
   @date       May 2002
   @author     Gabrielle Allen
   @desc 
      Check parameters for constraint violating initial data
   @enddesc 
   @calls     
@@*/

void IDConstraintViolate_ParamChecker(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  /* Do we know how to deal with this type of metric ? */
  if(! CCTK_Equals(metric_type, "physical"))
  {
    CCTK_VParamWarn(CCTK_THORNSTRING,
		    "Unknown metric_type (%s): Known types: \"physical\"",
		    metric_type);
  }

  if (CCTK_Equals(initial_data,"constraint violating gaussian"))
  {
    CCTK_INFO("Metric set to (constraint violating) gaussian");
    CCTK_INFO("Extrinsic curvature set to zero");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
