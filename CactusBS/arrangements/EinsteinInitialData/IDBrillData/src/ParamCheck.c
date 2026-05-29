 /*@@
   @file      ParamCheck.c
   @date      Sun May  5 19:38:51 CEST 2002
   @author    David Rideout
   @desc 
   Check the parameters for IDBrillData
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_IDBrillData_ParamCheck_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void IDBrillData_ParamChecker(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

void IDBrillData_ParamChecker(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  /* Do we know how to deal with this type of metric ? */
  if( ! CCTK_Equals(metric_type, "physical") && 
      ! CCTK_Equals(metric_type, "static conformal"))
  {
    CCTK_VPARAMWARN("Unknown ADMBase::metric_type (%s): "
                    "known types are \"physical\" and \"static conformal\"",
                    metric_type);
  }

  if (CCTK_Equals(metric_type, "static conformal"))
  {
    if (!CCTK_Equals(conformal_storage,"factor"))
    {
      CCTK_PARAMWARN("BrillData only sets the conformal factor, not its derivatives. (This could easily be changed ... please ask.)");
    }
  }

  CCTK_INFO("Setting up Brill data");
  if (CCTK_Equals(metric_type, "static conformal"))
  {
     CCTK_VInfo(CCTK_THORNSTRING,"  ... using trivial conformal %s",conformal_storage);
  }
  else if (CCTK_Equals(metric_type, "physical"))
  {
    CCTK_INFO("  ... using physical metric");
  }

  if (CCTK_Equals(initial_data,"brilldata"))
  {
    CCTK_INFO("  ... using full 3D solver");
  }
  else if (CCTK_Equals(initial_data,"brilldata2D"))
  {
    CCTK_INFO("  ... using reduced 2D solver");
    if (brill3d_d != 0.0)
    {
      CCTK_WARN(0,"But 3D parameter brill3d_d set for the 2D solver !?!");
    }
  }

  if (CCTK_Equals(q_function,"exp"))
  {
    CCTK_INFO("  ... Exponential brill function");
  }
  else if (CCTK_Equals(q_function,"eppley"))
  {
    CCTK_INFO("  ... generalised Eppley Brill function");
  }
  else if (CCTK_Equals(q_function,"exp"))
  {
    CCTK_INFO("  ... Gundlach-Holz Brill function");
  }

}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
