 /*@@
   @file      IVP_ParamCheck.F90
   @date      Sun Jul  7 15:29:20 2002
   @author    Ian Hawke
   @desc 
   Check some parameters.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "cctk_Functions.h"

subroutine IVP_ParamCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  if (.not.CCTK_EQUALS(metric_type, "physical")) then
    call CCTK_WARN(1, "Currently the IVP solver will only do physical metrics.")
  end if

  if (spatial_order /= 2) then
    call CCTK_WARN(0, "Currently the IVP solver will work properly on boundaries only with second order spatial differencing (in ADMMacros).")
  end if

end subroutine IVP_ParamCheck
