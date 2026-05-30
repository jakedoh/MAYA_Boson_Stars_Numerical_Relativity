 /*@@
   @file      TVD_ParamCheck.F90
   @date      Sat Feb  9 23:48:01 2002
   @author    
   @desc 
   Parameter checking routine.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    TVD_ParamCheck
   @date       Sat Feb  9 23:48:43 2002
   @author     Ian Hawke
   @desc 
   Checks the parameters.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine TVD_ParamCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (whisky_stencil < 2) then
    call CCTK_PARAMWARN("The stencil size must be at least 2 to use TVD reconstruction")
  end if
  
end subroutine TVD_ParamCheck
