 /*@@
   @file      ReconstructPPM_ParamCheck.F90
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
   @routine    ReconstructPPM_ParamCheck
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

subroutine ReconstructPPM_ParamCheck(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (CCTK_EQUALS(ppm_flatten,"stencil_3").and.(whisky_stencil < 3)) then
    call CCTK_PARAMWARN("The stencil size must be at least 3 to use PPM reconstruction with the stencil-3 variant of the flattening procedure")
  else if (CCTK_EQUALS(ppm_flatten,"stencil_4").and.(whisky_stencil < 4)) then
    call CCTK_PARAMWARN("The stencil size must be at least 4 to use PPM reconstruction with the stencil-4 (original) flattening procedure")
  end if

end subroutine ReconstructPPM_ParamCheck
