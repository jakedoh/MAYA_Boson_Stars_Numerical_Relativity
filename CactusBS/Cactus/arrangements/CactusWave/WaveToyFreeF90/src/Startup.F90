 /*@@
   @file      Startup.F90
   @date      
   @author    Gabrielle Allen
   @desc 
              Register banner 
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

integer function WaveToyFreeF90_Startup ()

  implicit none
   
  integer ierr
  
  call CCTK_RegisterBanner(ierr, "WaveToyFreeF90: Evolutions of a Scalar Field")

  WaveToyFreeF90_Startup = 0

end function WaveToyFreeF90_Startup
