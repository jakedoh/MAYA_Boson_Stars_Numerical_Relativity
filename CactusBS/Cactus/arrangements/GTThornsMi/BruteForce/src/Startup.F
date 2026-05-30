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

      integer function WaveToyF90_Startup ()

      implicit none
   
      integer ierr
      call CCTK_RegisterBanner(ierr, "WaveToyF90: Evolutions of a Scalar Field")
      WaveToyF90_Startup = 0

      end function WaveToyF90_Startup
