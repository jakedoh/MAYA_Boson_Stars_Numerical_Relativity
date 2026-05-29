 /*@@
   @file      Startup.F77
   @date      July 2 2000
   @author    Gabrielle Allen
   @desc 
              Register banner for WaveToy2DF77
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

      integer function WaveToy2DF77_Startup()

      implicit none
   
      integer ierr

      call CCTK_RegisterBanner(ierr, 
     &	"WaveToy2DF77: Evolution of a 2D Scalar Field")

      WaveToy2DF77_Startup = ierr

      end 
