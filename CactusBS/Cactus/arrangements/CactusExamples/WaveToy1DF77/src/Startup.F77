 /*@@
   @file      Startup.F77
   @date      June 13 2000
   @author    Gabrielle Allen
   @desc 
              Register banner for WaveToy1DF77
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"

      integer function WaveToy1DF77_Startup()

      implicit none
   
      integer ierr

      call CCTK_RegisterBanner(ierr, 
     &	"WaveToy1DF77: Evolution of a 1D Scalar Field")

      WaveToy1DF77_Startup = ierr

      end 
