 /*@@
   @file    InitSymBound.F77
   @date    
   @author  Gabrielle Allen
   @desc
            Sets the symmetries for Wave Toy
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    WaveToyF77_InitSymBound
   @date      
   @author     Gabrielle Allen
   @desc 
               Sets the symmetries for Wave Toy
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

      subroutine WaveToyF77_InitSymBound(CCTK_ARGUMENTS)

      implicit none
      
      DECLARE_CCTK_ARGUMENTS 
      
      INTEGER one
      PARAMETER (one = 1)
      INTEGER ierr
      INTEGER sym(3)

      sym(1) = one
      sym(2) = one
      sym(3) = one

      call SetCartSymVN(ierr, cctkGH, sym,'wavetoy::phi')

      return
      end 

