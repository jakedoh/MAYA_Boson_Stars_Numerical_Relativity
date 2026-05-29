 /*@@
   @file    InitSymBound.F
   @date    
   @author  Gabrielle Allen
   @desc
            Sets the symmetries across the coordinate axes
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    WaveToyF90_InitSymBound
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

      subroutine WaveToyF90_InitSymBound(CCTK_ARGUMENTS)

      implicit none
      
      DECLARE_CCTK_ARGUMENTS 

      INTEGER :: ierr
      INTEGER, DIMENSION(3) :: sym = 1

      call SetCartSymVN(ierr, cctkGH, sym, 'wavetoy::phi')

      end subroutine WaveToyF90_InitSymbound

