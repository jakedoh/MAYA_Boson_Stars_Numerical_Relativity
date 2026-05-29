 /*@@
   @file    InitSymBound.F90
   @date    
   @author  Gabrielle Allen
   @desc
            Sets the symmetries across the coordinate axes
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    WaveToyFreeF90_InitSymBound
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

subroutine WaveToyFreeF90_InitSymBound(CCTK_ARGUMENTS)
  
  implicit none
  
  DECLARE_CCTK_ARGUMENTS 

  INTEGER :: ierr  
  INTEGER, DIMENSION(3) :: sym = 1
  
  call SetCartSymVN(ierr, cctkGH, sym, 'wavetoy::phi')
  
end subroutine WaveToyFreeF90_InitSymbound
