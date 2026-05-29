 /*@@
   @file    InitSymBound.c
   @date    
   @author  Gabrielle Allen
   @desc
            Sets the symmetries for Wave Toy
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "Symmetry.h"

 /*@@
   @routine    WaveToyC_InitSymBound
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

extern "C" void WaveToyCXX_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
    
  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH, sym,"wavetoy::phi");

  return;
}
