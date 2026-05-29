/*
   @file      CheckParameters.F77
   @date      
   @author    Gabrielle Allen
   @desc 
              Check parameters for the wave equation initial data
   @enddesc 
 @@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 /*@@
   @routine    IDScalarWave_CheckParameters
   @date       
   @author     Gabrielle Allen
   @desc 
               Check parameters for the wave equation initial data
   @enddesc 
   @calls      
   @calledby   
   @history 
   @hdate Mon Oct 11 11:49:21 1999 @hauthor Tom Goodale
   @hdesc Converted to C++ 
   @endhistory 

@@*/

extern "C" void IDScalarWaveCXX_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(CCTK_Equals(initial_data,"box"))
  {
    if (!CCTK_Equals(type, "box"))
    {
      CCTK_PARAMWARN("Must have a box grid with box initial data");
    }

    if (kx == 0 || ky == 0 || kz == 0)
    {
      CCTK_PARAMWARN("Cannot have zero kx,ky,kz for box initial data");
    }
  }
  return;

  
  

}

