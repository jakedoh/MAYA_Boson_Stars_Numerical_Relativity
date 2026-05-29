/*
   @file      CheckParameters.c
   @date      
   @author    Gabrielle Allen
   @desc 
              Check parameters for the wave equation initial data
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusExamples_IDWaveMoL_CheckParameters_c)

void IDWaveMoL_CheckParameters(CCTK_ARGUMENTS);

 /*@@
   @routine    IDScalarWaveC_CheckParameters
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
   @hdate Thu Feb 17 09:20:41 2000 @hauthor Tom Goodale
   @hdesc Converted to C
   @endhistory 

@@*/

void IDWaveMoL_CheckParameters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_IDWaveMoL_CheckParameters;
  DECLARE_CCTK_PARAMETERS;

}

