 /*@@
   @file      NullExact_MoLRegister.c
   @date      Tue Jul 1 6:18:52:41 CET 2008
   @author    Maria Babiuc
   @desc 
   Routine to register the variables with the MoL thorn.  Based on Ian Hawke's WaveMoL
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void NullExact_MoLRegister(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  group = CCTK_GroupIndex("NullExact::News_MoL");
  rhs = CCTK_GroupIndex("NullExact::dotNews_MoL");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }


  if (ierr) CCTK_WARN(0,"Problems registering with MoL");

  /* return ierr; */
}
