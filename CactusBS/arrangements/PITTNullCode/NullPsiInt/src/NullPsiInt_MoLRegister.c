 /*@@
   @file      NullPsiInt_MoLRegister.c
   @date      Wed Oct 21 11:52:41 CET 2009
   @author    Maria Babiuc
   @desc 
   Routine to register the variables with the MoL thorn.  Based on Ian Hawke's WaveMoL
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void NullPsiInt_MoLRegister(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  group = CCTK_GroupIndex("NullPsiInt::re_PsiInt");
  rhs = CCTK_GroupIndex("NullPsiInt::re_dotNewsB");

  if (CCTK_IsFunctionAliased("MoLRegisterEvolvedGroup"))
  {
    ierr += MoLRegisterEvolvedGroup(group, rhs);
  }
  else
  {
    CCTK_WARN(0, "MoL function not aliased");
    ierr++;
  }

  group = CCTK_GroupIndex("NullPsiInt::im_PsiInt");
  rhs = CCTK_GroupIndex("NullPsiInt::im_dotNewsB");

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

  /* return ierr;*/
}
