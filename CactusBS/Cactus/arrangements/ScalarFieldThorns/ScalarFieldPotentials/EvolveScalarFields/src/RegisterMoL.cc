/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void EvolveScalarFields_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_EvolveScalarFields_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(EvolveScalarFields_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("EvolveScalarFields::phia"),  CCTK_VarIndex("EvolveScalarFields::phiarhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("EvolveScalarFields::phib"),  CCTK_VarIndex("EvolveScalarFields::phibrhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("EvolveScalarFields::pia"),  CCTK_VarIndex("EvolveScalarFields::piarhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("EvolveScalarFields::pib"),  CCTK_VarIndex("EvolveScalarFields::pibrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
