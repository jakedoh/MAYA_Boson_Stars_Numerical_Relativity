/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_WaveToy_CL_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_WaveToy_CL_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_WaveToy_CL_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToy_CL::rho"),  CCTK_VarIndex("ML_WaveToy_CL::rhorhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_WaveToy_CL::u"),  CCTK_VarIndex("ML_WaveToy_CL::urhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
