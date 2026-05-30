/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_CCZ4_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ML_CCZ4_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ML_CCZ4_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::At11"),  CCTK_VarIndex("ML_CCZ4::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::At12"),  CCTK_VarIndex("ML_CCZ4::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::At13"),  CCTK_VarIndex("ML_CCZ4::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::At22"),  CCTK_VarIndex("ML_CCZ4::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::At23"),  CCTK_VarIndex("ML_CCZ4::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::At33"),  CCTK_VarIndex("ML_CCZ4::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::A"),  CCTK_VarIndex("ML_CCZ4::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::B1"),  CCTK_VarIndex("ML_CCZ4::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::B2"),  CCTK_VarIndex("ML_CCZ4::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::B3"),  CCTK_VarIndex("ML_CCZ4::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::Xt1"),  CCTK_VarIndex("ML_CCZ4::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::Xt2"),  CCTK_VarIndex("ML_CCZ4::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::Xt3"),  CCTK_VarIndex("ML_CCZ4::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::alpha"),  CCTK_VarIndex("ML_CCZ4::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::phi"),  CCTK_VarIndex("ML_CCZ4::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::gt11"),  CCTK_VarIndex("ML_CCZ4::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::gt12"),  CCTK_VarIndex("ML_CCZ4::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::gt13"),  CCTK_VarIndex("ML_CCZ4::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::gt22"),  CCTK_VarIndex("ML_CCZ4::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::gt23"),  CCTK_VarIndex("ML_CCZ4::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::gt33"),  CCTK_VarIndex("ML_CCZ4::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::beta1"),  CCTK_VarIndex("ML_CCZ4::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::beta2"),  CCTK_VarIndex("ML_CCZ4::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::beta3"),  CCTK_VarIndex("ML_CCZ4::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::Theta"),  CCTK_VarIndex("ML_CCZ4::Thetarhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_CCZ4::trK"),  CCTK_VarIndex("ML_CCZ4::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
