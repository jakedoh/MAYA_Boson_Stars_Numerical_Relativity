/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void Kranc2BSSNChi_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::A11"),  CCTK_VarIndex("Kranc2BSSNChi::A11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::A21"),  CCTK_VarIndex("Kranc2BSSNChi::A21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::A31"),  CCTK_VarIndex("Kranc2BSSNChi::A31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::A22"),  CCTK_VarIndex("Kranc2BSSNChi::A22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::A32"),  CCTK_VarIndex("Kranc2BSSNChi::A32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::A33"),  CCTK_VarIndex("Kranc2BSSNChi::A33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::alpha"),  CCTK_VarIndex("Kranc2BSSNChi::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::beta1"),  CCTK_VarIndex("Kranc2BSSNChi::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::beta2"),  CCTK_VarIndex("Kranc2BSSNChi::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::beta3"),  CCTK_VarIndex("Kranc2BSSNChi::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::betat1"),  CCTK_VarIndex("Kranc2BSSNChi::betat1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::betat2"),  CCTK_VarIndex("Kranc2BSSNChi::betat2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::betat3"),  CCTK_VarIndex("Kranc2BSSNChi::betat3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::chi"),  CCTK_VarIndex("Kranc2BSSNChi::chirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::Gam1"),  CCTK_VarIndex("Kranc2BSSNChi::Gam1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::Gam2"),  CCTK_VarIndex("Kranc2BSSNChi::Gam2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::Gam3"),  CCTK_VarIndex("Kranc2BSSNChi::Gam3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::h11"),  CCTK_VarIndex("Kranc2BSSNChi::h11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::h21"),  CCTK_VarIndex("Kranc2BSSNChi::h21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::h31"),  CCTK_VarIndex("Kranc2BSSNChi::h31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::h22"),  CCTK_VarIndex("Kranc2BSSNChi::h22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::h32"),  CCTK_VarIndex("Kranc2BSSNChi::h32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::h33"),  CCTK_VarIndex("Kranc2BSSNChi::h33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChi::K"),  CCTK_VarIndex("Kranc2BSSNChi::Krhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
