/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void Kranc2BSSN_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::A11"),  CCTK_VarIndex("Kranc2BSSN::A11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::A21"),  CCTK_VarIndex("Kranc2BSSN::A21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::A31"),  CCTK_VarIndex("Kranc2BSSN::A31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::A22"),  CCTK_VarIndex("Kranc2BSSN::A22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::A32"),  CCTK_VarIndex("Kranc2BSSN::A32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::A33"),  CCTK_VarIndex("Kranc2BSSN::A33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::alpha"),  CCTK_VarIndex("Kranc2BSSN::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::beta1"),  CCTK_VarIndex("Kranc2BSSN::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::beta2"),  CCTK_VarIndex("Kranc2BSSN::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::beta3"),  CCTK_VarIndex("Kranc2BSSN::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::betat1"),  CCTK_VarIndex("Kranc2BSSN::betat1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::betat2"),  CCTK_VarIndex("Kranc2BSSN::betat2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::betat3"),  CCTK_VarIndex("Kranc2BSSN::betat3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::Gam1"),  CCTK_VarIndex("Kranc2BSSN::Gam1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::Gam2"),  CCTK_VarIndex("Kranc2BSSN::Gam2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::Gam3"),  CCTK_VarIndex("Kranc2BSSN::Gam3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::h11"),  CCTK_VarIndex("Kranc2BSSN::h11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::h21"),  CCTK_VarIndex("Kranc2BSSN::h21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::h31"),  CCTK_VarIndex("Kranc2BSSN::h31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::h22"),  CCTK_VarIndex("Kranc2BSSN::h22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::h32"),  CCTK_VarIndex("Kranc2BSSN::h32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::h33"),  CCTK_VarIndex("Kranc2BSSN::h33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::K"),  CCTK_VarIndex("Kranc2BSSN::Krhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSN::phi"),  CCTK_VarIndex("Kranc2BSSN::phirhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
