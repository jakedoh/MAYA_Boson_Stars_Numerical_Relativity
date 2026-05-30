/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void Kranc2BSSNChiMatter_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_Kranc2BSSNChiMatter_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(Kranc2BSSNChiMatter_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::A11"),  CCTK_VarIndex("Kranc2BSSNChiMatter::A11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::A21"),  CCTK_VarIndex("Kranc2BSSNChiMatter::A21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::A31"),  CCTK_VarIndex("Kranc2BSSNChiMatter::A31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::A22"),  CCTK_VarIndex("Kranc2BSSNChiMatter::A22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::A32"),  CCTK_VarIndex("Kranc2BSSNChiMatter::A32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::A33"),  CCTK_VarIndex("Kranc2BSSNChiMatter::A33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::alpha"),  CCTK_VarIndex("Kranc2BSSNChiMatter::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::beta1"),  CCTK_VarIndex("Kranc2BSSNChiMatter::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::beta2"),  CCTK_VarIndex("Kranc2BSSNChiMatter::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::beta3"),  CCTK_VarIndex("Kranc2BSSNChiMatter::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::betat1"),  CCTK_VarIndex("Kranc2BSSNChiMatter::betat1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::betat2"),  CCTK_VarIndex("Kranc2BSSNChiMatter::betat2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::betat3"),  CCTK_VarIndex("Kranc2BSSNChiMatter::betat3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::chi"),  CCTK_VarIndex("Kranc2BSSNChiMatter::chirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::Gam1"),  CCTK_VarIndex("Kranc2BSSNChiMatter::Gam1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::Gam2"),  CCTK_VarIndex("Kranc2BSSNChiMatter::Gam2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::Gam3"),  CCTK_VarIndex("Kranc2BSSNChiMatter::Gam3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::h11"),  CCTK_VarIndex("Kranc2BSSNChiMatter::h11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::h21"),  CCTK_VarIndex("Kranc2BSSNChiMatter::h21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::h31"),  CCTK_VarIndex("Kranc2BSSNChiMatter::h31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::h22"),  CCTK_VarIndex("Kranc2BSSNChiMatter::h22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::h32"),  CCTK_VarIndex("Kranc2BSSNChiMatter::h32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::h33"),  CCTK_VarIndex("Kranc2BSSNChiMatter::h33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNChiMatter::K"),  CCTK_VarIndex("Kranc2BSSNChiMatter::Krhs"));
  /* Register all the evolved Array functions with MoL */

  if (mol_register_tmunu && (CCTK_EQUALS(couple_matter, "yes") || CCTK_EQUALS(recalculate_constraints, "matter"))) {
    ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("TmunuBase::stress_energy_scalar"));
    ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("TmunuBase::stress_energy_vector"));
    ierr += MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex("TmunuBase::stress_energy_tensor"));
  }

  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ADMBase::metric"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ADMBase::curv"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ADMBase::lapse"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("ADMBase::shift"));
  return;
}
