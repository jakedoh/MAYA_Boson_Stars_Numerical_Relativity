/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void Kranc2BSSNMatter_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::A11"),  CCTK_VarIndex("Kranc2BSSNMatter::A11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::A21"),  CCTK_VarIndex("Kranc2BSSNMatter::A21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::A31"),  CCTK_VarIndex("Kranc2BSSNMatter::A31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::A22"),  CCTK_VarIndex("Kranc2BSSNMatter::A22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::A32"),  CCTK_VarIndex("Kranc2BSSNMatter::A32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::A33"),  CCTK_VarIndex("Kranc2BSSNMatter::A33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::alpha"),  CCTK_VarIndex("Kranc2BSSNMatter::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::beta1"),  CCTK_VarIndex("Kranc2BSSNMatter::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::beta2"),  CCTK_VarIndex("Kranc2BSSNMatter::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::beta3"),  CCTK_VarIndex("Kranc2BSSNMatter::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::betat1"),  CCTK_VarIndex("Kranc2BSSNMatter::betat1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::betat2"),  CCTK_VarIndex("Kranc2BSSNMatter::betat2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::betat3"),  CCTK_VarIndex("Kranc2BSSNMatter::betat3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::Gam1"),  CCTK_VarIndex("Kranc2BSSNMatter::Gam1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::Gam2"),  CCTK_VarIndex("Kranc2BSSNMatter::Gam2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::Gam3"),  CCTK_VarIndex("Kranc2BSSNMatter::Gam3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::h11"),  CCTK_VarIndex("Kranc2BSSNMatter::h11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::h21"),  CCTK_VarIndex("Kranc2BSSNMatter::h21rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::h31"),  CCTK_VarIndex("Kranc2BSSNMatter::h31rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::h22"),  CCTK_VarIndex("Kranc2BSSNMatter::h22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::h32"),  CCTK_VarIndex("Kranc2BSSNMatter::h32rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::h33"),  CCTK_VarIndex("Kranc2BSSNMatter::h33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::K"),  CCTK_VarIndex("Kranc2BSSNMatter::Krhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("Kranc2BSSNMatter::phi"),  CCTK_VarIndex("Kranc2BSSNMatter::phirhs"));
  /* Register all the evolved Array functions with MoL */

  if (CCTK_EQUALS(couple_matter, "yes") || CCTK_EQUALS(recalculate_constraints, "matter")) {
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
