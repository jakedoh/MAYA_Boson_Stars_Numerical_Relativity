#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



int
LWT_startup (void)
{
  int ierr;
  
  ierr = CCTK_RegisterBanner ("Llama Wave Toy");
  assert (! ierr);
  
  return 0;
}



void
LWT_register_MoL (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  int evolved_group, rhs_group;
  int ierr;
  
  evolved_group = CCTK_GroupIndex ("LlamaWaveToy::scalar");
  assert (evolved_group >= 0);
  rhs_group = CCTK_GroupIndex ("LlamaWaveToy::scalardot");
  assert (rhs_group >= 0);
  ierr = MoLRegisterEvolvedGroup (evolved_group, rhs_group);
  assert (! ierr);
  
  evolved_group = CCTK_GroupIndex ("LlamaWaveToy::density");
  assert (evolved_group >= 0);
  rhs_group = CCTK_GroupIndex ("LlamaWaveToy::densitydot");
  assert (rhs_group >= 0);
  ierr = MoLRegisterEvolvedGroup (evolved_group, rhs_group);
  assert (! ierr);

  evolved_group = CCTK_GroupIndex ("LlamaWaveToy::velocity");
  assert (evolved_group >= 0);
  rhs_group = CCTK_GroupIndex ("LlamaWaveToy::velocitydot");
  assert (rhs_group >= 0);
  ierr = MoLRegisterEvolvedGroup (evolved_group, rhs_group);
  assert (! ierr);

/*  evolved_group = CCTK_VarIndex ("LlamaWaveToy::vy");
  assert (evolved_group >= 0);
  rhs_group = CCTK_VarIndex ("LlamaWaveToy::vydot");
  assert (rhs_group >= 0);
  ierr = MoLRegisterEvolved (evolved_group, rhs_group);
  assert (! ierr);
  
  evolved_group = CCTK_VarIndex ("LlamaWaveToy::vy");
  assert (evolved_group >= 0);
  rhs_group = CCTK_VarIndex ("LlamaWaveToy::vydot");
  assert (rhs_group >= 0);
  ierr = MoLRegisterEvolved (evolved_group, rhs_group);
  assert (! ierr);
  
  evolved_group = CCTK_VarIndex ("LlamaWaveToy::vz");
  assert (evolved_group >= 0);
  rhs_group = CCTK_VarIndex ("LlamaWaveToy::vzdot");
  assert (rhs_group >= 0);
  ierr = MoLRegisterEvolved (evolved_group, rhs_group);
  assert (! ierr); */
}
