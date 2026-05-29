#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void TmunuBase_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT ierr = 0;
  
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_scalar"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_vector"));
  ierr += MoLRegisterConstrainedGroup(CCTK_GroupIndex("TmunuBase::stress_energy_tensor"));

  return;
}
