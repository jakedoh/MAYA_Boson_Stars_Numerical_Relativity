#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"

void TmunuBase_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  return;
}

void TmunuBase_ApplyBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_INT ierr = 0;


    /* these are no longer needed I believe and cause L3 warnings in Cactus */
    ierr = CartSymGN(cctkGH, "TmunuBase::stress_energy_scalar");
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to apply symmetry BC for TmunuBase::stress_energy_scalar!");
  
    ierr = CartSymGN(cctkGH, "TmunuBase::stress_energy_vector");
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to apply symmetry BC for TmunuBase::stress_energy_vector!");
  
    ierr = CartSymGN(cctkGH, "TmunuBase::stress_energy_tensor");
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to apply symmetry BC for TmunuBase::stress_energy_tensor!");

 
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "TmunuBase::stress_energy_scalar", "none");
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register A_group_bound BC for TmunuBase::stress_energy_scalar!");

    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "TmunuBase::stress_energy_vector", "none");
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register A_group_bound BC for TmunuBase::stress_energy_vector!");

    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "TmunuBase::stress_energy_tensor", "none");
    if (ierr < 0)
       CCTK_WARN(-1, "Failed to register A_group_bound BC for TmunuBase::stress_energy_tensor!");

  return;
}

