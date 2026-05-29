! NPScalars
! NPboundaries.F90: define symmetry boundaries for NP scalars
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine NPboundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1
  CCTK_INT width

  ! let's just for simplicity say that
  width = 1

  ! since these are just used for analysis, register all BCs as "flat"

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, width, -one,      &
       "NPScalars::NPPsi4R_group", "flat")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for NPScalars::NPPsi4R_group!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, width, -one,      &
       "NPScalars::NPPsi4I_group", "flat")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for NPScalars::NPPsi4I_group!")

end subroutine NPboundaries
