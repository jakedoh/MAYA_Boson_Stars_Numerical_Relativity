! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine proca_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::Ex" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::Ey" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::Ez" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::Ax" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::Ay" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::Az" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Zeta" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Aphi" )

end subroutine proca_symmetries
