#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine Proca_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaEvolve::rhs_Ex" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaEvolve::rhs_Ey" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaEvolve::rhs_Ez" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaEvolve::rhs_Ax" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaEvolve::rhs_Ay" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaEvolve::rhs_Az" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaEvolve::rhs_Zeta" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaEvolve::rhs_Aphi" )

end subroutine Proca_InitSymBound
