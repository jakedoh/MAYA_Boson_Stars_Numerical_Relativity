! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "gauge.h"


subroutine NullEvol_Diag(CCTK_ARGUMENTS)

  use NullDecomp_SpinDecomp, only: SpinDecompCoefs
  use NullEvol_DiagMod

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  logical, save :: first_time = .TRUE.
  logical :: truncate

  truncate = (IO_TruncateOutputFiles(cctkGH) .ne. 0) .and. first_time

  call NullEvol_DiagImArray(cctkGH, 'J_ev', truncate, 2_ik, jcn, jcs, diagtmp, null_xb, cctk_time)

  call NullEvol_DiagImArray(cctkGH, 'dxJ_ev', truncate, 2_ik, dxjcn, dxjcs, diagtmp, null_xb, cctk_time)

  call NullEvol_DiagImArray(cctkGH, 'Q_ev', truncate, 1_ik, qcn, qcs, diagtmp, null_xb, cctk_time)

  call NullEvol_DiagImArray(cctkGH, 'U_ev', truncate, 1_ik, ucn, ucs, diagtmp, null_xbh, cctk_time)

  if(first_order_scheme.ne.0) then

    call NullEvol_DiagImArray(cctkGH, 'CB_ev', truncate, 1_ik, cbcn, cbcs, diagtmp, null_xb, cctk_time)

    call NullEvol_DiagImArray(cctkGH, 'CK_ev', truncate, 1_ik, ckcn, ckcs, diagtmp, null_xb, cctk_time)

    call NullEvol_DiagImArray(cctkGH, 'NU_ev', truncate, 1_ik, nucn, nucs, diagtmp, null_xb, cctk_time)
  end if

  call NullEvol_DiagReArray(cctkGH, 'W_ev', truncate, wcn, wcs, diagtmp, null_xb, cctk_time)

  call NullEvol_DiagReArray(cctkGH, 'B_ev', truncate, bcn, bcs, diagtmp, null_xb, cctk_time)

  first_time = .FALSE.

end subroutine NullEvol_Diag
