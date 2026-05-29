#include "cctk.h"

subroutine IP2_calc_inv3 (g3, gu3)
  use matinv
  implicit none
  CCTK_REAL, intent(in)  :: g3(3,3)
  CCTK_REAL, intent(out) :: gu3(3,3)
  call calc_inv3 (g3, gu3)
end subroutine IP2_calc_inv3
