
#include "cctk.h"
#include "cctk_Arguments.h"

subroutine LeanBSSN_zero_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS

  rhs_conf_fac = 0

  rhs_hxx = 0
  rhs_hxy = 0
  rhs_hxz = 0
  rhs_hyy = 0
  rhs_hyz = 0
  rhs_hzz = 0

  rhs_tracek = 0

  rhs_axx = 0
  rhs_axy = 0
  rhs_axz = 0
  rhs_ayy = 0
  rhs_ayz = 0
  rhs_azz = 0

  rhs_gammatx = 0
  rhs_gammaty = 0
  rhs_gammatz = 0

  rhs_alp = 0

  rhs_betax = 0
  rhs_betay = 0
  rhs_betaz = 0

end subroutine LeanBSSN_zero_rhs
