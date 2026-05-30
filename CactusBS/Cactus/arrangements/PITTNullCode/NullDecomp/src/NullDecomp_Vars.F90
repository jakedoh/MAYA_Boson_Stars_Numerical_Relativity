#include "cctk.h"
module NullDecomp_Vars
  CCTK_INT :: MyProc = -1
  CCTK_INT :: NumProcs = -1
  CCTK_INT :: lmax = 9
  CCTK_INT :: lq, uq, lp, up
  CCTK_INT :: sum_handle
  CCTK_INT :: min_handle
  CCTK_INT :: max_handle

  CCTK_REAL, DIMENSION(:,:,:), pointer :: kern_p
  CCTK_REAL, DIMENSION(:,:), pointer :: area_p
end module NullDecomp_Vars
