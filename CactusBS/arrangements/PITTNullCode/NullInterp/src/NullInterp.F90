! vim: syntax=fortran
! $Header$
#include "cctk.h"

module NullInterp
  use cctk
  use NullInterp_Eth
  use NullInterp_Deriv
  use NullInterp_Interp
  implicit none
end module NullInterp
