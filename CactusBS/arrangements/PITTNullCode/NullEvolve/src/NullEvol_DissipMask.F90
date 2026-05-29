! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module NullEvol_DissipMask

contains

  subroutine NullEvol_set_dissipmask(dissip_mask)
    use cctk
    use NullGrid_Vars
    implicit none

    CCTK_REAL, dimension (lsh(1),lsh(2)), intent (out) :: dissip_mask
    CCTK_REAL :: rD0, rD1
    CCTK_REAL, parameter :: zero = 0.0
    CCTK_REAL, parameter :: one = 1.0

    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS

    if (CCTK_EQUALS(dissip_mask_type, "one")) then
      dissip_mask = one
      return
    end if

    if (CCTK_EQUALS(dissip_mask_type, "zero at eq, one at pole")) then
       rD0 = one
       rD1 = zero
    else if (CCTK_EQUALS(dissip_mask_type, "zero at rD0, one at pole")) then
       rD0 = one + N_dissip_zero_outside_eq * maxval(delta)
       rD1 = zero
    else if (CCTK_EQUALS(dissip_mask_type, "zero at rD0, one at eq")) then
       rD0 = one + N_dissip_zero_outside_eq * maxval(delta)
       rD1 = one
    else if (CCTK_EQUALS(dissip_mask_type, "zero at rD0, one at rD1")) then
       rD0 = one + N_dissip_zero_outside_eq * maxval(delta)
       rD1 = one + N_dissip_one_outside_eq * maxval(delta)
    else
      call CCTK_WARN(0, "unsupported dissipation mask type")
    end if

    dissip_mask = (rD0-dsqrt(ps**2+qs**2)) / (rD0-rD1)
    dissip_mask = max(zero, min(one, dissip_mask))

  end subroutine NullEvol_set_dissipmask

end module NullEvol_DissipMask
