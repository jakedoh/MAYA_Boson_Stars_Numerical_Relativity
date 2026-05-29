! Module with various common variable.
! $Header$
#include "cctk.h"

module EHFinder_mod

  use EHFinder_Constants

  CCTK_INT, dimension(0:5), parameter :: ll = (/ 1, 2, 4, 8, 16, 32 /)
  CCTK_INT :: ixl, ixr, jyl, jyr, kzl, kzr
  CCTK_REAL :: hfac
  CCTK_INT :: nx, ny, nz, ngx, ngy, ngz
  CCTK_INT :: max_handle, min_handle, sum_handle
  CCTK_INT :: last_called_at_iteration = -1
  CCTK_INT :: ierr
  CCTK_REAL :: delta
  CCTK_REAL :: h
  CCTK_REAL :: ex_value
  CCTK_INT, dimension(3) :: imin_glob, imax_glob
  CCTK_REAL, dimension(3) :: fimin_glob, fimax_glob
  CCTK_INT :: githeta, gjphi
  logical :: mask_first = .true.
  logical :: symx, symy, symz
  logical :: s_symx, s_symy, s_symz
  logical, dimension(:), allocatable :: re_init_this_level_set, &
                                        re_initialize_undone
end module EHFinder_mod
