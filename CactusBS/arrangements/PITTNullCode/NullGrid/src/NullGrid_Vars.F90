! vim: syntax=fortran
! $Header$
#include "cctk.h"

module NullGrid_Vars

  implicit none

  CCTK_REAL,    dimension(:,:), pointer, save :: qs, ps, pp
  CCTK_COMPLEX, dimension(:,:), pointer, save :: zz
  CCTK_REAL,    dimension(:),   pointer, save :: xb, xbh, rb,rbh
  CCTK_INT,     dimension(:,:), pointer, save :: guard, EG

  CCTK_INT,               save :: nx
  CCTK_INT, dimension(2), save :: gsh, lsh, lbnd, ubnd
  CCTK_REAL, save :: delta(2)
  CCTK_REAL, save :: rwt, dx, xbin

   !----------------------------------------------------------------------
   ! patch indices
   !----------------------------------------------------------------------
 
   CCTK_INT, parameter :: ip_n = 1, ip_s = 2

end module NullGrid_Vars

