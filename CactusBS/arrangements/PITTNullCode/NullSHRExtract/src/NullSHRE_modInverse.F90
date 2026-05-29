! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modInverse

  use cctk
  implicit none

contains

  subroutine wt_gup3(g, detg, gup)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: g
    type (gf2d),                  intent (inout) :: detg
    type (gf2d), dimension (4,4), intent (inout) :: gup

    detg%d = g(1,1)%d * g(2,2)%d * g(3,3)%d &
           - g(1,1)%d * g(2,3)%d * g(2,3)%d &
           - g(1,2)%d * g(1,2)%d * g(3,3)%d &
           + g(1,2)%d * g(1,3)%d * g(2,3)%d &
           + g(1,3)%d * g(1,2)%d * g(2,3)%d &
           - g(1,3)%d * g(1,3)%d * g(2,2)%d

    gup(1,1)%d = (  g(2,2)%d * g(3,3)%d - g(2,3)%d * g(2,3)%d ) / detg%d
    gup(1,2)%d = (- g(1,2)%d * g(3,3)%d + g(1,3)%d * g(2,3)%d ) / detg%d
    gup(1,3)%d = (  g(1,2)%d * g(2,3)%d - g(1,3)%d * g(2,2)%d ) / detg%d
    gup(2,2)%d = (  g(1,1)%d * g(3,3)%d - g(1,3)%d * g(1,3)%d ) / detg%d
    gup(2,3)%d = (- g(1,1)%d * g(2,3)%d + g(1,3)%d * g(1,2)%d ) / detg%d
    gup(3,3)%d = (  g(1,1)%d * g(2,2)%d - g(1,2)%d * g(1,2)%d ) / detg%d

  end subroutine wt_gup3


  subroutine wt_gup4(alpha, na, gup)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                  intent (in)    :: alpha
    type (gf2d), dimension (4),   intent (in)    :: na
    type (gf2d), dimension (4,4), intent (inout) :: gup

    integer :: i, j

    do i = 1, 3
       do j = i, 3
          gup(i,j)%d = gup(i,j)%d - na(i)%d * na(j)%d
       end do
    end do

    do i = 1, 4
       gup(i,4)%d = - na(i)%d / alpha%d
    end do

  end subroutine wt_gup4

end module NullSHRE_modInverse
