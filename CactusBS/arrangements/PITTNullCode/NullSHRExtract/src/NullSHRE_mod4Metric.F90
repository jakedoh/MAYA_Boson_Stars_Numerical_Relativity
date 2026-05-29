! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_mod4Metric

  use cctk
  implicit none

contains

! computes the git components of the metric, eq. (21) 

  subroutine wt_g(alpha, beta, g)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                  intent (in)    :: alpha
    type (gf2d), dimension (3),   intent (in)    :: beta
    type (gf2d), dimension (4,4), intent (inout) :: g

    g(1,4)%d = g(1,1)%d * beta(1)%d + g(1,2)%d * beta(2)%d + g(1,3)%d&
         & * beta(3)%d
    g(2,4)%d = g(2,1)%d * beta(1)%d + g(2,2)%d * beta(2)%d + g(2,3)%d&
         & * beta(3)%d
    g(3,4)%d = g(3,1)%d * beta(1)%d + g(3,2)%d * beta(2)%d + g(3,3)%d&
         & * beta(3)%d

    g(4,4)%d = - alpha%d ** 2 + g(1,4)%d * beta(1)%d + g(2,4)%d *&
         & beta(2)%d + g(3,4)%d * beta(3)%d

  end subroutine wt_g

  subroutine wt_dg (alpha, beta, dalpha, dbeta, g, dg)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                    intent (in)    :: alpha
    type (gf2d), dimension (3),     intent (in)    :: beta
    type (gf2d), dimension (4),     intent (in)    :: dalpha
    type (gf2d), dimension (3,4),   intent (in)    :: dbeta
    type (gf2d), dimension (4,4),   intent (in)    :: g
    type (gf2d), dimension (4,4,4), intent (inout) :: dg

    integer i, j, l

    ! derivatives of g_i4 - note that g_ij,t is already known

    do i = 1, 3
       do l = 1, 4
          dg(i,4,l)%d = 0.d0
          do j = 1, 3
             dg(i,4,l)%d = dg(i,4,l)%d + dg(i,j,l)%d * beta(j)%d +&
                  & g(i ,j)%d * dbeta(j,l)%d
          end do
       end do
    end do

    ! derivatives of g_44

    do l = 1, 4
       dg(4,4,l)%d = -2.d0 * alpha%d * dalpha(l)%d
       do i = 1, 3
          do j = 1, 3
             dg(4,4,l)%d = dg(4,4,l)%d + dg(i,j,l)%d * beta(i)%d *&
                  & beta(j)%d + 2.d0 * g(i,j)%d * beta(i)%d * dbeta(j,l)%d
          end do
       end do
    end do

  end subroutine wt_dg


end module NullSHRE_mod4Metric
