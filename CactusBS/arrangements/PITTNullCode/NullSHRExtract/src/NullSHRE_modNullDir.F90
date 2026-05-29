! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modNullDir 

  use cctk
  implicit none

contains

  subroutine wt_na(alpha, beta, na) ! eq. (22)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                intent (in)    :: alpha
    type (gf2d), dimension (3), intent (inout) :: beta
    type (gf2d), dimension (4), intent (inout) :: na

    CCTK_INT :: i

    ! unit (contravariant) normal to the slices at t=const.

    do i = 1, 3
       na(i)%d = - beta(i)%d / alpha%d
    end do
    na(4)%d = 1.d0 / alpha%d

  end subroutine wt_na


  subroutine wt_sa(gup, sigma, sigma2, sa) ! eq. (25)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: gup
    type (gf2d), dimension (3),   intent (in)    :: sigma
    type (gf2d),                  intent (inout) :: sigma2
    type (gf2d), dimension (4),   intent (inout) :: sa

    CCTK_INT :: i, j

    ! raise sigma_a, store in s^a

    do i = 1, 3
       sa(i)%d = 0.d0
       do j = 1, 3
          sa(i)%d = sa(i)%d + gup(i,j)%d * sigma(j)%d
       end do
    end do

    ! norm of sigma

    sigma2%d = 0.d0
    do i = 1, 3
       sigma2%d = sigma2%d + sa(i)%d * sigma(i)%d
    end do

    ! normalize sa^a

    do i = 1, 3
       sa(i)%d = sa(i)%d / sqrt(sigma2%d)
    end do
    sa(4)%d = 0.d0

  end subroutine wt_sa


  subroutine wt_ell(g, alpha, beta, sa, na, elld, ell, halt_on_negative_elld, elld_min) ! eq. (26)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: g
    type (gf2d),                  intent (in)    :: alpha
    type (gf2d), dimension (3),   intent (in)    :: beta
    type (gf2d), dimension (4),   intent (in)    :: sa, na
    type (gf2d),                  intent (inout) :: elld
    type (gf2d), dimension (4),   intent (inout) :: ell
    CCTK_INT,                     intent (in)    :: halt_on_negative_elld
    CCTK_REAL,                    intent (in)    :: elld_min

    character(200) message

    CCTK_INT :: i, j

    elld%d = alpha%d
    do i = 1, 3
       do j = 1, 3
          elld%d = elld%d - g(i,j)%d * beta(i)%d * sa(j)%d
       end do
    end do

    if (halt_on_negative_elld.ne.0) then
       if (any(elld%d.lt.0)) call CCTK_WARN(0, "negative value for null vector denominator")
    else
       if (any(elld%d.lt.elld_min)) then
          write (message,'("null vector denominator is less than minimum accepted value: ",g15.5," < ",g15.5)') minval(elld%d), elld_min
          call CCTK_INFO(message)
          call CCTK_INFO("adding artificial correction term")
          elld%d = sqrt( elld%d**2 + elld_min**2 )
       end if
    end if

    do i = 1, 4
       ell(i)%d = (na(i)%d + sa(i)%d) / elld%d
    end do

  end subroutine wt_ell


  subroutine wt_g1(dg, ell, j0inv, g1) ! eq. (44)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4,4), intent (in)  :: dg
    type (gf2d), dimension (4),     intent (in)  :: ell
    type (gf2d), dimension (4,4),   intent (in)  :: j0inv
    type (gf2d), dimension (4,4),   intent (out) :: g1

    CCTK_INT :: i, j, k, l

    ! get the o(lambda) part of the cauchy metric

    do i = 1, 4
       do j = 1, 4
          g1(i,j)%d = 0.d0 
          do k = 1, 4
             do l = 1, 4
                g1(i,j)%d = g1(i,j)%d + j0inv(l,k)%d*ell(k)%d *dg(i,j,l)%d
             end do
          end do
       end do
    end do

  end subroutine wt_g1


  subroutine wt_dna (alpha, beta, dalpha, dbeta, dna) ! eq. (35,36,37)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                    intent (in)    :: alpha
    type (gf2d), dimension (3),     intent (in)    :: beta
    type (gf2d), dimension (4),     intent (in)    :: dalpha
    type (gf2d), dimension (3,4),   intent (in)    :: dbeta
    type (gf2d), dimension (4,2:4), intent (inout) :: dna

    CCTK_INT :: i, a

    ! angular derivative of the spatial part $n^{i}$

    do i = 1, 3
       do a = 2, 3
          dna(i,a)%d = 0.d0
          dna(i,a)%d = dna(i,a)%d + (- dbeta(i,a)%d / alpha%d +&
                     & beta(i)%d * dalpha(a)%d / alpha%d ** 2)
       end do
    end do

    ! angular derivative of the $n^{4}$ component

    do a = 2, 3
       dna(4,a)%d = 0.d0
       dna(4,a)%d = dna(4,a)%d - dalpha(a)%d / alpha%d ** 2
    end do

    ! time derivative of $n^{a}$

    do i = 1, 3
       dna(i,4)%d = - dbeta(i,4)%d / alpha%d + beta(i)%d * dalpha(4)&
                  & %d / alpha%d ** 2
    end do
    dna(4,4)%d = - dalpha(4)%d / alpha%d ** 2

  end subroutine wt_dna


  subroutine wt_dsa (sigma2, dsigma, sa, gup, dg, dsa)  ! eq (41..43)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                    intent (in)    :: sigma2
    type (gf2d), dimension (3,2:3), intent (in)    :: dsigma
    type (gf2d), dimension (4),     intent (in)    :: sa
    type (gf2d), dimension (4,4),   intent (in)    :: gup
    type (gf2d), dimension (4,4,4), intent (in)    :: dg
    type (gf2d), dimension (4,2:4), intent (inout) :: dsa

    CCTK_INT :: i, a, k, l, mu

    ! angular derivatives of $s^{i}$

    do i = 1, 3
       do a = 2, 3
          dsa(i,a)%d = 0.d0
          do k = 1, 3
             do l = 1, 3
                dsa(i,a)%d = dsa(i,a)%d + (- gup(i,k)%d + 0.5d0 *&
                        & sa(i)%d * sa(k)%d) * sa(l)%d * dg(k,l,a)%d
             end do
          end do
          do k = 1, 3
             dsa(i,a)%d = dsa(i,a)%d + (gup(i,k)%d - sa(i)%d * sa(k)&
                  & %d) * dsigma(k,a)%d / sqrt(sigma2%d)
          end do
       end do
    end do

    ! time derivative of $s^{i}$ eq. (38) 

    do i = 1, 3
       dsa(i,4)%d = 0.d0
       do k = 1, 3
          do l = 1, 3
             dsa(i,4)%d = dsa(i,4)%d + dg(k,l,4)%d * sa(l)%d * (&
                  & -gup(i,k)%d + 0.5d0 * sa(i)%d * sa(k)%d)
          end do
       end do
    end do

    ! the zero components

    do mu = 2, 4
       dsa(4,mu)%d = 0.d0
    end do

  end subroutine wt_dsa

end module NullSHRE_modNullDir
