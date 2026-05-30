! vim: syntax=fortran
#include "cctk.h"

module NullSHRE_modEta

  use cctk
  implicit none

contains

  subroutine wt_j1(g, gup, dg, dalpha, beta, dbeta, na, sa, ell, dsa, &
                  dna, j0inv, elld, delld, j1)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4),   intent (in)    :: g
    type (gf2d), dimension (4,4),   intent (in)    :: gup
    type (gf2d), dimension (4,4,4), intent (in)    :: dg
    type (gf2d), dimension (4),     intent (in)    :: dalpha
    type (gf2d), dimension (3),     intent (in)    :: beta
    type (gf2d), dimension (3,4),   intent (in)    :: dbeta
    type (gf2d), dimension (4),     intent (in)    :: na
    type (gf2d), dimension (4),     intent (in)    :: sa
    type (gf2d), dimension (4),     intent (in)    :: ell
    type (gf2d), dimension (4,2:4), intent (in)    :: dsa
    type (gf2d), dimension (4,2:4), intent (in)    :: dna
    type (gf2d), dimension (4,4),   intent (in)    :: j0inv
    type (gf2d),                    intent (in)    :: elld
    type (gf2d), dimension (4),     intent (inout) :: delld
    type (gf2d), dimension (4,4),   intent (inout) :: j1

   CCTK_INT :: mu, nu, rho, tau, a, l, i, j

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! lambda derivative of $\ell^{mu}$
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do mu = 1, 4
       j1(mu,1)%d = 0.d0
       do nu = 1, 4
          do rho = 1, 4
             do tau = 1, 4
                do l = 1, 4
                   j1(mu,1)%d = j1(mu,1)%d &
                   + gup(mu,tau)%d * ell(nu)%d * ell(rho)%d &
                   * ( 0.5d0 * dg(nu,rho,l)%d * j0inv(l,tau)%d &
                   - dg(rho,tau,l)%d * j0inv(l,nu)%d)
                end do
             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! angular derivatives of $L = \alpha - g_{ij} \beta^{i} s^{j}$
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! $L_{,A} = \alpha_{,a} - g_{ij,a} \beta^{i} s^{j} 
    ! - g_{ij} \beta^{i}_{,a} s^{j} - g_{ij} \beta^{i} s^{j}_{,A}  $

    do a = 2, 3
       delld(a)%d = 0.d0
       do i = 1, 3
          do j = 1, 3
             delld(a)%d = delld(a)%d&
                        - g(i,j)%d * dbeta(i,a)%d * sa(j)%d&
                        - g(i,j)%d * beta(i)%d * dsa(j,a)%d&
                        - beta(i)%d * sa(j)%d *dg(i,j,a)%d
          end do
       end do
    end do

    ! $\ldots + \alpha_{,a} \ldots $

    do a = 2, 3
          delld(a)%d = delld(a)%d + dalpha(a)%d
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time derivative of $L = \alpha - g_{ij} \beta^{i} s^{j}$
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! L_{u} = \alpha_{,t} \ldots $

    delld(1)%d = 0.d0
    delld(4)%d = dalpha(4)%d

    ! $\ldots - g_{ij,t} \beta^{i} s^{j} - g_{ij} \beta^{i}_{,t} s^{j}
    ! - g_{ij} \beta^{i} s^{j}_{,t}

    do i = 1, 3
       do j = 1, 3
          delld(4)%d = delld(4)%d &
                     - g(i,j)%d * dbeta(i,4)%d * sa(j)%d&
                     - g(i,j)%d * beta(i)%d * dsa(j,4)%d&
                     - beta(i)%d * sa(j)%d *dg(i,j,4)%d 
       end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! angular and time derivatives of $\ell^{\mu}$
    !  ell(i)%d = (na(i)%d + sa(i)%d) / elld%d
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do mu = 1, 4
       do nu = 2, 4
          j1(mu,nu)%d = (dna(mu,nu)%d + dsa(mu,nu)%d) / elld%d &
          - (na(mu)%d + sa(mu)%d) * delld(nu)%d / elld%d ** 2
       end do
    end do

  end subroutine wt_j1

 
  subroutine wt_eta0(g, j0, eta0)
!null metric zero order coefficients eq. (46)
    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: g, j0
    type (gf2d), dimension (4,4), intent (inout) :: eta0

    CCTK_INT :: i, j, k, l
   
    do i = 1, 4
       do j = i, 4
          eta0(i,j)%d = 0.d0
          do k = 1, 4
             do l = 1, 4
                eta0(i,j)%d = eta0(i,j)%d + j0(k,i)%d * j0(l,j)%d * g(k,l)%d
             end do
          end do
       end do
    end do

    eta0(1,1)%d = 0.d0
    eta0(1,2)%d = 0.d0
    eta0(1,3)%d = 0.d0
    eta0(1,4)%d = -1.d0

  end subroutine wt_eta0
 

  subroutine wt_eta1(g, g1, j0, j1, eta1)
! lambda derivative of the null metric coefficients -- eq. (47)
    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: g, g1
    type (gf2d), dimension (4,4), intent (in)    :: j0, j1
    type (gf2d), dimension (4,4), intent (inout) :: eta1

    CCTK_INT :: i, j, k, l

    do i = 1, 4
       do j = i, 4
          eta1(i,j)%d = 0.d0
          do k = 1, 4
             do l = 1, 4
                eta1(i,j)%d = eta1(i,j)%d &
                + j0(k,i)%d * j0(l,j)%d * g1(k,l)%d &
                + (j0(k,i)%d * j1(l,j)%d + j1(k,i)%d * j0(l,j)%d) * g(k,l)%d
             end do
          end do
       end do
    end do

    eta1(1,1)%d = 0.d0
    eta1(1,2)%d = 0.d0
    eta1(1,3)%d = 0.d0
    eta1(1,4)%d = 0.d0

  end subroutine wt_eta1
 

  subroutine wt_etaup0(eta0, temp, etaup0)
! contravariant null metric eq. (50)-(51)

    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: eta0
    type (gf2d),                  intent (inout) :: temp
    type (gf2d), dimension (4,4), intent (inout) :: etaup0

    ! (1) invert eta0_{ab}
    temp%d = eta0(2,2)%d * eta0(3,3)%d - eta0(2,3)%d * eta0(2,3)%d

    etaup0(2,2)%d =   eta0(3,3)%d / temp%d
    etaup0(2,3)%d = - eta0(2,3)%d / temp%d
    etaup0(3,3)%d =   eta0(2,2)%d / temp%d

    ! (3) eta^{1a}

    etaup0(1,2)%d = etaup0(2,2)%d * eta0(2,4)%d &
                  + etaup0(2,3)%d * eta0(3,4)%d
    etaup0(1,3)%d = etaup0(3,2)%d * eta0(2,4)%d &
                  + etaup0(3,3)%d * eta0(3,4)%d

    ! (4) eta^{11}

    etaup0(1,1)%d = - eta0(4,4)%d &
                  + etaup0(1,2)%d * eta0(2,4)%d &
                  + etaup0(1,3)%d * eta0(3,4)%d

    etaup0(4,1)%d = -1.d0
    etaup0(4,2)%d = 0.d0
    etaup0(4,3)%d = 0.d0
    etaup0(4,4)%d = 0.d0

  end subroutine wt_etaup0
 

  subroutine wt_etaup1(eta1, etaup0, etaup1)
! lambda derivative of contravariant null metric eq. (52)
    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4), intent (in)    :: eta1
    type (gf2d), dimension (4,4), intent (in)    :: etaup0
    type (gf2d), dimension (4,4), intent (inout) :: etaup1

    CCTK_INT :: i, j, k, l

    do i = 1, 4
       do j = i, 4
          etaup1(i,j)%d = 0.d0
          do k = 1, 4
             do l = 1, 4
                etaup1(i,j)%d = etaup1(i,j)%d &
                - etaup0(k,i)%d * etaup0(l,j)%d * eta1(k,l)%d
             end do
          end do
       end do
    end do

    etaup1(4,1)%d = 0.d0
    etaup1(4,2)%d = 0.d0
    etaup1(4,3)%d = 0.d0
    etaup1(4,4)%d = 0.d0

  end subroutine wt_etaup1


end module NullSHRE_modEta
