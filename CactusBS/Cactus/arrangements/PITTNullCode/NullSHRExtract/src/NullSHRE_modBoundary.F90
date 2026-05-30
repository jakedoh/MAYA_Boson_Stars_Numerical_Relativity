! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"

module NullSHRE_modBoundary

  use cctk
  implicit none

contains

  subroutine wt_null_boundary_mask (nn1, nn2, dx,&
             rwt, xin, r0, x_wt, boundary_masks, buffer_width)

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                       intent (in)    :: nn1, nn2 
    CCTK_REAL,                      intent (in)    :: dx, rwt, xin, buffer_width
    type (gf2d),                    intent (in)    :: r0
    type (gf2d),                    intent (inout) :: x_wt
    CCTK_INT, dimension (nn1,nn2),  intent (inout) :: boundary_masks

    DECLARE_CCTK_FUNCTIONS

    ! set the mask employed by the hypersurface and evolution routines
    ! note that we do keep a minimum distance of dx/4 between the WT
    ! and the first point to be updated by the evolution code

    x_wt%d = r0%d / (rwt + r0%d)
    boundary_masks(1:nn1,1:nn2) = int (buffer_width + (x_wt%d-xin)/dx) + 1
    ! int((x_wt%d - 0.5d0) / dx) + 1

    if ((minval(boundary_masks)).lt.3) then
       call CCTK_INFO("WARNING: Extraction world tube collapsing too close to innermost Bondi gridpoint")
    end if

  end subroutine wt_null_boundary_mask


  subroutine wt_jb_simple (nn1, nn2, pp, eta0, r0, qa, j_wt)
! the metric on the sphere -- eqs. (62), (63), (64)
    use NullSHRE_modGFdef
    implicit none
 
    CCTK_INT,                       intent (in) :: nn1, nn2
    CCTK_REAL, dimension (nn1, nn2), intent (in) :: pp
    type (gf2d),  dimension (4,4), intent (in)    :: eta0
    type (gf2d),                   intent (in)    :: r0

    type (gf2dc), dimension (2:3), intent (inout) :: qa
    type (gf2dc),                  intent (inout) :: j_wt

    CCTK_INT :: a, b
    CCTK_COMPLEX ii
    ii = (0.0d0, 1.0d0)

    ! q^{a}

    qa(2)%d = 0.5d0 * pp
    qa(3)%d = 0.5d0 * pp * ii

    ! this is the metric component j

    j_wt%d = (0., 0.d0)
    do a = 2, 3
       do b = 2, 3
          j_wt%d = j_wt%d + 0.5d0 * qa(a)%d * qa(b)%d * eta0(a,b)%d / r0%d ** 2
       end do
    end do

  end subroutine wt_jb_simple 



  subroutine wt_jb (nn1, nn2, pp, eta0, eta1, r0, dr0, qa, j_wt, j_l)
! the metric on the sphere -- eqs. (62), (63), (64)
    use NullSHRE_modGFdef
    implicit none
 
    CCTK_INT,                       intent (in) :: nn1, nn2
    CCTK_REAL, dimension (nn1, nn2), intent (in) :: pp
    type (gf2d),  dimension (4,4), intent (in)    :: eta0,&
         & eta1
    type (gf2d),                   intent (in)    :: r0

    type (gf2d),  dimension (4),   intent (in)    :: dr0
    type (gf2dc), dimension (2:3), intent (inout) :: qa
    type (gf2dc),                  intent (inout) :: j_wt,&
         & j_l

    CCTK_INT :: a, b
    CCTK_COMPLEX ii
    ii = (0.0d0, 1.0d0)

    ! q^{a}

    qa(2)%d = 0.5d0 * pp
    qa(3)%d = 0.5d0 * pp * ii

    ! this is the metric component j

    j_wt%d = (0., 0.d0)
    do a = 2, 3
       do b = 2, 3
          j_wt%d = j_wt%d + 0.5d0 * qa(a)%d * qa(b)%d * eta0(a,b)%d / r0%d ** 2
       end do
    end do

    ! j_{,\lambda}

    j_l%d = - 2.d0 * dr0(1)%d / r0%d * j_wt%d
    do a = 2, 3
       do b = 2, 3
          j_l%d = j_l%d + 0.5d0 * qa(a)%d * qa(b)%d * eta1(a,b)%d / r0 %d ** 2
       end do
    end do

  end subroutine wt_jb


  subroutine wt_bb (r0, dr0, j_wt, j_l, beta_wt, beta_l)
! the expansion factor eqs. (67) - (68)
    use NullSHRE_modGFdef
    implicit none

    type (gf2d),                 intent (in)    :: r0
    type (gf2d),  dimension (4), intent (in)    :: dr0
    type (gf2dc),                intent (in)    :: j_wt, j_l
    type (gf2d),                 intent (inout) :: beta_wt, beta_l

    ! beta eq. (67)

    beta_wt%d = -0.5d0 * log(dr0(1)%d)

    ! beta_{,\lambda} eq. (68)

    beta_l%d = r0%d / (8.d0 * dr0(1)%d) * (j_l%d * conjg(j_l%d) -&
         & dble(conjg(j_wt%d) * j_l%d) ** 2 / (1.d0 + dble(j_wt%d) ** 2&
         & + dimag(j_wt%d) ** 2))

  end subroutine wt_bb


  subroutine wt_ub (qs, ps, nn1, nn2, etaup0, etaup1,&
             dr0, dr1, beta_l, qa, u_wt, u_l)
! the shift U
    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                     intent (in) :: nn1, nn2
    CCTK_REAL,   dimension (nn1, nn2), intent (in) :: qs, ps
    type (gf2d), dimension (4,4), intent (in)    :: etaup0, etaup1
    type (gf2d), dimension (4),   intent (in)    :: dr0, dr1
    type (gf2d),                   intent (in)    :: beta_l
    type (gf2dc), dimension (2:3), intent (inout) :: qa
    type (gf2dc),                  intent (inout) :: u_wt, u_l

    CCTK_INT a, b
    CCTK_COMPLEX ii
    ii = (0.0d0, 1.0d0)

    ! q_{a}

    qa(2)%d = 2.d0 / (1.d0 + qs ** 2 + ps ** 2)
    qa(3)%d = 2.d0 / (1.d0 + qs ** 2 + ps ** 2) * ii

    ! u eq. (74)

    u_wt%d = (0., 0.d0)

    do a = 2, 3
       u_wt%d = u_wt%d - etaup0(1,a)%d * qa(a)%d
    end do

    do a = 2, 3
       do b = 2, 3
          u_wt%d = u_wt%d - dr0(b)%d / dr0(1)%d * etaup0(a,b)%d * qa(a)%d
       end do
    end do

    ! u_l eq. (75)

    u_l%d = 2.d0 * beta_l%d * (u_wt%d + etaup0(1,2)%d * qa(2)%d +&
         & etaup0(1,3)%d * qa(3)%d)

    u_l%d = u_l%d  - etaup1(1,2)%d * qa(2)%d - etaup1(1,3)%d * qa(3)%d

    do a = 2, 3
       do b = 2, 3
          u_l%d = u_l%d - dr1(b)%d * etaup0(a,b)%d * qa(a)%d / dr0(1)%d
       end do
    end do

    do a = 2, 3
       do b = 2, 3
          u_l%d = u_l%d - dr0(b)%d * etaup1(a,b)%d * qa(a)%d / dr0(1)%d
       end do
    end do

  end subroutine wt_ub


  subroutine wt_qb(rwt, r0, j_wt, u_l, q_wt)
! the auxiliary variable of the U equation
    use NullSHRE_modGFdef
    implicit none

    CCTK_REAL,                      intent (in)    :: rwt
    type (gf2d),                    intent (in)    :: r0
    type (gf2dc),                   intent (in)    :: j_wt, u_l
    type (gf2dc),                   intent (inout) :: q_wt

    q_wt%d = r0%d**2 * ( j_wt%d * conjg(u_l%d) &
           + ( sqrt(1. + j_wt%d * conjg(j_wt%d)) ) * u_l%d )

  end subroutine wt_qb


  subroutine wt_wb (etaup0, etaup1, r0, dr0, dr1, beta_l, w_wt, w_l)
! the mas aspect W
    use NullSHRE_modGFdef
    implicit none

    type (gf2d), dimension (4,4),   intent (in)    :: etaup0, etaup1
    type (gf2d),                    intent (in)    :: r0
    type (gf2d), dimension (4),     intent (in)    :: dr0
    type (gf2d), dimension (4),     intent (inout) :: dr1
    type (gf2d),                    intent (in)    :: beta_l
    type (gf2d),                    intent (inout) :: w_wt, w_l
    CCTK_INT :: a, b

    !------------------------------------------------------------------
    ! w eq. (78)
    !------------------------------------------------------------------

    w_wt%d = (dr0(1)%d * etaup0(1,1)%d &
         - 2.d0 * dr0(4)%d &
         + 2.d0 * dr0(2)%d * etaup0(1,2)%d &
         + 2.d0 * dr0(3)%d * etaup0(1,3)%d &
         + (  dr0(2)%d ** 2 * etaup0(2,2)%d &
         + 2.d0 * dr0(2)%d * dr0(3)%d * etaup0(2,3)%d &
         + dr0(3)%d ** 2 * etaup0(3,3)%d &
         ) / dr0(1)%d - 1.d0 &
         ) / r0%d

    !------------------------------------------------------------------
    ! w_l eq. (79)
    ! radial terms
    !---------------------------------------------------------------

    w_l%d = - dr0(1)%d / r0%d * w_wt%d &
         + ( - 2.d0 * dr1(4)%d + dr0(1)%d * etaup1(1,1)%d &
         - 2.d0 * dr0(1)%d * beta_l%d * etaup0(1,1)%d ) / r0%d

    !---------------------------------------------------------------
    ! angular terms
    !---------------------------------------------------------------

    do a = 2, 3
       w_l%d = w_l%d + 2.d0 * ( dr0(a)%d * etaup1(1,a)%d  &
            + dr1(a)%d * etaup0(1,a)%d ) / r0%d
    end do

    do a = 2, 3
       do b = 2, 3
          w_l%d = w_l%d &
               + ( dr0(a)%d * dr0(b)%d / dr0(1)%d &
               * ( 2.d0 * beta_l%d * etaup0(a,b)%d + etaup1(a,b)%d ) &
               + 2.d0 * dr1(a)%d * dr0(b)%d / dr0(1)%d * etaup0(a,b)%d ) &
               / r0%d
       end do
    end do

  end subroutine wt_wb


  subroutine wt_cb(nn1, nn2, pp, dr0, dr1, qa, beta_l, cb_wt) 

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                           intent (in)    :: nn1, nn2
    CCTK_REAL,    dimension (nn1, nn2), intent (in)    :: pp
    type (gf2d),                        intent (in)    :: beta_l
    type (gf2d),  dimension (4),        intent (in)    :: dr0, dr1
    type (gf2dc), dimension (2:3),      intent (inout) :: qa
    type (gf2dc),                       intent (inout) :: cb_wt

    CCTK_INT :: a
    CCTK_COMPLEX ii
    ii = (0.0d0, 1.0d0)

    qa(2)%d = 0.5d0 * pp
    qa(3)%d = 0.5d0 * pp * ii

    cb_wt%d = (0., 0.d0)

    do a =  2, 3
       cb_wt%d = cb_wt%d -0.5d0* qa(a)%d * dr1(a)%d / dr0(1)%d
    end do

    cb_wt%d = cb_wt%d - 0.5d0*pp*(dr0(2)%d + ii * dr0(3)%d)*beta_l%d/dr0(1)%d

  end subroutine wt_cb


  subroutine wt_nu_ck(nn1, nn2, pp, qs, ps, temp, eta0, deta0,&
             r0, dr0, qa, j_wt, j_l, nu_wt, ck_wt) 

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                             intent (in) :: nn1, nn2
    CCTK_REAL,   dimension (nn1, nn2),    intent (in) :: pp, qs, ps
    type (gf2d), dimension (4,4),         intent (in) :: eta0
    type (gf2d),                          intent (inout) :: temp
    type (gf2d), dimension (2:3,2:3,2:3), intent (inout) :: deta0
    type (gf2d),                          intent (in)    :: r0
    type (gf2d),  dimension (4),          intent (in)    :: dr0 
    type (gf2dc), dimension (2:3),        intent (inout) :: qa
    type (gf2dc),                         intent (in) :: j_wt, j_l
    type (gf2dc),                         intent (inout) :: nu_wt, ck_wt
    CCTK_INT :: a, b, c
    CCTK_COMPLEX ii
    ii = (0.0d0, 1.0d0)

    qa(2)%d = 0.5d0 * pp
    qa(3)%d = 0.5d0 * pp * ii

    nu_wt%d = (0., 0.d0)
    ck_wt%d = (0., 0.d0)

    do c =  2, 3
       do b = 2, 3
          do a = 2, 3

             temp%d = (deta0(a,b,c)%d - 2 * eta0(a,b)%d * dr0(c)%d / r0%d) / r0%d ** 2

             nu_wt%d = nu_wt%d + 0.5d0 * qa(a)%d *qa(b)%d * conjg(qa(c)%d) * temp%d
             ck_wt%d = ck_wt%d + 0.5d0 * qa(a)%d *conjg(qa(b)%d) * qa(c)%d * temp%d

          end do
       end do
    end do

    temp%d = sqrt(1+j_wt%d*conjg(j_wt%d)) ! i.e., store K in temp

    ck_wt%d = ck_wt%d + 2 * dcmplx(qs, ps) * temp%d

    nu_wt%d = nu_wt%d - 0.5d0*pp*(dr0(2)%d + ii * dr0(3)%d)/dr0(1)%d*j_l%d
    ck_wt%d = ck_wt%d - 0.5d0*pp*(dr0(2)%d + ii * dr0(3)%d)/dr0(1)%d &
             *dble(conjg(j_wt%d)*j_l%d)/temp%d

  end subroutine wt_nu_ck


  subroutine wt_eth_expand(cctkGH, gi_max, nx, nn1, nn2, tmp_cgfn, tmp_cgfs,&
                           bcn, bcs, jcn, jcs, cbcn, cbcs, nucn, nucs, ckcn, ckcs)

    use NullInterp
    implicit none

    CCTK_INT, parameter :: izero = 0
    integer, parameter :: ik = kind(izero)

    CCTK_INT :: i, gi_max
    CCTK_INT,                             intent (in)   :: nx, nn1, nn2
    CCTK_COMPLEX, dimension(nn1, nn2),    intent(inout) :: tmp_cgfn, tmp_cgfs

    CCTK_COMPLEX, dimension(nn1, nn2, nx), intent(in) :: jcn, jcs
    CCTK_REAL,    dimension(nn1, nn2, nx), intent(in) :: bcn, bcs

    CCTK_COMPLEX, dimension(nn1, nn2, nx), intent(inout) :: &
         nucn, nucs, ckcn, ckcs, cbcn, cbcs 

    CCTK_POINTER, intent (in) :: cctkGH

    CCTK_COMPLEX, dimension(nn1, nn2, 2) :: F, eth_F

    do i= 1, gi_max

       ! eth beta

       F(:,:,1) = bcn(:,:,i)
       F(:,:,2) = bcs(:,:,i)

       call NullInterp_eth1(cctkGH, tmp_cgfn, tmp_cgfs, eth_F, F, 0_ik, 1_ik)

       cbcn(:,:, i) = eth_F(:,:,1)
       cbcs(:,:, i) = eth_F(:,:,2)

       ! ethb J

       F(:,:,1) = jcn(:,:,i)
       F(:,:,2) = jcs(:,:,i)
       call NullInterp_eth1(cctkGH, tmp_cgfn, tmp_cgfs, eth_F, F, 2_ik, -1_ik)
       nucn(:,:, i) = eth_F(:,:,1)
       nucs(:,:, i) = eth_F(:,:,2)

       ! eth K

       F(:,:,1) = sqrt(1+jcn(:,:,i)*conjg(jcn(:,:,i)))
       F(:,:,2) = sqrt(1+jcs(:,:,i)*conjg(jcs(:,:,i)))
       call NullInterp_eth1(cctkGH, tmp_cgfn, tmp_cgfs, eth_F, F, 0_ik, 1_ik)

       ckcn(:,:, i) = eth_F(:,:,1)
       ckcs(:,:, i) = eth_F(:,:,2)

    end do

  end subroutine  wt_eth_expand

  subroutine apply_cbdata (gi_max, nx, nn1, nn2,&
             Fns, F_wt, F_l,  dr0, rwt, x_wt, x)

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                     intent (in)    :: nx, nn1, nn2
    CCTK_COMPLEX, dimension (nn1,nn2,nx), intent (inout) :: Fns
    type (gf2dc),                  intent (inout) :: F_wt, F_l
    type (gf2d),  dimension (4),  intent (in)    :: dr0
    CCTK_REAL,                    intent (in)    :: rwt
    type (gf2d),                  intent (inout) :: x_wt
    CCTK_REAL,    dimension (nx), intent (in)    :: x

    CCTK_INT :: i, gi_max

    do i=1, gi_max
       Fns(1:nn1,1:nn2,i) = F_wt%d + F_l%d / dr0(1)%d * rwt / (1.d0 - x_wt%d)** 2 * (x(i) - x_wt%d)
    end do

  end subroutine apply_cbdata


  subroutine apply_rbdata (gi_max, nx, nn1, nn2,&
             Fns, F_wt, F_l, dr0, rwt, x_wt, x)

    use NullSHRE_modGFdef
    implicit none

    CCTK_INT,                     intent (in)    :: nx, nn1, nn2
    CCTK_REAL, dimension (nn1,nn2,nx), intent (inout) :: Fns
    type (gf2d),                  intent (inout) :: F_wt, F_l
    type (gf2d),  dimension (4),  intent (in)    :: dr0
    CCTK_REAL,                    intent (in)    :: rwt
    type (gf2d),                  intent (inout) :: x_wt
    CCTK_REAL,    dimension (nx), intent (in)    :: x

    CCTK_INT :: i, gi_max

    do i=1, gi_max

       Fns(1:nn1,1:nn2,i) = F_wt%d + F_l%d / dr0(1)%d * rwt / (1.d0 - x_wt%d)** 2 * (x(i) - x_wt%d)

    end do

  end subroutine apply_rbdata


end module NullSHRE_modBoundary
