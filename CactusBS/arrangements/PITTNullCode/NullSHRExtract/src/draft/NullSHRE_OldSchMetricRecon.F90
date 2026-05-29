! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine NullSHRE_SchMetricRecon(CCTK_ARGUMENTS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!reads the numerical expansion coefficients Clm
!and reconstructs the metric on the worldtube 
!with the spherical harmonics in (p,q) coord.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use cctk
  use NullGrid_Vars
  implicit none

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_FUNCTIONS
    CCTK_INT :: n, l, m, k, retval
    CCTK_REAL :: Nlm

    SHRE_gij = 0.d0
    SHRE_git = 0.d0
    SHRE_drgij = 0.d0
    SHRE_dqgij = 0.d0
    SHRE_dpgij = 0.d0
    SHRE_dtgij = 0.d0
    SHRE_beta = 0.d0
    SHRE_drbeta = 0.d0
    SHRE_dqbeta = 0.d0
    SHRE_dpbeta = 0.d0
    SHRE_dtbeta = 0.d0
    SHRE_alpha = 0.d0
    SHRE_dralpha = 0.d0
    SHRE_dqalpha = 0.d0
    SHRE_dpalpha = 0.d0
    SHRE_dtalpha = 0.d0
    SHRE_drgit = 0.d0
    SHRE_dqgit = 0.d0
    SHRE_dpgit = 0.d0
    SHRE_dtgit = 0.d0

    do l = 0, lmax
      do m = -l, l

        n = l*l + l + m + 1
        Nlm = dsqrt(dble(l*(l + l)))

        do k = 0, 5

!gxx=1(1/2), gxy=2(3/4), gxz=3(5/6), gyy=4(7/8), gyz=5(9,10), gzz=6(11/12)
          SHRE_gij(:,:,2*k+1) = SHRE_gij(:,:,2*k+1)&
                           + (reTC(1+k,n) * reTYN(:,:,n) - imTC(1+k,n) * imTYN(:,:,n))
          SHRE_gij(:,:,2*k+2) = SHRE_gij(:,:,2*k+2)&
                           + (reTC(1+k,n) * reTYS(:,:,n) - imTC(1+k,n) * imTYS(:,:,n))

!dr_gxx, dr_gxy, dr_gxz, dr_gyy, dr_gyz, dr_gzz
          SHRE_drgij(:,:,2*k+1) = SHRE_drgij(:,:,2*k+1)&
                           + (reTCr(1+k,n) * reTYN(:,:,n) - imTCr(1+k,n) * imTYN(:,:,n))
          SHRE_drgij(:,:,2*k+2) = SHRE_drgij(:,:,2*k+2)&
                           + (reTCr(1+k,n) * reTYS(:,:,n) - imTCr(1+k,n) * imTYS(:,:,n))

!dq_gxx, dq_gxy, dq_gxz, dq_gyy, dq_gyz, dq_gzz
          SHRE_dqgij(:,:,2*k+1) = SHRE_dqgij(:,:,2*k+1) + 2.d0 * Nlm / pp&
                           * (reTC(1+k,n) * re1TYN(:,:,n) - imTC(1+k,n) * im1TYN(:,:,n))
          SHRE_dqgij(:,:,2*k+2) = SHRE_dqgij(:,:,2*k+2) + 2.d0 * Nlm / pp&
                           * (reTC(1+k,n) * re1TYS(:,:,n) - imTC(1+k,n) * im1TYS(:,:,n))

!dp_gxx, dp_gxy, dp_gxz, dp_gyy, dp_gyz, dp_gzz
          SHRE_dpgij(:,:,2*k+1) = SHRE_dpgij(:,:,2*k+1) + 2.d0 * Nlm / pp&
                           * (reTC(1+k,n) * im1TYN(:,:,n) + imTC(1+k,n) * re1TYN(:,:,n))
          SHRE_dpgij(:,:,2*k+2) = SHRE_dpgij(:,:,2*k+2) + 2.d0 * Nlm / pp&
                           * (reTC(1+k,n) * im1TYS(:,:,n) + imTC(1+k,n) * re1TYS(:,:,n))

!dt_gxx, dt_gxy, dt_gxz, dt_gyy, dt_gyz, dt_gzz
          SHRE_dtgij(:,:,2*k+1) = SHRE_dtgij(:,:,2*k+1)&
                           + (reTCt(1+k,n) * reTYN(:,:,n) - imTCt(1+k,n) * imTYN(:,:,n))
          SHRE_dtgij(:,:,2*k+2) = SHRE_dtgij(:,:,2*k+2)&
                           + (reTCt(1+k,n) * reTYS(:,:,n) - imTCt(1+k,n) * imTYS(:,:,n))

        end do


        do k = 0, 2

!betax=7(1/2), betay=8(3/4), betaz=9(5/6) on both patches
          SHRE_beta(:,:,2*k+1) = SHRE_beta(:,:,2*k+1)&
                           + (reTC(7+k,n) * reTYN(:,:,n) - imTC(7+k,n) * imTYN(:,:,n))
          SHRE_beta(:,:,2*k+2) = SHRE_beta(:,:,2*k+2)&
                           + (reTC(7+k,n) * reTYS(:,:,n) - imTC(7+k,n) * imTYS(:,:,n))

!dr_betax, dr_betay, dr_betaz on both patches
          SHRE_drbeta(:,:,2*k+1) = SHRE_drbeta(:,:,2*k+1)&
                           + (reTCr(7+k,n) * reTYN(:,:,n) - imTCr(7+k,n) * imTYN(:,:,n))
          SHRE_drbeta(:,:,2*k+2) = SHRE_drbeta(:,:,2*k+2)&
                           + (reTCr(7+k,n) * reTYS(:,:,n) - imTCr(7+k,n) * imTYS(:,:,n))

!dq_betax, dq_betay, dq_betaz on both patches
          SHRE_dqbeta(:,:,2*k+1) = SHRE_dqbeta(:,:,2*k+1) + 2.d0 * Nlm / pp&
                           * (reTC(7+k,n) * re1TYN(:,:,n) - imTC(7+k,n) * im1TYN(:,:,n))
          SHRE_dqbeta(:,:,2*k+2) = SHRE_dqbeta(:,:,2*k+2) + 2.d0 * Nlm / pp&
                           * (reTC(7+k,n) * re1TYS(:,:,n) - imTC(7+k,n) * im1TYS(:,:,n))

!dp_betax, dp_betay, dp_betaz on both patches
          SHRE_dpbeta(:,:,2*k+1) = SHRE_dpbeta(:,:,2*k+1) + 2.d0 * Nlm / pp&
                           * (reTC(7+k,n) * im1TYN(:,:,n) + imTC(7+k,n) * re1TYN(:,:,n))
          SHRE_dpbeta(:,:,2*k+2) = SHRE_dpbeta(:,:,2*k+2) + 2.d0 * Nlm / pp&
                           * (reTC(7+k,n) * im1TYS(:,:,n) + imTC(7+k,n) * re1TYS(:,:,n))

!dt_betax, dt_betay, dt_betaz on both patches
          SHRE_dtbeta(:,:,2*k+1) = SHRE_dtbeta(:,:,2*k+1)&
                           + (reTCt(7+k,n) * reTYN(:,:,n) - imTCt(7+k,n) * imTYN(:,:,n))
          SHRE_dtbeta(:,:,2*k+2) = SHRE_dtbeta(:,:,2*k+2)&
                           + (reTCt(7+k,n) * reTYS(:,:,n) - imTCt(7+k,n) * imTYS(:,:,n))

        end do

!alpha on both patches          
        SHRE_alpha(:,:,1) = SHRE_alpha(:,:,1)&
                          + (reTC(10,n) * reTYN(:,:,n) - imTC(10,n) * imTYN(:,:,n))
        SHRE_alpha(:,:,2) = SHRE_alpha(:,:,2)&
                          + (reTC(10,n) * reTYS(:,:,n) - imTC(10,n) * imTYS(:,:,n))

!dr_alpha on both patches          
        SHRE_dralpha(:,:,1) = SHRE_dralpha(:,:,1)&
                            + (reTCr(10,n) * reTYN(:,:,n) - imTCr(10,n) * imTYN(:,:,n))
        SHRE_dralpha(:,:,2) = SHRE_dtalpha(:,:,2)&
                            + (reTCr(10,n) * reTYS(:,:,n) - imTCr(10,n) * imTYS(:,:,n))

!dq_alpha on both patches          
        SHRE_dqalpha(:,:,1) = SHRE_dqalpha(:,:,1)  + 2.d0 * Nlm / pp&
                            * (reTC(10,n) * re1TYN(:,:,n) - imTC(10,n) * im1TYN(:,:,n))
        SHRE_dqalpha(:,:,2) = SHRE_dqalpha(:,:,2)&
                            * (reTC(10,n) * re1TYS(:,:,n) - imTC(10,n) * im1TYS(:,:,n))

!dp_alpha on both patches          
        SHRE_dpalpha(:,:,1) = SHRE_dpalpha(:,:,1) + 2.d0 * Nlm / pp&
                            * (reTC(10,n) * im1TYN(:,:,n) + imTC(10,n) * re1TYN(:,:,n))
        SHRE_dpalpha(:,:,2) = SHRE_dpalpha(:,:,2)&
                            * (reTC(10,n) * im1TYS(:,:,n) + imTC(10,n) * re1TYS(:,:,n))

!dt_alpha on both patches          
        SHRE_dtalpha(:,:,1) = SHRE_dtalpha(:,:,1)&
                            + (reTCt(10,n) * reTYN(:,:,n) - imTCt(10,n) * imTYN(:,:,n))
        SHRE_dtalpha(:,:,2) = SHRE_dtalpha(:,:,2)&
                            + (reTCt(10,n) * reTYS(:,:,n) - imTCt(10,n) * imTYS(:,:,n))

      end do
    end do

 end subroutine NullSHRE_SchMetricRecon

