! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine NullSHRE_FullMetricRecon(CCTK_ARGUMENTS)

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

!initialize arrays to zero, to enable the summation 
    SHRE_gij     = 0.d0
    SHRE_git     = 0.d0
    SHRE_drgij   = 0.d0
    SHRE_dqgij   = 0.d0
    SHRE_dpgij   = 0.d0
    SHRE_dtgij   = 0.d0
    SHRE_drgit   = 0.d0
    SHRE_dqgit   = 0.d0
    SHRE_dpgit   = 0.d0
    SHRE_dtgit   = 0.d0
    SHRE_beta    = 0.d0
    SHRE_drbeta  = 0.d0
    SHRE_dqbeta  = 0.d0
    SHRE_dpbeta  = 0.d0
    SHRE_dtbeta  = 0.d0
    SHRE_alpha   = 0.d0
    SHRE_dralpha = 0.d0
    SHRE_dqalpha = 0.d0
    SHRE_dpalpha = 0.d0
    SHRE_dtalpha = 0.d0

    do l = 0, l_max
      do m = -l, l

      retval = GetCurrentExtractionCoefs(l, m, RC, IC, RCr, ICr, RCt, ICt)

        n = l*l + l + m + 1
        Nlm = dsqrt(dble(l*(l + 1)))

        do k = 0, 5

!gxx=1(1/2), gxy=2(3/4), gxz=3(5/6), gyy=4(7/8), gyz=5(9,10), gzz=6(11/12)
          SHRE_gij(:,:,2*k+1) = SHRE_gij(:,:,2*k+1)&
                           + (RC(1+k) * reYN(:,:,n) - IC(1+k) * imYN(:,:,n))
          SHRE_gij(:,:,2*k+2) = SHRE_gij(:,:,2*k+2)&
                           + (RC(1+k) * reYS(:,:,n) - IC(1+k) * imYS(:,:,n))

!dr_gxx, dr_gxy, dr_gxz, dr_gyy, dr_gyz, dr_gzz
          SHRE_drgij(:,:,2*k+1) = SHRE_drgij(:,:,2*k+1)&
                           + (RCr(1+k) * reYN(:,:,n) - ICr(1+k) * imYN(:,:,n))
          SHRE_drgij(:,:,2*k+2) = SHRE_drgij(:,:,2*k+2)&
                           + (RCr(1+k) * reYS(:,:,n) - ICr(1+k) * imYS(:,:,n))

!dq_gxx, dq_gxy, dq_gxz, dq_gyy, dq_gyz, dq_gzz
          SHRE_dqgij(:,:,2*k+1) = SHRE_dqgij(:,:,2*k+1)&
                           + 2.d0 * Nlm / pp&
                           * (RC(1+k) * re1YN(:,:,n) - IC(1+k) * im1YN(:,:,n))
          SHRE_dqgij(:,:,2*k+2) = SHRE_dqgij(:,:,2*k+2)&
                           + 2.d0 * Nlm / pp&
                           * (RC(1+k) * re1YS(:,:,n) - IC(1+k) * im1YS(:,:,n))

!dp_gxx, dp_gxy, dp_gxz, dp_gyy, dp_gyz, dp_gzz
          SHRE_dpgij(:,:,2*k+1) = SHRE_dpgij(:,:,2*k+1)&
                           + 2.d0 * Nlm / pp&
                           * (RC(1+k) * im1YN(:,:,n) + IC(1+k) * re1YN(:,:,n))
          SHRE_dpgij(:,:,2*k+2) = SHRE_dpgij(:,:,2*k+2)&
                           + 2.d0 * Nlm / pp&
                           * (RC(1+k) * im1YS(:,:,n) + IC(1+k) * re1YS(:,:,n))

!dt_gxx, dt_gxy, dt_gxz, dt_gyy, dt_gyz, dt_gzz
          SHRE_dtgij(:,:,2*k+1) = SHRE_dtgij(:,:,2*k+1)&
                           + (RCt(1+k) * reYN(:,:,n) - ICt(1+k) * imYN(:,:,n))
          SHRE_dtgij(:,:,2*k+2) = SHRE_dtgij(:,:,2*k+2)&
                           + (RCt(1+k) * reYS(:,:,n) - ICt(1+k) * imYS(:,:,n))

        end do


        do k = 0, 2

!betax=7(1/2), betay=8(3/4), betaz=9(5/6) on both patches
          SHRE_beta(:,:,2*k+1) = SHRE_beta(:,:,2*k+1)&
                           + (RC(7+k) * reYN(:,:,n) - IC(7+k) * imYN(:,:,n))
          SHRE_beta(:,:,2*k+2) = SHRE_beta(:,:,2*k+2)&
                           + (RC(7+k) * reYS(:,:,n) - IC(7+k) * imYS(:,:,n))

!dr_betax, dr_betay, dr_betaz on both patches
          SHRE_drbeta(:,:,2*k+1) = SHRE_drbeta(:,:,2*k+1)&
                           + (RCr(7+k) * reYN(:,:,n) - ICr(7+k) * imYN(:,:,n))
          SHRE_drbeta(:,:,2*k+2) = SHRE_drbeta(:,:,2*k+2)&
                           + (RCr(7+k) * reYS(:,:,n) - ICr(7+k) * imYS(:,:,n))

!dq_betax, dq_betay, dq_betaz on both patches
          SHRE_dqbeta(:,:,2*k+1) = SHRE_dqbeta(:,:,2*k+1)&
                           + 2.d0 * Nlm / pp&
                           * (RC(7+k) * re1YN(:,:,n) - IC(7+k) * im1YN(:,:,n))
          SHRE_dqbeta(:,:,2*k+2) = SHRE_dqbeta(:,:,2*k+2)&
                           + 2.d0 * Nlm / pp&
                           * (RC(7+k) * re1YS(:,:,n) - IC(7+k) * im1YS(:,:,n))

!dp_betax, dp_betay, dp_betaz on both patches
          SHRE_dpbeta(:,:,2*k+1) = SHRE_dpbeta(:,:,2*k+1)&
                           + 2.d0 * Nlm / pp&
                           * (RC(7+k) * im1YN(:,:,n) + IC(7+k) * re1YN(:,:,n))
          SHRE_dpbeta(:,:,2*k+2) = SHRE_dpbeta(:,:,2*k+2)&
                           + 2.d0 * Nlm / pp&
                           * (RC(7+k) * im1YS(:,:,n) + IC(7+k) * re1YS(:,:,n))

!dt_betax, dt_betay, dt_betaz on both patches
          SHRE_dtbeta(:,:,2*k+1) = SHRE_dtbeta(:,:,2*k+1)&
                           + (RCt(7+k) * reYN(:,:,n) - ICt(7+k) * imYN(:,:,n))
          SHRE_dtbeta(:,:,2*k+2) = SHRE_dtbeta(:,:,2*k+2)&
                           + (RCt(7+k) * reYS(:,:,n) - ICt(7+k) * imYS(:,:,n))

        end do

!alpha on both patches          
        SHRE_alpha(:,:,1) = SHRE_alpha(:,:,1)&
                          + (RC(10) * reYN(:,:,n) - IC(10) * imYN(:,:,n))
        SHRE_alpha(:,:,2) = SHRE_alpha(:,:,2)&
                          + (RC(10) * reYS(:,:,n) - IC(10) * imYS(:,:,n))

!dr_alpha on both patches          
        SHRE_dralpha(:,:,1) = SHRE_dralpha(:,:,1)&
                            + (RCr(10) * reYN(:,:,n) - ICr(10) * imYN(:,:,n))
        SHRE_dralpha(:,:,2) = SHRE_dralpha(:,:,2)&
                            + (RCr(10) * reYS(:,:,n) - ICr(10) * imYS(:,:,n))

!dq_alpha on both patches          
        SHRE_dqalpha(:,:,1) = SHRE_dqalpha(:,:,1)&
                            + 2.d0 * Nlm / pp&
                            * (RC(10) * re1YN(:,:,n) - IC(10) * im1YN(:,:,n))
        SHRE_dqalpha(:,:,2) = SHRE_dqalpha(:,:,2)&
                            + 2.d0 * Nlm / pp&
                            * (RC(10) * re1YS(:,:,n) - IC(10) * im1YS(:,:,n))

!dp_alpha on both patches          
        SHRE_dpalpha(:,:,1) = SHRE_dpalpha(:,:,1)&
                            + 2.d0 * Nlm / pp&
                            * (RC(10) * im1YN(:,:,n) + IC(10) * re1YN(:,:,n))
        SHRE_dpalpha(:,:,2) = SHRE_dpalpha(:,:,2)&
                            + 2.d0 * Nlm / pp&
                            * (RC(10) * im1YS(:,:,n) + IC(10) * re1YS(:,:,n))

!dt_alpha on both patches          
        SHRE_dtalpha(:,:,1) = SHRE_dtalpha(:,:,1)&
                            + (RCt(10) * reYN(:,:,n) - ICt(10) * imYN(:,:,n))
        SHRE_dtalpha(:,:,2) = SHRE_dtalpha(:,:,2)&
                            + (RCt(10) * reYS(:,:,n) - ICt(10) * imYS(:,:,n))

      end do
    end do

 end subroutine NullSHRE_FullMetricRecon
