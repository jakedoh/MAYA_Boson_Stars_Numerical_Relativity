C     This routine sets the lapse and/or shift by calling a routine
C     that does it pointwise. Note that it could be easily modified
C     to set the Bona-Masso variables B_xx etc. 
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine Exact__gauge(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS_EXACT__GAUGE
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k
      integer nx,ny,nz
      logical set_lapse, set_dtlapse, set_shift, set_dtshift

      CCTK_REAL tt, xx, yy, zz
      CCTK_REAL gxxtmp, gyytmp, gzztmp, 
     $     gxytmp, gyztmp, gxztmp,
     $     hxxtmp, hyytmp, hzztmp, 
     $     hxytmp, hyztmp, hxztmp,
     $     dxgxxtmp, dxgyytmp, dxgzztmp, 
     $     dxgxytmp, dxgyztmp, dxgxztmp,
     $     dygxxtmp, dygyytmp, dygzztmp, 
     $     dygxytmp, dygyztmp, dygxztmp,
     $     dzgxxtmp, dzgyytmp, dzgzztmp, 
     $     dzgxytmp, dzgyztmp, dzgxztmp,
     $     alptmp, dtalptmp, axtmp, aytmp, aztmp, 
     $     betaxtmp, betaytmp, betaztmp,
     $     dtbetaxtmp, dtbetaytmp, dtbetaztmp,
     $     bxxtmp, bxytmp, bxztmp, 
     $     byxtmp, byytmp, byztmp, 
     $     bzxtmp, bzytmp, bzztmp
      CCTK_REAL
     $     exact_psi,
     $     exact_psix, exact_psiy, exact_psiz,
     $     exact_psixx, exact_psiyy, exact_psizz,
     $     exact_psixy, exact_psiyz, exact_psixz
      LOGICAL is_initial_slice, is_later_slice

C     are we on the initial slice or some later slice?
C     n.b. the logical expressions later in this function involving
C          these flags below would be *so* much nicer if Fortran
C          grokked C-style conditional expressions... :) :)
      is_initial_slice = cctk_iteration .eq. 0
      is_later_slice   = cctk_iteration .ne. 0

C     Grid parameters.
      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

C     This code used to set  t = time + dt/2  to get 2nd order accuracy,
C     but this leads to the initial data being set at the wrong time. :(
C     In the context of MoL, we want to set variables at the standard Cactus
C     time (cctk_time), because MoL takes care of calling us at each MoL
C     iteration, and updating the field variables appropriately.
C
C     Alas, setting at cctk_time probably gives O(dt) errors for non-MoL
C     evoutions where Exact is used to set stuff at each time step.
C     Fixing this [unless we just declare all non-MoL stuff obselete :) ]
C     probably requires cleaning up our (++messy)  schedule.ccl , which
C     is why this remains a bug for now... :( :(

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C     Set lapse and/or shift?
C
      if ( is_initial_slice ) then
         set_lapse   = CCTK_Equals(initial_lapse,   "exact").ne.0
         set_shift   = CCTK_Equals(initial_shift,   "exact").ne.0
         set_dtlapse = CCTK_Equals(initial_dtlapse, "exact").ne.0
         set_dtshift = CCTK_Equals(initial_dtshift, "exact").ne.0
      end if
      if ( is_later_slice ) then
         set_lapse   = CCTK_Equals(lapse_evolution_method, "exact").ne.0
         set_shift   = CCTK_Equals(shift_evolution_method, "exact").ne.0
         set_dtlapse = CCTK_Equals(dtlapse_evolution_method, "exact").ne.0
         set_dtshift = CCTK_Equals(dtshift_evolution_method, "exact").ne.0
      end if

      if ( set_lapse .or. set_shift .or. set_dtlapse .or. set_dtshift) then

C$omp    parallel do private (
C$omp$     i, j, k,
C$omp$     tt, xx, yy, zz,
C$omp$     alptmp, dtalptmp, axtmp, aytmp, aztmp, 
C$omp$     betaxtmp, betaytmp, betaztmp,
C$omp$     dtbetaxtmp, dtbetaytmp, dtbetaztmp,
C$omp$     bxxtmp, bxytmp, bxztmp, 
C$omp$     byxtmp, byytmp, byztmp, 
C$omp$     bzxtmp, bzytmp, bzztmp,
C$omp$     dxgxxtmp, dxgyytmp, dxgzztmp, 
C$omp$     dxgxytmp, dxgyztmp, dxgxztmp,
C$omp$     dygxxtmp, dygyytmp, dygzztmp, 
C$omp$     dygxytmp, dygyztmp, dygxztmp,
C$omp$     dzgxxtmp, dzgyytmp, dzgzztmp, 
C$omp$     dzgxytmp, dzgyztmp, dzgxztmp,
C$omp$     exact_psi,
C$omp$     exact_psix, exact_psiy, exact_psiz,
C$omp$     exact_psixx, exact_psiyy, exact_psizz,
C$omp$     exact_psixy, exact_psiyz, exact_psixz)
         do k=1,nz
            do j=1,ny
               do i=1,nx

                  tt = cctk_time
                  xx = x(i,j,k) - cctk_time * shift_add_x
                  yy = y(i,j,k) - cctk_time * shift_add_y
                  zz = z(i,j,k) - cctk_time * shift_add_z

C     Initialize the psi of exact
C     (also to tell the models about the conformal_state)
                  if (conformal_state .ne. 0) then
                    exact_psi = 1.0D0
                  else
                    exact_psi = 0.0D0
                  end if
                  exact_psix = 0.0D0
                  exact_psiy = 0.0D0
                  exact_psiz = 0.0D0
                  exact_psixx = 0.0D0
                  exact_psiyy = 0.0D0
                  exact_psizz = 0.0D0
                  exact_psixy = 0.0D0
                  exact_psiyz = 0.0D0
                  exact_psixz = 0.0D0

                  call Exact__Bona_Masso_data(
     $                 decoded_exact_model,
     $                 xx, yy, zz, tt,
     $                 gxxtmp, gyytmp, gzztmp, 
     $                 gxytmp, gyztmp, gxztmp,
     $                 hxxtmp, hyytmp, hzztmp, 
     $                 hxytmp, hyztmp, hxztmp,
     $                 exact_psi,
     $                 exact_psix, exact_psiy, exact_psiz,
     $                 exact_psixx, exact_psiyy, exact_psizz,
     $                 exact_psixy, exact_psiyz, exact_psixz,
     $                 dxgxxtmp, dxgyytmp, dxgzztmp, 
     $                 dxgxytmp, dxgyztmp, dxgxztmp,
     $                 dygxxtmp, dygyytmp, dygzztmp, 
     $                 dygxytmp, dygyztmp, dygxztmp,
     $                 dzgxxtmp, dzgyytmp, dzgzztmp, 
     $                 dzgxytmp, dzgyztmp, dzgxztmp,
     $                 alptmp, dtalptmp, axtmp, aytmp, aztmp, 
     $                 betaxtmp, betaytmp, betaztmp,
     $                 dtbetaxtmp, dtbetaytmp, dtbetaztmp,
     $                 bxxtmp, bxytmp, bxztmp, 
     $                 byxtmp, byytmp, byztmp, 
     $                 bzxtmp, bzytmp, bzztmp)

                  if ( set_lapse ) then
                     alp(i,j,k) = alptmp
                  end if
                  if ( set_shift ) then
                     betax(i,j,k) = betaxtmp + shift_add_x
                     betay(i,j,k) = betaytmp + shift_add_y
                     betaz(i,j,k) = betaztmp + shift_add_z
                  end if
                  if ( set_dtlapse ) then
                     dtalp(i,j,k) = dtalptmp
                  end if
                  if ( set_dtshift ) then
                     dtbetax(i,j,k) = dtbetaxtmp
                     dtbetay(i,j,k) = dtbetaytmp
                     dtbetaz(i,j,k) = dtbetaztmp
                  end if
               end do
            end do
         end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      else
         call CCTK_WARN(1,'Exact__gauge has been called without doing anything')
      end if

      return
      end
