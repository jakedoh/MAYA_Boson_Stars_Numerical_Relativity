C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine Exact__blended_boundary(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      logical doKij, doGij, doLapse, doShift

      integer i,j,k
      integer nx,ny,nz

      CCTK_REAL router, rinner, frac, onemfrac
      CCTK_REAL gxxe, gyye, gzze, gxye, gyze, gxze
      CCTK_REAL kxxe, kyye, kzze, kxye, kyze, kxze
      CCTK_REAL dxgxxe, dxgyye, dxgzze, dxgxye, dxgyze, dxgxze
      CCTK_REAL dygxxe, dygyye, dygzze, dygxye, dygyze, dygxze
      CCTK_REAL dzgxxe, dzgyye, dzgzze, dzgxye, dzgyze, dzgxze

      CCTK_REAL
     $     exact_psi,
     $     exact_psix, exact_psiy, exact_psiz,
     $     exact_psixx, exact_psiyy, exact_psizz,
     $     exact_psixy, exact_psiyz, exact_psixz

      CCTK_REAL alpe, dtalpe, axe, aye, aze
      CCTK_REAL betaxe,betaye,betaze, dtbetaxe,dtbetaye,dtbetaze
      CCTK_REAL bxxe,bxye,bxze,byxe,byye,byze,bzxe,bzye,bzze

      CCTK_REAL det, uxx, uxy, uxz, uyy, uyz, uzz

      CCTK_REAL dx,dy,dz,time
      CCTK_REAL xmx,xmn,ymx,ymn,zmx,zmn,rmx
      integer ierr

C     Grid parameters.

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      dx = cctk_delta_space(1)
      dy = cctk_delta_space(2)
      dz = cctk_delta_space(3)

      time = cctk_time

C     Other parameters.

      doKij = (exblend_Ks.eq.1)
      doGij = (exblend_gs.eq.1)

      doLapse  = ((exblend_gauge.eq.1).and.
     $            (CCTK_Equals(lapse_evolution_method,"exact").ne.0))
      doShift  = ((exblend_gauge.eq.1).and.
     $            (CCTK_Equals(shift_evolution_method,"exact").ne.0))

      call CCTK_CoordRange(ierr,cctkGH,xmn,xmx,-1,"x","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,ymn,ymx,-1,"y","cart3d")
      call CCTK_CoordRange(ierr,cctkGH,zmn,zmx,-1,"z","cart3d")

      rmx = min(xmx,ymx,zmx)

      if (exblend_rout.lt.0) then 
         router = rmx - 2.0d0*dx
      else
         router = exblend_rout
      endif

      if (exblend_width.lt.0) then 
         rinner = router + exblend_width*dx
      else
         rinner = router - exblend_width
      endif

      do k=1,nz
         do j=1,ny
            do i=1,nx

c     We only do anything if r >= rinner so only evaluate exact data 
c     there.

               if (r(i,j,k) .ge. rinner) then

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
     $                 x(i,j,k), y(i,j,k), z(i,j,k), time,
     $                 gxxe, gyye, gzze, gxye, gyze, gxze,
     $                 kxxe, kyye, kzze, kxye, kyze, kxze,
     $                 exact_psi,
     $                 exact_psix, exact_psiy, exact_psiz,
     $                 exact_psixx, exact_psiyy, exact_psizz,
     $                 exact_psixy, exact_psiyz, exact_psixz,
     $                 dxgxxe, dxgyye, dxgzze, dxgxye, dxgyze, dxgxze,
     $                 dygxxe, dygyye, dygzze, dygxye, dygyze, dygxze,
     $                 dzgxxe, dzgyye, dzgzze, dzgxye, dzgyze, dzgxze,
     $                 alpe, dtalpe, axe, aye, aze,
     $                 betaxe, betaye, betaze, dtbetaxe, dtbetaye, dtbetaze,
     $                 bxxe, bxye, bxze, byxe, 
     $                 byye, byze, bzxe, bzye, bzze)

c     This sucks, but we want the exact vs so we can blend them also.

                  det = -(gxze**2*gyye) 
     &                 + 2.d0*gxye*gxze*gyze 
     &                 - gxxe*gyze**2 
     &                 - gxye**2*gzze 
     &                 + gxxe*gyye*gzze
    
                  uxx=(-gyze**2 + gyye*gzze)/det
                  uxy=(gxze*gyze - gxye*gzze)/det
                  uyy=(-gxze**2 + gxxe*gzze)/det
                  uxz=(-gxze*gyye + gxye*gyze)/det
                  uyz=(gxye*gxze - gxxe*gyze)/det
                  uzz=(-gxye**2 + gxxe*gyye)/det

c     Outside of router we want to place exact data on our grid

                  if (r(i,j,k) .gt. router) then

c     This is one of those things I will invariably screw up if I type
c     it in so let the computer do it

#define exassign(q) q(i,j,k) = q e
#define exassign_grp(p) \
        exassign(p xx) &&\
        exassign(p xy) &&\
        exassign(p xz) &&\
        exassign(p yy) &&\
        exassign(p yz) &&\
        exassign(p zz)

c     Note this plays on the nasty trick that fortran doesnt give a
c     damn about spaces so gxx e is the same as gxxe for the parser...
c     Grody but effective!

                    if (doGij) then 
                       exassign_grp(g)
                    endif
                    if (doKij) then 
                       exassign_grp(k)
                    endif

                    if (doLapse) then 
                       exassign(alp)
                    endif
                    
                    if (doShift.and.(shift_state.ne.0)) then 
                       exassign(betax)
                       exassign(betay)
                       exassign(betaz)
                    endif

c     OK so we dont want to blend so use a goto to jump.
                 else

c     Evaluate the linear weighting fraction. Obvious...

                    frac = (r(i,j,k) - rinner) / (router - rinner)
                    onemfrac = 1.0D0 - frac

c     Once again some c-preprocessor tricks based on the whole fortran
c     space thing...

#define INTPOINT(f,v) f(i,j,k) = frac * v + onemfrac * f(i,j,k)
#define intone(f) INTPOINT(f, f e)
#define int_grp(p) \
        intone(p xx) &&\
        intone(p xy) &&\
        intone(p xz) &&\
        intone(p yy) &&\
        intone(p yz) &&\
        intone(p zz)

                    if (doGij) then 
                       int_grp(g)
                    endif
                    if (doKij) then
                       int_grp(k)
                    endif

                    if (doLapse) then 
                       intone(alp)
                    endif

                    if (doShift.and.(shift_state.ne.0)) then 
                       intone(betax)
                       intone(betay)
                       intone(betaz)
                    endif

                 endif          ! r > router else
              endif             ! r > rinner

            enddo
         enddo
      enddo

      return
      end 
