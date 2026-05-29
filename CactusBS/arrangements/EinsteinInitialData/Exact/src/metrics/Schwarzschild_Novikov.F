c     The metric given here corresponds to the novikov solution
c     in isotropic coordinates, as presented first in Bruegman96 
c     then in correct form in Cactus paper 1. This code is the code
c     which was used for the comparisons in cactus paper 1, and is written
c     by PW with input from BB.
C
C Author: unknown
C Copyright/License: unknown
C
C $Header$


#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Schwarzschild_Novikov(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x, y, z, t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL eps, mass
      CCTK_REAL r,c,psi4

      CCTK_REAL nov_dr_drmax, nov_rmax, nov_r
      CCTK_REAL grr, gqq, detg

      CCTK_REAL psi4_o_r2

c constants
      CCTK_REAL zero,one,two
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution.

      eps = Schwarzschild_Novikov__epsilon
      mass= Schwarzschild_Novikov__mass

      r = max(sqrt(x**2 + y**2 + z**2), eps)

c     Find r.
      r = sqrt(x**2 + y**2 + z**2)

c     Find conformal factor.
      c = mass/(two*r)
      psi4 = (one + c)**4

c     Evaluate novikov stuff. Note abs(t) since the data is time
c     symmetric (the metric is, at least...)
      grr = nov_dr_drmax(abs(t),abs(r))
      gqq = nov_r(abs(t),abs(r))

      grr = grr **2
      gqq = gqq**2 / nov_rmax(abs(r))**2


c     Find metric components.
      psi4_o_r2 = psi4 / r**2

      gdtt = - 1.0D0

      gdtx = zero
      gdty = zero
      gdtz = zero

c     This is just straightforward spherical -> cartesian I hope... ;-)
c     Note at t=0 (grr = gqq = 1) this gives the expected result
c     (namely diagonal psi^4, since psi4_o_r2 = psi^4 / r^2)
      gdxx = (grr * x**2 + gqq * (y**2 + z**2)) * psi4_o_r2
      gdyy = (grr * y**2 + gqq * (x**2 + z**2)) * psi4_o_r2
      gdzz = (grr * z**2 + gqq * (x**2 + y**2)) * psi4_o_r2
 
      gdxy = (grr - gqq) * x * y * psi4_o_r2
      gdzx = (grr - gqq) * x * z * psi4_o_r2
      gdyz = (grr - gqq) * y * z * psi4_o_r2

c     Find inverse metric.
      gutt = one/gdtt
      gutx = zero
      guty = zero
      gutz = zero
      detg = gdtt*gdxx*gdyy*gdzz-gdtt*gdxx*gdyz**2/4.D0-gdtt*gdxy**2*gdz
     $z/4.D0+gdtt*gdxy*gdzx*gdyz/4.D0-gdtt*gdzx**2*gdyy/4.D0
      guxx = -gdtt*(-4.D0*gdyy*gdzz+gdyz**2)/(4.D0 * detg)
      guxy = -gdtt*(2.D0*gdxy*gdzz-gdzx*gdyz)/(4.D0 * detg)
      guzx = gdtt*(gdxy*gdyz-2.D0*gdzx*gdyy)/(4.D0 * detg)
      guyy = gdtt*(4.D0*gdxx*gdzz-gdzx**2)/(4.D0 * detg)
      guyz = -gdtt*(2.D0*gdxx*gdyz-gdxy*gdzx)/(4.D0 * detg)
      guzz = gdtt*(4.D0*gdxx*gdyy-gdxy**2)/(4.D0 * detg)

      guxx = one/psi4
      guyy = one/psi4
      guzz = one/psi4

      guxy = zero
      guyz = zero
      guzx = zero

      return
      end

c     These are functions which evaluate the novikov stuff.

c     dr/drmax
      CCTK_REAL function nov_dr_drmax(tauin,rbarin)

      implicit none

      CCTK_REAL rbarin, tauin
      CCTK_REAL rt, nov_r, rmaxt, nov_rmax

      rt = nov_r(tauin, rbarin)
      rmaxt = nov_rmax(rbarin)
      nov_dr_drmax = 1.5D0 - rt / (2.0D0 * rmaxt) + 
     $     1.5D0 * sqrt(rmaxt / rt - 1.0D0) *
     $     acos(sqrt(rt/rmaxt))

      return
      end



c
c     Bisection to invert the function below. This is pretty crappy
c     but it works.
c
      CCTK_REAL function nov_r(tauin, rbarin)
      implicit none
c     input
      CCTK_REAL tauin, rbarin

c     funtions
      CCTK_REAL nov_rmax, nov_tau

c     temps
      CCTK_REAL rg, drg, delt, ttmp, rmt
      CCTK_REAL eps
      integer nit
      nit = 0
      delt = 1000.0D0
      rmt = nov_rmax(rbarin)
      rg  = rmt
      drg = rg / 2.0D0
      eps = 1.d-6 * rmt
      do while (delt .gt. eps .and. nit .lt. 100)
         ttmp = nov_tau(rg, rmt)
         delt = abs(tauin - ttmp)
         if (delt .gt. eps) then 
            if (ttmp .gt. tauin .or. rg .lt. drg) then 
               rg = rg + drg
c     Enforce upper bound
               if (rg .gt. rmt) rg = rmt
               drg = drg / 2.0D0
            else
               rg = rg - drg
            endif
         endif
c         write (*,*) rg, ttmp, tauin
         nit = nit + 1
      enddo
      if (nit .ge. 100) then
         write (*,*) "Novikov: inversion did not converge"
      endif
      nov_r = rg
      return
      end

c     Evaluate tau as a function of r and rmax
      CCTK_REAL function nov_tau(r, rmax)
      implicit none
      CCTK_REAL r, rmax

      nov_tau=  rmax * sqrt(0.5D0 * r * (1.0D0 - r / rmax)) +
     $     2.0D0 * (rmax / 2)**(3.0/2.0) *
     $     acos (sqrt(r/rmax))

      return
      end

c     Evaluate rmax as a function of rbar
      CCTK_REAL function nov_rmax(rbar)
      implicit none
      CCTK_REAL rbar
      nov_rmax = (1.0D0 + 2.0D0*rbar)**2 / (4.0D0 * rbar)
      return
      end
