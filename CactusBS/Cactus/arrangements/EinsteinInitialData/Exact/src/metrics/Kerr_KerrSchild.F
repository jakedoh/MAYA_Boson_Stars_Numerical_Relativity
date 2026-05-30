C     Kerr-Schild form of boosted rotating black hole.
C     Program g_ab = eta_ab + H l_a l_b, g^ab = eta^ab - H l^a l^b.
C     Here eta_ab is Minkowski in Cartesian coordinates, H is a scalar,
C     and l is a null vector.
C
C Author: unknown
C Copyright/License: unknown
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

      subroutine Exact__Kerr_KerrSchild(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

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
      CCTK_REAL boostv, eps, m, a
      integer   power

c local variables
      CCTK_REAL gamma, t0, z0, x0, y0, rho02, r02, r0, costheta0, 
     $     lt0, lx0, ly0, lz0, hh, lt, lx, ly, lz

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution,
C     and convert from parameter file spin parameter J/m^2
C     to the J/m definition used in the code here.

      boostv = Kerr_KerrSchild__boost_v
      eps    = Kerr_KerrSchild__epsilon
      power  = Kerr_KerrSchild__power
      m      = Kerr_KerrSchild__mass
      a      = m*Kerr_KerrSchild__spin

C     Boost factor.

      gamma = 1.d0 / sqrt(1.d0 - boostv**2)

C     Lorentz transform t,x,y,z -> t0,x0,y0,z0. 
C     t0 is never used, but is here for illustration, and we introduce
C     x0 and y0 also only for clarity.
C     Note that z0 = 0 means z = vt for the BH.

      t0 = gamma * ((t - Kerr_KerrSchild__t) - boostv * (z - Kerr_KerrSchild__z))
      z0 = gamma * ((z - Kerr_KerrSchild__z) - boostv * (t - Kerr_KerrSchild__t))
      x0 = x - Kerr_KerrSchild__x
      y0 = y - Kerr_KerrSchild__y

C     Coordinate distance to center of black hole. Note it moves!

      rho02 = x0**2 + y0**2 + z0**2

C     Spherical auxiliary coordinate r and angle theta in BH rest frame.

      r02 = 0.5d0 * (rho02 - a**2) 
     $     + sqrt(0.25d0 * (rho02 - a**2)**2 + a**2 * z0**2)
      r0 = sqrt(max(0.0d0,r02))
      if (Kerr_KerrSchild__parabolic .eq. 0) then
C        Use a power law to avoid the singularity
         r0 = (r0**power + eps**power)**(1.0d0/power)
      else
         if (r0 .lt. eps) then
            if (power .eq. 0) then
               r0 = eps
            else if (power .eq. 2) then
               r0 = eps/2 + r0**2 * 1/(2*eps)
            else if (power .eq. 4) then
               r0 = 3*eps/8 + r0**2 * (3/(4*eps) - r0**2 * 1/(8*eps**3))
            else if (power .eq. 6) then
               r0 = 5*eps/16 + r0**2 * (15/(16*eps) + r0**2 * (-5/(16*eps**3) + r0**2 * 1/(16*eps**5)))
            else if (power .eq. 8) then
               r0 = 35*eps/128 + r0**2 * (35/(32*eps) + r0**2 * (-35/(64*eps**3) + r0**2 * (7/(32*eps**5) - r0**2 * 5/(128*eps**7))))
            else
               call CCTK_WARN (CCTK_WARN_ABORT, "Unsupported value of parameter Kerr_KerrSchild__power")
            end if
         end if
      end if
C     Another idea:
C        r0 = r0 + eps * exp(-x/eps)

      costheta0 = z0 / r0
      
C     Coefficient H. Note this transforms as a scalar, so does not carry
C     the suffix 0.
      hh = m * r0 / (r0**2 + a**2 * costheta0**2)
      
C     Components of l_a in rest frame. Note indices down.
      lt0 = 1.d0
      lx0 = (r0 * x0 + a * y0) / (r0**2 + a**2)
      ly0 = (r0 * y0 - a * x0) / (r0**2 + a**2)
      lz0 = z0 / r0

C     Now boost it to coordinates x, y, z, t.
C     This is the reverse Lorentz transformation, but applied
C     to a one-form, so the sign of boostv is the same as the forward
C     Lorentz transformation applied to the coordinates.

      lt = gamma * (lt0 - boostv * lz0) 
      lz = gamma * (lz0 - boostv * lt0) 
      lx = lx0
      ly = ly0
      
C     Down metric. g_ab = flat_ab + H l_a l_b

      gdtt = - 1.d0 + 2.d0 * hh * lt * lt
      gdtx = 2.d0 * hh * lt * lx
      gdty = 2.d0 * hh * lt * ly
      gdtz = 2.d0 * hh * lt * lz
      gdxx = 1.d0 + 2.d0 * hh * lx * lx
      gdyy = 1.d0 + 2.d0 * hh * ly * ly
      gdzz = 1.d0 + 2.d0 * hh * lz * lz
      gdxy = 2.d0 * hh * lx * ly
      gdyz = 2.d0 * hh * ly * lz
      gdzx = 2.d0 * hh * lz * lx

C     Up metric. g^ab = flat^ab - H l^a l^b.
C     Notice that g^ab = g_ab and l^i = l_i and l^0 = - l_0 in flat spacetime.
      gutt = - 1.d0 - 2.d0 * hh * lt * lt 
      gutx = 2.d0 * hh * lt * lx
      guty = 2.d0 * hh * lt * ly
      gutz = 2.d0 * hh * lt * lz
      guxx = 1.d0 - 2.d0 * hh * lx * lx
      guyy = 1.d0 - 2.d0 * hh * ly * ly
      guzz = 1.d0 - 2.d0 * hh * lz * lz
      guxy = - 2.d0 * hh * lx * ly
      guyz = - 2.d0 * hh * ly * lz
      guzx = - 2.d0 * hh * lz * lx

      return
      end
