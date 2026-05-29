c     The metric given here corresponds to that of flat spacetime
c     but with non-trivial spatial coordinates.  Basically, I take
c     the flat metric in spherical coordinates, define a new
c     radial coordinate such that:
c
c     r  =  r    ( 1  -  a f(r   ) )
c            new              new
c
c     where f(r) is a gaussian centered at r=0 with amplitude 1.
c     Finally, I transform back to cartesian coordinates.
c     For  0 <= a < 1, the transformation above is monotonic.
C
C Author: unknown
C Copyright/License: unknown
C
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Minkowski_funny(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_DECLARE(CCTK_REAL, x,)
      CCTK_DECLARE(CCTK_REAL, y,)
      CCTK_DECLARE(CCTK_REAL, z,)
      CCTK_DECLARE(CCTK_REAL, t,)

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $          gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $          gutt, gutx, guty, gutz, 
     $          guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL a,s
      CCTK_REAL r,r2,x2,y2,z2
      CCTK_REAL f,fp,g11,g22
      CCTK_REAL det

c constants
      CCTK_REAL zero, one, two
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)

c This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

c     Read parameters

      a = Minkowski_funny__amplitude
      s = Minkowski_funny__sigma

c     Find transformation function.

      x2 = x**2
      y2 = y**2
      z2 = z**2

      r2 = x2 + y2 + z2
      r  = sqrt(r2)

      f  = exp(-r2/s**2)
      fp = - two*r/s**2*f

c     Give metric components.

      g11 = (one - a*(f + r*fp))**2
      g22 = (one - a*f)**2

      gdtt = - one
      gdtx = zero
      gdty = zero
      gdtz = zero

      gdxx = (x2*g11 + (y2 + z2)*g22)/r2
      gdyy = (y2*g11 + (x2 + z2)*g22)/r2
      gdzz = (z2*g11 + (x2 + y2)*g22)/r2

      gdxy = x*y*(g11 - g22)/r2
      gdzx = x*z*(g11 - g22)/r2
      gdyz = y*z*(g11 - g22)/r2

c     Find inverse metric.

      gutt = - one
      gutx = zero
      guty = zero
      gutz = zero

      det = gdxx*gdyy*gdzz + two*gdxy*gdzx*gdyz
     .    - gdxx*gdyz**2 - gdyy*gdzx**2 - gdzz*gdxy**2

      guxx = (gdyy*gdzz - gdyz**2)/det
      guyy = (gdxx*gdzz - gdzx**2)/det
      guzz = (gdxx*gdyy - gdxy**2)/det

      guxy = (gdzx*gdyz - gdzz*gdxy)/det
      guyz = (gdxy*gdzx - gdxx*gdyz)/det
      guzx = (gdxy*gdyz - gdyy*gdzx)/det

      return
      end
