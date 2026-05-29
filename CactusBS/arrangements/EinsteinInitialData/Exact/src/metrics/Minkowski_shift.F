c     The metric given here corresponds to that of flat spacetime
c     but with non-trivial slicing and a shift vector such that
c     the resulting metric is still time independent.  I take
c     the flat metric in spherical coordinates and define a new
c     time coordinate as:
c
c     t    =   t   -  a f(r)
c      new
c
c     where f(r) is a gaussian centered at r=0 with amplitude 1.
c     Finally, I transform back to cartesian coordinates.
c     For  -1 < fp < 1, the transformation above results in spatial
c     slices.
C
C Author: unknown
C Copyright/License: unknown
C
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Minkowski_shift(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x,y,z
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
      CCTK_REAL f,fp,fpr,fpr2

c constants
      CCTK_REAL zero,one,two
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)

c This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

c     Read parameters

      a = Minkowski_shift__amplitude
      s = Minkowski_shift__sigma

c     Find transformation function.

      x2 = x**2
      y2 = y**2
      z2 = z**2

      r2 = x2 + y2 + z2
      r  = sqrt(r2)

      f    = a*exp(-r2/s**2)
      fp   = - two*f*r/s**2
      fpr  = fp/r
      fpr2 = fpr**2

c     Give metric components.

      gdtt = - one

      gdtx = - x*fpr
      gdty = - y*fpr
      gdtz = - z*fpr

      gdxx = one - x2*fpr2
      gdyy = one - y2*fpr2
      gdzz = one - z2*fpr2

      gdxy = - x*y*fpr2
      gdzx = - x*z*fpr2
      gdyz = - y*z*fpr2

c     Inverse metric. And yes, it is this simple, simpler
c     than the metric itself.  I tripled checked!

      gutt = - one + fp**2

      gutx = - x*fpr
      guty = - y*fpr
      gutz = - z*fpr

      guxx = one
      guyy = one
      guzz = one

      guxy = zero
      guzx = zero
      guyz = zero

      return
      end
