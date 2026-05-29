c Schwarzschild spacetime in Painleve [e-acute on last e]-Gullstrand
c coordinates.  These have a *flat* 3-metric, and are described in
c detail in Martel and Poisson, gr-qc/0001069
C
C Author: unknown
C Copyright/License: unknown
C
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Schwarzschild_PG(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x, y, z
      CCTK_DECLARE(CCTK_REAL, t,)

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL eps, m
      CCTK_REAL r, bx, by, bz, b2

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution.

      eps = Schwarzschild_PG__epsilon
      m   = Schwarzschild_PG__mass

      r = max(sqrt(x**2 + y**2 + z**2), eps)
      bx = sqrt(2.d0 * m / r) * x / r
      by = sqrt(2.d0 * m / r) * y / r
      bz = sqrt(2.d0 * m / r) * z / r
      b2 = 2.d0 * m / r

      gdtt = - 1.d0 + b2
      gdtx = bx
      gdty = by
      gdtz = bz
      gdxx = 1.d0 
      gdyy = 1.d0 
      gdzz = 1.d0 
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = - 1.d0
      gutx = bx
      guty = by
      gutz = bz
      guxx = 1.d0 - bx**2
      guyy = 1.d0 - by**2
      guzz = 1.d0 - bz**2
      guxy = - bx * by
      guyz = - by * bz
      guzx = - bz * bx

      return
      end
