c Schwarzschild metric in Eddington-Finkelstein coordinates,
c as per MTW box 31.2
C
C Author: unknown
C Copyright/License: unknown
C
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Schwarzschild_EF(
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

c local static variables
      CCTK_REAL eps, m

c local variables
      CCTK_REAL r

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution.

      eps = Schwarzschild_EF__epsilon
      m   = Schwarzschild_EF__mass

      r = max(sqrt(x**2 + y**2 + z**2), eps)

      gdtt = - (1.d0 - 2.d0 * m / r)
      gdtx = 2.d0 * m * x / r**2
      gdty = 2.d0 * m * y / r**2
      gdtz = 2.d0 * m * z / r**2
      gdxx = 1.d0 + 2.d0 * m * x**2 / r**3
      gdyy = 1.d0 + 2.d0 * m * y**2 / r**3
      gdzz = 1.d0 + 2.d0 * m * z**2 / r**3
      gdxy = 2.d0 * m * x * y / r**3
      gdyz = 2.d0 * m * y * z / r**3
      gdzx = 2.d0 * m * z * x / r**3

      gutt = - (1.d0 + 2.d0 * m / r)
      gutx = 2.d0 * m * x / r**2
      guty = 2.d0 * m * y / r**2
      gutz = 2.d0 * m * z / r**2
      guxx = 1.d0 - 2.d0 * m * x**2 / r**3
      guyy = 1.d0 - 2.d0 * m * y**2 / r**3
      guzz = 1.d0 - 2.d0 * m * z**2 / r**3
      guxy = - 2.d0 * m * x * y / r**3
      guyz = - 2.d0 * m * y * z / r**3
      guzx = - 2.d0 * m * z * x / r**3

      return
      end
