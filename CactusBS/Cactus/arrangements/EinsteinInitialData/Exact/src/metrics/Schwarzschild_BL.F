c Schwarzschild spacetime in Brill-Lindquist coordinates.
C
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Schwarzschild_BL(
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

c local variables
      CCTK_REAL r, psi4

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution.

      eps = Schwarzschild_BL__epsilon
      m   = Schwarzschild_BL__mass

      r = ((x**2 + y**2 + z**2)**2 + eps**4) ** 0.25d0
      psi4 = (1 + m / (2 * r)) ** 4

      gdtt = -1
      gdtx = 0
      gdty = 0
      gdtz = 0
      gdxx = psi4
      gdyy = psi4
      gdzz = psi4
      gdxy = 0
      gdyz = 0
      gdzx = 0

      gutt = -1
      gutx = 0
      guty = 0
      gutz = 0
      guxx = 1 / psi4
      guyy = 1 / psi4
      guzz = 1 / psi4
      guxy = 0
      guyz = 0
      guzx = 0

      return
      end
