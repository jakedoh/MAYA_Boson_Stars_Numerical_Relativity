C  Schwarzschild-Lemaitre metric
c  (Schwarzschild black hole with cosmological constant)
C  Proposed by Lemaitre in 1932 
C
C Author : D. Vulcanov (Timisoara, Romania)
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Schwarzschild_Lemaitre(
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
      CCTK_REAL lam, mas
      CCTK_REAL r2, ppp, unu, doi

C This model has a cosmological constant
      Tmunu_flag = .true.

      lam = Schwarzschild_Lemaitre__Lambda
      mas = Schwarzschild_Lemaitre__mass
 
      r2 =x*x+y*y+z*z

      ppp=1.0D0 -2.0D0*mas/sqrt(r2) -r2*lam/3.0D0

      unu=(1.0D0-ppp)/ppp/r2      
 
      doi=(ppp - 1.0D0)/r2

      gdtt = -ppp
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = 1.0D0 + x*x*unu
      gdyy = 1.0D0 + y*y*unu
      gdzz = 1.0D0 + z*z*unu
      gdxy =         x*y*unu
      gdyz =         y*z*unu
      gdzx =         z*x*unu


      gutt = -1.0D0/ppp
      gutx = 0.0D0
      guty = 0.0D0
      gutz = 0.0D0
      guxx = 1.0D0 + x*x*doi
      guyy = 1.0D0 + y*y*doi
      guzz = 1.0D0 + z*z*doi
      guxy =         x*y*doi
      guyz =         y*z*doi
      guzx =         x*z*doi


      return
      end
