C Goedel spacetime !!!!
C See: S. Hawking, G.F.R. Ellis, The Large
C Scale Structure of space-time, Cambridge, 1973
C
C Author : D. Vulcanov (Timsioara, Romania)
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Goedel(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x
      CCTK_DECLARE(CCTK_REAL, y,)
      CCTK_DECLARE(CCTK_REAL, z,)
      CCTK_DECLARE(CCTK_REAL, t,)

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL arad

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

      arad = Goedel__scale
     
      gdtt = -arad*arad
      gdtx = 0.d0
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = 0.5d0*arad*arad*exp(x)*exp(x)
      gdyy = -arad*arad
      gdzz =  arad*arad*exp(x)
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = arad*arad*exp(x)

      gutt = -1d0/(arad*arad)
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = -2d0/(arad*arad*exp(x)*exp(x))
      guyy = -1d0/(arad*arad)
      guzz = -1d0/(arad*arad*exp(x))
      guxy = 0.d0
      guyz = 0.d0 
      guzx = 2d0/(arad*arad*exp(x))

      return
      end
