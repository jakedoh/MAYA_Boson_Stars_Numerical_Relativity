C Axisymmetric Kasner solution !
C
C Author : D. Vulcanov
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Kasner_axisymmetric(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL t
      CCTK_DECLARE(CCTK_REAL, x,)
      CCTK_DECLARE(CCTK_REAL, y,)
      CCTK_DECLARE(CCTK_REAL, z,)

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

      gdtt = -t**(-0.5d0) 
      gdtx = 0.d0 
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = t**(-0.5d0)
      gdyy = t
      gdzz = t
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -t**(0.5d0)
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = t**(0.5d0)
      guyy = 1.d0/t
      guzz = 1.d0/t
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0
 


      return
      end
