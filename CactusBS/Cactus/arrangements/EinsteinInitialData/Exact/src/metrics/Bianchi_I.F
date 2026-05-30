C Bianchi-I spacetime !!!! Fake one...
C It is not a real Bianchi I spacetime, it just emulates
C one, taking the two BianchiI functions as harmonic ones...
C
C Author : D. Vulcanov (Timisoara, Romania)
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Bianchi_I(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x, t
      CCTK_DECLARE(CCTK_REAL, y,)
      CCTK_DECLARE(CCTK_REAL, z,)
     

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL arad
      CCTK_REAL bx, by

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

      arad = Bianchi_I__scale

      bx = arad*sin(x+t)
      by = arad*cos(x+t)
     
      gdtt = - 1.d0
      gdtx = 0.d0
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = bx**2 
      gdyy = by**2 
      gdzz = by**2
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = - 1.d0
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = bx**(-2)
      guyy = by**(-2)
      guzz = by**(-2)
      guxy = 0.d0
      guyz = 0.d0 
      guzx = 0.d0

      return
      end
