C Milne spacetime metric ?!?!?
C Suggested by Matteo Rossi and E. Onofri (Univ. di Parma, Italy)
C They inted to use thsi metric for simulating an Pre-Big-Bang
C Cosmology, as proposed by Veneziano some year ago 
C
C Author : D. Vulcanov (Timsoara, Romania)
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"

      subroutine Exact__Milne(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

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
      CCTK_REAL coef, x2,y2,z2,t2

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

      x2= x*x
      y2= y*y
      z2= z*z
      t2= t*t
      coef= t2/(1.d0+ x2 +y2+ z2)

      gdtt = -1.d0 
      gdtx = 0.d0  
      gdty = 0.d0 
      gdtz = 0.d0 
      gdxx =  coef*(1.d0+y2+z2) 
      gdyy =  coef*(1.d0+x2+z2) 
      gdzz =  coef*(1.d0+x2+y2)
      gdxy = -coef*x*y 
      gdyz = -coef*y*z 
      gdzx = -coef*x*z

      gutt = -1.d0 
      gutx = 0.d0 
      guty = 0.d0 
      gutz = 0.d0 
      guxx = (1.d0+x2)/(t2)
      guyy = (1.d0+y2)/(t2)
      guzz = (1.d0+z2)/(t2)
      guxy = x*y/(t2)
      guyz = y*z/(t2) 
      guzx = x*z/(t2)


      return
      end
