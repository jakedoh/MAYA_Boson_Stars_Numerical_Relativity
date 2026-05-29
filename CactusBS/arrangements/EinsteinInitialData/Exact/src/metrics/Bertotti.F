C Bertotti spacetime with cosmological constant
C
C Author : D. Vulcanov (Timisoara, Romania)
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Bertotti(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x, z
      CCTK_DECLARE(CCTK_REAL, y,)
      CCTK_DECLARE(CCTK_REAL, t,)

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL baza
      CCTK_REAL unu, doi
 
C this model has a cosmological constant
C ==> it sets the stress-energy tensor in the "CalcTmunu" code
      Tmunu_flag = .true.

      baza = Bertotti__Lambda

      unu=exp(2.0D0*sqrt(-baza)*x)
      doi=exp(2.0D0*sqrt(-baza)*z)   

       
      gdtt = -unu
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = 1.0D0
      gdyy = doi
      gdzz = 1.0D0
      gdxy = 0.0D0
      gdyz = 0.0D0
      gdzx = 0.0D0

      gutt = -1.0D0/unu
      gutx = 0.0D0
      guty = 0.0D0
      gutz = 0.0D0
      guxx = 1.0D0
      guyy = 1.0D0/doi
      guzz = 1.0D0
      guxy = 0.0D0
      guyz = 0.0D0
      guzx = 0.0D0


      return
      end
