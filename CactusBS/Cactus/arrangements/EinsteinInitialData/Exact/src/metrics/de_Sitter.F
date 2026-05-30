C Einstein-DeSitter metric spacetime !!!!
C It emulates the Robertson-Walker universe
C near t=0, with zero pressure, and k=0
C See :J.N. Islam, An Introduction to 
C Mathematical Cosmology, Cambridge, 1992 and 
C S. Hawking, G.F.R. Ellis, The Large Scale
C Structure of space-time, Cambridge, 1973
C
C Author : D. Vulcanov (Timisoara, Romania)
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__de_Sitter(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_DECLARE(CCTK_REAL, x,)
      CCTK_DECLARE(CCTK_REAL, y,)
      CCTK_DECLARE(CCTK_REAL, z,)
      CCTK_REAL t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL arad
      CCTK_REAL am

C This model has a cosmological constant
C ==> it sets the stress-energy tensor in the "CalcTmunu" code
      Tmunu_flag = .true.

      arad = de_Sitter__scale
 
      am=arad*t**(4.0D0/3.0D0)     

      gdtt = -1.0D0 
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = am
      gdyy = am
      gdzz = am
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -1.d0
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = 1.d0/am
      guyy = 1.d0/am
      guzz = 1.d0/am
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0


      return
      end
