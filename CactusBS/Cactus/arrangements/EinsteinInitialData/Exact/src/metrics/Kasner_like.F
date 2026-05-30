C Kasner-like metric spacetime !!!!
C See: L. Pimentel, Int. Journ. of Theor. Physics,
C, vol. 32, No. 6, 1993, p. 979 (and the references
C cited here), and S. Gotlober, et. al. Early Evolution
C of the Universe and Formation Structure, Akad. Verlag, 1990
C
C Author : D. Vulcanov
C see ../../README for copyright & licensing info
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Kasner_like(
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

c local variables
      CCTK_REAL qq,a1,a3
 
c this model sets the stress-energy tensor in the "CalcTmunu" code
      Tmunu_flag = .true.

      qq = Kasner_like__q

      a1= t**(2.d0*qq)
      a3= t**(2.d0-4.d0*qq)    


      gdtt = -1.d0 
      gdtx = 0.d0 
      gdty = 0.d0
      gdtz = 0.d0
      gdxx = a1 
      gdyy = a1
      gdzz = a3
      gdxy = 0.d0
      gdyz = 0.d0
      gdzx = 0.d0

      gutt = -1.d0
      gutx = 0.d0
      guty = 0.d0
      gutz = 0.d0
      guxx = 1.d0/a1
      guyy = 1.d0/a1
      guzz = 1.d0/a3
      guxy = 0.d0
      guyz = 0.d0
      guzx = 0.d0
 

      return
      end
