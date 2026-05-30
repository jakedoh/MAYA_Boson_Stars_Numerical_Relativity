C Robertson-Walker universe 
C See: J.N. Islam, An Introduction to 
C Mathematical Cosmology, Cambridge, 1992
C
C Author : D. Vulcanov (Timisoara, Romania)
C see ../../README for copyright & licensing info
C
C $Header$
C
C FIXME:
C This metric doesn't work.  The argument  rama  is the R(t) in the
C Robertson-Walker metric, and if this is passed correctly then this
C subroutine computes the correct metric.  But the rest of this thorn
C doesn't know to pass this value. :( :(  See Mitica's Cosmo thorn
C for a better way to get this metric.

#include "cctk.h"
#include "cctk_Parameters.h"


      subroutine Exact__Robertson_Walker(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag,
     $     rama)

      implicit none

      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS      

c input arguments
      CCTK_REAL x, y, z, t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_REAL psi
      LOGICAL   Tmunu_flag
      CCTK_REAL rama

c local static variables
      logical firstcall
      CCTK_REAL kapa
      data firstcall /.true./
      save firstcall, kapa

c local variables
      CCTK_REAL r2,x2,y2,z2,am,ag

C This model may set the stress-energy tensor
      Tmunu_flag = .true.

      if (firstcall) then
          kapa = Robertson_Walker__k
          firstcall = .false.
      end if

      x2=x*x
      y2=y*y
      z2=z*z
      r2 = x2+y2+z2

      ag = kapa/(1.0D0-kapa*r2)    
       
      am = rama*rama
       


      gdtt = -1.0D0 
      gdtx = 0.0D0 
      gdty = 0.0D0
      gdtz = 0.0D0
      gdxx = am*(1.0D0 + ag*x2) 
      gdyy = am*(1.0D0 + ag*y2)
      gdzz = am*(1.0D0 + ag*z2)
      gdxy = am*ag*x*y
      gdyz = am*ag*y*z
      gdzx = am*ag*z*x

      gutt = -1.0D0
      gutx = 0.0D0
      guty = 0.0D0
      gutz = 0.0D0
      guxx = (1.0D0-kapa*x2)/am 
      guyy = (1.0D0-kapa*y2)/am
      guzz = (1.0D0-kapa*z2)/am
      guxy = -kapa*x*y/am
      guyz = -kapa*y*z/am
      guzx = -kapa*z*x/am
 


      return
      end
