c     The metric given here corresponds to a constant
c     density star, also known as a "Schwarzschild" star.
c     There is corresponding code in
c        include/Scalar_CalcTmunu.inc
c     to set up the matter variables.
C
c     Author: Mitica Vulcanov 
C     see ../../README for copyright & licensing info
C
c $Header$

c
c     The metric is given as a conformally flat metric.
c     Turns out that in the original areal radius, the
c     metric variables have a kink at the surface of the
c     star, but they are smooth in the conformal form.
c
c     Thanks to Philippos Papadopoulos for suggesting
c     the use of this metric.

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__constant_density_star(
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
      CCTK_REAL mass,radius
      CCTK_REAL r,c,psi4

c constants
      CCTK_REAL zero,one,two
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)

c This model has matter
c ==> it sets the stress-energy tensor in the "CalcTmunu" code
      Tmunu_flag = .true.

c     Get parameters of the metric.

      mass   = constant_density_star__mass
      radius = constant_density_star__radius

c     Find r.

      r = sqrt(x**2 + y**2 + z**2)

c     Find conformal factor.

      if (r.le.radius) then

         c = mass/(two*radius)

         psi4 = (one + c)**6/(one + c*(r/radius)**2)**2

      else

         c = mass/(two*r)

         psi4 = (one + c)**4

      end if

c     Find metric components.

      gdtt = - psi4

      gdtx = zero
      gdty = zero
      gdtz = zero

      gdxx = psi4
      gdyy = psi4
      gdzz = psi4

      gdxy = zero
      gdyz = zero
      gdzx = zero

c     Find inverse metric.

      gutt = -one/psi4

      gutx = zero
      guty = zero
      gutz = zero

      guxx = one/psi4
      guyy = one/psi4
      guzz = one/psi4

      guxy = zero
      guyz = zero
      guzx = zero

      return
      end
