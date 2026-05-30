C     fakebinary.F
C     Bernd Bruegmann, 6/98
C
C     Compute Thorne four-metric that, although not a solution to the 
C     Einstein equations, has several characteristic features of a binary
C     star system. See gr-qc/9808024.
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

      subroutine Exact__Thorne_fakebinary(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guxz,
     $     psi, Tmunu_flag)

      implicit none
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

c input arguments
      CCTK_REAL x, y, z, t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdxz,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guxz
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      logical firstcall
      CCTK_REAL eps, m, a0, Omega0, bround, atype, aretarded
      data firstcall /.true./
      save firstcall, eps, m, a0, Omega0, bround, atype, aretarded

C     temps
      CCTK_REAL a, Omega, tau, f
      CCTK_REAL c, c0, c1, c2, c3 
      CCTK_REAL rho, r, sinp, cosp, phi, sint, cost, tx, ty, tz, px, py, pz
      CCTK_REAL a2, b2, bx, by, bz, detgd

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     get parameters of the exact solution.

      if (firstcall) then
         firstcall = .false.

         eps = Thorne_fakebinary__epsilon
         m   = Thorne_fakebinary__mass
         a0  = Thorne_fakebinary__separation

         Omega0 = Thorne_fakebinary__Omega0
         bround = Thorne_fakebinary__smoothing

         bround = max(bround, eps)

         if     (CCTK_Equals(Thorne_fakebinary__atype, "constant").ne.0) then
            atype = 0.d0
         elseif (CCTK_Equals(Thorne_fakebinary__atype,"quadrupole").ne.0) then
            atype = 1.d0
         else
            call CCTK_Warn(0,
     $           "Unknown value of parameter Thorne_fakebinary__atype")
         endif

         if (Thorne_fakebinary__retarded.ne.0) then
            aretarded = 1.d0
         else
            aretarded = 0.d0
         endif
      end if

C     spherical coordinates

      rho = max(sqrt(x**2 + y**2), eps)
      r = sqrt(rho**2 + z**2)
      sinp = y / rho
      cosp = x / rho
      phi = acos(cosp)
      sint = rho / r
      cost = z / r
      tx = cost*cosp
      ty = cost*sinp
      tz = sint
      px = - sinp
      py = cosp
      pz = 0

C     distance function a(T-R)

      tau = 5.d0/128.d0 * a0**4 / m**3
      a = a0 * (1.d0 - atype * 4.d0*(t - aretarded*r)/tau)**(0.25d0)

C     orbital frequency Omega(T-R)

      Omega = 0.5d0*(m/a**3)**2

C     1/r type potential f

      c = y**2 + z**2 + bround**2;
      f = ((x-a)**2 + c)**(-0.5d0) + ((x+a)**2 + c)**(-0.5d0)

C     the three metric, tt part

      c3 = 2.d0*(phi + Omega*r)
      c0 =  - 4.d0 * m * a**2 * Omega**3 * (Omega*r)**4 
     .     / (1 + (Omega*r)**2)**(2.5d0)
      c1 = (1 + cost**2) * cos(c3) * c0
      c2 = - 2.d0 * cost * sin(c3) * c0
      gdxx = c1 * (tx*tx - px*px) +  c2 * (tx*px + px*tx)
      gdxy = c1 * (tx*ty - px*py) +  c2 * (tx*py + px*ty)
      gdxz = c1 * (tx*tz - px*pz) +  c2 * (tx*pz + px*tz)
      gdyy = c1 * (ty*ty - py*py) +  c2 * (ty*py + py*ty)
      gdyz = c1 * (ty*tz - py*pz) +  c2 * (ty*pz + py*tz)
      gdzz = c1 * (tz*tz - pz*pz) +  c2 * (tz*pz + pz*tz)

C     the three metric, add conformally flat part

      c = (1.d0 + m * f)**2
      gdxx = gdxx + c
      gdyy = gdyy + c
      gdzz = gdzz + c

C     the shift vector and covector

      c = (1.d0 -  2*m*a**2/(r**2+a**2) * f) * Omega * rho
      bx = c * px
      by = c * py
      bz = c * pz
      gdtx = gdxx*bx + gdxy*by + gdxz*bz
      gdty = gdxy*bx + gdyy*by + gdyz*bz
      gdtz = gdxz*bx + gdyz*by + gdzz*bz
      b2 = gdtx*bx + gdty*by + gdtz*bz

C     lapse squard and time-time component of the four metric

      a2 = (1.d0 - m * f)**2
      gdtt = b2 - a2

C     done with metric, find its inverse
C     inverse three metric

      detgd = -(gdxz**2*gdyy) + 2*gdxy*gdxz*gdyz - gdxx*gdyz**2
     .        - gdxy**2*gdzz -  gdxx*gdyy*gdzz
      guxx = (-gdyz**2 + gdyy*gdzz) / detgd
      guxy = (gdxz*gdyz - gdxy*gdzz) / detgd
      guxz = (-(gdxz*gdyy) + gdxy*gdyz) / detgd
      guyy = (-gdxz**2 + gdxx*gdzz) / detgd 
      guyz = (gdxy*gdxz - gdxx*gdyz) / detgd 
      guzz = (-gdxy**2 + gdxx*gdyy) / detgd 

C     inverse four metric

      gutt = - 1.d0/a2
      gutx = bx/a2
      guty = by/a2
      gutz = bz/a2
      guxx = guxx - bx*bx/a2
      guxy = guxy - bx*by/a2
      guxz = guxz - bx*bz/a2
      guyy = guyy - by*by/a2
      guyz = guyz - by*bz/a2
      guzz = guzz - bz*bz/a2

C     done!
      return
      end
