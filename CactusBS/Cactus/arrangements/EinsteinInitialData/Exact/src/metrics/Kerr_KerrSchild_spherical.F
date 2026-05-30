C     Kerr-Schild form of boosted rotating black hole.
C     Program g_ab = eta_ab + H l_a l_b, g^ab = eta^ab - H l^a l^b.
C     Here eta_ab is Minkowski in Cartesian coordinates, H is a scalar,
C     and l is a null vector.
C     
C     The coordinates are distorted, such that the event horizon is
C     a coordinate sphere.
C
C     Author: Erik Schnetter <schnetter@cct.lsu.edu>
C     This formulation was invented by Nils Dorband <dorband@cct.lsu.edu>
C
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

      subroutine Exact__Kerr_KerrSchild_spherical (
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

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
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local variables
      CCTK_REAL boostv, eps, m, a

c local variables
      CCTK_REAL t1, x1, y1, z1, rho1,
     $     gamma, t0, z0, x0, y0, r0, costheta0, 
     $     lt0, lx0, ly0, lz0, hh, lt, lx, ly, lz

      CCTK_REAL gd(3,3), gdt(3,3), det, jac(3,3)

      integer i, j, k, l

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution,
C     and convert from parameter file spin parameter J/m^2
C     to the J/m definition used in the code here.

      boostv = Kerr_KerrSchild__boost_v
      eps    = Kerr_KerrSchild__epsilon
      m      = Kerr_KerrSchild__mass
      a      = m*Kerr_KerrSchild__spin

C     Distort the coordinates such that the event horizon is a
C     coordinate sphere
      rho1 = sqrt (x**2 + y**2 + z**2)
      rho1 = (rho1**4 + eps**4) ** 0.25d0
      t1 = t
      x1 = x - a * y / rho1
      y1 = y + a * x / rho1
      z1 = z

C     Boost factor.

      gamma = 1 / sqrt(1 - boostv**2)

C     Lorentz transform t,x,y,z -> t0,x0,y0,z0. 
C     t0 is never used, but is here for illustration, and we introduce
C     x0 and y0 also only for clarity.
C     Note that z0 = 0 means z = vt for the BH.

      t0 = gamma * ((t1 - Kerr_KerrSchild__t) - boostv * (z1 - Kerr_KerrSchild__z))
      z0 = gamma * ((z1 - Kerr_KerrSchild__z) - boostv * (t1 - Kerr_KerrSchild__t))
      x0 = x1 - Kerr_KerrSchild__x
      y0 = y1 - Kerr_KerrSchild__y

C     Spherical auxiliary coordinate r and angle theta in BH rest frame.

      r0 = rho1
      costheta0 = z0 / r0
      
C     Coefficient H.  Note this transforms as a scalar, so does not carry
C     the suffix 0.
      hh = m * r0 / (r0**2 + a**2 * costheta0**2)
      
C     Components of l_a in rest frame. Note indices down.
      lt0 = 1.d0
      lx0 = (r0 * x0 + a * y0) / (r0**2 + a**2)
      ly0 = (r0 * y0 - a * x0) / (r0**2 + a**2)
      lz0 = z0 / r0

C     Now boost it to coordinates x, y, z, t.
C     This is the reverse Lorentz transformation, but applied
C     to a one-form, so the sign of boostv is the same as the forward
C     Lorentz transformation applied to the coordinates.

      lt = gamma * (lt0 - boostv * lz0) 
      lz = gamma * (lz0 - boostv * lt0) 
      lx = lx0
      ly = ly0

C     Down metric.  g_ab = flat_ab + H l_a l_b

      gdtt = - 1.d0 + 2.d0 * hh * lt * lt
      gdtx = 2.d0 * hh * lt * lx
      gdty = 2.d0 * hh * lt * ly
      gdtz = 2.d0 * hh * lt * lz

      gd(1,1) = 1.d0 + 2.d0 * hh * lx * lx
      gd(2,2) = 1.d0 + 2.d0 * hh * ly * ly
      gd(3,3) = 1.d0 + 2.d0 * hh * lz * lz
      gd(1,2) = 2.d0 * hh * lx * ly
      gd(2,3) = 2.d0 * hh * ly * lz
      gd(3,1) = 2.d0 * hh * lz * lx
      gd(2,1) = gd(1,2)
      gd(3,2) = gd(2,3)
      gd(1,3) = gd(3,1)

C     Transform the tensor basis back
      jac(1,1) = 1 + (a*x*y) / (rho1**3)
      jac(1,2) = - a * (x**2+z**2) / (rho1**3)
      jac(1,3) = a*y*z / (rho1**3)
      
      jac(2,1) = a * (y**2+z**2) / (rho1**3)
      jac(2,2) = 1 - a*x*y / (rho1**3)
      jac(2,3) = - a*x*z / (rho1**3)
      
      jac(3,1) = 0
      jac(3,2) = 0
      jac(3,3) = 1

      do i = 1, 3
         do j = 1, 3
            gdt(i,j) = 0
            do k = 1, 3
               do l = 1, 3
                  gdt(i,j) = gdt(i,j) + gd(k,l) * jac(k,i) * jac(l,j)
               end do
            end do
         end do
      end do

      gdxx = gdt(1,1)
      gdyy = gdt(2,2)
      gdzz = gdt(3,3)
      gdxy = gdt(1,2)
      gdyz = gdt(2,3)
      gdzx = gdt(3,1)

C     Up metric.  g^ab = flat^ab - H l^a l^b.
C     Notice that g^ab = g_ab and l^i = l_i and l^0 = - l_0 in flat spacetime.
      gutt = - 1.d0 - 2.d0 * hh * lt * lt 
      gutx = 2.d0 * hh * lt * lx
      guty = 2.d0 * hh * lt * ly
      gutz = 2.d0 * hh * lt * lz

      det = gdxx*gdyy*gdzz + 2*gdxy*gdzx*gdyz
     .    - gdxx*gdyz**2 - gdyy*gdzx**2 - gdzz*gdxy**2

      guxx = (gdyy*gdzz - gdyz**2)/det
      guyy = (gdxx*gdzz - gdzx**2)/det
      guzz = (gdxx*gdyy - gdxy**2)/det

      guxy = (gdzx*gdyz - gdzz*gdxy)/det
      guyz = (gdxy*gdzx - gdxx*gdyz)/det
      guzx = (gdxy*gdyz - gdyy*gdzx)/det

      end
