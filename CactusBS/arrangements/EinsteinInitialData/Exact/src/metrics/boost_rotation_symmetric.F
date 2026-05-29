c Boost-Rotation symmetric metric
c see Pravda and Pravdova [a-acute accent on last a], gr-qc/0003067
C
C Author: unknown
C Copyright/License: unknown
C
c $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__boost_rotation_symmetric(
     $     x, y, z, t,
     $     gdtt, gdtx, gdty, gdtz, 
     $     gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $     gutt, gutx, guty, gutz, 
     $     guxx, guyy, guzz, guxy, guyz, guzx,
     $     psi, Tmunu_flag)

      implicit none

      DECLARE_CCTK_PARAMETERS

c input arguments
      CCTK_REAL x, y, z, t

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz,
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz,
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c functions local to this file
      CCTK_REAL gfunc

c local variables
      CCTK_REAL h, d, numlim
      CCTK_REAL a, b, mu0, mu1, lam1, mu2, lam2,
     $    lam3, mu4, lam4, mu5, lam5, num, div, f,
     $    elam, emu0, delta, tmp

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution.

      h      = boost_rotation_symmetric__scale
      d      = boost_rotation_symmetric__amp
      numlim = boost_rotation_symmetric__min_d

C     Intermediate quantities.

      a = x**2 + y**2
      b = z**2 - t**2

      num = (0.5d0 * (a + b) - h)**2 + 2.d0 * h * a

C     Make sure we are not sitting on one of the two source wordlines,
C     given by x = y = 0, z = +/- sqrt(h^2 + t^2)

      if (num / h**4 .le. numlim) then
      	 call CCTK_WARN (0, "too close to source wordline")
      end if

      div = 1.d0 / sqrt(num**3)
      f = d**2  * ((0.25d0 * (a + b)**2 - h**2)**2 
     $     - 0.5d0 * h**2 * a * b) / num**4

      mu0 = - d * div * (0.5d0 * a**2 + h * a)
      mu1 = - d * div * (0.5d0 * b + a - h)
      lam1 = d * div * (0.5d0 * b - h) - a * f
      mu2 = gfunc(b, mu1)
      lam2 =  gfunc(b, lam1)

      lam3 = d * div * (0.5d0 * b**2 - h * b)
      lam4 = - d * div * (0.5d0 * a + h) - b * f
      mu4 = - d * div * (0.5d0 * a + b + h) 
      mu5 = gfunc(a, - mu4)
      lam5 = gfunc(a, lam4)

      elam = exp(lam3 + a * lam4)
      emu0 = exp(mu0)
      delta = exp(lam3) * (mu5 - lam5)

C     All nonvanishing metric coefficients (downstairs).

      gdxx = elam + y**2 * Delta
      gdyy = elam + x**2 * Delta
      gdxy = - x * y * Delta
      gdzz = emu0 * (1.d0 + lam2 * z**2 - mu2 * t**2)
      gdtz = - emu0 * z * t * (lam2 - mu2)
      gdtt = - emu0 * (1.d0 + mu2 * z**2 - lam2 * t**2)

C     Others.

      gdzx = 0.d0
      gdyz = 0.d0
      gdtx = 0.d0
      gdty = 0.d0

C     It is clear that the 3-metric is always spacelike in the xy plane. So
C     we only need to check that gdzz is positive.

      if (gdzz .le. 0.d0) then
         write(*,*) 'WARNING 3-metric not spacelike in boostrot at'
         write(*,*) 't =', t, 'z =', z
         write(*,*) 'x =', x, 'y =', y
	 call CCTK_WARN (0, "aborting")
      end if

C     Calculate inverse metric. That is not too difficult as it is
c     in block-diagonal form.

      tmp = gdtt * gdzz - gdtz**2

      if (tmp .eq. 0.d0) then
	 call CCTK_WARN (0, "boostrot metric inverse failed in tz plane")
      end if

      gutt = gdzz / tmp
      guzz = gdtt / tmp
      gutz = - gdtz / tmp
      
      tmp = gdxx * gdyy - gdxy**2

      if (tmp .eq. 0.d0) then
         call CCTK_WARN (0, "boostrot metric inverse failed in xy plane")
      end if

      guxx = gdyy / tmp
      guyy = gdxx / tmp
      guxy = - gdxy / tmp

      guzx = 0.d0
      guyz = 0.d0
      gutx = 0.d0
      guty = 0.d0

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Calculates g = [exp (x f) - 1] / x as a power series for small x, 
C     so that the expression is regular at x = 0.

      CCTK_REAL function gfunc(x, f)
    
      implicit none

      integer n

      CCTK_REAL x, f
      CCTK_REAL sum, tmp

      if (abs(x*f) .ge. 1.d-6) then
         gfunc = (exp(x*f) - 1.d0) / x
      else
         sum = 0.d0
         tmp = f
         do n=1,10
            tmp = tmp / dble(n) 
            sum = sum + tmp
            tmp = tmp * x * f
         end do
         gfunc = sum
      end if

      return
      end
