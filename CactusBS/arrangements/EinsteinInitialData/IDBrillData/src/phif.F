/*@@
  @file      phif.F
  @date      October 1998
  @author    Miguel Alcubierre
  @desc
             Function to find angle phi from coordinates {x,y}.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h"

      CCTK_REAL function phif(x,y)

      implicit none

      CCTK_REAL x,y
      CCTK_REAL zero,half,one,two,pi

c     Numbers.

      zero = 0.0D0
      half = 0.5D0
      one  = 1.0D0
      two  = 2.0D0

      pi = acos(-one)

c     Find angle between 0 and pi/2 such that tan(phi) = |y/x|.

      if (abs(x).gt.abs(y)) then
         phif = atan(abs(y/x))
      else if (abs(x).lt.abs(y)) then
         phif = half*pi - atan(abs(x/y))
      else
         phif = 0.25D0*pi
      end if

c     Use signs of {x,y} to move to correct quadrant.

      if ((x.eq.zero).and.(y.eq.zero)) then
         phif = zero
      else if ((x.le.zero).and.(y.ge.zero)) then
         phif = pi - phif
      else if ((x.le.zero).and.(y.le.zero)) then
         phif = pi + phif
      else if ((x.ge.zero).and.(y.le.zero)) then
         phif = two*pi - phif
      end if

      return
      end
