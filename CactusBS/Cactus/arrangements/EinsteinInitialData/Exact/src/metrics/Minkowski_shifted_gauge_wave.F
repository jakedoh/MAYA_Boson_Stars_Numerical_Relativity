C This subroutine sets up Minkowski spacetime with a gague wave.
C
c $Header$
C
C Author: unknown
C Copyright/License: unknown
C

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#define Pi (4 * atan(1.d0))


      subroutine Exact__Minkowski_shifted_gauge_wave(
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
      CCTK_REAL x, y, t
      CCTK_DECLARE(CCTK_REAL, z,)

c output arguments
      CCTK_REAL gdtt, gdtx, gdty, gdtz, 
     $       gdxx, gdyy, gdzz, gdxy, gdyz, gdzx,
     $       gutt, gutx, guty, gutz, 
     $       guxx, guyy, guzz, guxy, guyz, guzx
      CCTK_DECLARE(CCTK_REAL, psi,)
      LOGICAL   Tmunu_flag

c local parameter copies
      CCTK_REAL a, o, d, fs

c local variables
      CCTK_REAL H
      character*100 warn_buffer

c constants
      CCTK_REAL zero,half,one
      parameter (zero = 0.0d0, half=0.5d0, one=1.0d0)

C This is a vacuum spacetime with no cosmological constant
      Tmunu_flag = .false.

C     Get parameters of the exact solution.
      a  = Minkowski_gauge_wave__amplitude
      o  = Minkowski_gauge_wave__omega
      d  = Minkowski_gauge_wave__lambda
      fs = Minkowski_gauge_wave__phase

C     How should the wave look like.
      if (CCTK_EQUALS(Minkowski_gauge_wave__what_fn,"sin")) then
         d = Minkowski_gauge_wave__lambda * half / Pi
         if (Minkowski_gauge_wave__diagonal.ne.0) then
            H = a * sin((x-y)/d - o*t/d - fs)
         else
            H = a * sin((x-o*t)/d - fs)
         end if
      else
        write (warn_buffer, '(a,a,a)')
     $        'Unknown Minkowski_gauge_wave__what_fn = "',
     $         Minkowski_gauge_wave__what_fn, '"'
C       silence compiler warning about unused variable
        H = one
        call CCTK_WARN(0, warn_buffer)
      end if

C     write metric.

      if (Minkowski_gauge_wave__diagonal.ne.0) then

        gdxx = half * H + one
        gdxy = - half * H
        gdyy = half * H + one

        guxx = - half * H + one
        guxy = half * H
        guyy = - half * H + one

      else

        gdxx = H + one
        gdxy = zero
        gdyy = one

        guxx = - H + one
        guxy = zero
        guyy = one

      end if
      
      gdtt = H - one
      gdtx = - H
      gdty = zero
      gdtz = zero

      gdzx = zero
      gdyz = zero
      gdzz = one

C     and upper metric.

      gutt = - H - one
      gutx = - H
      guty = zero
      gutz = zero

      guyz = zero
      guzx = zero
      guzz = one

      end
