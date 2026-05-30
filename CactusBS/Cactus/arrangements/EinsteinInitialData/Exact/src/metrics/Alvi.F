C The Alvi metric. Full description of this metric is given
C in gr-qc/9912113 . 
C Author: Nina Jansen (jansen@aei-potsdam.mpg.de)
C $Header$

#include "cctk.h"
#include "cctk_Parameters.h"

      subroutine Exact__Alvi(
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
      CCTK_REAL psi
      LOGICAL   Tmunu_flag

c locals
      CCTK_REAL m1,m2,b
      CCTK_REAL rin1,rin2,rout,x1,x2, r1, r2, r3, results(10)

C this is a vacuum solution with no cosmological constant
C ==> it does not set the stress-energy tensor
      Tmunu_flag = .false.

      m1 = Alvi__mass1
      m2 = Alvi__mass2
      b  = Alvi__separation

      rin1 = dsqrt(m1*b)
      rin2 = dsqrt(m2*b)
      rout = b*dsqrt(b)/(2.0d0*dsqrt(m1+m2))

      x1 = b

      x2 = -b

      r1 = dsqrt((x-x1)**2.0D0 + y**2.0D0 + z**2.0D0)


      r2 = dsqrt((x-x2)**2.0D0 + y**2.0D0 + z**2.0D0)


      r3 = dsqrt(x**2.0D0 + y**2.0D0 + z**2.0D0)

      if (r1 .le. rin1) then

      call Alvireg1(x,y,z,m1,m2,b,results)
      
      gdtt = results(1)
      gdtx = results(2)
      gdty = results(3)
      gdtz = results(4)
      gdxx = results(5)
      gdyy = results(8)
      gdzz = results(10)
      gdxy = results(9)
      gdyz = results(9)
      gdzx = results(7)
      
      else if (r2 .le. rin2) then

         call Alvireg2(x,y,z,m1,m2,b,results)

         gdtt = results(1)
         gdtx = results(2)
         gdty = results(3)
         gdtz = results(4)
         gdxx = results(5)
         gdyy = results(8)
         gdzz = results(10)
         gdxy = results(9)
         gdyz = results(9)
         gdzx = results(7)
         
      else if ((r3 .le. rout) .and. (r2 .gt. rin2) .and. (r1 .gt. rin1)) then
         call Alvireg3(x,y,z,m1,m2,b,results)

         gdtt = results(1)
         gdtx = results(2)
         gdty = results(3)
         gdtz = results(4)
         gdxx = results(5)
         gdyy = results(8)
         gdzz = results(10)
         gdxy = results(9)
         gdyz = results(9)
         gdzx = results(7)
         
      else if (r3 .gt. rout) then
         call Alvireg4(x,y,z,m1,m2,b,results)

         gdtt = results(1)
         gdtx = results(2)
         gdty = results(3)
         gdtz = results(4)
         gdxx = results(5)
         gdyy = results(8)
         gdzz = results(10)
         gdxy = results(9)
         gdyz = results(9)
         gdzx = results(7)

      else 

         print *,'problem!'

      end if         
      
      gutt = 
     $ (gdzx**2*gdyy - 2*gdxy*gdzx*gdyz + gdxy**2*gdzz + gdxx*(gdyz**2
     $ - gdyy*gdzz))/ (gdtt*gdzx**2*gdyy + gdtz**2*(-gdxy**2 +
     $ gdxx*gdyy) - 2*gdtt*gdxy*gdzx*gdyz - gdtx**2*gdyz**2 +
     $ gdtt*gdxx*gdyz**2 + 2*gdtz*(gdty*gdxy*gdzx - gdtx*gdzx*gdyy -
     $ gdty*gdxx*gdyz + gdtx*gdxy*gdyz) + gdtt*gdxy**2*gdzz +
     $ gdtx**2*gdyy*gdzz - gdtt*gdxx*gdyy*gdzz + gdty**2*(-gdzx**2 +
     $ gdxx*gdzz) + 2*gdtx*gdty*(gdzx*gdyz - gdxy*gdzz))

      gutx = 
     $ (-(gdtz*gdzx*gdyy) + gdtz*gdxy*gdyz + gdty*gdzx*gdyz -
     $ gdtx*gdyz**2 - gdty*gdxy*gdzz + gdtx*gdyy*gdzz)/
     $ (gdtt*gdzx**2*gdyy + gdtz**2*(-gdxy**2 + gdxx*gdyy) -
     $ 2*gdtt*gdxy*gdzx*gdyz - gdtx**2*gdyz**2 + gdtt*gdxx*gdyz**2 +
     $ 2*gdtz*(gdty*gdxy*gdzx - gdtx*gdzx*gdyy - gdty*gdxx*gdyz +
     $ gdtx*gdxy*gdyz) + gdtt*gdxy**2*gdzz + gdtx**2*gdyy*gdzz -
     $ gdtt*gdxx*gdyy*gdzz + gdty**2*(-gdzx**2 + gdxx*gdzz) +
     $ 2*gdtx*gdty*(gdzx*gdyz - gdxy*gdzz))


      guty = 
     $ (-(gdtz*gdxy*gdzx) + gdty*gdzx**2 + gdtz*gdxx*gdyz -
     $ gdtx*gdzx*gdyz - gdty*gdxx*gdzz + gdtx*gdxy*gdzz)/
     $ (-(gdtt*gdzx**2*gdyy) + gdtz**2*(gdxy**2 - gdxx*gdyy) +
     $ 2*gdtt*gdxy*gdzx*gdyz + gdtx**2*gdyz**2 - gdtt*gdxx*gdyz**2 -
     $ 2*gdtz*(gdty*gdxy*gdzx - gdtx*gdzx*gdyy - gdty*gdxx*gdyz +
     $ gdtx*gdxy*gdyz) - gdtt*gdxy**2*gdzz - gdtx**2*gdyy*gdzz +
     $ gdtt*gdxx*gdyy*gdzz + gdty**2*(gdzx**2 - gdxx*gdzz) +
     $ gdty*(-2*gdtx*gdzx*gdyz + 2*gdtx*gdxy*gdzz))





      gutz = 
     $ (-(gdty*gdxy*gdzx) + gdtx*gdzx*gdyy + gdtz*(gdxy**2 - gdxx*gdyy)
     $ + gdty*gdxx*gdyz - gdtx*gdxy*gdyz)/ (-(gdtt*gdzx**2*gdyy) +
     $ gdtz**2*(gdxy**2 - gdxx*gdyy) + 2*gdtt*gdxy*gdzx*gdyz +
     $ gdtx**2*gdyz**2 - gdtt*gdxx*gdyz**2 - 2*gdtz*(gdty*gdxy*gdzx -
     $ gdtx*gdzx*gdyy - gdty*gdxx*gdyz + gdtx*gdxy*gdyz) -
     $ gdtt*gdxy**2*gdzz - gdtx**2*gdyy*gdzz + gdtt*gdxx*gdyy*gdzz +
     $ gdty**2*(gdzx**2 - gdxx*gdzz) + gdty*(-2*gdtx*gdzx*gdyz +
     $ 2*gdtx*gdxy*gdzz)) 
      guxx =
     $ (gdtz**2*gdyy - 2*gdty*gdtz*gdyz + gdty**2*gdzz + gdtt*(gdyz**2
     $ - gdyy*gdzz))/ (gdtt*gdzx**2*gdyy + gdtz**2*(-gdxy**2 +
     $ gdxx*gdyy) - 2*gdtt*gdxy*gdzx*gdyz - gdtx**2*gdyz**2 +
     $ gdtt*gdxx*gdyz**2 + 2*gdtz*(gdty*gdxy*gdzx - gdtx*gdzx*gdyy -
     $ gdty*gdxx*gdyz + gdtx*gdxy*gdyz) + gdtt*gdxy**2*gdzz +
     $ gdtx**2*gdyy*gdzz - gdtt*gdxx*gdyy*gdzz + gdty**2*(-gdzx**2 +
     $ gdxx*gdzz) + 2*gdtx*gdty*(gdzx*gdyz - gdxy*gdzz)) 
      guyy = 
     $ (gdtz**2*gdxx - 2*gdtx*gdtz*gdzx + gdtx**2*gdzz + gdtt*(gdzx**2
     $ - gdxx*gdzz))/ (gdtt*gdzx**2*gdyy + gdtz**2*(-gdxy**2 +
     $ gdxx*gdyy) - 2*gdtt*gdxy*gdzx*gdyz - gdtx**2*gdyz**2 +
     $ gdtt*gdxx*gdyz**2 + 2*gdtz*(gdty*gdxy*gdzx - gdtx*gdzx*gdyy -
     $ gdty*gdxx*gdyz + gdtx*gdxy*gdyz) + gdtt*gdxy**2*gdzz +
     $ gdtx**2*gdyy*gdzz - gdtt*gdxx*gdyy*gdzz + gdty**2*(-gdzx**2 +
     $ gdxx*gdzz) + 2*gdtx*gdty*(gdzx*gdyz - gdxy*gdzz)) 
      guzz = 
     $ (gdty**2*gdxx - 2*gdtx*gdty*gdxy + gdtx**2*gdyy + gdtt*(gdxy**2
     $ - gdxx*gdyy))/ (gdtt*gdzx**2*gdyy + gdtz**2*(-gdxy**2 +
     $ gdxx*gdyy) - 2*gdtt*gdxy*gdzx*gdyz - gdtx**2*gdyz**2 +
     $ gdtt*gdxx*gdyz**2 + 2*gdtz*(gdty*gdxy*gdzx - gdtx*gdzx*gdyy -
     $ gdty*gdxx*gdyz + gdtx*gdxy*gdyz) + gdtt*gdxy**2*gdzz +
     $ gdtx**2*gdyy*gdzz - gdtt*gdxx*gdyy*gdzz + gdty**2*(-gdzx**2 +
     $ gdxx*gdzz) + 2*gdtx*gdty*(gdzx*gdyz - gdxy*gdzz)) 
      guxy = 
     $ (gdtz**2*gdxy + gdtt*gdzx*gdyz - gdtz*(gdty*gdzx + gdtx*gdyz) +
     $ gdtx*gdty*gdzz - gdtt*gdxy*gdzz)/ (-(gdtt*gdzx**2*gdyy) +
     $ gdtz**2*(gdxy**2 - gdxx*gdyy) + 2*gdtt*gdxy*gdzx*gdyz +
     $ gdtx**2*gdyz**2 - gdtt*gdxx*gdyz**2 - 2*gdtz*(gdty*gdxy*gdzx -
     $ gdtx*gdzx*gdyy - gdty*gdxx*gdyz + gdtx*gdxy*gdyz) -
     $ gdtt*gdxy**2*gdzz - gdtx**2*gdyy*gdzz + gdtt*gdxx*gdyy*gdzz +
     $ gdty**2*(gdzx**2 - gdxx*gdzz) + gdty*(-2*gdtx*gdzx*gdyz +
     $ 2*gdtx*gdxy*gdzz)) 
      guyz = 
     $ (-(gdty*gdtz*gdxx) + gdtx*gdtz*gdxy + gdtx*gdty*gdzx -
     $ gdtt*gdxy*gdzx - gdtx**2*gdyz + gdtt*gdxx*gdyz)/
     $ (gdtt*gdzx**2*gdyy + gdtz**2*(-gdxy**2 + gdxx*gdyy) -
     $ 2*gdtt*gdxy*gdzx*gdyz - gdtx**2*gdyz**2 + gdtt*gdxx*gdyz**2 +
     $ 2*gdtz*(gdty*gdxy*gdzx - gdtx*gdzx*gdyy - gdty*gdxx*gdyz +
     $ gdtx*gdxy*gdyz) + gdtt*gdxy**2*gdzz + gdtx**2*gdyy*gdzz -
     $ gdtt*gdxx*gdyy*gdzz + gdty**2*(-gdzx**2 + gdxx*gdzz) +
     $ 2*gdtx*gdty*(gdzx*gdyz - gdxy*gdzz)) 
      guzx = 
     $ (gdty**2*gdzx + gdtx*gdtz*gdyy - gdtt*gdzx*gdyy + gdtt*gdxy*gdyz
     $ - gdty*(gdtz*gdxy + gdtx*gdyz))/ (-(gdtt*gdzx**2*gdyy) +
     $ gdtz**2*(gdxy**2 - gdxx*gdyy) + 2*gdtt*gdxy*gdzx*gdyz +
     $ gdtx**2*gdyz**2 - gdtt*gdxx*gdyz**2 - 2*gdtz*(gdty*gdxy*gdzx -
     $ gdtx*gdzx*gdyy - gdty*gdxx*gdyz + gdtx*gdxy*gdyz) -
     $ gdtt*gdxy**2*gdzz - gdtx**2*gdyy*gdzz + gdtt*gdxx*gdyy*gdzz +
     $ gdty**2*(gdzx**2 - gdxx*gdzz) + gdty*(-2*gdtx*gdzx*gdyz +
     $ 2*gdtx*gdxy*gdzz)) 
      
      return

      end
      
      
      
      
