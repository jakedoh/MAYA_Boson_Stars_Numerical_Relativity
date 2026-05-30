/*@@
  @file      finishbrilldata.F
  @date
  @author    Carsten Gundlach (Cactus 4, Miguel Alcubierre)
  @desc
             Reconstruct metric from conformal factor.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine IDBrillData_Finish(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k

      CCTK_REAL x1,y1,z1,rho1,rho2
      CCTK_REAL phi,psi4,e2q
      CCTK_REAL zero,one

      CCTK_REAL brillq,phif

c     Numbers.

      zero = 0.0
      one  = 1.0

c     Replace flat metric left over from elliptic solve by
c     Brill metric calculated from q and Psi.

      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho2 = x1*x1 + y1*y1
               rho1 = sqrt(rho2)

               phi = phif(x1,y1)

               e2q  = exp(2.0*brillq(rho1,z1,phi))

c              Fudge division by rho^2 on axis. (Physically, y^/rho^2,
c              x^2/rho^2 and xy/rho^2 are of course regular.)
c              Transform Brills form of the physical metric to Cartesian
c              coordinates via
c
c              e^2q (drho^2 + dz^2) + rho^2 dphi^2 =
c              e^2q (dx^2 + dy^2 + dz^2) + (1-e^2q) (xdy-ydx)^2/rho^2
c
c              The individual coefficients can be read off as

               if (rho1.gt.rhofudge) then

                  gxx(i,j,k) = (e2q + (one - e2q)*y1*y1/rho2)
                  gyy(i,j,k) = (e2q + (one - e2q)*x1*x1/rho2)
                  gzz(i,j,k) = e2q
                  gxy(i,j,k) = - (one - e2q)*x1*y1/rho2

               else

c                 This fudge assumes that q = O(rho^2) near the axis. Which
c                 it should be, or the data will be singular.

                  gxx(i,j,k) = one
                  gyy(i,j,k) = one
                  gzz(i,j,k) = one
                  gxy(i,j,k) = zero

               end if

               gxz(i,j,k) = zero
               gyz(i,j,k) = zero

               kxx(i,j,k) = zero
               kyy(i,j,k) = zero
               kzz(i,j,k) = zero
               kxy(i,j,k) = zero
               kxz(i,j,k) = zero
               kyz(i,j,k) = zero

            end do
         end do
      end do

      if (CCTK_EQUALS(metric_type,"static conformal")) then

         conformal_state = 1

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  psi(i,j,k) = brillpsi(i,j,k)
                  
               end do
            end do
         end do

      else

         conformal_state = 0

         do k=1,cctk_lsh(3)
            do j=1,cctk_lsh(2)
               do i=1,cctk_lsh(1)

                  psi4  = brillpsi(i,j,k)**4

                  gxx(i,j,k) = psi4*gxx(i,j,k)
                  gyy(i,j,k) = psi4*gyy(i,j,k)
                  gzz(i,j,k) = psi4*gzz(i,j,k)
                  gxy(i,j,k) = psi4*gxy(i,j,k)
                  
               end do
            end do
         end do

      end if

      return
      end
