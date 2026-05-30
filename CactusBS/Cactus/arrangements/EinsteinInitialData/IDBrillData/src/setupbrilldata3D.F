/*@@
  @file      setupbrilldata3D.F
  @date      October 1998
  @author    Miguel Alcubierre
  @desc
             Set up non-axisymmetric Brill data for elliptic solve.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine setupbrilldata3D(CCTK_ARGUMENTS)

c     Set up 3D Brill data for elliptic solve.  The elliptic
c     equation we need to solve is:
c
c     __
c     \/  psi  -  psi R  /  8  =  0
c       c              c
c
c     where:
c     __
c     \/  =  Laplacian operator for conformal metric.
c       c
c
c     R   =  Ricci scalar for conformal metric. 
c      c
c
c     The Ricci scalar for the conformal metric turns out to be:
c
c               /  -2q    2      2             -2            2       2      \
c     R  =  - 2 | e    ( d q  + d   q )  +  rho   ( 3 (d   q)  +  2 d   q ) |
c      c        \         z      rho                    phi          phi    /


      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer i,j,k

      CCTK_REAL x1,y1,z1,rho1,rho2
      CCTK_REAL phi,phif,e2q
      CCTK_REAL brillq,eps
      CCTK_REAL zero,one,two,three

c     Define numbers

      zero  = 0.0D0
      one   = 1.0D0
      two   = 2.0D0
      three = 3.0D0

c     Epsilon for finite differencing.

      eps = cctk_delta_space(1)

      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho2 = x1*x1 + y1*y1
               rho1 = sqrt(rho2)
               phi = phif(x1,y1)

               e2q  = exp(two*brillq(rho1,z1,phi))

c              Initialise psi

               brillpsi(i,j,k) = one

c              Set up metric and coefficient arrays for elliptic solve.
c              Notice that the Cactus conventions are:
c              __
c              \/ psi +  Mlinear*psi  +  Nsource  =  0

c              Find M using centered differences over a small
c              interval.

c              Here we assume that for very small rho, the
c              phi derivatives are essentially zero.  This
c              must always be true otherwise the function
c              is not regular on the axis.

               if (rho1.gt.rhofudge) then

                  gxx(i,j,k) = e2q + (one - e2q)*y1*y1/rho2
                  gyy(i,j,k) = e2q + (one - e2q)*x1*x1/rho2
                  gzz(i,j,k) = e2q
                  gxy(i,j,k) = - (one - e2q)*x1*y1/rho2

                  brillMlinear(i,j,k) = 0.25/e2q
     .                 *(brillq(rho1,z1+eps,phi) 
     .                 + brillq(rho1,z1-eps,phi) 
     .                 + brillq(rho1+eps,z1,phi) 
     .                 + brillq(rho1-eps,z1,phi) 
     .                 - 4.0*brillq(rho1,z1,phi))
     .                 / eps**2

                  brillMlinear(i,j,k) = brillMlinear(i,j,k) + 0.25/rho2
     .               *(three*0.25*(brillq(rho1,z1,phi+eps)
     .               - brillq(rho1,z1,phi-eps))**2
     .               + two*(brillq(rho1,z1,phi+eps)
     .               - two*brillq(rho1,z1,phi)
     .               + brillq(rho1,z1,phi-eps)))/eps**2

               else

                  gxx(i,j,k) = one 
                  gyy(i,j,k) = one 
                  gzz(i,j,k) = one
                  gxy(i,j,k) = zero

                  brillMlinear(i,j,k) = 0.25/e2q
     .                 *(brillq(rho1,z1+eps,phi) 
     .                 + brillq(rho1,z1-eps,phi) 
     .                 + two*brillq(eps,z1,phi) 
     .                 - two*brillq(rho1,z1,phi))
     .                 / eps**2

               end if

               gxz(i,j,k) = zero
               gyz(i,j,k) = zero

c              Set coefficient Nsource = 0
               brillNsource(i,j,k) = zero

            end do
         end do
      end do


c     Synchronization and boundaries.

c      call CCTK_SyncGroup(ierr,cctkGH,'idbrilldata::brillelliptic')
c      call CartSymGN(ierr,cctkGH,'idbrilldata::brillelliptic')

      return
      end
