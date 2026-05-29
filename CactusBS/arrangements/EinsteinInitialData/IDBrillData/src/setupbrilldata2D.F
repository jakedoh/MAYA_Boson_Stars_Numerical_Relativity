/*@@
  @file      setupbrilldata2D.F
  @date
  @author    Carsten Gundlach (Cactus 4, Miguel Alcubierre)
  @desc
             Set up axisymmetric Brill data for elliptic solve.
  @enddesc
  @version   $Header$
@@*/

#include "cctk.h" 
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

      subroutine setupbrilldata2D(CCTK_ARGUMENTS)

c     Set up axisymmetric Brill data for elliptic solve. The elliptic
c     equation we need to solve is:
c
c     __                 2      2
c     \/  psi  +  psi ( d q  + d   q ) / 4  =  0
c       f                z      rho
c
c     where:
c
c     __
c     \/  =  Flat space Laplacian.
c       f

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k

      CCTK_REAL x1,y1,z1,rho1
      CCTK_REAL brillq,eps
      CCTK_REAL zp,zm,rhop,rhom
      CCTK_REAL zero,one

      external brillq

c     Numbers.

      zero = 0.0D0
      one  = 1.0D0

c     Epsilon for finite differencing.

      eps = cctk_delta_space(1)

c     Set up coefficient arrays for elliptic solve.
c     Notice that the Cactus conventions are:
c     __
c     \/ psi +  Mlinear*psi  +  Nsource  =  0

      do k=1,cctk_lsh(3)
         do j=1,cctk_lsh(2)
            do i=1,cctk_lsh(1)

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho1 = sqrt(x1*x1 + y1*y1)

               brillpsi(i,j,k) = one

               gxx(i,j,k) = one
               gyy(i,j,k) = one
               gzz(i,j,k) = one

               gxy(i,j,k) = zero
               gxz(i,j,k) = zero
               gyz(i,j,k) = zero

c              Centered derivatives. Note that here we may be calling brillq
c              with a small negative rho, but that should be ok as long as
c              brillq is even in rho - physically it must be, or the data
c              will not be regular on the axis.

               zp = z1 + eps
               zm = z1 - eps

               rhop = rho1 + eps
               rhom = rho1 - eps

               brillMlinear(i,j,k) = 0.25
     .              *(brillq(rho1,zp,zero)
     .              + brillq(rho1,zm,zero)
     .              + brillq(rhop,z1,zero)
     .              + brillq(rhom,z1,zero)
     .              - 4.0*brillq(rho1,z1,zero))/eps**2

               brillNsource(i,j,k) = zero

            end do
         end do
      end do

c     Synchronization and boundaries.

c      call CCTK_SyncGroup(ierr,cctkGH,'idbrilldata::brillelliptic')
c      call CartSymGN(ierr,cctkGH,'idbrilldata::brillelliptic')

      return
      end
