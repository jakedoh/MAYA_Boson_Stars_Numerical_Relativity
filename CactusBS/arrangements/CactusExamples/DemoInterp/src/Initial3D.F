#include "cctk.h"
#include "cctk_Arguments.h"

      subroutine SetUp3DGFs(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer nx,ny,nz
      integer i,j,k
      CCTK_REAL xp,yp,zp,rp
      CCTK_REAL dx,dy,dz

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)
      nz = cctk_lsh(3)

      dx = 1.0d0/(cctk_gsh(1)-1)
      dy = 1.0d0/(cctk_gsh(2)-1)
      dz = 1.0d0/(cctk_gsh(3)-1)

      do k=1,nz
         do j=1,ny
            do i=1,nx
               xp = (cctk_lbnd(1)+i-1)*dx-0.5d0
               yp = (cctk_lbnd(2)+j-1)*dy-0.5d0
               zp = (cctk_lbnd(3)+k-1)*dz-0.5d0
               rp = sqrt( xp**2+yp**2+zp**2 )
               realgf3(i,j,k) = rp
               compgf3(i,j,k) = cmplx(rp,cos(rp))
             end do
          end do
       end do

       cctk_delta_space(1) = dx
       cctk_delta_space(2) = dy
       cctk_delta_space(3) = dz

       return
       end
