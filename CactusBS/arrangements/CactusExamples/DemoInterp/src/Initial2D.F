#include "cctk.h"
#include "cctk_Arguments.h"

      subroutine SetUp2DGFs(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer nx,ny
      integer i,j
      CCTK_REAL xp,yp,rp
      CCTK_REAL dx,dy

      nx = cctk_lsh(1)
      ny = cctk_lsh(2)

      dx = cctk_delta_space(1)
      dy = cctk_delta_space(2)

      do j=1,ny
         do i=1,nx
            xp = (cctk_lbnd(1)+i-1)*dx-0.5d0
            yp = (cctk_lbnd(2)+j-1)*dy-0.5d0
            rp = sqrt( xp**2+yp**2 )
            realgf2(i,j) = rp
            compgf2(i,j) = cmplx(rp,cos(rp))
         end do
      end do

      return
      end
