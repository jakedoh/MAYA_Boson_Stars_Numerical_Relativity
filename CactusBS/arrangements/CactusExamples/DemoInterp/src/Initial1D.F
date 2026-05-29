#include "cctk.h"
#include "cctk_Arguments.h"

      subroutine SetUp1DGFs(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS

      integer nx,ny
      integer i,j
      CCTK_REAL xp,rp
      CCTK_REAL dx

      nx = cctk_lsh(1)

      dx = cctk_delta_space(1)

      do i=1,nx
            xp = (cctk_lbnd(1)+i-1)*dx-0.5d0
            rp = xp
            realgf1(i) = rp
            compgf1(i) = cmplx(rp,rp)
      end do

      return
      end
