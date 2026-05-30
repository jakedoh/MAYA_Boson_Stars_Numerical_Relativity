#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#define Bvx(i,j,kk) Bvec(i,j,kk,1)
#define Bvy(i,j,kk) Bvec(i,j,kk,2)
#define Bvz(i,j,kk) Bvec(i,j,kk,3)

! Second order f.d.

#define DIFF_X_2(q) 0.5d0 * (q(i+1,j,kk) - q(i-1,j,kk)) * idx
#define DIFF_Y_2(q) 0.5d0 * (q(i,j+1,kk) - q(i,j-1,kk)) * idy
#define DIFF_Z_2(q) 0.5d0 * (q(i,j,kk+1) - q(i,j,kk-1)) * idz

! Fourth order f.d.

#define DIFF_X_4(q) (-q(i+2,j,kk) + 8.d0 * q(i+1,j,kk) - 8.d0 * q(i-1,j,kk) + \
                      q(i-2,j,kk)) / 12.d0 * idx
#define DIFF_Y_4(q) (-q(i,j+2,kk) + 8.d0 * q(i,j+1,kk) - 8.d0 * q(i,j-1,kk) + \
                      q(i,j-2,kk)) / 12.d0 * idy
#define DIFF_Z_4(q) (-q(i,j,kk+2) + 8.d0 * q(i,j,kk+1) - 8.d0 * q(i,j,kk-1) + \
                      q(i,j,kk-2)) / 12.d0 * idz

 /*@@
   @routine    whisky_recalc_DivB
   @author     Tanja Bode
   @desc 
   Recalculate div(B).
   @enddesc 
@@*/
 
subroutine whisky_recalc_DivB( CCTK_ARGUMENTS )

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i,j,kk, nx,ny,nz
  CCTK_REAL :: dxBx, dyBy, dzBz, detg
  CCTK_REAL :: dx, dy, dz, idx, idy, idz

  !! Currently uses straight-up finite differencing at the cell center
  !!  XXX: Upgrade to do simple reconstruction of Bvec
  !!       to cell faces to monitor cell average of divB
  !!       instead?

  divB = 0.d0
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  dx = CCTK_DELTA_SPACE(1)
  dy = CCTK_DELTA_SPACE(2)
  dz = CCTK_DELTA_SPACE(3)
  idx = 1.d0/dx
  idy = 1.d0/dy
  idz = 1.d0/dz

  !! Loop over interior
  do kk = whisky_stencil+1, nz-whisky_stencil
     do j = whisky_stencil+1, ny-whisky_stencil
        do i = whisky_stencil+1, nx-whisky_stencil

           !! if ( spatial_order.eq.2 ) then

              dxBx = DIFF_X_2(Bnx)
              dyBy = DIFF_Y_2(Bny)
              dzBz = DIFF_Z_2(Bnz)

           !!else

           !!   dxBx = DIFF_X_4(Bnx)
           !!   dyBy = DIFF_Y_4(Bny)
           !!   dzBz = DIFF_Z_4(Bnz)

           !!end if

           detg = - gxz(i,j,kk)**2 * gyy(i,j,kk) &
             - gxx(i,j,kk)*gyz(i,j,kk)**2 - gxy(i,j,kk)**2*gzz(i,j,kk) &
             + 2.d0*gxy(i,j,kk)*gxz(i,j,kk)*gyz(i,j,kk) &
             + gxx(i,j,kk)*gyy(i,j,kk)*gzz(i,j,kk)

           divB(i,j,kk) = (dxBx+dyBy+dzBz) !!/sqrt(detg)

        end do
     end do
  end do

end subroutine whisky_recalc_DivB
