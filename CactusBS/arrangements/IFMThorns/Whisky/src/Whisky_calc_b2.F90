#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#define Bvx(i,j,kk) Bvec(i,j,kk,1)
#define Bvy(i,j,kk) Bvec(i,j,kk,2)
#define Bvz(i,j,kk) Bvec(i,j,kk,3)

 /*@@
   @routine    whisky_b2comoving
   @author     Tanja Bode
   @desc 
   Calculate and store b2 comoving 
   @enddesc 
@@*/
 
subroutine Whisky_b2comoving( CCTK_ARGUMENTS )

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT :: i,j,kk, nx,ny,nz
  CCTK_REAL :: vlowx, vlowy, vlowz, Bvec2, vdotB

  CCTK_REAL :: lgxx, lgxy, lgxz, lgyy, lgyz, lgzz
  CCTK_REAL :: lBx, lBy, lBz

  !! Loop over everywhere
  ! kill isolated cells or cells whose density is too low
  ! NB: this will leave a one cell thick zone of alive cells, which are required
  ! to have material flow into grid level
  !$OMP PARALLEL DO PRIVATE (i,j,kk,inv_w2,vlowx,vlowy,vlowz,Bvec2,&
  !$OMP     lgxx,lgxy,lgxz,lgyy,lgyz,lgzz, lBx,lBy,lBz, vdotB)
  do kk = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           !! Local copies
           lgxx = gxx(i,j,kk) 
           lgxy = gxy(i,j,kk) 
           lgxz = gxz(i,j,kk) 
           lgyy = gyy(i,j,kk) 
           lgyz = gyz(i,j,kk) 
           lgzz = gzz(i,j,kk) 

           lBx = Bvx(i,j,kk)
           lBy = Bvy(i,j,kk)
           lBz = Bvz(i,j,kk)

           !! Projections
           vlowx = lgxx*velx(i,j,kk) + lgxy*vely(i,j,kk) + lgxz*velz(i,j,kk)
           vlowy = lgxy*velx(i,j,kk) + lgyy*vely(i,j,kk) + lgyz*velz(i,j,kk)
           vlowz = lgxz*velx(i,j,kk) + lgyz*vely(i,j,kk) + lgzz*velz(i,j,kk)

           Bvec2 = lgxx*lBx**2 + lgyy*lBy**2 + lgzz*lBz**2 + 2.d0*lgxy*lBx*lBy &
                   + 2.d0*lgxz*lBx*lBz + 2.d0*lgyz*lBy*lBz
           vdotB = vlowx*lBx + vlowy*lBy + vlowz*lBz

           !! Our goal
           b2_comoving(i,j,kk) = Bvec2/(w_lorentz(i,j,kk)**2) + vdotB**2 

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine Whisky_b2comoving
