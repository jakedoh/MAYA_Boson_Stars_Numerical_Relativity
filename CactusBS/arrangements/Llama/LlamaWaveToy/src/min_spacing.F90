#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine LWT_min_spacing (CCTK_ARGUMENTS)
  use lapack
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL :: gg(3,3)
  CCTK_REAL :: ev(3)
  
  integer, parameter :: lwork = 100
  CCTK_REAL :: work(lwork)
  integer   :: info
  
  integer   :: i, j, k
  
  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
           
           gg(1,1) = gxx(i,j,k)
           gg(1,2) = gxy(i,j,k)
           gg(1,3) = gxz(i,j,k)
           gg(2,2) = gyy(i,j,k)
           gg(2,3) = gyz(i,j,k)
           gg(3,3) = gzz(i,j,k)
           gg(2,1) = gg(1,2)
           gg(3,1) = gg(1,3)
           gg(3,2) = gg(2,3)
           
           call dsyev ('N', 'U', 3, gg, 3, ev, work, lwork, info)
           if (info /= 0) then
              call CCTK_WARN (1, "error in call to SYEV")
              ev(:) = 0
           end if
           
           min_spacing(i,j,k) = sqrt(ev(1))
           
        end do
     end do
  end do
  
end subroutine LWT_min_spacing
