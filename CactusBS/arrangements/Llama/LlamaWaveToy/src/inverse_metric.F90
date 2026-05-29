#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"



subroutine LWT_calc_inverse_metric (CCTK_ARGUMENTS)
  use cctk
  use tensor
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  
  CCTK_REAL :: gg(3,3), dtg, gu(3,3)
  integer   :: i, j, k
  
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           
           gg(1,1) = gxx(i,j,k)
           gg(1,2) = gxy(i,j,k)
           gg(1,3) = gxz(i,j,k)
           gg(2,2) = gyy(i,j,k)
           gg(2,3) = gyz(i,j,k)
           gg(3,3) = gzz(i,j,k)
           gg(2,1) = gg(1,2)
           gg(3,1) = gg(1,3)
           gg(3,2) = gg(2,3)
           
           call calc_det (gg, dtg)
           call calc_inv (gg, dtg, gu)
           
           epsilon(i,j,k) = sqrt(dtg)
           
           guxx(i,j,k) = gu(1,1)
           guxy(i,j,k) = gu(1,2)
           guxz(i,j,k) = gu(1,3)
           guyy(i,j,k) = gu(2,2)
           guyz(i,j,k) = gu(2,3)
           guzz(i,j,k) = gu(3,3)
           
        end do
     end do
  end do
  
end subroutine LWT_calc_inverse_metric
