! vim: syntax=fortran
#include "cctk.h"

module NullEvol_Mask
 contains
subroutine NullEvol_remask (k, boundary_mask, evolution_mask)
use NullGrid_Vars
   implicit none

   CCTK_INT,                  intent (in)  ::  k
   CCTK_INT, dimension (lsh(1),lsh(2)), intent (in)  :: boundary_mask
   CCTK_INT, dimension (lsh(1),lsh(2)), intent (out) :: evolution_mask

   integer n1, n2, i, j

   n1 = size(boundary_mask,1)
   n2 = size(boundary_mask,2)

   ! we could use a where, but the performance of the where sucks...

   do j = 1, n2
      do i = 1, n1
         if (k .gt. boundary_mask(i,j)) then
            evolution_mask(i,j) = 1
         else
            evolution_mask(i,j) = 0
         end if
      end do
   end do

end subroutine NullEvol_remask
end module NullEvol_Mask
