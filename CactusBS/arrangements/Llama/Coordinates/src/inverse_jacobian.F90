#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"


subroutine Coordinates_SetInverseJacobian (CCTK_ARGUMENTS)

  use matdet
  use matinv
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: i, j, k
  integer :: nx, ny, nz

  CCTK_REAL                 :: detJ
  CCTK_REAL, dimension(3,3) :: Jac, iJac

  volume_form_state      = store_volume_form
  inverse_Jacobian_state = store_inverse_jacobian

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

  !$omp parallel do private(i,j,k, detJ, Jac, iJac)
  do k = 1, nz
     do j = 1, ny
        do i = 1, nx

           Jac(1,1) = J11(i,j,k)
           Jac(1,2) = J12(i,j,k)
           Jac(1,3) = J13(i,j,k)
           Jac(2,1) = J21(i,j,k)
           Jac(2,2) = J22(i,j,k)
           Jac(2,3) = J23(i,j,k)
           Jac(3,1) = J31(i,j,k)
           Jac(3,2) = J32(i,j,k)
           Jac(3,3) = J33(i,j,k)

           if (store_volume_form /= 0) then
              call calc_det3 (Jac, detJ)

              ! TODO: The volume form d^3x is defined as
              !       d^3x = da db dc * det(dx / da) * M,
              !       where M is the mask storing the nominal cell volume
              volume_form(i,j,k) = detJ
           end if

           if (store_inverse_jacobian /= 0) then
              call calc_inv3 (Jac, iJac)

              iJ11(i,j,k) = iJac(1,1)
              iJ12(i,j,k) = iJac(1,2)
              iJ13(i,j,k) = iJac(1,3)
              iJ21(i,j,k) = iJac(2,1)
              iJ22(i,j,k) = iJac(2,2)
              iJ23(i,j,k) = iJac(2,3)
              iJ31(i,j,k) = iJac(3,1)
              iJ32(i,j,k) = iJac(3,2)
              iJ33(i,j,k) = iJac(3,3)
           end if

        end do
     end do
  end do

end subroutine Coordinates_SetInverseJacobian
