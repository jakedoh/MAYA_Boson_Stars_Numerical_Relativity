#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"



subroutine LWT_outerboundary (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS  
  
  integer :: i, j, k
  
!$OMP PARALLEL DO private(i,j,k)
  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)
           if (Sn(i,j,k) == -2) then
              
              if (CCTK_EQUALS(outer_bound,"solution")) then
                 call exactsolution (u(i,j,k), rho(i,j,k), x(i,j,k), y(i,j,k), z(i,j,k), cctk_time)
              else if (CCTK_EQUALS(outer_bound, "dirichlet")) then
                 u(i,j,k) = 0
                 rho(i,j,k) = 0
              else if (CCTK_EQUALS(outer_bound, "radiative") .OR. CCTK_EQUALS(outer_bound, "none")) then
                 ! Radiative outer boundary is applied to the RHS
                 continue;
              end if
              
           end if
        end do
     end do
  end do
!$OMP END PARALLEL DO
  
end subroutine LWT_outerboundary


subroutine LWT_RHS_outerboundary (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  integer :: i, j, k

  if (CCTK_EQUALS(outer_bound, "radiative")) then
    if ( (CCTK_Equals(coordinate_system, "Thornburg04").EQ.0) .AND. (CCTK_Equals(coordinate_system, "Thornburg04nc").EQ.0) ) then
      call CCTK_WARN (0, "Radiative boundary condition only supports Thornburg04 coordinates.");
    end if
!$OMP PARALLEL DO private(i,j,k)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        if (Sn(i,j,cctk_lsh(3)) == -2) then
          rhodot(i,j,cctk_lsh(3)) = -(rho(i,j,cctk_lsh(3)) - rho(i,j,cctk_lsh(3)-1)) / CCTK_DELTA_SPACE(3) - rho(i,j,cctk_lsh(3)) / r(i,j,cctk_lsh(3));
        end if
      end do
    end do
!$OMP END PARALLEL DO
   end if
end subroutine LWT_RHS_outerboundary
