! vim: syntax=fortran
#include "cctk.h"
!#include "cctk_Functions.h"
!#include "cctk_Arguments.h"
!#include "cctk_Parameters.h"


module NullNews_ScriUtil

  implicit none

contains

  subroutine NullNews_cscrival (lsh, nx, jnn, jns, J)

    implicit none

    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_INT,               intent(in) :: nx

    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx) :: jnn, jns
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),2)  :: J

    J(:,:,1) = jnn(:,:,nx)
    J(:,:,2) = jns(:,:,nx)

  end subroutine NullNews_cscrival

  subroutine NullNews_rscrival (lsh, nx, bnn, bns, beta)

    implicit none

    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_INT,               intent(in) :: nx 

    CCTK_REAL, dimension (lsh(1),lsh(2),nx) :: bnn, bns
    CCTK_REAL, dimension (lsh(1),lsh(2),2)  :: beta

    beta(:,:,1) = bnn(:,:,nx)
    beta(:,:,2) = bns(:,:,nx)

  end subroutine NullNews_rscrival

  subroutine NullNews_cscrivalh (lsh, nx, unn, uns, U)

    implicit none

    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_INT,               intent(in) :: nx

    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx) :: unn, uns
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),2)  :: U

    U(:,:,1) = 0.5d0 * (unn(:,:,nx) + unn(:,:,nx-1))
    U(:,:,2) = 0.5d0 * (uns(:,:,nx) + uns(:,:,nx-1))

  end subroutine NullNews_cscrivalh

  subroutine NullNews_cscridbydl (l_deriv, lsh, nx, dx, rwt, jnn, jns, J_l)
    implicit none

    CCTK_INT, intent(in) :: l_deriv
    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_INT,               intent(in) :: nx
    CCTK_REAL,              intent(in) :: rwt, dx

    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx) :: jnn, jns
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),2)  :: J_l
    CCTK_REAL :: A, B, C, D, E
   
    select case (l_deriv)
    case (1)
      A = 1.d0; B=-1.d0
      J_l(:,:,1)=-(A*jnn(:,:,nx) + B*jnn(:,:,nx-1))/(dx)*rwt
      J_l(:,:,2)=-(A*jns(:,:,nx) + B*jns(:,:,nx-1))/(dx)*rwt
    case (2)
      A = 3.d0/2.d0; B = -2.d0; C = 1/2.d0
      J_l(:,:,1) = -(A*jnn(:,:,nx) + B*jnn(:,:,nx-1) + C*jnn(:,:,nx-2))/(dx)*rwt
      J_l(:,:,2) = -(A*jns(:,:,nx) + B*jns(:,:,nx-1) + C*jns(:,:,nx-2))/(dx)*rwt

!    J_l(:,:,1) = - 0.5d0 * ( 3.0d0 * jnn(:,:,nx) - 4.0d0 * jnn(:,:,nx-1) &
!         + jnn(:,:,nx-2) ) / (dx ) * rwt
!
!   J_l(:,:,2) = - 0.5d0 * ( 3.0d0 * jns(:,:,nx) - 4.0d0 * jns(:,:,nx-1) &
!         + jns(:,:,nx-2) ) / (dx ) * rwt
    case (3)
      A = 11.d0/6.d0; B = -3.d0; C = 3.d0/2.d0; D = -1/3.d0
      J_l(:,:,1) = -(A*jnn(:,:,nx) + B*jnn(:,:,nx-1) &
                 + C*jnn(:,:,nx-2) + D*jnn(:,:,nx-3))/(dx)*rwt
      J_l(:,:,2) = -(A*jns(:,:,nx) + B*jns(:,:,nx-1) &
                 + C*jns(:,:,nx-2) + D*jns(:,:,nx-3))/(dx)*rwt
    case (4)
      A = 25.d0/12.d0; B = -4.d0; C = 3.d0; D = -4.d0/3.d0; E = 1/4.d0 
      J_l(:,:,1) = -(A*jnn(:,:,nx) + B*jnn(:,:,nx-1) + C*jnn(:,:,nx-2)&
                 + D*jnn(:,:,nx-3) + E*jnn(:,:,nx-4))/(dx)*rwt
      J_l(:,:,2) = -(A*jns(:,:,nx) + B*jns(:,:,nx-1) + C*jns(:,:,nx-2) &
                 + D*jns(:,:,nx-3) + E*jns(:,:,nx-4))/(dx)*rwt
    case default 
         call CCTK_WARN(1, "unrecognized value for l_deriv")
    end select

  end subroutine NullNews_cscridbydl

  subroutine NullNews_cscridbydlh (lsh, nx, dx, rwt, unn, uns, U_l)
    implicit none

    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_INT,               intent(in) :: nx
    CCTK_REAL,              intent(in) :: rwt, dx

    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx) :: unn, uns
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),2)  :: U_l

    U_l(:,:,1) = - 0.5d0 * (unn(:,:,nx) - unn(:,:,nx-1) ) / (dx ) * rwt

    U_l(:,:,2) = - 0.5d0 * (uns(:,:,nx) - uns(:,:,nx-1) ) / (dx ) * rwt

  end subroutine NullNews_cscridbydlh

  subroutine NullNews_cscridbydl2h (lsh, nx, dx, rwt, unn, uns, U_l_l)
    implicit none

    CCTK_INT, dimension(2), intent(in) :: lsh
    CCTK_INT,               intent(in) :: nx
    CCTK_REAL,              intent(in) :: rwt, dx

    CCTK_COMPLEX, dimension (lsh(1),lsh(2),nx) :: unn, uns
    CCTK_COMPLEX, dimension (lsh(1),lsh(2),2)  :: U_l_l

    U_l_l(:,:,1) = rwt**2 * ( (unn(:,:,nx) - unn(:,:,nx-1) ) / dx  &
         + (3.0d0 * unn(:,:,nx) - 7.0d0 * unn(:,:,nx-1) + &
         5.0d0* unn(:,:,nx-2) -unn(:,:,nx-3)) /2/dx**2 )

    U_l_l(:,:,2) = rwt**2 * ( (uns(:,:,nx) - uns(:,:,nx-1) ) / dx  &
         + (3.0d0 * unn(:,:,nx) - 7.0d0 * uns(:,:,nx-1) + &
         5.0d0* uns(:,:,nx-2) -uns(:,:,nx-3)) /2/dx**2 )

  end subroutine NullNews_cscridbydl2h

  subroutine NullNews_ResetInactive(lsh, fscri)
    use NullGrid_Vars, only: EG
    implicit none

    CCTK_INT,     dimension(2),               intent(in)    :: lsh
    CCTK_COMPLEX, dimension(lsh(1),lsh(2),2), intent(inout) :: fscri

    integer :: i, j

    do j = 1, lsh(2)
       do i = 1, lsh(1)
          if (EG(i,j).lt.0.1) then
             fscri(i,j,:) = 0
          end if
       end do
    end do
!   fscri(:,:,1) = fscri(:,:,1) * dble(EG)
!   fscri(:,:,2) = fscri(:,:,2) * dble(EG)

  end subroutine NullNews_ResetInactive

  subroutine NullNews_ResetInactiveRe(lsh, fscri, reset_value)
    use NullGrid_Vars, only: EG
    implicit none

    CCTK_INT,  dimension(2),               intent(in)    :: lsh
    CCTK_REAL, dimension(lsh(1),lsh(2),2), intent(inout) :: fscri
    CCTK_REAL,                             intent(in)    :: reset_value

    integer :: i, j

    do j = 1, lsh(2)
       do i = 1, lsh(1)
          if (EG(i,j).lt.0.1) then
             fscri(i,j,:) = reset_value
          end if
       end do
    end do
!   fscri(:,:,1) = fscri(:,:,1) * dble(EG)
!   fscri(:,:,2) = fscri(:,:,2) * dble(EG)

  end subroutine NullNews_ResetInactiveRe

end module NullNews_ScriUtil

