
  module arrays

! --->   MODULE FOR ARRAYS   <---

  implicit none
!  real(kind=8), allocatable, dimension (:,:,:) :: xx,yy,zz
  real(kind=8), allocatable, dimension (:) :: x_my
!  real(kind=8), allocatable, dimension (:,:) :: u

  real(kind=8), allocatable, dimension (:) :: a, alpha, m

  real(kind=8), allocatable, dimension (:) :: myvel,WW,p,rho0
  real(kind=8), allocatable, dimension (:) :: hh,epsilon,sound

!  real(kind=8), allocatable, dimension (:,:,:) :: r1,r2,rho_3d

! --->   END   <---

  end module arrays
