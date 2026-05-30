#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


      
subroutine LWT_CalcEnergy (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  CCTK_REAL :: dx(3), dv
  integer :: nx, ny, nz
  CCTK_REAL, dimension(:,:,:), allocatable &
       :: hxx, hxy, hxz, hyy, hyz, hzz, vv2
  
  nx = cctk_lsh(1); ny = cctk_lsh(2); nz = cctk_lsh(3)

  allocate ( hxx(nx,ny,nz), hxy(nx,ny,nz), hxz(nx,ny,nz), &
             hyy(nx,ny,nz), hyz(nx,ny,nz), hzz(nx,ny,nz), vv2(nx,ny,nz) )

  dx(:) = CCTK_DELTA_SPACE(:)
  dv = product(dx)
  
  ! H^ij = alpha^2 g^ij - beta^i beta^j
  
  hxx = guxx * alpha**2 - betax * betax
  hxy = guxy * alpha**2 - betax * betay
  hxz = guxz * alpha**2 - betax * betaz
  hyy = guyy * alpha**2 - betay * betay
  hyz = guyz * alpha**2 - betay * betaz
  hzz = guzz * alpha**2 - betaz * betaz
  
  ! E = 1/2 mask (rho^2 + H^ij v_i v_j) epsilon / alpha dV
  
  vv2 =        hxx * vx**2   &
       & + 2 * hxy * vx * vy &
       & + 2 * hxz * vx * vz &
       & +     hyy * vy**2   &
       & + 2 * hyz * vy * vz &
       & +     hzz * vz**2
  
  energy = 0.5d0 * nmask * (rho**2 + vv2) * epsilon / alpha * dv
  
  deallocate ( hxx, hxy, hxz, hyy, hyz, hzz, vv2 )

end subroutine LWT_CalcEnergy
