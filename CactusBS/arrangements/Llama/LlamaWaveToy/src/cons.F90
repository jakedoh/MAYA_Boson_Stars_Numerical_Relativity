#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine LWT_calc_constraints (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: nx, ny, nz
  CCTK_REAL, dimension(:,:,:), allocatable &
       :: dyvx, dzvx, dxvy, dzvy, dxvz, dyvz, dxu, dyu, dzu
  
  nx = cctk_lsh(1); ny = cctk_lsh(2); nz = cctk_lsh(3)

  allocate ( dyvx(nx,ny,nz), dzvx(nx,ny,nz), dxvy(nx,ny,nz), &
             dzvy(nx,ny,nz), dxvz(nx,ny,nz), dyvz(nx,ny,nz) )

  call Diff_gv (cctkGH, 1, vx, dyvx, -1)
  call Diff_gv (cctkGH, 2, vx, dzvx, -1)
  call Diff_gv (cctkGH, 0, vy, dxvy, -1)
  call Diff_gv (cctkGH, 2, vy, dzvy, -1)
  call Diff_gv (cctkGH, 0, vz, dxvz, -1)
  call Diff_gv (cctkGH, 1, vz, dyvz, -1)
  
  ! w = curl v
  wx = dyvz - dzvy
  wy = dzvx - dxvz
  wz = dxvy - dyvx
  
  deallocate ( dyvx, dzvx, dxvy, dzvy, dxvz, dyvz)

  allocate ( dxu(nx,ny,nz), dyu(nx,ny,nz), dzu(nx,ny,nz) )

  call Diff_gv (cctkGH, 0, u , dxu , -1)
  call Diff_gv (cctkGH, 1, u , dyu , -1)
  call Diff_gv (cctkGH, 2, u , dzu , -1)

  ! diff_v = v - grad u
  diff_vx = vx - dxu
  diff_vy = vy - dyu
  diff_vz = vz - dzu
  
  deallocate ( dxu, dyu, dzu )

  ! v2 = g^ij v_i v_j
  v2 =         guxx * vx**2   &
       & + 2 * guxy * vx * vy &
       & + 2 * guxz * vx * vz &
       & +     guyy * vy**2   &
       & + 2 * guyz * vy * vz &
       & +     guzz * vz**2
  
end subroutine LWT_calc_constraints



subroutine LWT_constraint_boundaries (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: ierr
  
  ierr = Boundary_SelectGroupForBC &
       (cctkGH, CCTK_ALL_FACES, +1, -1, &
       "LlamaWaveToy::constraints", "scalar")
  if (ierr/=0) call CCTK_WARN (0, "internal error")
  
  ierr = Boundary_SelectGroupForBC &
       (cctkGH, CCTK_ALL_FACES, +1, -1, &
       "LlamaWaveToy::difference_v", "scalar")
  if (ierr/=0) call CCTK_WARN (0, "internal error")
  
end subroutine LWT_constraint_boundaries
