#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine LWT_error (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  integer :: nx, ny, nz
  CCTK_REAL :: timetmp
  CCTK_REAL, dimension(:,:,:), allocatable &
       :: utmp, rhotmp

  nx = cctk_lsh(1); ny = cctk_lsh(2); nz = cctk_lsh(3)

  allocate ( utmp(nx,ny,nz), rhotmp(nx,ny,nz) )
  
  !!!!!!!!!!!!! Compute the error with respect to the analytic solution !!!!!!!!!!!!!!!!  
  utmp   = u
  rhotmp = rho
  
  call LWT_init (CCTK_PASS_FTOF)
!  call LWT_transform (CCTK_PASS_FTOF)
  
  exact = u
  exact_rho = rho

  error = utmp - u
  error_rho = rhotmp - rho

  u   = utmp
  rho = rhotmp
  
 !!!!!!!!!!!! Compute the error with respect to a solution that is only known at some
 !!!!!!!!!!!! periodic point in time !!!!!!!!!!!!
  utmp   = u
  rhotmp = rho
  timetmp = cctk_time
  cctk_time = 0 

  call LWT_init (CCTK_PASS_FTOF)

  errorperiodic = utmp - u
  errorperiodic_rho = rhotmp - rho

  u   = utmp
  rho = rhotmp
  cctk_time = timetmp
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  deallocate ( utmp, rhotmp )
end subroutine LWT_error
