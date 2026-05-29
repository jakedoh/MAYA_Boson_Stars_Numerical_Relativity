! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


integer function NullGrid_RegisterCoords()
  implicit none

  integer :: retval,ierr
  integer ::  one = 1
  integer ::  two = 2
  retval = 0

  call CCTK_CoordRegisterSystem(ierr,two,"stereo")

  call CCTK_CoordRegisterData(ierr, one, "NullGrid::stereo_q", "q", "stereo")

  if (ierr .lt. 0) then
     call CCTK_WARN(one,"Problem with registering coordinate q")
     retval = -1
  end if
  call CCTK_CoordRegisterData(ierr, two, "NullGrid::stereo_p", "p", "stereo")
  if (ierr .lt. 0) then
     call CCTK_WARN(one,"Problem with registering coordinate p")
     retval = -1
  end if

  NullGrid_RegisterCoords = retval

  return
end function NullGrid_RegisterCoords



subroutine NullGrid_Coord(CCTK_ARGUMENTS)

  use cctk
  implicit none

  INTEGER   ::   i,j,ierr
  integer   :: two = 2
  integer   :: tmp(2)

  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  !   Need to do this because of mixing gf dims
  call CCTK_GroupgshGN(ierr,cctkGH,two,tmp,"NullGrid::StCrd")
  null_gsh = tmp
  if (ierr .ne. 0 ) then 
     call CCTK_WARN(0,"Could not get gsh")
  endif
  call CCTK_GrouplshGN(ierr,cctkGH,two,tmp,"NullGrid::StCrd")
  null_lsh = tmp
  if (ierr .ne. 0 ) then 
     call CCTK_WARN(0,"Could not get lsh")
  endif
  call CCTK_GrouplbndGN(ierr,cctkGH,two,tmp,"NullGrid::StCrd")
  null_lbnd = tmp
  if (ierr .ne. 0 ) then 
     call CCTK_WARN(0,"Could not get lbnd")
  endif
  call CCTK_GroupubndGN(ierr,cctkGH,two,tmp,"NullGrid::StCrd")
  null_ubnd = tmp
  if (ierr .ne. 0 ) then 
     call CCTK_WARN(0,"Could not get ubnd")
  endif

  if (null_gsh(1).ne.null_gsh(2)) call CCTK_WARN(0, "grid setup error");

  null_delta = 2.0d0 / (N_ang_pts_inside_eq-1)
  qsize = 1 + null_delta(1) * ( N_ang_stencil_size + N_ang_ev_outside_eq )
! write (*,*) 'qsize = ', qsize

  do j=1,null_lsh(2)
     do i=1,null_lsh(1)
        stereo_q(i,j)  = -qsize + null_delta(1)*(i-1+null_lbnd(1))
        stereo_p(i,j)  = -qsize + null_delta(2)*(j-1+null_lbnd(2))
     end do
  end do

! write (*,*) 'local qrange:', minval(stereo_q), maxval(stereo_q)
! write (*,*) 'local prange:', minval(stereo_p), maxval(stereo_p)

  zeta      = dcmplx(stereo_q,stereo_p)
  stereo_pp = 1. + stereo_q**2 + stereo_p**2

  null_dx = (1.0d0 - null_xin) / dble(N_radial_pts - 1)

  do i = 1, N_radial_pts
     null_xb(i)  = null_xin + (i - 1  ) * null_dx
     null_xbh(i) = null_xin + (i - 0.5) * null_dx
  end do

  do i = 1, N_radial_pts - 1
     null_rb(i)  = null_rwt * null_xb(i) / (1. - null_xb(i))
     null_rbh(i) = null_rwt * null_xbh(i) / (1. - null_xbh(i))
  end do
  null_rb(N_radial_pts) = null_rb(N_radial_pts-1)
  null_rbh(N_radial_pts) = null_rbh(N_radial_pts-1)


end subroutine  NullGrid_Coord
