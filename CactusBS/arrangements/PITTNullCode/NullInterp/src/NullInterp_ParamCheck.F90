! vim: syntax=fortran
! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NullInterp_ParamCheck(CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (2*N_ang_ghost_pts.lt.interpolation_order+1) then
     call CCTK_PARAMWARN("Must set NullGrid::N_ang_ghost_pts >= ( NullGrid::interpolation_order+1 ) / 2")
  end if

  if (N_ang_ghost_pts.lt.N_ang_stencil_size) then
     call CCTK_PARAMWARN("Must set NullGrid::N_ang_ghost_pts >= NullGrid::N_ang_stencil_size")
  end if

  if (deriv_accuracy.eq.4 .and. N_ang_ghost_pts.lt.4) then
     call CCTK_PARAMWARN("If deriv_accuracy==4, then N_ang_ghost_pts >=4 is required")
  end if

  if (deriv_accuracy.eq.4 .and. N_ang_stencil_size.lt.4) then
     call CCTK_PARAMWARN("If deriv_accuracy==4, then N_stencil_size >=4 is required")
  end if


! derivative_accuracy = deriv_accuracy

! if (deriv_accuracy.eq.4) then
!    interp_width = 2
!    select case (interpolation_order)
!    case (1,2)
!       call CCTK_PARAMWARN("too low interpolation accuracy for this derivative accuracy")
!    case (3)
!       call CCTK_WARN(1, "should use interpolation_order=3 for this derivative accuracy")
!    end select
! else 
!    interp_width = 1
! end if

end subroutine NullInterp_ParamCheck
