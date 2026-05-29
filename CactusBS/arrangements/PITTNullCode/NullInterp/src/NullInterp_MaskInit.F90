! vim: syntax=fortran
#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

subroutine NullInterp_MaskInit(CCTK_ARGUMENTS)

  use cctk
  use NullGrid_Vars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  if (CCTK_EQUALS(stereo_patch_type,"square")) then
     call CCTK_INFO("will set up square-shaped guard point mask");
     call setup_square_masks
  else if (CCTK_EQUALS(stereo_patch_type,"circle")) then
     call CCTK_INFO("will set up circle-shaped guard point mask");
     call setup_circle_masks
  else
     call CCTK_WARN(0, "unsupported stereographic patch type");
  end if

  if (extended_guard_zone.ne.0) call extend_guard_zone

  ! check that guard zone is far enough from equator

  call test_interpolation

contains

  subroutine extend_guard_zone

    guard_mask = guard_mask + (1-EG_mask)
    EG_mask = 1
    if (any((guard_mask.ne.0).and.(guard_mask.ne.1))) stop 'inconsistent guard mask'

  end subroutine extend_guard_zone

  subroutine setup_square_masks

    CCTK_INT :: i, j

    EG_mask = 1
    guard_mask = 0

    if (lbnd(1).eq.0) guard_mask(1:N_ang_stencil_size,:) = 1
    if (ubnd(1).eq.gsh(1)-1) guard_mask(lsh(1)-N_ang_stencil_size+1:,:) = 1

    if (lbnd(2).eq.0) guard_mask(:,1:N_ang_stencil_size) = 1
    if (ubnd(2).eq.gsh(2)-1) guard_mask(:,lsh(2)-N_ang_stencil_size+1:) = 1


    ! mark equator points and evolution (finite difference) points

    EQ_mask = 1
    EV_mask = 1

    do j = 1, null_lsh(2)
       do i = 1, null_lsh(1)
          if (dsqrt(qs(i,j)**2 + ps(i,j)**2) > (1 - 1.0e-10)) EQ_mask(i,j) = 0
       end do
    end do

    EV_mask = EV_mask - guard_mask

  end subroutine setup_square_masks

  subroutine setup_circle_masks

    CCTK_REAL :: r2_in
    CCTK_INT :: i, j, imin, imax, jmin, jmax
    character(len=500) :: message

    r2_in  = (1 + N_ang_ev_outside_eq * maxval(null_delta))**2 - 1.0e-10
    evolution_radius = 1 + N_ang_ev_outside_eq * maxval(null_delta)

    write (message,*) 'setting up a circular evolution mask of radius', evolution_radius
    call CCTK_INFO(message)

    write (message,*) 'the guard point shell is set for a max. stencil size of', N_ang_stencil_size
    call CCTK_INFO(message)

    ! mark evolution points

    EG_mask = 1

    do j = 1, null_lsh(2)
       do i = 1, null_lsh(1)
          if (qs(i,j)**2 + ps(i,j)**2 > r2_in) EG_mask(i,j) = 0
       end do
    end do


    ! check that there is enough room left
    ! for guard points near the patch boundaries

    if (((lbnd(1).eq.0).and.any(EG_mask(1:N_ang_stencil_size,:).ne.0)) .or. &
         ((ubnd(1).eq.gsh(1)-1).and.any(EG_mask(lsh(1)-N_ang_stencil_size+1:,:).ne.0)) .or. &
         ((lbnd(2).eq.0).and.any(EG_mask(:,1:N_ang_stencil_size).ne.0)) .or. &
         ((ubnd(2).eq.gsh(2)-1).and.any(EG_mask(:,lsh(2)-N_ang_stencil_size+1:).ne.0))) &
         call CCTK_WARN(0, "mask setup error");

    ! mark guard points

    ! identify union of evolution and guard points

    guard_mask = EG_mask

    do j = 1, null_lsh(2)
       do i = 1, null_lsh(1)
          if (EG_mask(i,j).eq.1) then

             imin = max(1, i-N_ang_stencil_size)
             jmin = max(1, j-N_ang_stencil_size)

             imax = min(null_lsh(1), i+N_ang_stencil_size)
             jmax = min(null_lsh(2), j+N_ang_stencil_size)

             guard_mask(imin:imax,jmin:jmax) = 1

          end if
       end do
    end do

    ! substract set of evolution points

    guard_mask = guard_mask - EG_mask

    if (any(guard_mask.lt.0).or.any(guard_mask.gt.1))&
         call CCTK_WARN(0, "error in setting up guard mask");

    ! EG mask is union of evolution and guard points

    EG_mask = EG_mask + guard_mask

    if (any(EG_mask.lt.0).or.any(EG_mask.gt.1))&
         call CCTK_WARN(0, "error in setting up EG_mask");

    ! mark equator points and evolution (finite difference) points

    EQ_mask = 1

    do j = 1, null_lsh(2)
       do i = 1, null_lsh(1)
          if (dsqrt(qs(i,j)**2 + ps(i,j)**2) > (1 - 1.0e-10)) EQ_mask(i,j) = 0
       end do
    end do

    EV_mask = EG_mask - guard_mask

  end subroutine setup_circle_masks

  subroutine test_interpolation
    use NullInterp_Interp

    integer :: i, retval, reduction_handle
    CCTK_REAL :: global_value, interp_err, minimum_r_guard
    character(len=500) :: message

    tmp_maskn = (1-(EG_mask-guard_mask)) * poison_value
    tmp_masks = (1-(EG_mask-guard_mask)) * poison_value

    call NullInterp_rinterp(cctkGH, tmp_rgfn, tmp_rgfs, tmp_maskn, tmp_masks)

    interp_err = max(maxval(abs(tmp_maskn), guard_mask.eq.1),&
         maxval(abs(tmp_masks), guard_mask.eq.1))

    call CCTK_ReductionArrayHandle(reduction_handle, "maximum");
    if (reduction_handle .lt. 0 ) then
       call CCTK_WARN(0,"Could not get reduction handle")
    end if

    call CCTK_ReduceLocScalar(retval, cctkGH, -1, reduction_handle, interp_err, global_value, CCTK_VARIABLE_REAL)
    if (retval .ne. 0 ) then
       call CCTK_WARN(0,"Error in obtaining outer evolution point radius")
    endif

    if (global_value .lt. mask_testinterp_tolerance) return;

    tmp_maskn = abs(tmp_maskn)+abs(tmp_masks);

    minimum_r_guard = dsqrt(minval( qs**2+ps**2, guard_mask.eq.1 .and. tmp_maskn.lt.mask_testinterp_tolerance ))

    call CCTK_ReductionArrayHandle(reduction_handle, "minimum");
    if (reduction_handle .lt. 0 ) then
       call CCTK_WARN(0,"Could not get reduction handle")
    end if

    call CCTK_ReduceLocScalar(retval, cctkGH, -1, reduction_handle, minimum_r_guard, global_value, CCTK_VARIABLE_REAL)
    if (retval .ne. 0 ) then
       call CCTK_WARN(0,"Error in obtaining outer evolution point radius")
    endif

    i = int((global_value - 1)/null_delta(1))+1

    write (message,*) 'guard point shell is too close to equator.  Try setting NullGrid::N_ang_ev_outside_eq >=', i;
    call CCTK_WARN(0, trim(message))

  end subroutine test_interpolation

end subroutine NullInterp_MaskInit
