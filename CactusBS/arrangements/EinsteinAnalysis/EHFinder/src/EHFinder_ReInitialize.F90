! Re-Initialization of the level set function
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! Control routine for re-initialization of the level set function.

subroutine EHFinder_ReInitializeControl(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
 
  CCTK_REAL, dimension(eh_number_level_sets) :: fmin, fmax
  CCTK_INT :: i, j, l !, re_init_control
  character(len=128) :: info_message, variable_name
  character(11) :: form(4) = (/ '(a12,i1,a1)', '(a12,i2,a1)', &
                                '(a12,i3,a1)', '(a12,i4,a1)' /)

! Initialize the control variables. The default means no re-initialization.
  re_init_control = 0

! Check if it is time for a pde re-initialization and set the control
! variable to 1.
  if ( re_initialize_every > 0 ) then
    if ( mod(cctk_iteration,re_initialize_every) == 0 ) then
      re_init_control = 1
    end if 
  end if 
    
  if ( re_init_control == 1 ) then

  ! Get reduction handles for minimum and maximum reductions.
    call CCTK_ReductionHandle ( max_handle, 'maximum' )
    if ( max_handle < 0 ) then
      call CCTK_WARN(0,'Could not obtain a handle for maximum reduction')
    end if
    call CCTK_ReductionHandle ( min_handle, 'minimum' )
    if ( min_handle < 0 ) then
      call CCTK_WARN(0,'Could not obtain a handle for minimum reduction')
    end if

    do l = 1, eh_number_level_sets
      select case (l)
      case (1:9)
        j = 1
      case(10:99)
        j = 2
      case(100:999)
        j = 3
      case(1000:9999)
        j = 4
      case(10000:)
        call CCTK_WARN(0,'To many simultaneous level set functions')
      end select
      write (variable_name,form(j)) 'ehfinder::f[',l-1,']'
      call CCTK_VarIndex ( i, variable_name )
      if ( i < 0 ) then
        write (info_message,'(a9,a14,a16)') 'Variable ',variable_name, &
                                            ' does not exist'
        call CCTK_WARN(0,info_message)
      end if
      call CCTK_Reduce ( ierr, cctkGH, -1, max_handle, 1, &
                         CCTK_VARIABLE_REAL, fmax(l), 1, i) 
      if ( ierr .ne. 0 ) then
        call CCTK_WARN(0,'Reduction of fmax failed')
      end if
      call CCTK_Reduce ( ierr, cctkGH, -1, min_handle, 1, &
                         CCTK_VARIABLE_REAL, fmin(l), 1, i) 
      if ( ierr .ne. 0 ) then
        call CCTK_WARN(0,'Reduction of fmin failed')
      end if
    end do

    re_init_this_level_set = .true.

  ! If fmin and fmax have the same sign, there is no surface present and
  ! re-initialization will not be performed for the given surface.
    do l = 1, eh_number_level_sets
      if ( fmin(l) * fmax(l) > zero ) then
        re_init_this_level_set(l) = .false.
      end if
    end do

  ! If none of the surfaces should be re-initialized set the control variables
  ! to zero.
    if ( all ( .not. re_init_this_level_set ) ) then
      re_init_control = 0
      info_message = 'No zero-points of the level set functions.'
      info_message = trim(info_message)//' No re-initialization performed.'
      call CCTK_INFO ( trim(info_message) )
    end if

  end if

  if ( re_init_control == 1 ) then

!   If we are called more than once at the same iteration number (for example
!   during the 3-level data initialization routine with Carpet) only do
!   something the first time.
    if ( cctk_iteration == last_called_at_iteration ) then
      re_init_control = 0
    else
      last_called_at_iteration = cctk_iteration
    end if
  end if 

  if ( re_init_control == 1 ) then

  ! Set the courant factor for the 'Euler' re-initialization scheme.
    if ( CCTK_EQUALS ( re_init_int_method, 'Euler' ) ) then
      hfac = one / four
    end if

  ! Set the courant factor for the 'rk2' re-initialization scheme.
    if ( CCTK_EQUALS ( re_init_int_method, 'rk2' ) ) then
      hfac = half
    end if

    niter_reinit = 0

    call CCTK_INFO ('Re-Initialization started')
  end if

end subroutine EHFinder_ReInitializeControl

subroutine EHFinder_ReInitializeInitialize(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS


  if ( re_init_control == 1 ) then

  ! Set up counter for the number of iterations. Copy f into ftmp and
  ! fbak and eh_mask into eh_mask_bak (in case the re-initialization has
  ! to be undone) and find the min and max values of f.
    ftmp = f
    fbak = f
    eh_mask_bak = eh_mask

  end if
end subroutine EHFinder_ReInitializeInitialize

subroutine EHFinder_ReInitializeRK2_1(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  CCTK_REAL :: idx, idy, idz
  CCTK_REAL :: al, ar, bl, br, cl, cr
  CCTK_REAL :: alminus, alplus, blminus, blplus, clminus, clplus
  CCTK_REAL :: arminus, arplus, brminus, brplus, crminus, crplus

! Find the physical part of the computational domain.
#include "include/physical_part.h"

! Set the step size from the courant factor and the minimum of the
! grid spacing.
  h = hfac*min(CCTK_DELTA_SPACE(1),CCTK_DELTA_SPACE(2),CCTK_DELTA_SPACE(3))

! If centered finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'centered' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = half / CCTK_DELTA_SPACE(1)
    idy = half / CCTK_DELTA_SPACE(2)
    idz = half / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/centered_second2.h"
          end do
        end do
      end do
    end do
  end if

! If first order upwinded finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'upwind' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = one / CCTK_DELTA_SPACE(1)
    idy = one / CCTK_DELTA_SPACE(2)
    idz = one / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
# include "include/upwind_first.h"
          end do
        end do
      end do
    end do
  end if

! If second order upwinded finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'upwind2' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = half / CCTK_DELTA_SPACE(1)
    idy = half / CCTK_DELTA_SPACE(2)
    idz = half / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/upwind_second2.h"
          end do
        end do
      end do
    end do
  end if

! For all active grid points take half an euler step.
  where ( eh_mask >= 0 )
    sftmp =  sqrt( dfx**2 + dfy**2 + dfz**2 )
    sftmp = - half * h * f / sqrt(f**2+one) * ( sftmp - one )
    f = f + sftmp
  elsewhere
    sftmp = zero
  end where

end subroutine EHFinder_ReInitializeRK2_1


subroutine EHFinder_ReInitializeRK2_2(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  CCTK_REAL :: idx, idy, idz
  CCTK_REAL, dimension(eh_number_level_sets) :: maxdfloc, maxdf, &
                                                mindfloc, mindf
  CCTK_REAL :: al, ar, bl, br, cl, cr
  CCTK_REAL :: alminus, alplus, blminus, blplus, clminus, clplus
  CCTK_REAL :: arminus, arplus, brminus, brplus, crminus, crplus
  character(len=128) :: info_message

#include "include/physical_part.h"

! If centered finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'centered' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = half / CCTK_DELTA_SPACE(1)
    idy = half / CCTK_DELTA_SPACE(2)
    idz = half / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/centered_second2.h"
          end do
        end do
      end do
    end do
  end if

! If first order upwinded finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'upwind' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = one / CCTK_DELTA_SPACE(1)
    idy = one / CCTK_DELTA_SPACE(2)
    idz = one / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/upwind_first.h"
          end do
        end do
      end do
    end do
  end if

! If second order upwinded finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'upwind2' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = half / CCTK_DELTA_SPACE(1)
    idy = half / CCTK_DELTA_SPACE(2)
    idz = half / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/upwind_second2.h"
          end do
        end do
      end do
    end do
  end if

! For all active grid points take a full euler step using the right hand
! side calculated at the half step.
  where ( eh_mask >= 0 )
    sftmp =  sqrt( dfx**2 + dfy**2 + dfz**2 )
    sftmp = - h * f / sqrt(f**2+one) * ( sftmp - one )
    f = ftmp + sftmp
  elsewhere
    sftmp = zero
  end where

! Find the maximum and minimum of the changes to the level set functions
! done in this step for all active grid points on the local processor
  do l = 1, eh_number_level_sets
    maxdfloc(l) = maxval ( abs(sftmp(ixl:ixr,jyl:jyr,kzl:kzr,l)), &
                             mask = eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) >= 0)
    mindfloc(l) = minval ( abs(sftmp(ixl:ixr,jyl:jyr,kzl:kzr,l)), &
                             mask = eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) >= 0)
  end do

! Find the maximum of the changes to the level set functions done on all
! processors
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, max_handle, &
                                    maxdfloc, maxdf, eh_number_level_sets, &
                                    CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0,'Reduction of maxdf failed')
  end if

! Find the minimum of the changes to the level set functions done on all
! processors
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, min_handle, &
                                    mindfloc, mindf, eh_number_level_sets, &
                                    CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0,'Reduction of mindf failed')
  end if

! If the maximum change on all processors are small enough signal that
! we are done with re-initialization by setting re_init_control = 0.
  if ( all ( maxdf < h*min(CCTK_DELTA_SPACE(1),CCTK_DELTA_SPACE(2), &
                           CCTK_DELTA_SPACE(3))**2 ) ) then
!    pugh_re_init_control = 0
    re_init_control = 0
  end if

! Update the number of calls counter.
  niter_reinit = niter_reinit + 1

! If we are finished write an info message.
!  if ( pugh_re_init_control == 0 ) then
  if ( re_init_control == 0 ) then
    write(info_message,'(a35,i5,a12)') &
        'PDE re-initialization complete in ',niter_reinit,' iterations.'
    call CCTK_INFO( trim(info_message) )
  end if

! If we ran out of iterations without converging, print a level 1 warning
! message.
  if ( niter_reinit > re_init_max_iter ) then
    call CCTK_WARN(1,'Re-initialization failed to converge')
!    pugh_re_init_control = 0
    re_init_control = 0
  end if

! Copy the new values for the level set function into the temporary variable,
! so we are ready for the next step.
  ftmp = f

end subroutine EHFinder_ReInitializeRK2_2


subroutine EHFinder_ReInitializeEuler(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  CCTK_REAL :: idx, idy, idz
  CCTK_REAL :: al, ar, bl, br, cl, cr
  CCTK_REAL :: alminus, alplus, blminus, blplus, clminus, clplus
  CCTK_REAL :: arminus, arplus, brminus, brplus, crminus, crplus

! Find the physical part of the computational domain.
#include "include/physical_part.h"

! Set the step size from the courant factor and the minimum of the
! grid spacing.
  h = hfac*min(CCTK_DELTA_SPACE(1),CCTK_DELTA_SPACE(2),CCTK_DELTA_SPACE(3))

! If centered finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'centered' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = half / CCTK_DELTA_SPACE(1)
    idy = half / CCTK_DELTA_SPACE(2)
    idz = half / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/centered_second2.h"
          end do
        end do
      end do
    end do
  end if

! If first order upwinded finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'upwind' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = one / CCTK_DELTA_SPACE(1)
    idy = one / CCTK_DELTA_SPACE(2)
    idz = one / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/upwind_first.h"
          end do
        end do
      end do
    end do
  end if

! If second order upwinded finite differences are required...
  if ( CCTK_EQUALS ( pde_differences, 'upwind2' ) ) then

!   Caclulate the factor used in the finite differences used in the
!   include file.
    idx = half / CCTK_DELTA_SPACE(1)
    idy = half / CCTK_DELTA_SPACE(2)
    idz = half / CCTK_DELTA_SPACE(3)

!   For all level sets and all physical points...
    do l = 1, eh_number_level_sets
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
# include "include/upwind_second2.h"
          end do
        end do
      end do
    end do
  end if

! For all active grid points take a full euler step using.
  where ( eh_mask >= 0 )
    sftmp =  sqrt( dfx**2 + dfy**2 + dfz**2 )
    sftmp =  - h * f / sqrt(f**2+one) * ( sftmp - one )
    f = f + sftmp
  elsewhere
    sftmp = zero
  end where
  
! Store the absulute value of the change to be used in the exit criteria,
! scaled by an appropriate factor.

  sftmp = abs ( sftmp ) / ( h*min(CCTK_DELTA_SPACE(1),CCTK_DELTA_SPACE(2), &
                                  CCTK_DELTA_SPACE(3))**2 )

end subroutine EHFinder_ReInitializeEuler

subroutine EHFinder_ReInitializePostStep(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, dimension(eh_number_level_sets) :: maxdf
  CCTK_INT :: i,j,l
  character(len=128) :: info_message, variable_name
  character(15) :: form(4) = (/ '(a16,i1,a1)', '(a16,i2,a1)', &
                                '(a16,i3,a1)', '(a16,i4,a1)' /)

! Get reduction handles for maximum reduction.
  call CCTK_ReductionHandle ( max_handle, 'maximum' )
  if ( max_handle < 0 ) then
    call CCTK_WARN(0,'Could not obtain a handle for maximum reduction')
  end if

  do l = 1, eh_number_level_sets
    select case (l)
    case (1:9)
      j = 1
    case(10:99)
      j = 2
    case(100:999)
      j = 3
    case(1000:9999)
      j = 4
    case(10000:)
      call CCTK_WARN(0,'To many simultaneous level set functions')
    end select

    write (variable_name,form(j)) 'ehfinder::sftmp[',l-1,']'

    call CCTK_VarIndex ( i, variable_name )

    if ( i < 0 ) then
      write (info_message,'(a9,a18,a16)') 'Variable ',variable_name, &
                                          ' does not exist'
      call CCTK_WARN(0,info_message)
    end if

    call CCTK_Reduce ( ierr, cctkGH, -1, max_handle, 1, &
                       CCTK_VARIABLE_REAL, maxdf(l), 1, i)

    if ( ierr .ne. 0 ) then
      call CCTK_WARN(0,'Reduction of fmax failed')
    end if

  end do

  if ( re_init_verbose /= 0 ) then
    write(info_message,'(a13,i5,a9,es12.5)') 'At iteration ', niter_reinit, &
                                            ' maxdf = ', maxval(maxdf)
    call CCTK_INFO ( trim(info_message) )
  end if

! If the maximum change on all processors are small enough signal that
! we are done with re-initialization by setting *_re_init_control = 0.
  if ( all ( maxdf < one ) ) then
      re_init_control = 0
  endif

  if ( CCTK_IsThornACtive ( 'PUGH' ) /= 0 ) then
!   Update the number of calls counter.
    niter_reinit = niter_reinit + 1
  end if

! If we are finished write an info message.
!  if ( pugh_re_init_control == 0  .and. carpet_re_init_control == 0 ) then
  if ( re_init_control == 0 ) then
    write(info_message,'(a35,i5,a12)') &
        'PDE re-initialization complete in ',niter_reinit,' iterations.'
    call CCTK_INFO( trim(info_message) )
  end if

! If we ran out of iterations without converging, print a level 1 warning
! message.
  if ( niter_reinit > re_init_max_iter ) then
    call CCTK_WARN(1,'Re-initialization failed to converge')
    re_init_control = 0
  end if

end subroutine EHFinder_ReInitializePostStep
