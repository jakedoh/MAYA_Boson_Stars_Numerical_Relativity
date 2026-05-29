! Check to see if the re-initialization has to be undone.
! $Header$
!option opt(o(0)))

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_ReInitialize_Check(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  CCTK_REAL :: fmin_extremum, fmin_extremum_loc, fmax_loc
  character(len=128) :: info_message
  logical :: check_re_init

! Initializa flags.
  check_re_init = .false.
  re_initialize_undone = .false.

! If re-initialization has just be done we need to check and see if it
! needs to be undone.
  if ( re_initialize_every .gt. 0 ) then
    if ( mod(cctk_iteration,re_initialize_every) .eq. 0 ) then
      check_re_init = .true.
    end if 
  end if 
    
#include "include/physical_part2.h"

  if ( check_re_init ) then

!   Loop over the different level sets.
    do l = 1, eh_number_level_sets
    
!     Initialize the variables holding the processor local maximum
!     values in the extremum.
      fmax_loc = zero
  
!     On the local processor find the maximum value of f and the maximum value
!     of f in any local extrema where f is negative.
      fmin_extremum_loc = ex_value
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
            fmax_loc = max ( f(i,j,k,l), fmax_loc )
            if ( f(i,j,k,l) .lt. zero ) then
              if ( ( (f(i,j,k,l) - f(i-1,j,k,l)) * &
                     (f(i+1,j,k,l) - f(i,j,k,l)) .le. zero ) .and. &
                   ( (f(i,j,k,l) - f(i,j-1,k,l)) * &
                     (f(i,j+1,k,l) - f(i,j,k,l)) .le. zero ) .and. &
                   ( (f(i,j,k,l) - f(i,j,k-1,l)) * &
                     (f(i,j,k+1,l) - f(i,j,k,l)) .le. zero ) ) then
                fmin_extremum_loc = max ( f(i,j,k,l), fmin_extremum_loc )
              end if
            end if
          end do
        end do
      end do
  
!     Maximum reduction over all processors.
      call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
                         fmin_extremum_loc, fmin_extremum, CCTK_VARIABLE_REAL )
  
      if ( ierr .ne. 0 ) then
        call CCTK_WARN(0,'Reduction of fmax_extremum failed')
      end if
  
!     Write an info message about the local extremum.
      write(info_message,'(a39,i3,a6,f12.5)') &
          'Minimum f at the extrema for level set ', l, ' is : ',fmin_extremum
      call CCTK_INFO(trim(info_message))
  
!     If the local extremum indicate and imminent pinch off undo the
!     re-initialization if this is requested by the user.
      if ( re_init_undo .gt. 0 ) then
        if ( fmin_extremum .gt. -two * min ( CCTK_DELTA_SPACE(1), &
                                             CCTK_DELTA_SPACE(2), &
                                             CCTK_DELTA_SPACE(3) ) ) then
          write(info_message,'(a45,i3,a33)') &
               'Re-initialization of level set ', l, &
               ' undone due to imminent pinch-off'
          call CCTK_INFO(trim(info_message))
          f(:,:,:,l) = fbak(:,:,:,l)
          eh_mask(:,:,:,l) = eh_mask_bak(:,:,:,l)
          re_initialize_undone(l) = .true.
        end if
      end if
    end do
  end if

end subroutine EHFinder_ReInitialize_Check
