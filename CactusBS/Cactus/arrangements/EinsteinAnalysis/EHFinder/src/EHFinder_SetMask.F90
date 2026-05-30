! Mask modification functions.
! The grid function eh_mask is used to encode excision and boundary
! information. An excised point has the value -1, while an active point
! is 0 or positive. If it is zero the point has neighbours in all 
! directions while if it is positive the value encodes which neighbours are
! missing. The mask uses 6 bits of an integer to encode this information.
! If the neighbour in the -x-direction is missing the mask is set 
! to 1 (b000001). If the neighbour in the +x-direction is missing the mask
! is set to 2 (b000010). For the y-directions the values are 4 (b000100) and 
! 8 (b001000). For the z-direction they are 16 (b010000) and 32 (b100000).
! For example a point with no active neighbours in the +x, -y and +z direction
! will have a mask value of 38 (b100110).
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

! This routine is called only once to initialise the mask at the physical
! outer boundaries. The value of the mask in these points should never
! be changed again.
subroutine EHFinder_MaskInit(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

! Initialise the mask to 0.
  eh_mask = 0

! Figure out the local dimensions of the grid excluding ghostzones and
! symmetry zones.
# include "include/physical_part.h"

! Check if the point is located at a physical outer boundary and set the
! eh_mask appropriately. The array ll is defined in EHFinder_mod.F90
! and contains the values ( 1, 2, 4, 8, 16, 32 ).
  if ( ( cctk_bbox(1) .eq. 1 ) .and. ( ixl .eq. 1 ) ) then
    eh_mask(1,:,:,:) = eh_mask(1,:,:,:) + ll(0)
  end if
  if ( ( cctk_bbox(2) .eq. 1 ) .and. ( ixr .eq. nx ) ) then
    eh_mask(nx,:,:,:) = eh_mask(nx,:,:,:) + ll(1)
  end if
  if ( ( cctk_bbox(3) .eq. 1 ) .and. ( jyl .eq. 1 ) ) then
    eh_mask(:,1,:,:) = eh_mask(:,1,:,:) + ll(2)
  end if
  if ( ( cctk_bbox(4) .eq. 1 ) .and. ( jyr .eq. ny ) ) then
    eh_mask(:,ny,:,:) = eh_mask(:,ny,:,:) + ll(3)
  end if
  if ( ( cctk_bbox(5) .eq. 1 ) .and. ( kzl .eq. 1 ) ) then
    eh_mask(:,:,1,:) = eh_mask(:,:,1,:) + ll(4)
  end if
  if ( ( cctk_bbox(6) .eq. 1 ) .and. ( kzr .eq. nz ) ) then
    eh_mask(:,:,nz,:) = eh_mask(:,:,nz,:) + ll(5)
  end if

  return
end subroutine EHFinder_MaskInit


! This routine will excise points and re-activate excised points
! as needed. EHFinder_Setmask2 will then find the boundary of the
! excsion region and make sure that the mask is correct there.
subroutine EHFinder_SetMask1(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l, i1, i2, j1, j2, k1, k2
  logical :: active

  CCTK_INT, dimension(3) :: imin_loc, imax_loc, imin_n, imax_n
  CCTK_REAL, dimension(3) :: fimin_loc, fimax_loc

  active = .false.

! If the mask has not been set before, we should do it now so 
! we set active=.true.
! mask_first is initialized to .true. in EHFinder_mod.F90
  if ( mask_first ) active = .true.

! If re-parametrization has just been done, we should reset the mask so we
! set active=.true.
  if ( re_initialize_every .gt. 0 ) then
    if ( mod(cctk_iteration,re_initialize_every) .eq. 0 ) active = .true.
  end if

! If the reparametrization was not undone...
  if ( active .and. .not. all(re_initialize_undone) ) then

!   Get the minimum and maximum index excluding ghost and symmetry cells.
# include "include/physical_part.h"

    loop_over_l: do l = 1, eh_number_level_sets

!     Store the current mask and level set function
      tm_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) = eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l)
      ftmp(ixl:ixr,jyl:jyr,kzl:kzr,l) = f(ixl:ixr,jyl:jyr,kzl:kzr,l)

      not_undone: if ( .not. re_initialize_undone(l) ) then

!     Next we try to locate the minimum and maximum of the global indeces of
!     the currently active cells giving us the smallest rectangular box
!     containing all active cells.

        if ( use_outer_excision .gt. 0 ) then
!         First initialize some variables.
!         fimin_loc, and fimax_loc are 3 element arrays that will contain the 
!         minimum value of f on the boundaries of this rectangular box. They
!         will be used to decide if the box should be changed. They are
!         initialized to the negative of the value used in excised cells.

!         The 3 element arrays imin_loc and imax_loc will contain the min and
!         max global indices of the ractangular box. They are initialised to
!         the maximum global index and 0, respectively.

          fimin_loc = -ex_value; fimax_loc = -ex_value
          imin_loc = cctk_gsh; imax_loc = 0

!         Find all cells with f<shell_width*delta and find their minimum
!         and maximum global cell index.
          do k = kzl, kzr
            do j = jyl, jyr
              do i = ixl, ixr
                if ( eh_mask(i,j,k,l) .ge. 0 ) then
                  if ( f(i,j,k,l) .lt. shell_width * delta ) then
                    if ( i+cctk_lbnd(1) .le. imin_loc(1) ) then
                      imin_loc(1) = i+cctk_lbnd(1)
                    end if
                    if ( i+cctk_lbnd(1) .ge. imax_loc(1) ) then
                      imax_loc(1) = i+cctk_lbnd(1)
                    end if
                    if ( j+cctk_lbnd(2) .le. imin_loc(2) ) then
                      imin_loc(2) = j+cctk_lbnd(2)
                    end if
                    if ( j+cctk_lbnd(2) .ge. imax_loc(2) ) then
                      imax_loc(2) = j+cctk_lbnd(2)
                    end if
                    if ( k+cctk_lbnd(3) .le. imin_loc(3) ) then
                      imin_loc(3) = k+cctk_lbnd(3)
                    end if
                    if ( k+cctk_lbnd(3) .ge. imax_loc(3) ) then
                      imax_loc(3) = k+cctk_lbnd(3)
                    end if
                  end if
                end if
              end do
            end do
          end do

!         Reduce over all processors to get the global indeces of the
!         rectangular box containing all cells with f<shell_width*delta
          call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, min_handle, &
                   CCTK_PointerTo( imin_loc ), CCTK_PointerTo( imin_glob ), &
                   3, CCTK_VARIABLE_INT )
          if ( ierr .ne. zero ) then
            call CCTK_WARN ( 0, 'Reduction of array imin_loc failed' )
          end if

          call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, max_handle, &
                   CCTK_PointerTo( imax_loc ), CCTK_PointerTo( imax_glob ), &
                   3, CCTK_VARIABLE_INT )
          if ( ierr .ne. zero ) then
            call CCTK_WARN ( 0, 'Reduction of array imax_loc failed' )
          end if

!         Convert into local indeces. Note that these might be less than 1 or
!         larger than the local grid size if the box does not overlap with the
!         local grid.
          i1 = imin_glob(1)-cctk_lbnd(1)
          i2 = imax_glob(1)-cctk_lbnd(1)
          j1 = imin_glob(2)-cctk_lbnd(2)
          j2 = imax_glob(2)-cctk_lbnd(2)
          k1 = imin_glob(3)-cctk_lbnd(3)
          k2 = imax_glob(3)-cctk_lbnd(3)

!         Find the minimum value of f on the various faces of the rectangular
!         box if part of the face is present on the current grid.
          if ( ( 1 .le. i1 ) .and. ( i1 .le. nx ) ) then
            fimin_loc(1) = minval ( f(i1,jyl:jyr,kzl:kzr,l) )
          end if
          if ( ( 1 .le. i2 ) .and. ( i2 .le. nx ) ) then
            fimax_loc(1) = minval ( f(i2,jyl:jyr,kzl:kzr,l) )
          end if
          if ( ( 1 .le. j1 ) .and. ( j1 .le. ny ) ) then
            fimin_loc(2) = minval ( f(ixl:ixr,j1,kzl:kzr,l) )
          end if
          if ( ( 1 .le. j2 ) .and. ( j2 .le. ny ) ) then
            fimax_loc(2) = minval ( f(ixl:ixr,j2,kzl:kzr,l) )
          end if
          if ( ( 1 .le. k1 ) .and. ( k1 .le. nz ) ) then
            fimin_loc(3) = minval ( f(ixl:ixr,jyl:jyr,k1,l) )
          end if
          if ( ( 1 .le. k2 ) .and. ( k2 .le. nz ) ) then
            fimax_loc(3) = minval ( f(ixl:ixr,jyl:jyr,k2,l) )
          end if

!         Reduce to get the minimum values of f on the faces of the box
          call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, min_handle, &
                 CCTK_PointerTo( fimin_loc ), CCTK_PointerTo( fimin_glob ), &
                 3, CCTK_VARIABLE_REAL )
          if ( ierr .ne. zero ) then
            call CCTK_WARN ( 0, 'Reduction of array fimin_loc failed' )
          end if

          call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, min_handle, &
                 CCTK_PointerTo( fimax_loc ), CCTK_PointerTo( fimax_glob ), &
                 3, CCTK_VARIABLE_REAL )
          if ( ierr .ne. zero ) then
            call CCTK_WARN ( 0, 'Reduction of array fimax_loc failed' )
          end if

        end if
           
!       Now check and see if any interior excised points need to be activated.
        if ( use_inner_excision .gt. 0 ) then
          do k = kzl, kzr
            do j = jyl, jyr
              do i = ixl, ixr

!               If an interior excised point has a non-excised neighbour where
!               f-delta>-shell_width*delta then activate the point by setting
!               the temporary mask to zero. The new value of f will be the value
!               if f in its neighbour point - delta. Do this for all directions.

                if ( ( eh_mask(i,j,k,l) .eq. -1 ) .and. &
                     ( f(i,j,k,l) .lt. 0 ) ) then
                  if ( eh_mask(i-1,j,k,l) .ge. 0 ) then
                    if ( f(i-1,j,k,l) - delta .gt. -shell_width * delta ) then
                      tm_mask(i,j,k,l) = 0
                      ftmp(i,j,k,l) = f(i-1,j,k,l) - delta
                    end if
                  end if

                  if ( eh_mask(i+1,j,k,l) .ge. 0 ) then
                    if ( f(i+1,j,k,l) - delta .gt. -shell_width * delta ) then
                      tm_mask(i,j,k,l) = 0
                      ftmp(i,j,k,l) = f(i+1,j,k,l) - delta
                    end if
                  end if

                  if ( eh_mask(i,j-1,k,l) .ge. 0 ) then
                    if ( f(i,j-1,k,l) - delta .gt. -shell_width * delta ) then
                      tm_mask(i,j,k,l) = 0
                      ftmp(i,j,k,l) = f(i,j-1,k,l) - delta
                    end if
                  end if

                  if ( eh_mask(i,j+1,k,l) .ge. 0 ) then
                    if ( f(i,j+1,k,l) - delta .gt. -shell_width * delta ) then
                      tm_mask(i,j,k,l) = 0
                      ftmp(i,j,k,l) = f(i,j+1,k,l) - delta
                    end if
                  end if

                  if ( eh_mask(i,j,k-1,l) .ge. 0 ) then
                    if ( f(i,j,k-1,l) - delta .gt. -shell_width * delta ) then
                      tm_mask(i,j,k,l) = 0
                      ftmp(i,j,k,l) = f(i,j,k-1,l) - delta
                    end if
                  end if

                  if ( eh_mask(i,j,k+1,l) .ge. 0 ) then
                    if ( f(i,j,k+1,l) - delta .gt. -shell_width * delta ) then
                      tm_mask(i,j,k,l) = 0
                      ftmp(i,j,k,l) = f(i,j,k+1,l) - delta
                    end if
                  end if

                end if
              end do
            end do
          end do
        end if

!       Check and see if the boundary of the exterior excision region should
!       be changed and if so find the indices describing the new excision
!       region.
        if ( use_outer_excision .gt. 0 ) then
          do i = 1, 3
            imin_n(i) = imin_glob(i)
            imax_n(i) = imax_glob(i)
!           If the minimum value of f on a face on the box plus delta is less
!           than shell_width * delta then the active region has to be increased
!           in size in the corresponding region, i.e. inactive cells has to be
!           activated.
            if ( fimin_glob(i) + delta .lt. shell_width * delta ) then
              if ( imin_glob(i) .gt. 1 ) then
                imin_n(i) = imin_glob(i) - 1
              endif
            end if
            if ( fimax_glob(i) + delta .lt. shell_width * delta ) then
              if ( imax_glob(i) .lt. cctk_gsh(i) ) then
                imax_n(i) = imax_glob(i) + 1
              endif
            end if
          end do

!         Use the new indices to actually activate points, taking care to
!         activate points that are actually on the local grid.  First do the
!         faces of the rectangular box.
          if ( imin_n(1) .ne. imin_glob(1) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) ) then
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( ( j1 .le. j2 ) .and. ( k1 .le. k2 ) ) then
                tm_mask(i1,j1:j2,k1:k2,l) = 0
                ftmp(i1,j1:j2,k1:k2,l) = ftmp(i1+1,j1:j2,k1:k2,l) + delta
              end if
            end if
          end if
          if ( imax_n(1) .ne. imax_glob(1) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            if ( ( ixl .lt. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) ) then
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( ( j1 .le. j2 ) .and. ( k1 .le. k2 ) ) then
                tm_mask(i2,j1:j2,k1:k2,l) = 0
                ftmp(i2,j1:j2,k1:k2,l) = ftmp(i2-1,j1:j2,k1:k2,l) + delta
              end if
            end if
          end if
          if ( imin_n(2) .ne. imin_glob(2) ) then
            j1 = imin_n(2) - cctk_lbnd(2)
            if ( ( max(jyl,2) .le. j1 ) .and. ( j1 .lt. jyr ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( ( i1 .le. i2 ) .and. ( k1 .le. k2 ) ) then
                tm_mask(i1:i2,j1,k1:k2,l) = 0
                ftmp(i1:i2,j1,k1:k2,l) = ftmp(i1:i2,j1+1,k1:k2,l) + delta
              end if
            end if
          end if
          if ( imax_n(2) .ne. imax_glob(2) ) then
            j2 = imax_n(2) - cctk_lbnd(2)
            if ( ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( ( i1 .le. i2 ) .and. ( k1 .le. k2 ) ) then
                tm_mask(i1:i2,j2,k1:k2,l) = 0
                ftmp(i1:i2,j2,k1:k2,l) = f(i1:i2,j2-1,k1:k2,l) + delta
              end if
            end if
          end if
          if ( imin_n(3) .ne. imin_glob(3) ) then
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              if ( ( i1 .le. i2 ) .and. ( j1 .le. j2 ) ) then
                tm_mask(i1:i2,j1:j2,k1,l) = 0
                ftmp(i1:i2,j1:j2,k1,l) = ftmp(i1:i2,j1:j2,k1+1,l) + delta
              end if
            end if
          end if
          if ( imax_n(3) .ne. imax_glob(3) ) then
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( kzl .lt. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              if ( ( i1 .le. i2 ) .and. ( j1 .le. j2 ) ) then
                tm_mask(i1:i2,j1:j2,k2,l) = 0
                ftmp(i1:i2,j1:j2,k2,l) = ftmp(i1:i2,j1:j2,k2-1,l) + delta
              end if
            end if
          end if

!         Then do the edges of the box.
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imin_n(2) .ne. imin_glob(2) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            j1 = imin_n(2) - cctk_lbnd(2)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) ) then
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( k1 .le. k2 ) then
                tm_mask(i1,j1,k1:k2,l) = 0
                ftmp(i1,j1,k1:k2,l) = f(i1+1,j1+1,k1:k2,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imax_n(2) .ne. imax_glob(2) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            j2 = imax_n(2) - cctk_lbnd(2)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) ) then
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( k1 .le. k2 ) then
                tm_mask(i1,j2,k1:k2,l) = 0
                ftmp(i1,j2,k1:k2,l) = f(i1+1,j2-1,k1:k2,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              if ( j1 .le. j2 ) then
                tm_mask(i1,j1:j2,k1,l) = 0
                ftmp(i1,j1:j2,k1,l) = f(i1+1,j1:j2,k1+1,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              if ( j1 .le. j2 ) then
                tm_mask(i1,j1:j2,k2,l) = 0
                ftmp(i1,j1:j2,k2,l) = f(i1+1,j1:j2,k2-1,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imin_n(2) .ne. imin_glob(2) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            j1 = imin_n(2) - cctk_lbnd(2)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) ) then
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( k1 .le. k2 ) then
                tm_mask(i2,j1,k1:k2,l) = 0
                ftmp(i2,j1,k1:k2,l) = f(i2-1,j1+1,k1:k2,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imax_n(2) .ne. imax_glob(2) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            j2 = imax_n(2) - cctk_lbnd(2)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) ) then
              k1 = max ( max ( imin_n(3), imin_glob(3) ) - cctk_lbnd(3), kzl )
              k2 = min ( min ( imax_n(3), imax_glob(3) ) - cctk_lbnd(3), kzr )
              if ( k1 .le. k2 ) then
                tm_mask(i2,j2,k1:k2,l) = 0
                ftmp(i2,j2,k1:k2,l) = f(i2-1,j2-1,k1:k2,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              if ( j1 .le. j2 ) then
                tm_mask(i2,j1:j2,k1,l) = 0
                ftmp(i2,j1:j2,k1,l) = f(i2-1,j1:j2,k1+1,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              j1 = max ( max ( imin_n(2), imin_glob(2) ) - cctk_lbnd(2), jyl )
              j2 = min ( min ( imax_n(2), imax_glob(2) ) - cctk_lbnd(2), jyr )
              if ( j1 .le. j2 ) then
                tm_mask(i2,j1:j2,k2,l) = 0
                ftmp(i2,j1:j2,k2,l) = f(i2-1,j1:j2,k2-1,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imin_n(2) .ne. imin_glob(2) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            j1 = imin_n(2) - cctk_lbnd(2)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              if ( i1 .le. i2 ) then
                tm_mask(i1:i2,j1,k1,l) = 0
                ftmp(i1:i2,j1,k1,l) = f(i1:i2,j1+1,k1+1,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imin_n(2) .ne. imin_glob(2) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            j1 = imin_n(2) - cctk_lbnd(2)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              if ( i1 .le. i2 ) then
                tm_mask(i1:i2,j1,k2,l) = 0
                ftmp(i1:i2,j1,k2,l) = f(i1:i2,j1+1,k2-1,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imax_n(2) .ne. imax_glob(2) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            j2 = imax_n(2) - cctk_lbnd(2)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              if ( i1 .le. i2 ) then
                tm_mask(i1:i2,j2,k1,l) = 0
                ftmp(i1:i2,j2,k1,l) = f(i1:i2,j2-1,k1+1,l) + sqrt(two)*delta
              end if
            end if
          end if
          if ( ( imax_n(2) .ne. imax_glob(2) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            j2 = imax_n(2) - cctk_lbnd(2)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              i1 = max ( max ( imin_n(1), imin_glob(1) ) - cctk_lbnd(1), ixl )
              i2 = min ( min ( imax_n(1), imax_glob(1) ) - cctk_lbnd(1), ixr )
              if ( i1 .le. i2 ) then
                tm_mask(i1:i2,j2,k2,l) = 0
                ftmp(i1:i2,j2,k2,l) = f(i1:i2,j2-1,k2-1,l) + sqrt(two)*delta
              end if
            end if
          end if

!         And finally do the corners.
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imin_n(2) .ne. imin_glob(2) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            j1 = imin_n(2) - cctk_lbnd(2)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              tm_mask(i1,j1,k1,l) = 0
              ftmp(i1,j1,k1,l) = f(i1+1,j1+1,k1+1,l) + sqrt(three)*delta
            end if
          end if
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imin_n(2) .ne. imin_glob(2) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            j1 = imin_n(2) - cctk_lbnd(2)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              tm_mask(i1,j1,k2,l) = 0
              ftmp(i1,j1,k2,l) = f(i1+1,j1+1,k2-1,l) + sqrt(three)*delta
            end if
          end if
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imax_n(2) .ne. imax_glob(2) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            j2 = imax_n(2) - cctk_lbnd(2)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              tm_mask(i1,j2,k1,l) = 0
              ftmp(i1,j2,k1,l) = f(i1+1,j2-1,k1+1,l) + sqrt(three)*delta
            end if
          end if
          if ( ( imin_n(1) .ne. imin_glob(1) ) .and. &
               ( imax_n(2) .ne. imax_glob(2) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            i1 = imin_n(1) - cctk_lbnd(1)
            j2 = imax_n(2) - cctk_lbnd(2)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( max(ixl,2) .le. i1 ) .and. ( i1 .le. ixr ) .and. &
                 ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              tm_mask(i1,j2,k2,l) = 0
              ftmp(i1,j2,k2,l) = f(i1+1,j2-1,k2-1,l) + sqrt(three)*delta
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imin_n(2) .ne. imin_glob(2) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            j1 = imin_n(2) - cctk_lbnd(2)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              tm_mask(i2,j1,k1,l) = 0
              ftmp(i2,j1,k1,l) = f(i2-1,j1+1,k1+1,l) + sqrt(three)*delta
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imin_n(2) .ne. imin_glob(2) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            j1 = imin_n(2) - cctk_lbnd(2)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( max(jyl,2) .le. j1 ) .and. ( j1 .le. jyr ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              tm_mask(i2,j1,k2,l) = 0
              ftmp(i2,j1,k2,l) = f(i2-1,j1+1,k2-1,l) + sqrt(three)*delta
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imax_n(2) .ne. imax_glob(2) ) .and. &
               ( imin_n(3) .ne. imin_glob(3) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            j2 = imax_n(2) - cctk_lbnd(2)
            k1 = imin_n(3) - cctk_lbnd(3)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) .and. &
                 ( max(kzl,2) .le. k1 ) .and. ( k1 .le. kzr ) ) then
              tm_mask(i2,j2,k1,l) = 0
              ftmp(i2,j2,k1,l) = f(i2-1,j2-1,k1+1,l) + sqrt(three)*delta
            end if
          end if
          if ( ( imax_n(1) .ne. imax_glob(1) ) .and. &
               ( imax_n(2) .ne. imax_glob(2) ) .and. &
               ( imax_n(3) .ne. imax_glob(3) ) ) then
            i2 = imax_n(1) - cctk_lbnd(1)
            j2 = imax_n(2) - cctk_lbnd(2)
            k2 = imax_n(3) - cctk_lbnd(3)
            if ( ( ixl .le. i2 ) .and. ( i2 .le. min(ixr,nx-1) ) .and. &
                 ( jyl .le. j2 ) .and. ( j2 .le. min(jyr,ny-1) ) .and. &
                 ( kzl .le. k2 ) .and. ( k2 .le. min(kzr,nz-1) ) ) then
              tm_mask(i2,j2,k2,l) = 0
              ftmp(i2,j2,k2,l) = f(i2-1,j2-1,k2-1,l) + sqrt(three)*delta
            end if
          end if
        end if

!       Copy the modified mask and level set function into the proper place.
        eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) = tm_mask(ixl:ixr,jyl:jyr,kzl:kzr,l)
        f(ixl:ixr,jyl:jyr,kzl:kzr,l) = ftmp(ixl:ixr,jyl:jyr,kzl:kzr,l)

!       Now check if any interior points are far enough away from the f=0
!       surface and if so excise them.
        if ( use_inner_excision .gt. 0 ) then
          where ( ( eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) .ge. 0 ) .and. &
                  ( f(ixl:ixr,jyl:jyr,kzl:kzr,l) .lt. -shell_width*delta ) ) 
            eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) = -1
            f(ixl:ixr,jyl:jyr,kzl:kzr,l) = ex_value
          end where
        end if

!       Now check if any remaining active points where actually excised in the
!       numerical run generating the metric data.
  
        if ( use_mask .gt. 0 ) then
          where ( emask(ixl:ixr,jyl:jyr,kzl:kzr) .lt. three*quarter )
            eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) = -1
            f(ixl:ixr,jyl:jyr,kzl:kzr,l) = ex_value
          end where
        end if

!       Make sure to mark all points outside of the rectangular box as excised
!       points and set the value of f to -ex_value.
        if ( use_outer_excision .gt. 0 ) then
          do k = kzl, kzr
            do j = jyl, jyr
              do i = ixl, ixr
                if ( ( ( i + cctk_lbnd(1) .lt. imin_n(1) ) .or. &
                       ( i + cctk_lbnd(1) .gt. imax_n(1) ) .or. &
                       ( j + cctk_lbnd(2) .lt. imin_n(2) ) .or. &
                       ( j + cctk_lbnd(2) .gt. imax_n(2) ) .or. &
                       ( k + cctk_lbnd(3) .lt. imin_n(3) ) .or. &
                       ( k + cctk_lbnd(3) .gt. imax_n(3) ) ) .and. &
                       ( eh_mask(i,j,k,l) .ge. 0 ) ) then
                  eh_mask(i,j,k,l) = -1
                  f(i,j,k,l) = -ex_value
                end if
              end do
            end do
          end do
        end if
      end if not_undone
    end do loop_over_l
  end if
end subroutine EHFinder_SetMask1


! This routine finds the excision boundary after it has been changed and
! makes sure that it has the right mask value.
subroutine EHFinder_SetMask2(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  logical :: active

  active = .false.

! If re-parametrization has just been done, we are in the process of
! resetting the mask so we set active=.true.
  if ( mask_first ) active = .true.
  if ( re_initialize_every .gt. 0 ) then
    if ( mod(cctk_iteration,re_initialize_every) .eq. 0 ) active = .true.
  end if

! If the reparametrization was not undone...
  if ( active .and. .not. all(re_initialize_undone) ) then

! Get the minimum and maximum index excluding ghost and symmetry cells.
# include "include/physical_part.h"
    
!   Make sure we are not on the physical outer boundary.
    if ( ixl .eq. 1 ) ixl = 2
    if ( ixr .eq. nx ) ixr = nx - 1
    if ( jyl .eq. 1 ) jyl = 2
    if ( jyr .eq. ny ) jyr = ny - 1
    if ( kzl .eq. 1 ) kzl = 2
    if ( kzr .eq. nz ) kzr = nz - 1

    do l = 1, eh_number_level_sets
!     Initialize the temporary mask to zero.
      tm_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) = 0

      if ( .not. re_initialize_undone(l) ) then

!       We loop over all points and adjust the mask if the point is
!       on the excision boundary.
        do k = kzl, kzr
          do j = jyl, jyr
            do i = ixl, ixr

!             If the point is active.....
              if ( eh_mask(i,j,k,l) .ge. 0 ) then

!               If the neighbour in the -x-directions is excised....
                if ( eh_mask(i-1,j,k,l) .eq. -1 ) then
                  tm_mask(i,j,k,l) = tm_mask(i,j,k,l) + ll(0)
                end if

!               If the neighbour in the +x-directions is excised....
                if ( eh_mask(i+1,j,k,l) .eq. -1 ) then
                  tm_mask(i,j,k,l) = tm_mask(i,j,k,l) + ll(1)
                end if

!               If the neighbour in the -y-directions is excised....
                if ( eh_mask(i,j-1,k,l) .eq. -1 ) then
                  tm_mask(i,j,k,l) = tm_mask(i,j,k,l) + ll(2)
                end if

!               If the neighbour in the +y-directions is excised....
                if ( eh_mask(i,j+1,k,l) .eq. -1 ) then
                  tm_mask(i,j,k,l) = tm_mask(i,j,k,l) + ll(3)
                end if

!               If the neighbour in the -z-directions is excised....
                if ( eh_mask(i,j,k-1,l) .eq. -1 ) then
                  tm_mask(i,j,k,l) = tm_mask(i,j,k,l) + ll(4)
                end if

!               If the neighbour in the +z-directions is excised....
                if ( eh_mask(i,j,k+1,l) .eq. -1 ) then
                  tm_mask(i,j,k,l) = tm_mask(i,j,k,l) + ll(5)
                end if

              end if

            end do
          end do
        end do

!       Copy the changes back into eh_mask
        where ( eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) .ge. 0 )
          eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,l) = &
                        tm_mask(ixl:ixr,jyl:jyr,kzl:kzr,l)
        end where
      end if
    end do

!   Indicate that it is not anymore the first time the mask is set.
    mask_first = .false.
  end if

end subroutine EHFinder_SetMask2


! This routine checks if any of the excision regions are to close together
! so that some points do not have enough active neigbouring points to be
! able to calculate derivatives. If this is the case these points are 
! excised as well.
subroutine EHFinder_SetMask3(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  logical active, mask_modified

  active = .false.

! If re-parametrization has just been done, we are in the process of
! resetting the mask so we set active=.true.
  if ( mask_first ) active = .true.
  if ( re_initialize_every .gt. 0 ) then
    if ( mod(cctk_iteration,re_initialize_every) .eq. 0 ) active = .true.
  end if

! If the reparametrization was not undone...
  if ( active .and. .not. all(re_initialize_undone) ) then

    mask_modified = .false.
# include "include/physical_part.h"

!   Make sure we are not on the physical outer boundary.
    if ( ixl .eq. 1 ) ixl = 2
    if ( ixr .eq. nx ) ixr = nx - 1
    if ( jyl .eq. 1 ) jyl = 2
    if ( jyr .eq. ny ) jyr = ny - 1
    if ( kzl .eq. 1 ) kzl = 2
    if ( kzr .eq. nz ) kzr = nz - 1

    do l = 1, eh_number_level_sets
      tm_mask(:,:,:,l) = eh_mask(:,:,:,l)

      if ( .not. re_initialize_undone(l) ) then

        do k = kzl, kzr
          do j = jyl, jyr
            do i = ixl, ixr
              if ( eh_mask(i,j,k,l) .ge. 0 ) then

!               If we are on an excision boundary in the -x-direction...
                if ( btest(eh_mask(i,j,k,l),0) ) then

!                 If any of the two closest neighbours in the +x-direction is
!                 excised we have to excise this point
                  if ( ( eh_mask(i+1,j,k,l) .eq. -1 ) .or. &
                       ( eh_mask(i+2,j,k,l) .eq. -1 ) ) then
                    tm_mask(i,j,k,l) = -1
                    mask_modified = .true.
                  end if
                end if

!               If we are on an excision boundary in the +x-direction...
                if ( btest(eh_mask(i,j,k,l),1) ) then

!                 If any of the two closest neighbours in the -x-direction is
!                 excised we have to excise this point
                  if ( ( eh_mask(i-1,j,k,l) .eq. -1 ) .or. &
                       ( eh_mask(i-2,j,k,l) .eq. -1 ) ) then
                    tm_mask(i,j,k,l) = -1
                    mask_modified = .true.
                  end if
                end if

!               If we are on an excision boundary in the -y-direction...
                if ( btest(eh_mask(i,j,k,l),2) ) then

!                 If any of the two closest neighbours in the +y-direction is
!                 excised we have to excise this point
                  if ( ( eh_mask(i,j+1,k,l) .eq. -1 ) .or. &
                       ( eh_mask(i,j+2,k,l) .eq. -1 ) ) then
                    tm_mask(i,j,k,l) = -1
                    mask_modified = .true.
                  end if
                end if

!               If we are on an excision boundary in the +y-direction...
                if ( btest(eh_mask(i,j,k,l),3) ) then

!                 If any of the two closest neighbours in the -y-direction is
!                 excised we have to excise this point
                  if ( ( eh_mask(i,j-1,k,l) .eq. -1 ) .or. &
                       ( eh_mask(i,j-2,k,l) .eq. -1 ) ) then
                    tm_mask(i,j,k,l) = -1
                    mask_modified = .true.
                  end if
                end if

!               If we are on an excision boundary in the -z-direction...
                if ( btest(eh_mask(i,j,k,l),4) ) then

!                 If any of the two closest neighbours in the +z-direction is
!                 excised we have to excise this point
                  if ( ( eh_mask(i,j,k+1,l) .eq. -1 ) .or. &
                       ( eh_mask(i,j,k+2,l) .eq. -1 ) ) then
                    tm_mask(i,j,k,l) = -1
                    mask_modified = .true.
                  end if
                end if

!               If we are on an excision boundary in the +z-direction...
                if ( btest(eh_mask(i,j,k,l),5) ) then

!                 If any of the two closest neighbours in the -z-direction is
!                 excised we have to excise this point
                  if ( ( eh_mask(i,j,k-1,l) .eq. -1 ) .or. &
                       ( eh_mask(i,j,k-2,l) .eq. -1 ) ) then
                    tm_mask(i,j,k,l) = -1
                    mask_modified = .true.
                  end if
                end if
              end if

            end do
          end do
        end do

        if ( mask_modified ) then
          call CCTK_WARN ( 1, 'Mask modified because it was not convex' )
        end if
      end if
    end do
    eh_mask = tm_mask
  end if

end subroutine EHFinder_SetMask3
