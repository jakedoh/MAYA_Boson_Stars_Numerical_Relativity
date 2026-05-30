! Finding the surface(s). Right now only in full mode but fully parallelised.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_FindSurface(CCTK_ARGUMENTS)

  use EHFinder_mod
  use EHFinder_fuzzy

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l, im, jm, sn
  CCTK_REAL, parameter :: eps = 1.0d-10
  CCTK_INT :: interp_handle, table_handle, status, coord_system_handle

  character(len=80) :: info_message
  character(len=128) :: warn_message
  character(len=200) :: surface_interp
  CCTK_INT :: surface_interp_len
  character(len=7) :: surface_order
  character(len=15) :: vname

  CCTK_INT, dimension(4) :: bbox
  CCTK_INT, dimension(2) :: gsh, lsh, lbnd, ubnd, nghost

  CCTK_REAL, dimension(3) :: center_loc
  CCTK_INT :: nc_loc, nc
  CCTK_REAL :: ncinv, maxdr_loc, maxdr, maxf_loc, maxf
  CCTK_REAL, dimension(3) :: crange_min_loc, crange_max_loc, &
                             crange_min_glob, crange_max_glob
  CCTK_REAL :: ltheta, utheta, lphi, uphi

  CCTK_REAL :: dtheta, dphi, dthetainv, dphiinv, min_delta, rloc
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(4) :: out_array
  CCTK_INT :: in_array
  CCTK_INT, dimension(4), parameter :: op_indices = (/ 0, 0, 0, 0 /), &
                                       op_codes = (/ 0, 1, 2, 3 /), &
                                       out_types = (/ CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL /)
  
! Find index ranges for points excluding ghost zones and symmetry points for
! the 3D grid functions.
#include "include/physical_part.h"

! Initialize find_surface_status to 0. If something goes wrong it will be
! set to -1.
  find_surface_status = 0
  sn = surface_counter - integrate_counter + 1

! Set the index to the current level set function.
  l = levelset_counter

! Write an info message about what is going on.
  write(info_message,'(a26,i3,a14,i3)') 'Finding points on surface ', sn, &
                                 ' in level set ',l
  call CCTK_INFO ( trim(info_message) )

! Find information about the local and global properties of the 2D
! Grid arrays.
  call CCTK_GroupbboxGN ( status, cctkGH, 4, bbox, 'ehfinder::surface_arrays' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get bounding box for surface arrays' )
  end if
  call CCTK_GroupgshGN ( status, cctkGH, 2, gsh, 'ehfinder::surface_arrays' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get global size for surface arrays' )
  end if
  call CCTK_GrouplbndGN ( status, cctkGH, 2, lbnd, 'ehfinder::surface_arrays' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get lower bounds for surface arrays' )
  end if
  call CCTK_GroupubndGN ( status, cctkGH, 2, ubnd, 'ehfinder::surface_arrays' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get upper bounds for surface arrays' )
  end if
  call CCTK_GrouplshGN ( status, cctkGH, 2, lsh, 'ehfinder::surface_arrays' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if
  call CCTK_GroupnghostzonesGN ( status, cctkGH, 2, nghost, &
                                         'ehfinder::surface_arrays' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if
  
! If the user defined center point should be used...
  if ( use_user_center .gt. 0 ) then
    center(1) = center_x
    center(2) = center_y
    center(3) = center_z
  else

!   If not find an approximation for it by averaging the coodinates that
!   are within one grid cell distance from the surface in question.
    nc_loc = 0; center_loc = zero
    min_delta = minval ( cctk_delta_space )
    do k = kzl, kzr
      do j = jyl, jyr
        do i = ixl, ixr
          if ( eh_mask(i,j,k,l) .eq. 0 ) then
            if ( ( ( f(i,j,k,l) * f(i-1,j,k,l) .lt. zero ) .or. &
                   ( f(i,j,k,l) * f(i+1,j,k,l) .lt. zero ) .or. &
                   ( f(i,j,k,l) * f(i,j-1,k,l) .lt. zero ) .or. &
                   ( f(i,j,k,l) * f(i,j+1,k,l) .lt. zero ) .or. &
                   ( f(i,j,k,l) * f(i,j,k-1,l) .lt. zero ) .or. &
                   ( f(i,j,k,l) * f(i,j,k+1,l) .lt. zero ) ) .and. &
                 ( fuzzy_eq ( sc(i,j,k), real(sn,wp), small ) .or. &
                   fuzzy_eq ( sc(i-1,j,k), real(sn,wp), small ) .or. & 
                   fuzzy_eq ( sc(i+1,j,k), real(sn,wp), small ) .or. &
                   fuzzy_eq ( sc(i,j-1,k), real(sn,wp), small ) .or. &
                   fuzzy_eq ( sc(i,j+1,k), real(sn,wp), small ) .or. &
                   fuzzy_eq ( sc(i,j,k-1), real(sn,wp), small ) .or. &
                   fuzzy_eq ( sc(i,j,k+1), real(sn,wp), small ) ) ) then
              nc_loc = nc_loc + 1
              center_loc(1) = center_loc(1) + x(i,j,k)
              center_loc(2) = center_loc(2) + y(i,j,k)
              center_loc(3) = center_loc(3) + z(i,j,k)
            end if
          end if
        end do
      end do
    end do

!   Reduce over the processors.
    call CCTK_ReduceLocScalar ( status, cctkGH, -1, sum_handle, &
                                CCTK_PointerTo( nc_loc ), CCTK_PointerTo( nc ), &
                                CCTK_VARIABLE_INT )
    call CCTK_ReduceLocArrayToArray1D ( status, cctkGH, -1, sum_handle, &
                                        CCTK_PointerTo( center_loc ), &
                                        CCTK_PointerTo( center ) , 3, &
                                        CCTK_VARIABLE_REAL )
    ncinv = one / nc
    center = ncinv * center
  end if

! Try to figure out the symmetries. First find the range in coordinates
! for the points inside the current surface.
  crange_min_loc(1) = minval ( x, mask = fuzzy_eq ( sc, real(sn,wp), small ) )
  crange_min_loc(2) = minval ( y, mask = fuzzy_eq ( sc, real(sn,wp), small ) )
  crange_min_loc(3) = minval ( z, mask = fuzzy_eq ( sc, real(sn,wp), small ) )
  crange_max_loc(1) = maxval ( x, mask = fuzzy_eq ( sc, real(sn,wp), small ) )
  crange_max_loc(2) = maxval ( y, mask = fuzzy_eq ( sc, real(sn,wp), small ) )
  crange_max_loc(3) = maxval ( z, mask = fuzzy_eq ( sc, real(sn,wp), small ) )

! Reduce over all processors.
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, min_handle, &
                                      CCTK_PointerTo( crange_min_loc ), &
                                      CCTK_PointerTo( crange_min_glob ), &
                                      3, CCTK_VARIABLE_REAL )
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, max_handle, &
                                      CCTK_PointerTo( crange_max_loc ), &
                                      CCTK_PointerTo( crange_max_glob ), &
                                      3, CCTK_VARIABLE_REAL )

! Full mode. Modify if necessary.
  ltheta = zero; utheta = pi
  lphi = zero; uphi = two * pi
  sym_factor = one
  s_symx = .false.; s_symy = .false.; s_symz = .false.

  if ( CCTK_EQUALS ( domain, 'bitant' ) ) then

!   Bitant with symmetry in the z-direction

    if ( CCTK_EQUALS ( bitant_plane, 'xy' ) ) then
!     If the minimum z-coordinate of the points inside the surface is less
!     than zero, then the surface must share the grid symmetry across the
!     xy-plane.
      if ( crange_min_glob(3) .lt. zero ) then
        ltheta = zero; utheta = half * pi
        lphi = zero; uphi = two * pi
        sym_factor = two
        center(3) = zero
        s_symz = .true.
      end if 
    end if

!   Bitant with symmetry in the y-direction
  
    if ( CCTK_EQUALS ( bitant_plane, 'xz' ) ) then
      if ( crange_min_glob(2) .lt. zero ) then
        ltheta = zero; utheta = pi
        lphi = zero; uphi = pi
        sym_factor = two
        center(2) = zero
        s_symy = .true.
      end if 
    end if

!   Bitant with symmetry in the x-direction

    if ( CCTK_EQUALS ( bitant_plane, 'yz' ) ) then
      if ( crange_min_glob(1) .lt. zero ) then
        ltheta = zero; utheta = pi
        lphi = -half * pi; uphi = half * pi
        sym_factor = two
        center(1) = zero
        s_symx = .true.
      end if 
    end if
  end if

  if ( CCTK_EQUALS ( domain, 'quadrant' ) ) then

!   Quadrant in x-direction

    if ( CCTK_EQUALS ( quadrant_direction, 'x' ) ) then
      if ( ( crange_min_glob(3) .lt. zero ) .and. &
           ( crange_min_glob(2) .lt. zero ) ) then
        ltheta = zero; utheta = half * pi
        lphi = zero; uphi = pi
        sym_factor = four
        center(3) = zero; center(2) = zero
        s_symz = .true.; s_symy = .true.
      else if ( crange_min_glob(3) .lt. zero ) then
        ltheta = zero; utheta = half * pi
        lphi = zero; uphi = two * pi
        sym_factor = two
        center(3) = zero
        s_symz = .true.
      else if ( crange_min_glob(2) .lt. zero ) then
        ltheta = zero; utheta = pi
        lphi = zero; uphi = pi
        sym_factor = two
        center(2) = zero
        s_symy = .true.
      end if
    end if
    
!   Quadrant in y-direction

    if ( CCTK_EQUALS ( quadrant_direction, 'y' ) ) then
      if ( ( crange_min_glob(3) .lt. zero ) .and. &
           ( crange_min_glob(1) .lt. zero ) ) then
        ltheta = zero; utheta = half * pi
        lphi = -half * pi; uphi = half * pi
        sym_factor = four
        center(3) = zero; center(1) = zero
        s_symz = .true.; s_symx = .true.
      else if ( crange_min_glob(3) .lt. zero ) then
        ltheta = zero; utheta = half * pi
        lphi = zero; uphi = two * pi
        sym_factor = two
        center(3) = zero
        s_symz = .true.
      else if ( crange_min_glob(1) .lt. zero ) then
        ltheta = zero; utheta = pi
        lphi = -half * pi; uphi = half * pi
        sym_factor = two
        center(1) = zero
        s_symx = .true.
      end if
    end if

!   Quadrant in z-direction

    if ( CCTK_EQUALS ( quadrant_direction,'z' ) ) then
      if ( ( crange_min_glob(2) .lt. zero ) .and. &
           ( crange_min_glob(1) .lt. zero ) ) then
        ltheta = zero; utheta = pi
        lphi = zero; uphi = half * pi
        sym_factor = four
        center(2) = zero; center(1) = zero
        s_symy = .true.; s_symx = .true.
      else if ( crange_min_glob(2) .lt. zero ) then
        ltheta = zero; utheta = pi
        lphi = zero; uphi = pi
        sym_factor = two
        center(2) = zero
        s_symy = .true.
      else if ( crange_min_glob(1) .lt. zero ) then
        ltheta = zero; utheta = pi
        lphi = -half * pi; uphi = half * pi
        sym_factor = two
        center(1) = zero
        s_symx = .true.
      end if
    end if
  end if
    
! Octant

  if ( CCTK_EQUALS ( domain, 'octant' ) ) then
    if ( ( crange_min_glob(3) .lt. zero ) .and. &
         ( crange_min_glob(2) .lt. zero ) .and. &
         ( crange_min_glob(1) .lt. zero ) ) then
      ltheta = zero; utheta = half * pi
      lphi = zero; uphi = half * pi
      sym_factor = eight
      center(3) = zero; center(2) = zero; center(1) = zero
      s_symz = .true.; s_symy = .true.; s_symx = .true.
    else if ( ( crange_min_glob(3) .lt. zero ) .and. &
              ( crange_min_glob(2) .lt. zero ) ) then
      ltheta = zero; utheta = half * pi
      lphi = zero; uphi = pi
      sym_factor = four
      center(3) = zero; center(2) = zero
      s_symz = .true.; s_symy = .true.
    else if ( ( crange_min_glob(3) .lt. zero ) .and. &
              ( crange_min_glob(1) .lt. zero ) ) then
      ltheta = zero; utheta = half * pi
      lphi = -half * pi; uphi = half * pi
      sym_factor = four
      center(3) = zero; center(1) = zero
      s_symz = .true.; s_symx = .true.
    else if ( ( crange_min_glob(2) .lt. zero ) .and. &
              ( crange_min_glob(1) .lt. zero ) ) then
      ltheta = zero; utheta = pi
      lphi = zero; uphi = half * pi
      sym_factor = four
      center(2) = zero; center(1) = zero
      s_symy = .true.; s_symx = .true.
    else if ( crange_min_glob(3) .lt. zero ) then 
      ltheta = zero; utheta = half * pi
      lphi = zero; uphi = two * pi
      sym_factor = two
      center(3) = zero
      s_symz = .true.
    else if ( crange_min_glob(2) .lt. zero ) then 
      ltheta = zero; utheta = pi
      lphi = zero; uphi = pi
      sym_factor = two
      center(2) = zero
      s_symy = .true.
    else if ( crange_min_glob(1) .lt. zero ) then 
      ltheta = zero; utheta = pi
      lphi = -half * pi; uphi = half * pi
      sym_factor = two
      center(1) = zero
      s_symx = .true.
    end if
  end if

! Set the theta index to be used for the equatorial circumference integration
! based on the symmetries.
  if ( .not. s_symz ) then
    githeta = ( ntheta - 1 ) / 2 + 1
  else
    githeta = ntheta
  end if
  
! Set the phi index to be used for the polar circumference integration based
! on the symmetries.
  if ( ( .not. s_symy ) .and. ( s_symx ) ) then
    gjphi = ( nphi - 1 ) / 2 + 1
  else
    gjphi = 1
  end if

! Set up the symmetry muliplication factors based on the range in the
! angular coordinates.
  theta_sym_factor = two * pi / ( utheta - ltheta )
  phi_sym_factor = two * pi / ( uphi - lphi )

! Find dtheta and dphi and initialise the theta and phi grid arrays.
! Here i + lbnd(1) - 1 is the global index for theta and
!      j + lbnd(2) - 1 is the global index for phi.

  dtheta = ( utheta - ltheta ) / ( ntheta - 1 )
  dphi = ( uphi - lphi ) / ( nphi - 1 )
  dthetainv = one / dtheta
  dphiinv = one / dphi
  do i = 1, lsh(1)
    ctheta(i,:) = ltheta + dtheta * ( i + lbnd(1) - 1 )
  end do
  do j = 1, lsh(2)
    cphi(:,j) = lphi + dphi * ( j +lbnd(2) - 1 )
  end do

! Calculate sines and cosines and store them in grid arrays since they
! are computationally expensive and memory for 2D grid arrays is cheap.
  sintheta = sin(ctheta)
  costheta = cos(ctheta)
  sinphi = sin(cphi)
  cosphi = cos(cphi)
  
! rsurf is the radial grid array (as function of theta and phi) and is
! initialised to zero. The grid array drsurf contains the changes to rsurf
! and is initialised to 2*Delta. The grid array n_since_last_reduction is
! a help in keeping track on how long ago drsurf for a given point has been
! reduced and is initialised to -1 to indicate that no reduction has taken
! place.
  drsurf = two * min_delta
  rsurf = zero
  n_since_last_reduction = -1

! Get and interpolation handle for required polynomial interpolation.
! If Lagrange polynomial interpolation does not work try Hermite
! polynomial interpolation, which avoids problems with non-continuous
! derivatives in the Newton iteration later on.

! First convert the string parameter to a Fortran string.
  call CCTK_FortranString ( surface_interp_len, surface_interpolator, &
                                                surface_interp )

! Then get the handle ...
  call CCTK_InterpHandle ( interp_handle, surface_interp(1:surface_interp_len) )

! The following is done because some C-preprocessors seem to have problems
! handling character strings that are continued on several lines.
  if ( interp_handle .lt. 0 ) then
    warn_message = 'Cannot get handle for interpolation. '
    warn_message = trim(warn_message)//'Forgot to activate an implementation '
    warn_message = trim(warn_message)//'providing interpolation operators?'
    call CCTK_WARN( 0, trim(warn_message) )
  end if

! Write the required interpolation order parameter into a fortran string.
  write(surface_order,'(a6,i1)') 'order=',surface_interpolation_order

! Then create a table from the string.
  call Util_TableCreateFromString ( table_handle, surface_order )
  if ( table_handle .lt. 0 ) then
    call CCTK_WARN( 0, 'Cannot create parameter table for interpolator' )
  end if
  
! Get the 3D coordinate system handle.
  call CCTK_CoordSystemHandle ( coord_system_handle, 'cart3d' )
  if ( coord_system_handle .lt. 0) then
    warn_message = 'Cannot get handle for cart3d coordinate system. '
    warn_message = trim(warn_message)//'Forgot to activate an implementation '
    warn_message = trim(warn_message)//'providing coordinates?'
    call CCTK_WARN( 0, trim(warn_message) )
  endif

! Get pointers to the grid arrays containing the interpolation points.
  interp_coords(1) = CCTK_PointerTo(interp_x)
  interp_coords(2) = CCTK_PointerTo(interp_y)
  interp_coords(3) = CCTK_PointerTo(interp_z)

! Get pointer to the interpolation output array
  out_array(1) = CCTK_PointerTo(f_interp)

! Get the variable index for f
  write(vname,'(a12,i1,a1)') 'ehfinder::f[', l-1, ']'
  call CCTK_VarIndex ( in_array, vname )

! Loop to find starting point in all directions at the same time. For
! simplicity all ghost zones are found on all processors
! The algorithm simply starts at the center and works it way outwards
! in all directions. When the function changes signs, it reduces the
! changes in radius until all points have negative f values and are
! within one grid cells difference from the surface.
  find_starting_points: do

    do j = 1, lsh(2)
      do i = 1, lsh(1)

!       The radius that are tested in this direction.
        rloc = rsurf(i,j) + drsurf(i,j)

!       Calculate the interpolation points.
        interp_x(i,j) = center(1) + rloc * sintheta(i,j) * cosphi(i,j)
        interp_y(i,j) = center(2) + rloc * sintheta(i,j) * sinphi(i,j)
        interp_z(i,j) = center(3) + rloc * costheta(i,j)

!       If the point has been reduced once, increment n_since_last_reduction.
        if ( n_since_last_reduction(i,j) .ge. 0 ) then 
          n_since_last_reduction(i,j) = n_since_last_reduction(i,j) + 1
        end if
      end do
    end do

!   Call the interpolator.
    call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
                                 table_handle, coord_system_handle, &
                                 lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                                 interp_coords, 1, in_array, &
                                 1, (/CCTK_VARIABLE_REAL/), out_array(1) )
    if ( status .lt. 0 ) then
       call CCTK_INFO ( 'Interpolation failed. Giving up on finding surfaces' )
       find_surface_status = -1
       integrate_counter = integrate_counter - 1
       return
    end if
                               
!   Find the maximal change in rsurf.
    maxdr_loc = maxval ( drsurf )
    call CCTK_ReduceLocScalar ( status, cctkGH, -1, max_handle, &
                                CCTK_PointerTo( maxdr_loc ), &
                                CCTK_PointerTo( maxdr ), CCTK_VARIABLE_REAL )

!   Find the maximal value of interpolated f values.
    maxf_loc = maxval ( abs ( f_interp ) )
    call CCTK_ReduceLocScalar ( status, cctkGH, -1, max_handle, &
                                CCTK_PointerTo( maxf_loc ), &
                                CCTK_PointerTo( maxf ), CCTK_VARIABLE_REAL )
    
!   If the maximal change in rsurf is small enough, we are satisfied that
!   the Newton iteration will be able to find the root.
    if ( maxdr .lt. min_delta ) exit find_starting_points

!   Check if we have crossed the zero point and if we have reduce the
!   change in rsurf we are going to try next. 
    do j = 1, lsh(2)
      do i = 1, lsh(1)
        if ( f_interp(i,j) .gt. zero ) then
          drsurf(i,j) = half * drsurf(i,j)
          n_since_last_reduction(i,j) = 0
        else
          rsurf(i,j) = rsurf(i,j) + drsurf(i,j)
          if ( n_since_last_reduction(i,j) .ge. 1 ) then
            drsurf(i,j) = half * drsurf(i,j)
          end if
        end if
      end do
    end do
    
  end do find_starting_points

! Now for the Newton iteration we also need the interpolates derivatives.
  out_array(2) = CCTK_PointerTo(dfdx_interp)
  out_array(3) = CCTK_PointerTo(dfdy_interp)
  out_array(4) = CCTK_PointerTo(dfdz_interp)

  call Util_TableSetIntArray ( status, table_handle, 4, &
                               op_indices, 'operand_indices' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'Cannot set operand indices array in parameter table' )
  endif
  call Util_TableSetIntArray ( status, table_handle, 4, &
                               op_codes, 'operation_codes' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'Cannot set operation codes array in parameter table' )
  endif

! find the surface points in all directions simultaneously.

  find_surface_points: do k = 1, 15

!   Calculate the interpolation points.
    do j = 1, lsh(2)
      do i = 1, lsh(1)
        interp_x(i,j) = center(1) + rsurf(i,j) * sintheta(i,j) * cosphi(i,j)
        interp_y(i,j) = center(2) + rsurf(i,j) * sintheta(i,j) * sinphi(i,j)
        interp_z(i,j) = center(3) + rsurf(i,j) * costheta(i,j)
      end do
    end do

!   Call the interpolator.
    call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
                                 table_handle, coord_system_handle, &
                                 lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                                 interp_coords, 1, in_array, &
                                 4, out_types, out_array )
    if ( status .lt. 0 ) then
       call CCTK_INFO ( 'Interpolation failed. Giving up on finding surfaces' )
       find_surface_status = -1
       integrate_counter = integrate_counter - 1
       return
    end if
 
!   Find the maximum value of interpolated f values.
    maxf_loc = maxval ( abs ( f_interp ) )
    call CCTK_ReduceLocScalar ( status, cctkGH, -1, max_handle, &
                                CCTK_PointerTo( maxf_loc ), &
                                CCTK_PointerTo( maxf ), CCTK_VARIABLE_REAL )
    
!   If the maximum interpolated value of f is to big, perform a Newton
!   iteration
    if ( maxf .gt. eps ) then
      do j = 1, lsh(2)
        do i = 1, lsh(1)
          rsurf(i,j) = rsurf(i,j) - f_interp(i,j) / &
                       ( dfdx_interp(i,j) * sintheta(i,j) * cosphi(i,j) + &
                         dfdy_interp(i,j) * sintheta(i,j) * sinphi(i,j) + &
                         dfdz_interp(i,j) * costheta(i,j) )
        end do
      end do
    else
!     Else we are done
      exit find_surface_points
    end if

  end do find_surface_points

! If we used to many iterations, the Newton iteration did not converge, so
! print a warning.  
  if ( k .gt. 15 ) then
    warn_message = 'Newton iteration did not converge! Area may be inaccurate'
    call CCTK_WARN(1, trim(warn_message) )
  end if

  drdtheta = zero
  drdphi = zero

! Calculate the derivatives of r.

! First for interior points
  do j = 2, lsh(2) - 1
    do i = 2, lsh(1) - 1
      drdtheta(i,j) = half * dthetainv * ( rsurf(i+1,j) - rsurf(i-1,j) )
      drdphi(i,j) = half * dphiinv * ( rsurf(i,j+1) - rsurf(i,j-1) )
    end do
  end do

! If boundary 1 is physical update those points.
  if ( bbox(1) .eq. 1 ) then
    do j = 2, lsh(2) - 1
      drdtheta(1,j) = half * dthetainv * ( - three * rsurf(1,j) + &
                                             four * rsurf(2,j) - &
                                             rsurf(3,j) )
      drdphi(1,j) = half * dphiinv * ( rsurf(1,j+1) - rsurf(1,j-1) )
    end do
  end if

! If boundary 2 is physical update those points.
  if ( bbox(2) .eq. 1 ) then
    im = lsh(1)
    do j = 2, lsh(2) - 1
      drdtheta(im,j) = half * dthetainv * ( three * rsurf(im,j) - &
                                            four * rsurf(im-1,j) + &
                                            rsurf(im-2,j) )
      drdphi(im,j) = half * dphiinv * ( rsurf(im,j+1) - rsurf(im,j-1) )
    end do
  end if

! If boundary 3 is physical update those points.
  if ( bbox(3) .eq. 1 ) then
    do i = 2, lsh(1) - 1
      drdtheta(i,1) = half * dthetainv * ( rsurf(i+1,1) - rsurf(i-1,1) )
      drdphi(i,1) = half * dphiinv * ( - three * rsurf(i,1) + &
                                         four * rsurf(i,2) - &
                                         rsurf(i,3) )
    end do
  end if

! If boundary 4 is physical update those points.
  if ( bbox(4) .eq. 1 ) then
    jm = lsh(2)
    do i = 2, lsh(1) - 1
      drdtheta(i,jm) = half * dthetainv * ( rsurf(i+1,jm) - rsurf(i-1,jm) )
      drdphi(i,jm) = half * dphiinv * ( three * rsurf(i,jm) - &
                                        four * rsurf(i,jm-1) + &
                                        rsurf(i,jm-2) )
    end do
  end if

! If corner 1 is physical update that point.
  if ( bbox(1) .eq. 1 .and. bbox(3) .eq. 1 ) then
    drdtheta(1,1) = half * dthetainv * ( - three * rsurf(1,1) + &
                                             four * rsurf(2,1) - &
                                             rsurf(3,1) )
    drdphi(1,1) = half * dphiinv * ( - three * rsurf(1,1) + &
                                       four * rsurf(1,2) - &
                                       rsurf(1,3) )
  end if

! If corner 2 is physical update that point.
  if ( bbox(1) .eq. 1 .and. bbox(4) .eq. 1 ) then
    jm = lsh(2)
    drdtheta(1,jm) = half * dthetainv * ( - three * rsurf(1,jm) + &
                                             four * rsurf(2,jm) - &
                                             rsurf(3,jm) )
    drdphi(1,jm) = half * dphiinv * ( three * rsurf(1,jm) - &
                                      four * rsurf(1,jm-1) + &
                                      rsurf(1,jm-2) )
  end if
 
! If corner 3 is physical update that point.
  if ( bbox(2) .eq. 1 .and. bbox(3) .eq. 1 ) then
    im = lsh(1)
    drdtheta(im,1) = half * dthetainv * ( three * rsurf(im,1) - &
                                          four * rsurf(im-1,1) + &
                                          rsurf(im-2,1) )
    drdphi(im,1) = half * dphiinv * ( - three * rsurf(im,1) + &
                                        four * rsurf(im,2) - &
                                        rsurf(im,3) )
  end if
 
! If corner 4 is physical update that point.
  if ( bbox(2) .eq. 1 .and. bbox(4) .eq. 1 ) then
    im = lsh(1)
    jm = lsh(2)
    drdtheta(im,jm) = half * dthetainv * ( three * rsurf(im,jm) - &
                                           four * rsurf(im-1,jm) + &
                                           rsurf(im-2,jm) )
    drdphi(im,jm) = half * dphiinv * ( three * rsurf(im,jm) - &
                                       four * rsurf(im,jm-1) + &
                                       rsurf(im,jm-2) )
  end if

! I do not remember why this is done. I will have to investigate.
  integrate_counter = integrate_counter - 1

! Get rid of the interpolation table.
  call Util_TableDestroy ( status, table_handle )
  
end subroutine EHFinder_FindSurface


subroutine EHFinder_LevelSetLoopInit(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  more_levelsets = 1
  levelset_counter = 1

end subroutine EHFinder_LevelSetLoopInit


subroutine EHFinder_UpdateCounter(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  if ( levelset_counter .lt. eh_number_level_sets ) then
    levelset_counter = levelset_counter + 1
    more_levelsets = 1
  else
    more_levelsets = 0
  end if

end subroutine EHFinder_UpdateCounter


subroutine EHFinder_CountSurfacesInit(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  more_surfaces = 1
  surface_counter = 0
  sc = zero

end subroutine EHFinder_CountSurfacesInit
  

subroutine EHFinder_CountSurfaces(CCTK_ARGUMENTS)

  use EHFinder_mod
  use EHFinder_fuzzy

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, dimension(3) :: minl
  CCTK_REAL :: minf_loc, minf_glob
  CCTK_INT :: procpos_loc, procpos_glob

! Find index ranges for points excluding ghost zones and symmetry points for
! the 3D grid functions.
#include "include/physical_part.h"

! Find the location of the minimum among active points not already found
! to be inside a surface.
  minl = minloc ( f(ixl:ixr,jyl:jyr,kzl:kzr,levelset_counter), &
                  fuzzy_eq ( sc(ixl:ixr,jyl:jyr,kzl:kzr), zero, small ) .and. &
                  ( eh_mask(ixl:ixr,jyl:jyr,kzl:kzr,levelset_counter) .ge. 0 ) )

! Find the value of f at that location.
  minf_loc = f(minl(1)+ixl-1,minl(2)+jyl-1,minl(3)+kzl-1,levelset_counter)

! Find the minimum over all processors.
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, min_handle, &
                              CCTK_PointerTo( minf_loc ), &
                              CCTK_PointerTo( minf_glob ), CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0,'Reduction of minf_glob failed')
  end if

! To figure out on which processor the minimum is located I set procpos_loc
! to the local processor number if the local minimum is equal to the global
! minimum otherwise I set it to the total number of processors. When I then
! find the minimum of procpos_loc over all processors I get the processor
! containing the minimum.
  if ( minf_loc .eq. minf_glob ) then
    procpos_loc = CCTK_MyProc ( cctkGH )
  else
    procpos_loc = CCTK_nProcs ( cctkGH )
  end if
  
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, min_handle, &
                              CCTK_PointerTo( procpos_loc ), &
                              CCTK_PointerTo( procpos_glob ), CCTK_VARIABLE_INT )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0,'Reduction of minf_glob failed')
  end if

! If minf_glob<0 there must exist another surface so lets set up the variables
! to start marking and counting these points.
  if ( minf_glob .lt. zero ) then
    surface_counter = surface_counter + 1
    points_counter = 0
    more_points = 1
    if ( procpos_loc .eq. procpos_glob ) then
      sc(minl(1)+ixl-1,minl(2)+jyl-1,minl(3)+kzl-1) = real(surface_counter,wp)
    end if
  else
    more_surfaces = 0
    more_points = 0
  end if

end subroutine EHFinder_CountSurfaces    


! This routine marks all points inside the current surface.
subroutine EHFinder_MarkSurfaces(CCTK_ARGUMENTS)

  use EHFinder_mod
  use EHFinder_fuzzy

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: n_loc, n_glob, i, j, k
  CCTK_INT :: i1, i2, j1, j2, k1, k2
  logical :: marked

! Find index ranges for points excluding ghost zones and symmetry points for
! the 3D grid functions.
#include "include/physical_part.h"
  
! Initialise the local point counter.
  n_loc = 0

! Check if any points have a marked point as a neighbour. Only look at
! active points with f<0 that are not already marked. Take care at the
! boundaries (physical and excision) not to look at inactive or not
! exisiting points.
  do k = kzl, kzr
    do j = jyl, jyr
      do i = ixl, ixr
        marked = .false.
        if ( ( eh_mask(i,j,k,levelset_counter) .ge. 0 ) .and. &
             ( f(i,j,k,levelset_counter) .lt. 0 ) .and. &
             fuzzy_ne ( sc(i,j,k), real(surface_counter,wp), small ) ) then

!         Find out if we are on the boundary
          if ( i + cctk_lbnd(1) .eq. 1 ) then
            i1 = 0
          else
            i1 = -1
          end if
          if ( i + cctk_lbnd(1) .eq. cctk_gsh(1) ) then
            i2 = 0
          else
            i2 = 1
          end if
          if ( j + cctk_lbnd(2) .eq. 1 ) then
            j1 = 0
          else
            j1 = -1
          end if
          if ( j + cctk_lbnd(2) .eq. cctk_gsh(2) ) then
            j2 = 0
          else
            j2 = 1
          end if
          if ( k + cctk_lbnd(3) .eq. 1 ) then
            k1 = 0
          else
            k1 = -1
          end if
          if ( k + cctk_lbnd(3) .eq. cctk_gsh(3) ) then
            k2 = 0
          else
            k2 = 1
          end if

!         Then check if any valid neighbours are marked
          if ( any ( fuzzy_eq ( sc(i+i1:i+i2,j+j1:j+j2,k), &
                     real(surface_counter,wp), small ) ) ) then
            marked = .true.
          end if
          if ( any ( fuzzy_eq ( sc(i+i1:i+i2,j,k+k1:k+k2), &
                     real(surface_counter,wp), small ) ) ) then
            marked = .true.
          end if
          if ( any ( fuzzy_eq ( sc(i,j+j1:j+j2,k+k1:k+k2), &
                     real(surface_counter,wp), small ) ) ) then
            marked = .true.
          end if

          if ( marked ) then
            sc(i,j,k) = real(surface_counter,wp)
            n_loc = n_loc + 1
          end if

        end if    
      end do
    end do
  end do

! Reduce the number of newly marked points over all processors
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, sum_handle, &
                              CCTK_PointerTo( n_loc ), &
                              CCTK_PointerTo( n_glob ), CCTK_VARIABLE_INT )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0,'Reduction of n_glob failed')
  end if

! If the total number of newly marked points is 0, set more_points to 0
! to indicate that we are finished. Otherwise continue and add the number
! of points to the points_counter.
  if ( n_glob .eq. 0 ) then
    more_points = 0
  end if
  points_counter = points_counter + n_glob

end subroutine EHFinder_MarkSurfaces


! This routine prints out information about the found surfaces.
subroutine EHFinder_InfoSurfaces(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  character(len=80) :: info_message

! Write out how many surfaces were found.
  if ( surface_counter .eq. 1 ) then
    write(info_message,'(a9,i3,a22,i3)') 'There is ',surface_counter, &
                                      ' surface in level set ',levelset_counter
  else
    write(info_message,'(a10,i3,a23,i3)') 'There are ',surface_counter, &
                                      ' surfaces in level set ',levelset_counter
  end if
  call CCTK_INFO( trim(info_message) )

  integrate_counter = surface_counter

end subroutine EHFinder_InfoSurfaces
