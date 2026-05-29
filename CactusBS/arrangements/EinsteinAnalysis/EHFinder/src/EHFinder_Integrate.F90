! Integrating over the surface(s). Right now only in full mode.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_FindSurfaceElement(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, l
  CCTK_INT :: interp_handle, table_handle, coord_system_handle

  character(len=200) :: area_interp
  character(len=128) :: warn_message
  CCTK_INT :: area_interp_len
  character(len=7) :: area_order

  CCTK_INT, dimension(4) :: bbox
  CCTK_INT, dimension(2) :: gsh, lsh, lbnd, ubnd, nghost

  CCTK_REAL :: dtheta, dphi, dthetainv, dphiinv
  CCTK_REAL :: dxdth, dxdph, dydth, dydph, dzdth, dzdph
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(7) :: out_array
  CCTK_INT, dimension(7) :: in_array
  CCTK_INT, dimension(7), parameter :: out_types = (/ CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL, &
                                                      CCTK_VARIABLE_REAL /)
  
! If finding of surface failed do not try to integrate but exit.
  if ( find_surface_status .lt. 0 ) then
    return
  endif

! Store the levelset_counter in a shorter named variable for convenience.
  l = levelset_counter

! Write an Info message about what we are doing.
  call CCTK_INFO ( 'Finding surface element' )

! Obtain the dimensions and characteristics of the grid arrays.
  call CCTK_GroupbboxGN ( ierr, cctkGH, 4, bbox, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get bounding box for surface arrays' )
  end if
  call CCTK_GroupgshGN ( ierr, cctkGH, 2, gsh, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get global size for surface arrays' )
  end if
  call CCTK_GrouplbndGN ( ierr, cctkGH, 2, lbnd, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get lower bounds for surface arrays' )
  end if
  call CCTK_GroupubndGN ( ierr, cctkGH, 2, ubnd, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get upper bounds for surface arrays' )
  end if
  call CCTK_GrouplshGN ( ierr, cctkGH, 2, lsh, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if
  call CCTK_GroupnghostzonesGN ( ierr, cctkGH, 2, nghost, &
                                           'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if

! Setup the theta and phi spacings and their inverse.
  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)
  dthetainv = one / dtheta
  dphiinv = one / dphi

! Convert the interpolator parameter into a fortran string.
  call CCTK_FortranString ( area_interp_len, area_interpolator, &
                                             area_interp )

! Get an interpolation handle.
  call CCTK_InterpHandle ( interp_handle, area_interp(1:area_interp_len) )

  if ( interp_handle .lt. 0 ) then
    warn_message = 'Cannot get handle for interpolation. '
    warn_message = trim(warn_message)//'Forgot to activate an implementation '
    warn_message = trim(warn_message)//'providing interpolation operators?'
    call CCTK_WARN( 0, trim(warn_message) )
  end if

! Write the interpolation order parameter into a string
  write(area_order,'(a6,i1)') 'order=',area_interpolation_order

! Create a table from the interpolation order string.
  call Util_TableCreateFromString ( table_handle, area_order )
  if ( table_handle .lt. 0 ) then
    call CCTK_WARN( 0, 'Cannot create parameter table for interpolator' )
  end if

! Get a handle for the coordinate system.
  call CCTK_CoordSystemHandle ( coord_system_handle, 'cart3d' )
  if ( coord_system_handle .lt. 0) then
    warn_message = 'Cannot get handle for cart3d coordinate system. '
    warn_message = trim(warn_message)//'Forgot to activate an implementation '
    warn_message = trim(warn_message)//'providing coordinates?'
    call CCTK_WARN( 0, trim(warn_message) )
  endif

! Get the pointers to the interpolation points.
  interp_coords(1) = CCTK_PointerTo(interp_x)
  interp_coords(2) = CCTK_PointerTo(interp_y)
  interp_coords(3) = CCTK_PointerTo(interp_z)

! Get the pointers to the interpolation return arrays.
  out_array(1) = CCTK_PointerTo(gxxi)
  out_array(2) = CCTK_PointerTo(gxyi)
  out_array(3) = CCTK_PointerTo(gxzi)
  out_array(4) = CCTK_PointerTo(gyyi)
  out_array(5) = CCTK_PointerTo(gyzi)
  out_array(6) = CCTK_PointerTo(gzzi)
  out_array(7) = CCTK_PointerTo(psii)

! find the cartesian coordinates for the interpolation points
  do j = 1, lsh(2)
    do i = 1, lsh(1)
      interp_x(i,j) = center(1) + rsurf(i,j) * sintheta(i,j) * cosphi(i,j)
      interp_y(i,j) = center(2) + rsurf(i,j) * sintheta(i,j) * sinphi(i,j)
      interp_z(i,j) = center(3) + rsurf(i,j) * costheta(i,j)
    end do
  end do

! Get the indices for the grid functions to be interpolated.
  call CCTK_VarIndex ( in_array(1), 'admbase::gxx' )
  call CCTK_VarIndex ( in_array(2), 'admbase::gxy' )
  call CCTK_VarIndex ( in_array(3), 'admbase::gxz' )
  call CCTK_VarIndex ( in_array(4), 'admbase::gyy' )
  call CCTK_VarIndex ( in_array(5), 'admbase::gyz' )
  call CCTK_VarIndex ( in_array(6), 'admbase::gzz' )

! If the static conformal factor is used, we also need to interpolate it.
  if ( CCTK_EQUALS ( metric_type, 'static conformal' ) ) then
    call CCTK_VarIndex ( in_array(7), 'staticconformal::psi' )

!   Call the interpolator for the metric and the conformal factor.
    call CCTK_InterpGridArrays ( ierr, cctkGH, 3, interp_handle, &
                                 table_handle, coord_system_handle, &
                                 lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                                 interp_coords, 7, in_array, &
                                 7, out_types, out_array )

!   Get the physical metric.
    gxxi = psii**4 * gxxi
    gxyi = psii**4 * gxyi
    gxzi = psii**4 * gxzi
    gyyi = psii**4 * gyyi
    gyzi = psii**4 * gyzi
    gzzi = psii**4 * gzzi
  else
!   Else call the interpolater only for the metric.
    call CCTK_InterpGridArrays ( ierr, cctkGH, 3, interp_handle, &
                                 table_handle, coord_system_handle, &
                                 lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                                 interp_coords, 6, in_array(1:6), &
                                 6, out_types(1:6), out_array(1:6) )
  end if

! For all points on the surface...
  do j = 1, lsh(2)
    do i = 1, lsh(1)

!     Calculate the coordinate transformation between the 3-d cartesian
!     coordinates and the 2-d (theta,phi) coordinates on the surface.
      dxdth = ( drdtheta(i,j) * sintheta(i,j) + &
                rsurf(i,j) * costheta(i,j) ) * cosphi(i,j)
      dydth = ( drdtheta(i,j) * sintheta(i,j) + &
                rsurf(i,j) * costheta(i,j) ) * sinphi(i,j)
      dzdth = drdtheta(i,j) * costheta(i,j) - rsurf(i,j) * sintheta(i,j)
      dxdph = ( drdphi(i,j) * cosphi(i,j) - &
                rsurf(i,j) * sinphi(i,j) ) * sintheta(i,j)
      dydph = ( drdphi(i,j) * sinphi(i,j) + &
                rsurf(i,j) * cosphi(i,j) ) * sintheta(i,j)
      dzdph = drdphi(i,j) * costheta(i,j)

!     Using this coordinate transformation, calculate the induced 2-metric
!     on the surface.
      gtt(i,j) = dxdth**2 * gxxi(i,j) + dydth**2 * gyyi(i,j) + &
            dzdth**2 * gzzi(i,j) + &
            two * ( dxdth * dydth * gxyi(i,j) + &
                    dxdth * dzdth * gxzi(i,j) + &
                    dydth * dzdth * gyzi(i,j) )
      gtp(i,j) = dxdth * dxdph * gxxi(i,j) + &
            ( dxdth * dydph + dydth * dxdph ) * gxyi(i,j) + &
            ( dxdth * dzdph + dzdth * dxdph ) * gxzi(i,j) + &
            dydth * dydph * gyyi(i,j) + &
            ( dydth * dzdph + dzdth * dydph ) * gyzi(i,j) + &
            dzdth * dzdph * gzzi(i,j)
      gpp(i,j) = dxdph**2 * gxxi(i,j) + dydph**2 * gyyi(i,j) + &
            dzdph**2 * gzzi(i,j) + &
            two * ( dxdph * dydph * gxyi(i,j) + &
                    dxdph * dzdph * gxzi(i,j) + &
                    dydph * dzdph * gyzi(i,j) )

!     Calculate the area element as the square-root of the determinant
!     of the two metric multiplied by dtheta*dphi.
      da(i,j) = dtheta * dphi * sqrt ( gtt(i,j) * gpp(i,j) - gtp(i,j)**2 )

!     Calculate the arc length element for the polar circumference integration.
      dltheta(i,j) = dtheta * sqrt ( gtt(i,j) )

!     Calculate the arc length element for the equatorial circumference
!     integration.
      dlphi(i,j) = dphi * sqrt ( gpp(i,j) )

    end do
  end do
 
end subroutine EHFinder_FindSurfaceElement


subroutine EHFinder_IntegrateArea(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: int_var, sn, l
  CCTK_REAL :: eh_area_tmp
  character(len=30) :: info_message
  CCTK_INT, dimension(1) :: lbnd, ubnd, lsh

! If finding of surface failed do not try to integrate but exit.
  if ( find_surface_status .lt. 0 ) then
    return
  endif

! Store the levelset_counter in a shorter named variable for convenience.
  l = levelset_counter

! Write an info message, to indicate what we are doing.
  call CCTK_INFO ( 'Integrating area' )

  sn = surface_counter - integrate_counter

! Taking into account the sym_factor and the weights, the area is simply
! The sum of the elements of this array.
  int_tmp = sym_factor * weights * da

! Get the index to the area integration array.
  call CCTK_VarIndex ( int_var, 'ehfinder::int_tmp' )
  if ( int_var .lt. 0 ) then
    call CCTK_WARN ( 0, 'Could not get index to grid array int_tmp' )
  end if

! Do a sum reduce over all processors. The result is the area.
  call CCTK_Reduce ( ierr, cctkGH, -1, sum_handle, 1, CCTK_VARIABLE_REAL, &
                     eh_area_tmp, 1, int_var )

! Write an info message reporting the area found.
  write(info_message,'(a15,f10.5)') 'Horizon area = ',eh_area_tmp

  call CCTK_INFO( trim(info_message) )

! Figure out the bounds of the area storage grid array.
  call CCTK_GrouplbndGN ( ierr, cctkGH, 1, lbnd, 'ehfinder::eh_area2' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get lower bounds for area array' )
  end if
  call CCTK_GroupubndGN ( ierr, cctkGH, 1, ubnd, 'ehfinder::eh_area2' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get upper bounds for area array' )
  end if
  call CCTK_GrouplshGN ( ierr, cctkGH, 1, lsh, 'ehfinder::eh_area2' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for area array' )
  end if

! If this processor contains the right part of the area storage grid array
! save the area, otherwise do nothing.
  if ( ( sn .ge. lbnd(1) + 1 ) .and. ( sn .le. ubnd(1) + 1 ) ) then
    eh_area2(sn-lbnd(1),l) = eh_area_tmp
  end if

end subroutine EHFinder_IntegrateArea


subroutine EHFinder_IntegrateCentroid(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: int_var, sn, k, l
  character(len=50) :: info_message
  CCTK_INT, dimension(1) :: lbnd, ubnd, lsh
  CCTK_REAL :: centroid_x, centroid_y, centroid_z

! If finding of surface failed do not try to integrate but exit.
  if ( find_surface_status .lt. 0 ) then
    return
  endif

! Store the levelset_counter in a shorter named variable for convenience.
  l = levelset_counter

! Write an info message, to indicate what we are doing.
  call CCTK_INFO ( 'Integrating centroid' )

  sn = surface_counter - integrate_counter

! Calculate the (x,y,z) coordinates of the points on the surface.
  interp_x = center(1) + rsurf * sintheta * cosphi
  interp_y = center(2) + rsurf * sintheta * sinphi
  interp_z = center(3) + rsurf * costheta

! Get the index to the temporary integration array.
  call CCTK_VarIndex ( int_var, 'ehfinder::int_tmp' )
  if ( int_var .lt. 0 ) then
    call CCTK_WARN ( 0, 'Could not get index to grid array int_tmp' )
  end if

! Take into account the symmetry, weights, x-coordinate and area element
! to calculate the integral of x*da over the surface. This gives 
! centroid_x * Area, so later we divide by the area.
  int_tmp = sym_factor * interp_x * weights * da

! Sum up the all elements of int_tmp.
  call CCTK_Reduce ( ierr, cctkGH, -1, sum_handle, 1, CCTK_VARIABLE_REAL, &
                     centroid_x, 1, int_var )

! Take into account the symmetry, weights, y-coordinate and area element
! to calculate the integral of y*da over the surface. This gives
! centroid_y * Area, so later we divide by the area.
  int_tmp = sym_factor * interp_y * weights * da

! Sum up the all elements of int_tmp.
  call CCTK_Reduce ( ierr, cctkGH, -1, sum_handle, 1, CCTK_VARIABLE_REAL, &
                     centroid_y, 1, int_var )

! Take into account the symmetry, weights, z-coordinate and area element
! to calculate the integral of z*da over the surface. This gives
! centroid_z * Area, so later we divide by the area.
  int_tmp = sym_factor * interp_z * weights * da

! Sum up the all elements of int_tmp.
  call CCTK_Reduce ( ierr, cctkGH, -1, sum_handle, 1, CCTK_VARIABLE_REAL, &
                     centroid_z, 1, int_var )
! If this surface has symmetries set the corresponding centroids to zero.
  if ( s_symx ) centroid_x = zero
  if ( s_symy ) centroid_y = zero
  if ( s_symz ) centroid_z = zero

! Get the bounds for the centroid arrays. It should be enough to do
! it just for the x-array, since they are all defined in the same way.
  call CCTK_GrouplbndGN ( ierr, cctkGH, 1, lbnd, 'ehfinder::eh_centroid2_x' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get lower bounds for centroid array' )
  end if
  call CCTK_GroupubndGN ( ierr, cctkGH, 1, ubnd, 'ehfinder::eh_centroid2_x' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get upper bounds for centroid array' )
  end if
  call CCTK_GrouplshGN ( ierr, cctkGH, 1, lsh, 'ehfinder::eh_centroid2_x' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for centroid array' )
  end if

! If we are on the right processor, store the centroid, otherwise do nothing.
  if ( ( sn .ge. lbnd(1) + 1 ) .and. ( sn .le. ubnd(1) + 1 ) ) then
    k = sn - lbnd(1)
    eh_centroid2_x(k,l) = centroid_x / eh_area2(k,l)
    eh_centroid2_y(k,l) = centroid_y / eh_area2(k,l)
    eh_centroid2_z(k,l) = centroid_z / eh_area2(k,l)
  end if

! Write an info message with the calculated centroids.
  write(info_message,'(a19,3(f10.5))') 'Horizon centroid = ', &
                                       eh_centroid2_x(sn,l), &
                                       eh_centroid2_y(sn,l), &
                                       eh_centroid2_z(sn,l)
  call CCTK_INFO( trim(info_message) )

end subroutine EHFinder_IntegrateCentroid


subroutine EHFinder_IntegrateCircumference(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: sn, k, l
  character(len=50) :: info_message
  CCTK_INT, dimension(1) :: lbnd1, ubnd1, lsh1
  CCTK_INT, dimension(2) :: lbnd, lsh, nghost
  CCTK_INT, dimension(4) :: bbox
  CCTK_INT :: ithl, ithr, jphl, jphr, itha, jpha
  CCTK_REAL :: eq_circ_loc, eq_circ, pol_circ_loc, pol_circ

! If finding of surface failed do not try to integrate but exit.
  if ( find_surface_status .lt. 0 ) then
    return
  endif

! Store the levelset_counter in a shorter named variable for convenience.
  l = levelset_counter

  sn = surface_counter - integrate_counter

  ! Find out the lower bounds of the distributed integration grid arrays.
  call CCTK_GrouplbndGN ( ierr, cctkGH, 2, lbnd, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get lower bounds for surface arrays' )
  end if

  ! Find out the local size of the distributed integration grid arrays
  call CCTK_GrouplshGN ( ierr, cctkGH, 2, lsh, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if

  ! Find out the bounding box of the distributed integration grid arrays
  call CCTK_GroupbboxGN ( ierr, cctkGH, 4, bbox, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get bounding box for surface arrays' )
  end if

  ! Find out the ghost size of the distributed integration grid arrays
  call CCTK_GroupnghostzonesGN ( ierr, cctkGH, 2, nghost, &
                                       'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get ghost zones for surface arrays' )
  end if
 
! Initialize the lower and upper local index for the theta and phi arrays.
  ithl = 1; ithr = lsh(1)
  jphl = 1; jphr = lsh(2)
  
! If the boundaries are internal, adjust the lower and upper local index
! to take into account the ghost cells.
  if ( bbox(1) .eq. 0 ) ithl = ithl + nghost(1)
  if ( bbox(2) .eq. 0 ) ithr = ithr - nghost(1)
  if ( bbox(3) .eq. 0 ) jphl = jphl + nghost(2)
  if ( bbox(4) .eq. 0 ) jphr = jphr - nghost(2)

! Write and info message about what we are doing.
  call CCTK_INFO ( 'Integrating equatorial circumference' )

! Setup the integration array for the equatorial circumference
! calculation. Note this is still a 2d-array, so we have to
! do the sum only for the correct 1d-slice of the array.
  int_tmp = phi_sym_factor * phiweights * dlphi

! If we are on a processor that contains part of the correct 1-d slice
! sum up the appropriate elemets, otherwise this processor will
! contribute zero.
  if ( ( ithl+lbnd(1) .le. githeta ) .and. ( githeta .le. ithr+lbnd(1) ) ) then
    itha = githeta - lbnd(1)
    eq_circ_loc = sum ( int_tmp(itha,jphl:jphr) )
  else
    eq_circ_loc = zero
  end if

! Then reduce the results over all processors and send the result to
! all processors.
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, sum_handle, &
                              eq_circ_loc, eq_circ, CCTK_VARIABLE_REAL )

! Get the upper and lower bounds for the circumference arrays.
  call CCTK_GrouplbndGN ( ierr, cctkGH, 1, lbnd1, 'ehfinder::eh_circ_eq2' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get lower bounds for area array' )
  end if
  call CCTK_GroupubndGN ( ierr, cctkGH, 1, ubnd1, 'ehfinder::eh_circ_eq2' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get upper bounds for area array' )
  end if
  call CCTK_GrouplshGN ( ierr, cctkGH, 1, lsh1, 'ehfinder::eh_circ_eq2' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for area array' )
  end if

! If we are on the correct processor store the result.
  if ( ( sn .ge. lbnd1(1) + 1 ) .and. ( sn .le. ubnd1(1) + 1 ) ) then
    k = sn - lbnd1(1)
    eh_circ_eq2(k,l) = eq_circ
  end if

! Write an info message with the calculated circunference.
  write(info_message,'(a35,f10.5)') 'Horizon equatorial circumference = ', &
                                    eq_circ
  call CCTK_INFO( trim(info_message) )

! Now do the same for the polar circomference.
  call CCTK_INFO ( 'Integrating polar circumference' )

! Setup the integration array for the polar circumference
! calculation. Note this is still a 2d-array, so we have to
! do the sum only for the correct 1d-slice of the array.
  int_tmp = theta_sym_factor * thetaweights * dltheta

! If we are on a processor that contains part of the correct 1-d slice
! sum up the appropriate elemets, otherwise this processor will
! contribute zero.
  if ( ( jphl+lbnd(2) .le. gjphi ) .and. ( gjphi .le. jphr+lbnd(2) ) ) then
    jpha = gjphi - lbnd(2)
    pol_circ_loc = sum ( int_tmp(ithl:ithr,jpha) )
  else
    pol_circ_loc = zero
  end if

! Then reduce the results over all processors and send the result to
! all processors.
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, sum_handle, &
                              pol_circ_loc, pol_circ, CCTK_VARIABLE_REAL )


! If we are on the correct processor store the result.
  if ( ( sn .ge. lbnd1(1) + 1 ) .and. ( sn .le. ubnd1(1) + 1 ) ) then
    k = sn - lbnd1(1)
    eh_circ_pol2(k,l) = pol_circ
  end if

! Write an info message with the calculated circumference.
  write(info_message,'(a30,f10.5)') 'Horizon polar circumference = ',pol_circ
  call CCTK_INFO( trim(info_message) )
  
end subroutine EHFinder_IntegrateCircumference


! This routine sets up the weights for the Simpsons rule integration
! over the surface.
subroutine EHFinder_InitWeights(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j

  CCTK_INT, dimension(2) :: lsh, lbnd

! Find out the lower bounds of the distributed integration grid arrays.
  call CCTK_GrouplbndGN ( ierr, cctkGH, 2, lbnd, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get lower bounds for surface arrays' )
  end if

! Find out the local size of the distributed integration grid arrays
  call CCTK_GrouplshGN ( ierr, cctkGH, 2, lsh, 'ehfinder::surface_arrays' )
  if ( ierr .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if

  weights = one

! Initialise the weight grid array for the 2D Simpsons rule integration.
! To do this I need to figure out the global location of the given point.
! There are 3 cases in the one dimensional case. If the point is on the
!  boundary the weight is 1/3. If it is at an even position the weight
! is 4/3 and if it is at an odd position the weight is 2/3. 

  do j = 1, lsh(2)

!   This is first done in the phi direction. Meaning that all points with
!   the same theta coordinate are set to the same weight.
    if ( ( lbnd(2)+j .eq. 1 ) .or. ( lbnd(2)+j .eq. nphi ) ) then
      weights(:,j) = onethird
      phiweights(:,j) = onethird
    else if ( mod(lbnd(2)+j,2) .eq. 0 ) then
      weights(:,j) = fourthirds
      phiweights(:,j) = fourthirds
    else
      weights(:,j) = twothirds
      phiweights(:,j) = twothirds
    end if

!   Then it is done in the theta direction with the one-dimensional
!   weights beeing multiplied.
    do i = 1, lsh(1)
      if ( ( lbnd(1)+i .eq. 1 ) .or. ( lbnd(1)+i .eq. ntheta ) ) then
        weights(i,j) = onethird * weights(i,j)
        thetaweights(i,j) = onethird
      else if ( mod(lbnd(1)+i,2) .eq. 0 ) then
        weights(i,j) = fourthirds * weights(i,j)
        thetaweights(i,j) = fourthirds
      else
        weights(i,j) = twothirds * weights(i,j)
        thetaweights(i,j) = twothirds
      end if
    end do
  end do

! The end result is a 2D array with the weights in the following pattern.

!  1/9   4/9  2/9   4/9  2/9   4/9  1/9
!  4/9  16/9  8/9  16/9  8/9  16/9  4/9
!  2/9   8/9  4/9   8/9  4/9   8/9  2/9
!  4/9  16/9  8/9  16/9  8/9  16/9  4/9
!  2/9   8/9  4/9   8/9  4/9   8/9  2/9
!  4/9  16/9  8/9  16/9  8/9  16/9  4/9
!  1/9   4/9  2/9   4/9  2/9   4/9  1/9
end subroutine EHFinder_InitWeights


! This is the routine where the stored areas (in eh_area2) are copied
! into the trigger variable eh_area.
subroutine EHFinder_CopyArea(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

! If finding of surface failed set area to zero.
  if ( find_surface_status .lt. 0 ) then
    eh_area = zero
    return
  else
    call CCTK_INFO('Copying area data')
    eh_area = eh_area2
  endif
end subroutine EHFinder_CopyArea


! This is the routine where the stored centroids (in eh_centroid2_[xyz])
! are copied into the trigger variable eh_centroid_[xyz].
subroutine EHFinder_CopyCentroid(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

! If finding of surface failed set centroids to zero.
  if ( find_surface_status .lt. 0 ) then
    eh_centroid_x = zero
    eh_centroid_y = zero
    eh_centroid_z = zero
    return
  else
    call CCTK_INFO('Copying centroid data')
    eh_centroid_x = eh_centroid2_x
    eh_centroid_y = eh_centroid2_y
    eh_centroid_z = eh_centroid2_z
  end if

end subroutine EHFinder_CopyCentroid


! This is the routine where the stored circumferences (in eh_circ_[eq,pol]2)
! are copied into the trigger variable eh_circ_[eq,pol].
subroutine EHFinder_CopyCircumference(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  if ( find_surface_status .lt. 0 ) then
    eh_circ_eq = zero
    eh_circ_pol = zero
    return
  else
    call CCTK_INFO('Copying circumference data')
    eh_circ_eq = eh_circ_eq2
    eh_circ_pol = eh_circ_pol2
  end if

end subroutine EHFinder_CopyCircumference
