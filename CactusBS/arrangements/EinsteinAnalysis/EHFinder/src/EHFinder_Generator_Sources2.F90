! Calculation of the sources for the level set function.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_Generator_Sources2(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  CCTK_INT :: interp_handle, table_handle, status, coord_system_handle

  character(len=200) :: gen_interp
  character(len=128) :: warn_message
  CCTK_INT :: gen_interp_len
  character(len=7) :: gen_order

  CCTK_INT, dimension(1) :: lsh
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(3) :: out_arrays
  CCTK_INT, dimension(3) :: in_arrays
  CCTK_INT, dimension(3), parameter :: op_indices = (/ 0, 1, 2 /), &
                                        op_codes = (/ 0, 0, 0 /)
  CCTK_INT, dimension(3) :: out_types
  CCTK_REAL :: alp2, psi4, dfux, dfuy, dfuz, factor, ssign
  CCTK_REAL :: idetg, guxx, guxy, guxz, guyy, guyz, guzz

! Set the sign depending on the surface direction.
  if ( CCTK_EQUALS ( surface_direction, 'outward' ) ) ssign = one
  if ( CCTK_EQUALS ( surface_direction, 'inward' ) ) ssign = -one

  out_types = CCTK_VARIABLE_REAL

! Convert the generator_interpolator string parameter to a Fortran string.
  call CCTK_FortranString ( gen_interp_len, generator_interpolator, &
                                            gen_interp )

! Get the corresponding interpolator handle.
  call CCTK_InterpHandle ( interp_handle, gen_interp )

  if ( interp_handle .lt. 0 ) then
    warn_message = 'Cannot get handle for interpolation. '
    warn_message = trim(warn_message)//'Forgot to activate an implementation '
    warn_message = trim(warn_message)//'providing interpolation operators?'
    call CCTK_WARN( 0, trim(warn_message) )
  end if

! Convert the interpolation order parameter to a Fortran string to be placed
! in the interpolator table. Note that the order is assumed to contain only
! 1 digit.
  write(gen_order,'(a6,i1)') 'order=',generator_interpolation_order

! Create the table directly from the string.
  call Util_TableCreateFromString ( table_handle, gen_order )
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

#include "include/physical_part.h"

! Find out how many interpolation points are located on this processor.
  call CCTK_GrouplshGN ( status, cctkGH, 1, lsh, 'ehfinder::xg' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if

! Set the indices to the input grid functions.
  call CCTK_VarIndex ( in_arrays(1), 'ehfinder::xgf' )
  call CCTK_VarIndex ( in_arrays(2), 'ehfinder::ygf' )
  call CCTK_VarIndex ( in_arrays(3), 'ehfinder::zgf' )

! Set the operand indices table entry, corresponding
! to interpolation of ehfinder::generator_gf (3)
  call Util_TableSetIntArray ( status, table_handle, 3, &
                               op_indices, 'operand_indices' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'Cannot set operand indices array in parameter table' )
  endif

! Set the corresponding table entry for the operation codes.
  call Util_TableSetIntArray ( status, table_handle, 3, &
                               op_codes, 'operation_codes' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'Cannot set operation codes array in parameter table' )
  endif

! Loop over the level sets
  do l = 1, eh_number_level_sets

!   Set the pointers to the points to be interpolated to.
    interp_coords(1) = CCTK_PointerTo(xg(:,l))
    interp_coords(2) = CCTK_PointerTo(yg(:,l))
    interp_coords(3) = CCTK_PointerTo(zg(:,l))

!   Set the pointers to the output arrays.
    out_arrays(1) = CCTK_PointerTo(dxg(:,l))
    out_arrays(2) = CCTK_PointerTo(dyg(:,l))
    out_arrays(3) = CCTK_PointerTo(dzg(:,l))

!   Check the metric type. At present physical and static_conformal are
!   supported.
    if ( CCTK_EQUALS ( metric_type, 'physical' ) ) then

      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
            if ( eh_mask(i,j,k,l) .ge. 0 ) then

!             calculate the square of the lapse.
              alp2 = alp(i,j,k)**2

!             Calculate the inverse of the 3-metric.
              guxx = gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k)**2
              guxy = gxz(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gzz(i,j,k)
              guxz = gxy(i,j,k) * gyz(i,j,k) - gxz(i,j,k) * gyy(i,j,k)

              idetg = one / ( gxx(i,j,k) * guxx + &
                              gxy(i,j,k) * guxy + &
                              gxz(i,j,k) * guxz )

              guxx = idetg * guxx
              guxy = idetg * guxy
              guxz = idetg * guxz

              guyy = ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k)**2 ) * idetg
              guyz = ( gxy(i,j,k) * gxz(i,j,k) - &
                       gxx(i,j,k) * gyz(i,j,k) ) * idetg
              guzz = ( gxx(i,j,k) * gyy(i,j,k) - &
                       gxy(i,j,k)**2 ) * idetg
        
!             Raise the index of the partial derivatives of f.
              dfux = guxx * dfx(i,j,k,l) + guxy * dfy(i,j,k,l) + &
                                           guxz * dfz(i,j,k,l)
              dfuy = guxy * dfx(i,j,k,l) + guyy * dfy(i,j,k,l) + &
                                           guyz * dfz(i,j,k,l)
              dfuz = guxz * dfx(i,j,k,l) + guyz * dfy(i,j,k,l) + &
                                           guzz * dfz(i,j,k,l)

!             Calculate the overall multiplication factor.
              factor = alp2 / sqrt ( alp2 * ( dfux * dfx(i,j,k,l) + &
                                              dfuy * dfy(i,j,k,l) + &
                                              dfuz * dfz(i,j,k,l) ) )

!             Finally obtain dx^i/dt.
              xgf(i,j,k) = - betax(i,j,k) + ssign * factor * dfux
              ygf(i,j,k) = - betay(i,j,k) + ssign * factor * dfuy
              zgf(i,j,k) = - betaz(i,j,k) + ssign * factor * dfuz
            else
              xgf(i,j,k) = zero
              ygf(i,j,k) = zero
              zgf(i,j,k) = zero
            end if
          end do
        end do
      end do

    else if ( CCTK_EQUALS ( metric_type, 'static conformal' ) ) then

      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
            if ( eh_mask(i,j,k,l) .ge. 0 ) then
              alp2 = alp(i,j,k)**2

!             The inverse of psi^4
              psi4 = one / psi(i,j,k)**4

              guxx = gyy(i,j,k) * gzz(i,j,k) - gyz(i,j,k)**2
              guxy = gxz(i,j,k) * gyz(i,j,k) - gxy(i,j,k) * gzz(i,j,k)
              guxz = gxy(i,j,k) * gyz(i,j,k) - gxz(i,j,k) * gyy(i,j,k)

!             The inverse of the determinant divided by psi^4.
              idetg = psi4 / ( gxx(i,j,k) * guxx + &
                               gxy(i,j,k) * guxy + &
                               gxz(i,j,k) * guxz )

!             The inverse metric. Since the determinant is already divided
!             by psi^4, this gives the inverse of the physical metric.
              guxx = idetg * guxx
              guxy = idetg * guxy
              guxz = idetg * guxz

              guyy = ( gxx(i,j,k) * gzz(i,j,k) - gxz(i,j,k)**2 ) * idetg
              guyz = ( gxy(i,j,k) * gxz(i,j,k) - &
                       gxx(i,j,k) * gyz(i,j,k) ) * idetg
              guzz = ( gxx(i,j,k) * gyy(i,j,k) - gxy(i,j,k)**2 ) * idetg
      
              dfux = guxx * dfx(i,j,k,l) + guxy * dfy(i,j,k,l) + &
                                           guxz * dfz(i,j,k,l)
              dfuy = guxy * dfx(i,j,k,l) + guyy * dfy(i,j,k,l) + &
                                           guyz * dfz(i,j,k,l)
              dfuz = guxz * dfx(i,j,k,l) + guyz * dfy(i,j,k,l) + &
                                           guzz * dfz(i,j,k,l)
              factor = alp2 / sqrt ( alp2 * ( dfux * dfx(i,j,k,l) + &
                                              dfuy * dfy(i,j,k,l) + &
                                              dfuz * dfz(i,j,k,l) ) )
              xgf(i,j,k) = - betax(i,j,k) + ssign * factor * dfux
              ygf(i,j,k) = - betay(i,j,k) + ssign * factor * dfuy
              zgf(i,j,k) = - betaz(i,j,k) + ssign * factor * dfuz
            else
              xgf(i,j,k) = zero
              ygf(i,j,k) = zero
              zgf(i,j,k) = zero
            end if
          end do
        end do
      end do
    end if

!   Call the interpolator.
    call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
                                 table_handle, coord_system_handle, &
                                 lsh(1), CCTK_VARIABLE_REAL, &
                                 interp_coords, 3, in_arrays, &
                                 3, out_types, out_arrays )

    if ( status .lt. 0 ) then
      call CCTK_INFO ( 'Interpolation failed.' )
    end if

  end do

  return
end subroutine EHFinder_Generator_Sources2
