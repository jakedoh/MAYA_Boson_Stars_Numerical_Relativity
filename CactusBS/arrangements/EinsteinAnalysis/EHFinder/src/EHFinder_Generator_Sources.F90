! Calculation of the sources for the level set function.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_Generator_Sources(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, l
  CCTK_INT :: interp_handle, table_handle, status, coord_system_handle

  character(len=200) :: gen_interp
  character(len=128) :: warn_message
  CCTK_INT :: gen_interp_len
  character(len=7) :: gen_order
  character(len=15) :: vname

  CCTK_INT, dimension(1) :: lsh
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(14) :: out_arrays
  CCTK_INT, dimension(12) :: in_arrays
  CCTK_INT, dimension(14), parameter :: op_indices = (/ 0, 1, 2, 3, 4, &
                                                        5, 6, 7, 8, 9, &
                                                        10, 10, 10, 11 /), &
                                        op_codes = (/ 0, 0, 0, 0, 0, &
                                                      0, 0, 0, 0, 0, &
                                                      1, 2, 3, 0 /)
  CCTK_INT, dimension(14) :: out_types
  CCTK_REAL :: alp2, psi4, dfux, dfuy, dfuz, factor, ssign
  CCTK_REAL :: idetg, guxx, guxy, guxz, guyy, guyz, guzz

! Set the sign correctly depending on the surface direction.
  if ( CCTK_EQUALS ( surface_direction, 'outward' ) ) ssign = one
  if ( CCTK_EQUALS ( surface_direction, 'inward' ) ) ssign = -one

  out_types = CCTK_VARIABLE_REAL

! Convert the generator_interpolator string parameter to a Fortran string.
  call CCTK_FortranString ( gen_interp_len, generator_interpolator, &
                                            gen_interp )

! Get the corresponding interpolator handle.
  call CCTK_InterpHandle ( interp_handle, gen_interp )

  if ( interp_handle .lt. 0 ) then
    warn_message = 'Cannot get handle for interpolation.'
    warn_message = trim(warn_message)//' Forgot to activate an implementation'
    warn_message = trim(warn_message)//' providing interpolation operators?'
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
    warn_message = 'Cannot get handle for cart3d coordinate system.'
    warn_message = trim(warn_message)//' Forgot to activate an implementation'
    warn_message = trim(warn_message)//' providing coordinates?'
    call CCTK_WARN( 0, trim(warn_message) )
  endif

! Find out how many interpolation points are located on this processor.
  call CCTK_GrouplshGN ( status, cctkGH, 1, lsh, 'ehfinder::xg' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if

! Set the pointers to the output arrays.
  out_arrays(1) = CCTK_PointerTo(alpg)
  out_arrays(2) = CCTK_PointerTo(betaxg)
  out_arrays(3) = CCTK_PointerTo(betayg)
  out_arrays(4) = CCTK_PointerTo(betazg)
  out_arrays(5) = CCTK_PointerTo(gxxg)
  out_arrays(6) = CCTK_PointerTo(gxyg)
  out_arrays(7) = CCTK_PointerTo(gxzg)
  out_arrays(8) = CCTK_PointerTo(gyyg)
  out_arrays(9) = CCTK_PointerTo(gyzg)
  out_arrays(10) = CCTK_PointerTo(gzzg)
  out_arrays(11) = CCTK_PointerTo(dfxg)
  out_arrays(12) = CCTK_PointerTo(dfyg)
  out_arrays(13) = CCTK_PointerTo(dfzg)
  out_arrays(14) = CCTK_PointerTo(psig)

! Loop over the level sets
  do l = 1, eh_number_level_sets

!   Set the pointers to the points to be interpolated to.
    interp_coords(1) = CCTK_PointerTo(xg(:,l))
    interp_coords(2) = CCTK_PointerTo(yg(:,l))
    interp_coords(3) = CCTK_PointerTo(zg(:,l))

!   Set the indices to the input grid functions.
    call CCTK_VarIndex ( in_arrays(1), 'admbase::alp' )
    call CCTK_VarIndex ( in_arrays(2), 'admbase::betax' )
    call CCTK_VarIndex ( in_arrays(3), 'admbase::betay' )
    call CCTK_VarIndex ( in_arrays(4), 'admbase::betaz' )
    call CCTK_VarIndex ( in_arrays(5), 'admbase::gxx' )
    call CCTK_VarIndex ( in_arrays(6), 'admbase::gxy' )
    call CCTK_VarIndex ( in_arrays(7), 'admbase::gxz' )
    call CCTK_VarIndex ( in_arrays(8), 'admbase::gyy' )
    call CCTK_VarIndex ( in_arrays(9), 'admbase::gyz' )
    call CCTK_VarIndex ( in_arrays(10), 'admbase::gzz' )
    write(vname,'(a12,i1,a1)') 'ehfinder::f[', l-1,']'
    call CCTK_VarIndex ( in_arrays(11), vname )
    call CCTK_VarIndex ( in_arrays(12), 'staticconformal::psi' )

!   Check the metric type. At present physical and static_conformal are
!   supported.
    if ( CCTK_EQUALS ( metric_type, 'physical' ) ) then

!     Set the operand indices table entry, corresponding
!     to interpolation of admbase::alp (1), admbase::shift (3),
!     admbase::metric(6) and the first derivatives of ehfinder::f (3) for
!     a total of 13 output arrays.
      call Util_TableSetIntArray ( status, table_handle, 13, &
                                 op_indices(1:13), 'operand_indices' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operand indices array in parameter table' 
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

!     Set the corresponding table entry for the operation codes.
      call Util_TableSetIntArray ( status, table_handle, 13, &
                                 op_codes(1:13), 'operation_codes' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operation codes array in parameter table'
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

!     Call the interpolator.
      call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
                                   table_handle, coord_system_handle, &
                                   lsh(1), CCTK_VARIABLE_REAL, &
                                   interp_coords, 11, in_arrays(1:11), &
                                   13, out_types(1:13), out_arrays(1:13) )

      if ( status .lt. 0 ) then
        call CCTK_INFO ( 'Interpolation failed.' )
      end if

!     For each point on this processor calculate the right hand side of the
!     characteristic evolution equation.
      do i = 1, lsh(1)

!       calculate the square of the lapse.
        alp2 = alpg(i)**2

!       Calculate the inverse of the 3-metric.
        guxx = gyyg(i) * gzzg(i) - gyzg(i)**2
        guxy = gxzg(i) * gyzg(i) - gxyg(i) * gzzg(i)
        guxz = gxyg(i) * gyzg(i) - gxzg(i) * gyyg(i)

        idetg = one / ( gxxg(i) * guxx + gxyg(i) * guxy + gxzg(i) * guxz )

        guxx = idetg * guxx
        guxy = idetg * guxy
        guxz = idetg * guxz

        guyy = ( gxxg(i) * gzzg(i) - gxzg(i)**2 ) * idetg
        guyz = ( gxyg(i) * gxzg(i) - gxxg(i) * gyzg(i) ) * idetg
        guzz = ( gxxg(i) * gyyg(i) - gxyg(i)**2 ) * idetg
      
!       Raise the index of the partial derivatives of f.
        dfux = guxx * dfxg(i) + guxy * dfyg(i) + guxz * dfzg(i)
        dfuy = guxy * dfxg(i) + guyy * dfyg(i) + guyz * dfzg(i)
        dfuz = guxz * dfxg(i) + guyz * dfyg(i) + guzz * dfzg(i)

!       Calculate the overall multiplication factor.
        factor = alp2 / sqrt ( alp2 * ( dfux * dfxg(i) + &
                                        dfuy * dfyg(i) + &
                                        dfuz * dfzg(i) ) )

!       Finally obtain dx^i/dt.
        dxg(i,l) = - betaxg(i) + ssign * factor * dfux
        dyg(i,l) = - betayg(i) + ssign * factor * dfuy
        dzg(i,l) = - betazg(i) + ssign * factor * dfuz

      end do
    else if ( CCTK_EQUALS ( metric_type, 'static conformal' ) ) then
      call Util_TableSetIntArray ( status, table_handle, 14, &
                                 op_indices, 'operand_indices' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operand indices array in parameter table'
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

      call Util_TableSetIntArray ( status, table_handle, 14, &
                                 op_codes, 'operation_codes' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operation codes array in parameter table'
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

      call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
                                   table_handle, coord_system_handle, &
                                   lsh(1), CCTK_VARIABLE_REAL, &
                                   interp_coords, 12, in_arrays, &
                                   14, out_types, out_arrays )

      if ( status .lt. 0 ) then
        call CCTK_INFO ( 'Interpolation failed.' )
      end if

      do i = 1, lsh(1)
        alp2 = alpg(i)**2

!       The inverse of psi^4
        psi4 = one / psig(i)**4

        guxx = gyyg(i) * gzzg(i) - gyzg(i)**2
        guxy = gxzg(i) * gyzg(i) - gxyg(i) * gzzg(i)
        guxz = gxyg(i) * gyzg(i) - gxzg(i) * gyyg(i)

!       The inverse of the determinant divided by psi^4.
        idetg = psi4 / ( gxxg(i) * guxx + &
                         gxyg(i) * guxy + &
                         gxzg(i) * guxz )

!       The inverse metric. Since the determinant is already divided
!       by psi^4, this gives the inverse of the physical metric.
        guxx = idetg * guxx
        guxy = idetg * guxy
        guxz = idetg * guxz

        guyy = ( gxxg(i) * gzzg(i) - gxzg(i)**2 ) * idetg
        guyz = ( gxyg(i) * gxzg(i) - gxxg(i) * gyzg(i) ) * idetg
        guzz = ( gxxg(i) * gyyg(i) - gxyg(i)**2 ) * idetg
      
        dfux = guxx * dfxg(i) + guxy * dfyg(i) + guxz * dfzg(i)
        dfuy = guxy * dfxg(i) + guyy * dfyg(i) + guyz * dfzg(i)
        dfuz = guxz * dfxg(i) + guyz * dfyg(i) + guzz * dfzg(i)
        factor = alp2 / sqrt ( alp2 * ( dfux * dfxg(i) + &
                                        dfuy * dfyg(i) + &
                                        dfuz * dfzg(i) ) )
        dxg(i,l) = - betaxg(i) + ssign * factor * dfux
        dyg(i,l) = - betayg(i) + ssign * factor * dfuy
        dzg(i,l) = - betazg(i) + ssign * factor * dfuz
      end do
    end if
  end do

  return
end subroutine EHFinder_Generator_Sources


subroutine EHFinder_Generator_Sources_2D(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, l
  CCTK_INT :: interp_handle, table_handle, status, coord_system_handle

  character(len=200) :: gen_interp
  character(len=128) :: warn_message
  CCTK_INT :: gen_interp_len
  character(len=7) :: gen_order
  character(len=15) :: vname

  CCTK_INT, dimension(2) :: lsh
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(14) :: out_arrays
  CCTK_INT, dimension(12) :: in_arrays
  CCTK_INT, dimension(14), parameter :: op_indices = (/ 0, 1, 2, 3, 4, &
                                                        5, 6, 7, 8, 9, &
                                                        10, 10, 10, 11 /), &
                                        op_codes = (/ 0, 0, 0, 0, 0, &
                                                      0, 0, 0, 0, 0, &
                                                      1, 2, 3, 0 /)
  CCTK_INT, dimension(14) :: out_types
  CCTK_REAL :: alp2, psi4, dfux, dfuy, dfuz, factor, ssign
  CCTK_REAL :: idetg, guxx, guxy, guxz, guyy, guyz, guzz

! Set the sign correctly depending on the surface direction.
  if ( CCTK_EQUALS ( surface_direction, 'outward' ) ) ssign = one
  if ( CCTK_EQUALS ( surface_direction, 'inward' ) ) ssign = -one

  out_types = CCTK_VARIABLE_REAL

! Convert the generator_interpolator string parameter to a Fortran string.
  call CCTK_FortranString ( gen_interp_len, generator_interpolator, &
                                            gen_interp )

! Get the corresponding interpolator handle.
  call CCTK_InterpHandle ( interp_handle, gen_interp )

  if ( interp_handle .lt. 0 ) then
    warn_message = 'Cannot get handle for interpolation.'
    warn_message = trim(warn_message)//' Forgot to activate an implementation'
    warn_message = trim(warn_message)//' providing interpolation operators?'
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
    warn_message = 'Cannot get handle for cart3d coordinate system.'
    warn_message = trim(warn_message)//' Forgot to activate an implementation'
    warn_message = trim(warn_message)//' providing coordinates?'
    call CCTK_WARN( 0, trim(warn_message) )
  endif

! Find out how many interpolation points are located on this processor.
  call CCTK_GrouplshGN ( status, cctkGH, 2, lsh, 'ehfinder::xg2' )
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, 'cannot get local size for surface arrays' )
  end if

  out_arrays(1) = CCTK_PointerTo(alpg2)
  out_arrays(2) = CCTK_PointerTo(betaxg2)
  out_arrays(3) = CCTK_PointerTo(betayg2)
  out_arrays(4) = CCTK_PointerTo(betazg2)
  out_arrays(5) = CCTK_PointerTo(gxxg2)
  out_arrays(6) = CCTK_PointerTo(gxyg2)
  out_arrays(7) = CCTK_PointerTo(gxzg2)
  out_arrays(8) = CCTK_PointerTo(gyyg2)
  out_arrays(9) = CCTK_PointerTo(gyzg2)
  out_arrays(10) = CCTK_PointerTo(gzzg2)
  out_arrays(11) = CCTK_PointerTo(dfxg2)
  out_arrays(12) = CCTK_PointerTo(dfyg2)
  out_arrays(13) = CCTK_PointerTo(dfzg2)
  out_arrays(14) = CCTK_PointerTo(psig2)

! Set the pointers to the output arrays.

! Loop over the level sets
  do l = 1, eh_number_level_sets

!   Set the pointers to the points to be interpolated to.

    interp_coords(1) = CCTK_PointerTo(xg2(1,1,l))
    interp_coords(2) = CCTK_PointerTo(yg2(1,1,l))
    interp_coords(3) = CCTK_PointerTo(zg2(1,1,l))

!   Set the indices to the input grid functions.
    call CCTK_VarIndex ( in_arrays(1), 'admbase::alp' )
    call CCTK_VarIndex ( in_arrays(2), 'admbase::betax' )
    call CCTK_VarIndex ( in_arrays(3), 'admbase::betay' )
    call CCTK_VarIndex ( in_arrays(4), 'admbase::betaz' )
    call CCTK_VarIndex ( in_arrays(5), 'admbase::gxx' )
    call CCTK_VarIndex ( in_arrays(6), 'admbase::gxy' )
    call CCTK_VarIndex ( in_arrays(7), 'admbase::gxz' )
    call CCTK_VarIndex ( in_arrays(8), 'admbase::gyy' )
    call CCTK_VarIndex ( in_arrays(9), 'admbase::gyz' )
    call CCTK_VarIndex ( in_arrays(10), 'admbase::gzz' )
    write(vname,'(a12,i1,a1)') 'ehfinder::f[', l-1,']'
    call CCTK_VarIndex ( in_arrays(11), vname )
    call CCTK_VarIndex ( in_arrays(12), 'staticconformal::psi' )

!   Check the metric type. At present physical and static_conformal are
!   supported.
    if ( CCTK_EQUALS ( metric_type, 'physical' ) ) then

!     Set the operand indices table entry, corresponding
!     to interpolation of admbase::alp (1), admbase::shift (3),
!     admbase::metric(6) and the first derivatives of ehfinder::f (3) for
!     a total of 13 output arrays.
      call Util_TableSetIntArray ( status, table_handle, 13, &
                                 op_indices(1:13), 'operand_indices' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operand indices array in parameter table' 
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

!     Set the corresponding table entry for the operation codes.
      call Util_TableSetIntArray ( status, table_handle, 13, &
                                 op_codes(1:13), 'operation_codes' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operation codes array in parameter table'
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

!     Call the interpolator.
      call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
                                   table_handle, coord_system_handle, &
                                   lsh(1)*lsh(2), CCTK_VARIABLE_REAL, &
                                   interp_coords, 11, in_arrays(1:11), &
                                   13, out_types(1:13), out_arrays(1:13) )

      if ( status .lt. 0 ) then
        call CCTK_INFO ( 'Interpolation failed.' )
      end if

!     For each point on this processor calculate the right hand side of the
!     characteristic evolution equation.
      do j = 1, lsh(2)
        do i = 1, lsh(1)

!         calculate the square of the lapse.
          alp2 = alpg2(i,j)**2

!         Calculate the inverse of the 3-metric.
          guxx = gyyg2(i,j) * gzzg2(i,j) - gyzg2(i,j)**2
          guxy = gxzg2(i,j) * gyzg2(i,j) - gxyg2(i,j) * gzzg2(i,j)
          guxz = gxyg2(i,j) * gyzg2(i,j) - gxzg2(i,j) * gyyg2(i,j)

          idetg = one / ( gxxg2(i,j) * guxx + gxyg2(i,j) * guxy &
                                            + gxzg2(i,j) * guxz )

          guxx = idetg * guxx
          guxy = idetg * guxy
          guxz = idetg * guxz

          guyy = ( gxxg2(i,j) * gzzg2(i,j) - gxzg2(i,j)**2 ) * idetg
          guyz = ( gxyg2(i,j) * gxzg2(i,j) - gxxg2(i,j) * gyzg2(i,j) ) * idetg
          guzz = ( gxxg2(i,j) * gyyg2(i,j) - gxyg2(i,j)**2 ) * idetg
      
!         Raise the index of the partial derivatives of f.
          dfux = guxx * dfxg2(i,j) + guxy * dfyg2(i,j) + guxz * dfzg2(i,j)
          dfuy = guxy * dfxg2(i,j) + guyy * dfyg2(i,j) + guyz * dfzg2(i,j)
          dfuz = guxz * dfxg2(i,j) + guyz * dfyg2(i,j) + guzz * dfzg2(i,j)

!         Calculate the overall multiplication factor.
          factor = alp2 / sqrt ( alp2 * ( dfux * dfxg2(i,j) + &
                                          dfuy * dfyg2(i,j) + &
                                          dfuz * dfzg2(i,j) ) )

!         Finally obtain dx^i/dt.
          dxg2(i,j,l) = - betaxg2(i,j) + ssign * factor * dfux
          dyg2(i,j,l) = - betayg2(i,j) + ssign * factor * dfuy
          dzg2(i,j,l) = - betazg2(i,j) + ssign * factor * dfuz

        end do
      end do
    else if ( CCTK_EQUALS ( metric_type, 'static conformal' ) ) then
      call Util_TableSetIntArray ( status, table_handle, 14, &
                                 op_indices, 'operand_indices' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operand indices array in parameter table'
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

      call Util_TableSetIntArray ( status, table_handle, 14, &
                                 op_codes, 'operation_codes' )
      if ( status .lt. 0 ) then
        warn_message = 'Cannot set operation codes array in parameter table'
        call CCTK_WARN ( 0, trim(warn_message) )
      endif

      call CCTK_InterpGridArrays ( status, cctkGH, 3, interp_handle, &
                                   table_handle, coord_system_handle, &
                                   lsh(1)*lsh(2), CCTK_VARIABLE_REAL, &
                                   interp_coords, 12, in_arrays, &
                                   14, out_types, out_arrays )

      if ( status .lt. 0 ) then
        call CCTK_INFO ( 'Interpolation failed.' )
      end if

      do j = 1, lsh(2)
        do i = 1, lsh(1)
          alp2 = alpg2(i,j)**2

!         The inverse of psi^4
          psi4 = one / psig2(i,j)**4

          guxx = gyyg2(i,j) * gzzg2(i,j) - gyzg2(i,j)**2
          guxy = gxzg2(i,j) * gyzg2(i,j) - gxyg2(i,j) * gzzg2(i,j)
          guxz = gxyg2(i,j) * gyzg2(i,j) - gxzg2(i,j) * gyyg2(i,j)

!         The inverse of the determinant divided by psi^4.
          idetg = psi4 / ( gxxg2(i,j) * guxx + &
                           gxyg2(i,j) * guxy + &
                           gxzg2(i,j) * guxz )

!         The inverse metric. Since the determinant is already divided
!         by psi^4, this gives the inverse of the physical metric.
          guxx = idetg * guxx
          guxy = idetg * guxy
          guxz = idetg * guxz

          guyy = ( gxxg2(i,j) * gzzg2(i,j) - gxzg2(i,j)**2 ) * idetg
          guyz = ( gxyg2(i,j) * gxzg2(i,j) - gxxg2(i,j) * gyzg2(i,j) ) * idetg
          guzz = ( gxxg2(i,j) * gyyg2(i,j) - gxyg2(i,j)**2 ) * idetg

          dfux = guxx * dfxg2(i,j) + guxy * dfyg2(i,j) + guxz * dfzg2(i,j)
          dfuy = guxy * dfxg2(i,j) + guyy * dfyg2(i,j) + guyz * dfzg2(i,j)
          dfuz = guxz * dfxg2(i,j) + guyz * dfyg2(i,j) + guzz * dfzg2(i,j)
          factor = alp2 / sqrt ( alp2 * ( dfux * dfxg2(i,j) + &
                                          dfuy * dfyg2(i,j) + &
                                          dfuz * dfzg2(i,j) ) )
          dxg2(i,j,l) = - betaxg2(i,j) + ssign * factor * dfux
          dyg2(i,j,l) = - betayg2(i,j) + ssign * factor * dfuy
          dzg2(i,j,l) = - betazg2(i,j) + ssign * factor * dfuz

        end do
      end do
    end if
  end do

  return
end subroutine EHFinder_Generator_Sources_2D
