#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

      subroutine DemoInterp_Interp1D(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      integer i
      integer ierror
      CCTK_REAL dx

      integer interp_handle, coord_system_handle, param_table_handle
      character(30) options_string
      CCTK_POINTER, dimension(1) :: interp_coords
      CCTK_INT, dimension(2) :: in_array_indices
      CCTK_POINTER, dimension(2) :: out_arrays
      CCTK_INT, dimension(2) :: out_array_type_codes


c     Set up coordinates to interpolate at (along line)
      dx = (2.0*interp_radius)/real(arrsize)
      do i=1,arrsize
         xcoord(i) = -interp_radius+(i-1)*dx
      end do

c     First we get the handles to our interpolator operator
c     and the coordinate system
      param_table_handle = -1
      interp_handle = -1
      coord_system_handle = -1

      options_string = "order = " // char(ichar('0') + interpolation_order)
      call Util_TableCreateFromString (param_table_handle, options_string)
      if (param_table_handle .lt. 0) then
        call CCTK_WARN(0,"Cannot create parameter table for interpolator")
      endif

      call CCTK_InterpHandle (interp_handle, "uniform cartesian")
      if (interp_handle < 0) then
        call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
      end if

      call CCTK_CoordSystemHandle (coord_system_handle, "cart1d")
      if (coord_system_handle < 0) then
         call CCTK_WARN(0,"Coordinate system 'cart1d' not registered")
      end if

c     fill in the input/output arrays for the interpolator
      interp_coords(1) = CCTK_PointerTo(xcoord)

      call CCTK_VarIndex(i, "DemoInterp::realgf1")
      in_array_indices(1) = i
      call CCTK_VarIndex(i, "DemoInterp::compgf1")
      in_array_indices(2) = i

      out_array_type_codes(1) = CCTK_VARIABLE_REAL
      out_array_type_codes(2) = CCTK_VARIABLE_COMPLEX

      out_arrays(1) = CCTK_PointerTo(realinterp1)
      out_arrays(2) = CCTK_PointerTo(compinterp1)

c     ==================================================================
c     Interpolation of both real and complex GFs using general call
c     ==================================================================

      call CCTK_InterpGridArrays (ierror, cctkGH, 1, interp_handle,
     .                            param_table_handle, coord_system_handle,
     .                            arrsize, CCTK_VARIABLE_REAL, interp_coords,
     .                            2, in_array_indices,
     .                            2, out_array_type_codes, out_arrays)


      if (ierror<0) then
         call CCTK_WARN(1,"Interpolation error")
      endif

      print *
      print *,"Real 1D interpolation using CCTK_InterpGridArrays()"
      print *,"==================================================="
      print *,"Interpolated values should be : "
      do i=1,arrsize
         print *,xcoord(i)
      end do
      print *,"Actual values are : "
      do i=1,arrsize
         print *,realinterp1(i)
      end do

      print *
      print *,"Complex 1D interpolation using CCTK_InterpGridArrays()"
      print *,"======================================================"
      print *,"Interpolated values should be : "
      do i=1,arrsize
         print *,xcoord(i),xcoord(i)
      end do
      print *,"Actual values are : "
      do i=1,arrsize
         print *,compinterp1(i)
      end do

      return
      end
