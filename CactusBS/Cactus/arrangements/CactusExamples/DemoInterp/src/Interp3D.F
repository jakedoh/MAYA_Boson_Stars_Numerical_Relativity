#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

      subroutine DemoInterp_Interp3D(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      integer three,one
      integer i
      integer ierror
      CCTK_REAL dth,phi
      CCTK_REAL originx,originy,originz

      integer interp_handle, coord_system_handle, param_table_handle
      character(30) options_string
      CCTK_POINTER, dimension(3) :: interp_coords
      CCTK_INT, dimension(2) :: in_array_indices
      CCTK_POINTER, dimension(2) :: out_arrays
      CCTK_INT, dimension(2) :: out_array_type_codes


c     Set up coordinates to interpolate at (circle at z=0)
      dth = 2.0*atan(1.0)/real(arrsize)
      phi = 0.0
      do i=1,arrsize
         xcoord(i) = interp_radius*sin(i*dth)*cos(phi)
         ycoord(i) = interp_radius*sin(i*dth)*sin(phi)
         zcoord(i) = interp_radius*cos(i*dth)
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

      call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
      if (coord_system_handle < 0) then
         call CCTK_WARN(0,"Coordinate system 'cart3d' not registered")
      end if

c     fill in the input/output arrays for the interpolator
      interp_coords(1) = CCTK_PointerTo(xcoord)
      interp_coords(2) = CCTK_PointerTo(ycoord)
      interp_coords(3) = CCTK_PointerTo(zcoord)

      call CCTK_VarIndex(i, "DemoInterp::realgf3")
      in_array_indices(1) = i
      call CCTK_VarIndex(i, "DemoInterp::compgf3")
      in_array_indices(2) = i

      out_array_type_codes(1) = CCTK_VARIABLE_REAL
      out_array_type_codes(2) = CCTK_VARIABLE_COMPLEX

      out_arrays(1) = CCTK_PointerTo(realinterp3)
      out_arrays(2) = CCTK_PointerTo(compinterp3)

c     ==================================================================
c     Interpolation of both real and complex GFs using general call
c     ==================================================================

      call CCTK_InterpGridArrays (ierror, cctkGH, 3, interp_handle,
     .                            param_table_handle, coord_system_handle,
     .                            arrsize, CCTK_VARIABLE_REAL, interp_coords,
     .                            2, in_array_indices,
     .                            2, out_array_type_codes, out_arrays)

      if (ierror<0) then
         call CCTK_WARN(1,"Interpolation error")
      endif

      print *
      print *,"Real 3D interpolation using CCTK_InterpGridArrays()"
      print *,"==================================================="
      print *,"Interpolated values should be ",interp_radius
      print *,"Actual values are : "
      do i=1,arrsize
         print *,realinterp3(i)
      end do

      print *
      print *,"Complex 3D interpolation using CCTK_InterpGridArrays()"
      print *,"======================================================"
      print *,"Interpolated values should be (",interp_radius,",",
     &     cos(interp_radius),")"
      print *,"Actual values are : "
      do i=1,arrsize
         print *,compinterp3(i)
      end do

      return
      end
