#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
 
 /*@@
   @routine    Hydro_Analysis_FindSeparation
   @date       Thu May 20 12:35:20 2004
   @author     Ian Hawke
   @desc 
   Finds the separation (in coordinate and proper distance)
   between the NS. Equal mass symmetry is assumed so it is
   actually the distance between the origin and the location
   of maximum density. This is along a straight line, not a
   geodesic, so still not perfect.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Hydro_Analysis_FindSeparation(CCTK_ARGUMENTS)

#ifdef HAVE_CAPABILITY_Fortran
  use util_Table
#endif

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: separation_dx, separation_dy, separation_dz
  CCTK_REAL, dimension(:), allocatable :: separation_x, separation_y, &
       separation_z
  CCTK_REAL, dimension(:), allocatable :: s_gxx, s_gxy, s_gxz, &
       s_gyy, s_gyz, s_gzz
  integer :: ierr, i
  integer :: param_table_handle, interp_handle, coord_system_handle
  CCTK_INT :: vindex

  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_INT, dimension(6) :: in_array_indices
  CCTK_POINTER, dimension(6) :: out_arrays
  CCTK_INT, dimension(6) :: out_array_type_codes

  ! Fortran copies of the Cactus CCTK_STRING parameters controllin the
  ! interpolator
  integer, parameter:: max_string_length = 500
  integer :: string_length
  character(len=max_string_length) :: Hydro_Analysis_interpolator_options_fstring
  character(len=max_string_length) :: Hydro_Analysis_interpolator_name_fstring
  character(len=max_string_length) :: Hydro_Analysis_interpolator_coordinates_fstring


!!$  Proper separation requires interpolation

  allocate(separation_x(Hydro_Analysis_rho_max_origin_distance_npoints), &
       separation_y(Hydro_Analysis_rho_max_origin_distance_npoints), &
       separation_z(Hydro_Analysis_rho_max_origin_distance_npoints), STAT=ierr)

  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Failed to allocate separation coordinate arrays")
  end if

  call CCTK_FortranString(string_length, Hydro_Analysis_interpolator_options, &
                          Hydro_Analysis_interpolator_options_fstring)
  if (string_length .gt. max_string_length) then
     call CCTK_WARN(CCTK_WARN_ALERT, "'Hydro_Analysis_interpolator_options' string too long!")
  end if
  call CCTK_FortranString(string_length, Hydro_Analysis_interpolator_name, &
                          Hydro_Analysis_interpolator_name_fstring)
  if (string_length .gt. max_string_length) then
     call CCTK_WARN(CCTK_WARN_ALERT, "'Hydro_Analysis_interpolator_name' string too long!")
  end if
  call CCTK_FortranString(string_length, Hydro_Analysis_interpolator_coordinates, &
                          Hydro_Analysis_interpolator_coordinates_fstring)
  if (string_length .gt. max_string_length) then
     call CCTK_WARN(CCTK_WARN_ALERT, "'Hydro_Analysis_interpolator_coordinates' string too long!")
  end if

  
  separation_dx = Hydro_Analysis_rho_max_loc(1) / dble(Hydro_Analysis_rho_max_origin_distance_npoints)
  separation_dy = Hydro_Analysis_rho_max_loc(2) / dble(Hydro_Analysis_rho_max_origin_distance_npoints)
  separation_dz = Hydro_Analysis_rho_max_loc(3) / dble(Hydro_Analysis_rho_max_origin_distance_npoints)

  do i = 1, Hydro_Analysis_rho_max_origin_distance_npoints
    
    separation_x(i) = i * separation_dx
    separation_y(i) = i * separation_dy
    separation_z(i) = i * separation_dz

  end do
  
  allocate(s_gxx(Hydro_Analysis_rho_max_origin_distance_npoints), &
       s_gxy(Hydro_Analysis_rho_max_origin_distance_npoints), &
       s_gxz(Hydro_Analysis_rho_max_origin_distance_npoints), &
       s_gyy(Hydro_Analysis_rho_max_origin_distance_npoints), &
       s_gyz(Hydro_Analysis_rho_max_origin_distance_npoints), &
       s_gzz(Hydro_Analysis_rho_max_origin_distance_npoints), &
       STAT=ierr)

  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Failed to allocate separation metric arrays")
  end if

  param_table_handle = -1
  interp_handle = -1
  coord_system_handle = -1

  call Util_TableCreateFromString (param_table_handle, &
                                   Hydro_Analysis_interpolator_options_fstring)
  if (param_table_handle .lt. 0) then
    call CCTK_WARN(0,"Cannot create parameter table for interpolator")
  endif

  call CCTK_InterpHandle (interp_handle, Hydro_Analysis_interpolator_name_fstring)
  if (interp_handle.lt.0) then
    call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators (e.g. LocalInterp)?")
  endif

  call CCTK_CoordSystemHandle (coord_system_handle, Hydro_Analysis_interpolator_coordinates_fstring)
  if (coord_system_handle .lt. 0) then
    call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates (e.g. LocalInterp)?")
  endif

!     fill in the input/output arrays for the interpolator
  interp_coords(1) = CCTK_PointerTo(separation_x)
  interp_coords(2) = CCTK_PointerTo(separation_y)
  interp_coords(3) = CCTK_PointerTo(separation_z)
  
  call CCTK_VarIndex (vindex, "admbase::gxx")
  in_array_indices(1) = vindex
  call CCTK_VarIndex (vindex, "admbase::gyy")
  in_array_indices(2) = vindex
  call CCTK_VarIndex (vindex, "admbase::gzz")
  in_array_indices(3) = vindex
  call CCTK_VarIndex (vindex, "admbase::gxy")
  in_array_indices(4) = vindex
  call CCTK_VarIndex (vindex, "admbase::gxz")
  in_array_indices(5) = vindex
  call CCTK_VarIndex (vindex, "admbase::gyz")
  in_array_indices(6) = vindex
  
  out_arrays(1) = CCTK_PointerTo(s_gxx)
  out_arrays(2) = CCTK_PointerTo(s_gyy)
  out_arrays(3) = CCTK_PointerTo(s_gzz)
  out_arrays(4) = CCTK_PointerTo(s_gxy)
  out_arrays(5) = CCTK_PointerTo(s_gxz)
  out_arrays(6) = CCTK_PointerTo(s_gyz)
  
  out_array_type_codes = CCTK_VARIABLE_REAL
  
!     Interpolation.

  call CCTK_InterpGridArrays (ierr, cctkGH, 3, interp_handle,&
       param_table_handle, coord_system_handle,&
       Hydro_Analysis_rho_max_origin_distance_npoints, CCTK_VARIABLE_REAL, interp_coords,&
       6, in_array_indices,&
       6, out_array_type_codes, out_arrays)
  if (ierr < 0) then
    call CCTK_WARN (1, "interpolator call returned an error code");
  endif
  
!     release parameter table
  call Util_TableDestroy (ierr, param_table_handle)
  
!        Integrate using trapezoidal rule.
  
  Hydro_Analysis_rho_max_origin_distance = 0.d0

  do i=2, Hydro_Analysis_rho_max_origin_distance_npoints
    
    Hydro_Analysis_rho_max_origin_distance = Hydro_Analysis_rho_max_origin_distance + 0.5d0 &
         *(sqrt(s_gxx(i-1)*separation_dx**2  + &
                s_gyy(i-1)*separation_dy**2  + &
                s_gzz(i-1)*separation_dz**2 &
        + 2.d0*(s_gxy(i-1)*separation_dx*separation_dy + &
                s_gxz(i-1)*separation_dx*separation_dz + &
                s_gyz(i-1)*separation_dy*separation_dz)) &
         + sqrt(s_gxx(i  )*separation_dx**2  + &
                s_gyy(i  )*separation_dy**2  + &
                s_gzz(i  )*separation_dz**2 &
        + 2.d0*(s_gxy(i  )*separation_dx*separation_dy + &
                s_gxz(i  )*separation_dx*separation_dz + &
                s_gyz(i  )*separation_dy*separation_dz)))
    
  end do
    
end subroutine Hydro_Analysis_FindSeparation

