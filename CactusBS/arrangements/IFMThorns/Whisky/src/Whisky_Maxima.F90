 /*@@
   @file      Whisky_Maxima.F90
   @date      Mon May 17 07:05:55 2004
   @author    Ian Hawke
   @desc 
   Find the maximum of the density over the whole (x>0) grid
   and use its location to set the DriftCorrection thorn.
   @enddesc 
 @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    Whisky_FindMaxima
   @date       Mon May 17 07:06:49 2004
   @author     Ian Hawke
   @desc 
   Find the maximum of the density over the whole (x>0) grid
   and use its location to set the DriftCorrection thorn.
   Some of this code is borrowed, with alterations, from
   Peter Dieners EHFinder.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_InitFindMaxima(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  whisky_maxima_iteration = -1

  maximum_density = -1.d0

end subroutine Whisky_InitFindMaxima


subroutine Whisky_FindMaxima(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: imin, imax, jmin, jmax, kmin, kmax 
  CCTK_INT :: ierr, max_handle, sum_handle
  CCTK_INT, dimension(3) :: maxrho_local_location, maxrho_location
  CCTK_REAL, dimension(3) :: maxrho_local_position, maxrho_position
  CCTK_REAL :: maxrho_local

  maxrho_local_location = 0 
  maxrho_location       = 0 
  maxrho_local_position = 0.d0
  maxrho_position       = 0.d0


  if (cctk_iteration > whisky_maxima_iteration) then
    whisky_maxima_iteration = cctk_iteration
    maximum_density = 0.d0
  end if
  ! remember that in fortran array indeces start from 1
  imin = whisky_stencil + 1
  jmin = whisky_stencil + 1
  kmin = whisky_stencil + 1
  imax = cctk_lsh(1) - whisky_stencil
  jmax = cctk_lsh(2) - whisky_stencil
  kmax = cctk_lsh(3) - whisky_stencil

  call CCTK_ReductionArrayHandle (max_handle, 'maximum')
  call CCTK_ReductionArrayHandle (sum_handle, 'sum')

  maxrho_local_location = maxloc ( rho(imin:imax,jmin:jmax,kmin:kmax) )

  maxrho_local = rho (&
       maxrho_local_location(1) + whisky_stencil, &
       maxrho_local_location(2) + whisky_stencil, &
       maxrho_local_location(3) + whisky_stencil  )

!!$  write(*,*)
!!$  write(*,*) "Itaration",cctk_iteration
!!$  write(*,*) "1) On reflevel",aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))
!!$  write(*,*) "Max rho on level",maxrho_local, maxval ( rho(imin:imax,jmin:jmax,kmin:kmax) )
!!$  write(*,*) "Mrho_loc_loc",maxrho_local_location
!!$
!!$
!!$  if (aint(log10(dble(cctk_levfac(1)))/log10(2.0d0)) == 3) then
!!$    do ierr=1,cctk_lsh(1)
!!$      write(*,*) "x(",ierr,")=",x(ierr,24,24),y(ierr,24,24),z(ierr,24,24)
!!$      write(*,*) "rho",rho(ierr,24,24)
!!$    end do
!!$  end if
  



!!$  Only check the x>0 part of the grid. This ensures that if there are 
!!$  two stars either side of x=0 only one will be picked every time.

  if ( Whisky_useFullGridForMaxima .eq. 0 ) then
    if ( x (&
         maxrho_local_location(1) + whisky_stencil, &
         maxrho_local_location(2) + whisky_stencil, &
         maxrho_local_location(3) + whisky_stencil  ) < 0.d0 ) then
      maxrho_local = 0.d0 ! this is how maxrho_local can be seen to be zero in some output
    end if
  end if
  
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, max_handle, &
       maxrho_local, maxrho_global, CCTK_VARIABLE_REAL)
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0, 'Reduction of maxrho failed')
  end if
  
  if ( maxrho_local .eq. maxrho_global ) then
    maxrho_local_position(1) = x (&
       maxrho_local_location(1) + whisky_stencil, &
       maxrho_local_location(2) + whisky_stencil, &
       maxrho_local_location(3) + whisky_stencil )
    maxrho_local_position(2) = y (&
       maxrho_local_location(1) + whisky_stencil, &
       maxrho_local_location(2) + whisky_stencil, &
       maxrho_local_location(3) + whisky_stencil )
    maxrho_local_position(3) = z (&
       maxrho_local_location(1) + whisky_stencil, &
       maxrho_local_location(2) + whisky_stencil, &
       maxrho_local_location(3) + whisky_stencil )
  else

    !if the processor does not have the global maximum, set maxrho_local_position to zero 
    !(for use in the below "sum" reduction)
    maxrho_local_position = 0.d0 !it is a vector


    !if the processor does not have the global maximum, set maxrho_local_location to zero 
    !(for use in the below "sum" reduction)
    maxrho_local_location = 0 !it is a vector

  end if

!!$  write(*,*)
!!$  write(*,*) "2) On reflevel",aint(log10(dble(cctk_levfac(1)))/log10(2.0d0))
!!$  write(*,*) "Max rho on level",maxrho_local
!!$  write(*,*) "Max rho this time step up to now",maximum_density
!!$  write(*,*) "Mrho_loc_loc",maxrho_local_location
!!$  write(*,*) "Mrho_loc_pos",maxrho_local_position

  !find maxrho_local_position on the processor that has the global maximum
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
       maxrho_local_position, maxrho_position, 3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0, 'Reduction of maxrho_position failed')
  end if

  !find maxrho_local_location on the processor that has the global maximum
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
       maxrho_local_location, maxrho_location, 3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0, 'Reduction of maxrho_position failed')
  end if

!!$  write(*,*) "Mrho_glob_loc",maxrho_location
!!$  write(*,*) "Mrho_glob_pos",maxrho_position


!!$  If the maximum on this iteration was larger on another level, 
!!$  do not set the position

  if (maxrho_global > maximum_density) then

    maximum_density = maxrho_global

    maxima_x = maxrho_position(1)
    maxima_y = maxrho_position(2)
    maxima_z = maxrho_position(3)

    maxima_i = maxrho_location(1)
    maxima_j = maxrho_location(2)
    maxima_k = maxrho_location(3)

  end if
  
!!$  write(*,*) "maxima location",maxima_i,maxima_j,maxima_k
!!$  write(*,*) "maxima position",maxima_x, maxima_y, maxima_z


end subroutine Whisky_FindMaxima

 /*@@
   @routine    Whisky_FindWeightedMaxima
   @date       Mon May 24 23:22:18 2004
   @author     Ian Hawke
   @desc 
   Find the volume weighted centre of mass, 
            /
           |       i
           |  rho x    dV
           |
  - i      /
  x  =     ----------------
           /             
           |             
           |  rho      dV 
           |              
           /              

   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_FindWeightedMaxima(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
  CCTK_INT :: ierr, sum_handle
  CCTK_REAL, dimension(3) :: weightedrho_local, weightedrho_global
  CCTK_REAL :: rho_local, rho_global

  imin = whisky_stencil
  jmin = whisky_stencil
  kmin = whisky_stencil
  imax = cctk_lsh(1) - whisky_stencil
  jmax = cctk_lsh(2) - whisky_stencil
  kmax = cctk_lsh(3) - whisky_stencil

  call CCTK_ReductionArrayHandle (sum_handle, 'sum')

  weightedrho_local = 0.d0
  rho_local = 0.d0

  !$OMP PARALLEL DO PRIVATE (i,j,k) &
  !$OMP REDUCTION (+ : weightedrho_local, rho_local)
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax

        weightedrho_local(1) = weightedrho_local(1) + &
             rho(i, j, k) * x(i, j, k)
        weightedrho_local(2) = weightedrho_local(2) + &
             rho(i, j, k) * y(i, j, k)
        weightedrho_local(3) = weightedrho_local(3) + &
!!$             rho(i, j, k) * z(i, j, k)
!!$  To get around problems with symmetry boundaries and because
!!$  the BNS do not have spin (yet) I just set \bar{z} to zero.             
             0.d0
        rho_local = rho_local + rho(i,j,k)

      end do
    end do
  end do
  !$OMP END PARALLEL DO
  
  call CCTK_ReduceLocArrayToArray1D ( ierr, cctkGH, -1, sum_handle, &
       weightedrho_local, weightedrho_global, 3, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0, 'Reduction of weightedrho failed')
  end if
  
  call CCTK_ReduceLocScalar ( ierr, cctkGH, -1, sum_handle, &
       rho_local, rho_global, CCTK_VARIABLE_REAL )
  if ( ierr .ne. 0 ) then
    call CCTK_WARN(0, 'Reduction of rho failed')
  end if

  maxima_x = weightedrho_global(1) / rho_global
  maxima_y = weightedrho_global(2) / rho_global
  maxima_z = weightedrho_global(3) / rho_global
  

end subroutine Whisky_FindWeightedMaxima

 /*@@
   @routine    Whisky_SetDCCentroid
   @date       Mon May 17 07:06:49 2004
   @author     Ian Hawke
   @desc 
   Find the maximum of the density over the whole (x>0) grid
   and use its location to set the DriftCorrection thorn.
   Some of this code is borrowed, with alterations, from
   Peter Dieners EHFinder.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine Whisky_SetDCCentroid(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: ierr

  if (set_dc_centroid .ne. 0) then

    call CCTK_IsFunctionAliased(ierr, "SetDriftCorrectPosition")
    if ( ierr .ne. 0 ) then
      call SetDriftCorrectPosition ( cctkGH, &
           maxima_x, maxima_y, maxima_z )
    end if
  
  end if
  
end subroutine Whisky_SetDCCentroid

 /*@@
   @routine    Whisky_FindSeparation
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

subroutine Whisky_FindSeparation(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL :: separation_dx, separation_dy, separation_dz
  CCTK_REAL, dimension(:), allocatable :: separation_x, separation_y, &
       separation_z
  CCTK_REAL, dimension(:), allocatable :: s_gxx, s_gxy, s_gxz, &
       s_gyy, s_gyz, s_gzz
  CCTK_INT :: ierr, i
  CCTK_INT :: param_table_handle, interp_handle, coord_system_handle
  CCTK_INT :: vindex

  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_INT, dimension(6) :: in_array_indices
  CCTK_POINTER, dimension(6) :: out_arrays
  CCTK_INT, dimension(6) :: out_array_type_codes

!!$  Coordinate separation is easy.

  whisky_separation = sqrt(maxima_x**2 + maxima_y**2 + maxima_z**2)

!!$  Proper separation requires interpolation

  allocate(separation_x(separation_npoints), &
       separation_y(separation_npoints), &
       separation_z(separation_npoints), STAT=ierr)

  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Failed to allocate separation coordinate arrays")
  end if
  
  separation_dx = maxima_x / dble(separation_npoints)
  separation_dy = maxima_y / dble(separation_npoints)
  separation_dz = maxima_z / dble(separation_npoints)

  do i = 1, separation_npoints
    
    separation_x(i) = i * separation_dx
    separation_y(i) = i * separation_dy
    separation_z(i) = i * separation_dz

  end do
  
  allocate(s_gxx(separation_npoints), &
       s_gxy(separation_npoints), &
       s_gxz(separation_npoints), &
       s_gyy(separation_npoints), &
       s_gyz(separation_npoints), &
       s_gzz(separation_npoints), &
       STAT=ierr)

  if (ierr .ne. 0) then
    call CCTK_WARN(0, "Failed to allocate separation metric arrays")
  end if

  param_table_handle = -1
  interp_handle = -1
  coord_system_handle = -1

  call Util_TableCreateFromString (param_table_handle, "order = 2")
  if (param_table_handle .lt. 0) then
    call CCTK_WARN(0,"Cannot create parameter table for interpolator")
  endif

  call CCTK_InterpHandle (interp_handle, "uniform cartesian")
  if (interp_handle.lt.0) then
    call CCTK_WARN(0,"Cannot get handle for interpolation ! Forgot to activate an implementation providing interpolation operators ??")
  endif

  call CCTK_CoordSystemHandle (coord_system_handle, "cart3d")
  if (coord_system_handle .lt. 0) then
    call CCTK_WARN(0,"Cannot get handle for cart3d coordinate system ! Forgot to activate an implementation providing coordinates ??")
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
       separation_npoints, CCTK_VARIABLE_REAL, interp_coords,&
       6, in_array_indices,&
       6, out_array_type_codes, out_arrays)
  if (ierr < 0) then
    call CCTK_WARN (1, "interpolator call returned an error code");
  endif
  
!     release parameter table
  call Util_TableDestroy (ierr, param_table_handle)
  
!        Integrate using trapezoidal rule.
  
  whisky_proper_separation = 0.d0

  do i=2, separation_npoints
    
    whisky_proper_separation = whisky_proper_separation + 0.5d0 &
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
    
end subroutine Whisky_FindSeparation
