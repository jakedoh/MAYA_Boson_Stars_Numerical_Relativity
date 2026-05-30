! $Id: SetupDetectors.F90,v 1.7 2005/07/26 12:53:57 herrmann Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine WavExtr_Setup_Detectors(CCTK_ARGUMENTS)

  use WavExtrConstants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL xmin,xmax
  CCTK_REAL ymin,ymax
  CCTK_REAL zmin,zmax

  CCTK_REAL dx,dy,dz

  CCTK_REAL rmax,rmin
  CCTK_REAL rmax_grid

  CCTK_REAL rstart,rend,dr,radius
  CCTK_INT idetector

  CCTK_INT ierror

  ! interpolator
  CCTK_INT :: interp_handle, table_handle, coord_system_handle
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_INT, dimension(1) :: in_array
  CCTK_POINTER, dimension(1) :: out_array
  character(len=200) :: interp
  CCTK_INT :: interp_len
  character(len=7) :: interp_order
  CCTK_INT, dimension(1) :: out_types
  integer it,ip
  CCTK_INT, dimension(2) :: lsh,lbnd

  integer rmax_found

  CCTK_REAL grid_acc

  integer :: i

  character(len=512) :: infoline

  logical :: firstcal

  data firstcal / .true. /
  save firstcal

! _________________________________________________________________


  if (verbose>2) then
    call CCTK_INFO("Setup of detectors")
  end if

  if (.not. firstcal) then
    return
  end if

  ! save firstcal var
  firstcal=.false.

  ! coarse grid information needed
  dx=cctk_delta_space(1)
  dy=cctk_delta_space(2)
  dz=cctk_delta_space(3)

  ! Check if dx=dy=dz
  if (cartoon.ne.0) then
    grid_acc=cartoon_grid_acc
  else
    grid_acc=1.d-5
  end if
  if ( (dx-dy)/dx.gt.grid_acc .or. (dx-dz)/dx.gt.grid_acc ) then
    call CCTK_WARN(1,"Strange grid: dx,dy,dz are not equal!")
    write(infoline,'(A7,G20.8,A6,G20.8,A6,G20.8)') '  dx = ',dx,' dy = ',dy,' dz = ',dz
    call CCTK_WARN(1,infoline)
    call CCTK_WARN(1,"No Extraction will be done.")
    do_nothing=1
    current_detector=0
    return
  end if


  ! Find the maximum radius we can use for our detectors. Coordinate Range:
  call CCTK_CoordRange(ierror,cctkGH,xmin,xmax,-1,"x","cart3d")
  call CCTK_CoordRange(ierror,cctkGH,ymin,ymax,-1,"y","cart3d")
  call CCTK_CoordRange(ierror,cctkGH,zmin,zmax,-1,"z","cart3d")

  rmax=max(xmax,ymax,zmax)
  rmin=min(xmin,ymin,zmin)
  rmax_grid=rmax

  if (interpolate_for_rmax>0) then

     ! call interpolator silently and check if a sphere with rmax would work.
     ! if not, reduce rmax by dx and try again. keep trying until it works.
     call CCTK_FortranString(interp_len, interpolation_operator, interp)
     
     call CCTK_InterpHandle(interp_handle,interp(1:interp_len))
     if ( interp_handle .lt. 0 ) then
        call CCTK_WARN( 0, "Cannot get handle for interpolation. Forgot to activate an implementation providing interpolation operators?" )
     end if

     write(interp_order,'(a6,i1)') 'order=',interpolation_order
     call Util_TableCreateFromString(table_handle,interp_order)
     if ( table_handle .lt. 0 ) then
        call CCTK_WARN( 0, "Cannot create parameter table for interpolator" )
     end if
     
     call Util_TableSetInt(ierror, table_handle, 0, "suppress_warnings")
     if (ierror .lt. 0) then
        call CCTK_WARN(1,"cannot suppress warnings of interpolator")
     end if

     call CCTK_CoordSystemHandle ( coord_system_handle, "cart3d" )
     if ( coord_system_handle .lt. 0) then
        call CCTK_WARN( 0, "Cannot get handle for observers coordinate system. Forgot to activate an implementation providing coordinates?" )
     end if

     ! suppress_warnings is needed to get rid of warning messages
     
     ! loop down rmax until interpolator does not return error anymore
     rmax_found=0
     do while (rmax_found .ne. 1)
        ! test xmax,ymax,zmax and diagonal element only
        interp_x = zero
        interp_y = zero
        interp_z = zero
        interp_x(1,1) = origin_x + rmax
        interp_y(1,1) = dx
        interp_z(1,1) = dx
        interp_x(2,1) = dx
        interp_y(2,1) = origin_y + rmax
        interp_z(2,1) = dx
        interp_x(3,1) = dx
        interp_y(3,1) = dx
        interp_z(3,1) = origin_z + rmax
        interp_x(4,1) = origin_x + 1.d0/sqrt(3.d0)*rmax
        interp_y(4,1) = origin_y + 1.d0/sqrt(3.d0)*rmax
        interp_z(4,1) = origin_z + 1.d0/sqrt(3.d0)*rmax
        lsh(1)=4
        lsh(2)=1
        
        interp_coords(1) = CCTK_PointerTo(interp_x)
        interp_coords(2) = CCTK_PointerTo(interp_y)
        interp_coords(3) = CCTK_PointerTo(interp_z)
        
        ! abuse the lapse to find interpolation range
        call CCTK_VarIndex(in_array(1), "ADMBase::alp")
        out_array(1) = CCTK_PointerTo(psii)
        out_types = CCTK_VARIABLE_REAL
        
        if (verbose > 3) then     
          print*,'test call to interpolator at rmax=',rmax
        end if
     
        call CCTK_InterpGridArrays(ierror, cctkGH, 3, interp_handle, &
             table_handle, coord_system_handle, &
             lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
             interp_coords, &
             1, &
             in_array, &
             1, &
             out_types, &
             out_array)
        if (ierror.ne.0) then
           call CCTK_WARN(2,"test interpolation to determine rmax failed")
           print*,'rmax is ',rmax
        else
           if (verbose>4) then
              print*,'rmax found and nailed at ',rmax
           end if
           rmax_found=1
        end if
        
        rmax=rmax-dx
     end do

     if (verbose>1) then
        write(infoline,'(A11,F6.2)') '  ... rmax=',rmax
        call CCTK_INFO(infoline)
     end if

  end if


  ! for minimum radius we just use 3*dx arbitrarily
  rmin = 3.d0*dx

  if (verbose>1) then
    write(infoline,'(A11,F6.2)') '  ... rmin=',rmin
    call CCTK_INFO(infoline)
  end if


  ! There are 2 different modes. In the normal mode the user specifies
  ! locations for a bunch of detectors. Alternatively the user can
  ! ask for a range in which detectors should be placed. In that case
  ! we put as many detectors as requested by number_of_detectors on
  ! the grid.

  ! Method 1: Detectors given
  if (CCTK_EQUALS(detection_mode,"specific detectors")) then
    ! Check each detector if it is out of range.
    do i=1,number_of_detectors
      if (surface_index(i)>=0) then
        continue
      end if
      if(detector_radius(i)>rmax) then
        write(infoline,'(A29,I2,A8,F6.2,A7,F6.2)') &
                    '  Out of RANGE! Detector No. ',&
                           i,' radius=',detector_radius(i),'> rmax=',rmax
        call CCTK_WARN(1,infoline)
        detector_radius(i)=zero
      end if
    end do

    ! Resort the detector array and reset number_of_detectors
    do i=1,number_of_detectors
      if (i>1) then
        if(detector_radius(i-1).eq.zero .and. detector_radius(i).ne.zero &
           .and. surface_index(i)<0 .and. surface_index(i-1)<0) then
          detector_radius(i-1)=detector_radius(i)
          detector_radius(i)=zero
        end if
      end if
      if(detector_radius(i).eq.zero .and. surface_index(i)<0) then
        number_of_detectors=number_of_detectors-1
      end if
    end do

    ! Check if we have no detectors left
    if(number_of_detectors .eq. 0) then
      do_nothing=1
      current_detector=0
      return
    end if

  ! Method 2: Cauchy extraction
  else if (CCTK_EQUALS(detection_mode,"Cauchy extraction")) then
    ! check start and end values
    if (Cauchy_radius_start_coord == -two .and. &
        Cauchy_radius_start_factor == -two .and. &
        Cauchy_radius_end_coord  == -two .and. &
        Cauchy_radius_end_factor == -two) Then
      write(infoline,'(A65)') &
        "You have specfied Cauchy extraction but deactivated all detectors"
      call CCTK_WARN(1,infoline)
      number_of_detectors=0
      do_nothing=1
      current_detector=0
      return
    end if

    ! set start and end radius
    ! 1) by coordinate location in par-file
    if (Cauchy_radius_start_coord+one .lt. 1.d-13) then
      rstart=zero
    else if (Cauchy_radius_start_coord > zero) then
      rstart=Cauchy_radius_start_coord
      print*,'set rstart'
    end if
    if (Cauchy_radius_end_coord == -one) then
      rend=rmax
    else if (Cauchy_radius_end_coord > -one) then
      rend=Cauchy_radius_end_coord
    end if

    print*,'rstart',rstart

    ! 2) by factor of grid size
    if (Cauchy_radius_start_factor .gt. zero) then
      rstart=rmax_grid*Cauchy_radius_start_factor
      print*,"set rstart again"
    end if
    if (Cauchy_radius_end_factor .gt. zero) then
      rend=rmax_grid*Cauchy_radius_end_factor
    end if

    ! work out dr
    ! user given or from number_of_detectors
    if (Cauchy_radius_dr .gt. zero) then
      dr=Cauchy_radius_dr
    else
      if (number_of_detectors>1) then
        dr=(rend - rstart)/(number_of_detectors -1)
      else
        call CCTK_WARN(1,"you need to set number_of_detectors>1")
        do_nothing=1
        current_detector=0
        return
      end if
      if (cauchy_radius_extreme.ne.0) then
        rstart=rmin
        dr=dx
      end if
    end if

    if (verbose>0) then
      write(infoline,'(A24,F10.4,A5,F10.4,A11,F10.4)') &
          'Cauchy detector setup: [',rstart,' --- ',rend,' ]  every: ',dr
      call CCTK_INFO(infoline)
    end if

    if (dr < dx) then
      write(infoline,'(A52,A5,F5.3,A7,F5.3)') &
       'dr used in extraction is smaller than grid spacing: ', &
       'dr = ',dr,' < dx= ',dx
      call CCTK_WARN(1,infoline)
      call CCTK_WARN(1,"This is not necessarily a problem ;-)")
    end if

    if (rstart < rmin) then
      call CCTK_WARN(1,"Don't specify the first detector < rmin. Increasing to rmin.")
      rstart = rmin
    end if

    if (rend>rmax) then
      call CCTK_WARN(1,"some detectors are out of the possible grid and will be dropped")
      call CCTK_WARN(1,"to get more detectors into the grid:")
      call CCTK_WARN(1,"     adjust interpolation stencil or increase resolution")
      do while (rend>rmax)
        write(infoline,'(A33,G20.8)') '  Dropping Detector at location: ',rend
        call CCTK_INFO(infoline)
        rend=rend-dr
      end do
      write(infoline,'(A33,G20.8)') '  Adjusted outer detector range: ',rend
      call CCTK_WARN(1,infoline)
    end if

    if (rstart >= rend) then
      call CCTK_WARN(1,"very strange detector setup, start range > end range")
      call CCTK_WARN(1,"we have no useable detectors - NO EXTRACTION")
      do_nothing=1
      current_detector=0
      return
    end if

    ! Compute the locations of the detector
    idetector=1
    radius=rstart
    number_of_detectors=1
    do while (radius<=rend)
      detector_radius(idetector)=radius

      radius=radius+dr
      idetector=idetector+1
      number_of_detectors=number_of_detectors+1
    end do
    number_of_detectors=number_of_detectors-1
  end if

  if (number_of_detectors.eq.0) then
    call CCTK_WARN(1,"no appropriate detectors found")
    do_nothing=1
    current_detector=0
    return
  end if

  do_nothing=0
  current_detector=0

  ! initialize the arrays for radiated energy and linear momentum
  Erad=0
  Pxrad=0
  Pyrad=0
  Pzrad=0

  grr_p=0
  grt_p=0
  grp_p=0
  gtt_p=0
  gtp_p=0
  gpp_p=0

  if (verbose > 1) then
    call CCTK_INFO("We are done with the detector setup.")
    write(infoline,'(A10,I5,A10)') '  We have ', &
          number_of_detectors,' detectors'
    call CCTK_INFO(infoline)
    idetector=1
    do idetector=1,number_of_detectors
      write(infoline,'(A12,I5,A10,F10.4)') '    detector', &
            idetector,' Radius ',detector_radius(idetector)
      call CCTK_INFO(infoline)
    end do
  end if

  end subroutine WavExtr_Setup_Detectors
