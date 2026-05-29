
! $Header: /numrelcvs/HerrmannCVS/Ylm_Decomp/src/Project_Sphere.F90,v 1.19 2005/01/11 15:16:14 herrmann Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine Ylm_Project_Sphere(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: im, jm, it, ip
  CCTK_INT :: interp_handle, table_handle, coord_system_handle

  character(len=200) :: interp
  CCTK_INT :: interp_len
  character(len=7) :: interp_order

  CCTK_REAL :: dthetainv, dphiinv
  CCTK_REAL :: dxdth, dxdph, dydth, dydph, dzdth, dzdph
!  CCTK_REAL :: gtt, gtp, gpp
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER out_array(nrdecompvars)
  CCTK_INT in_array(nrdecompvars)
  CCTK_INT out_types(nrdecompvars)
  CCTK_INT :: i, vind, myproc

! temp shortcuts
  CCTK_REAL :: tgxx, dgxx, &
               tgxy, dgxy, &
               tgxz, dgxz, &
               tgyy, dgyy, &
               tgyz, dgyz, &
               tgzz, dgzz

  CCTK_INT op_indices(nrdecompvars)
  CCTK_INT op_codes(nrdecompvars)

  CCTK_INT :: ierr,num_in_arrays,num_out_arrays,status,j
  CCTK_REAL :: rad,r2,ct,st,ct2,st2,cp,cp2,sp,sp2

  CCTK_REAL :: lapse,cor_angle
  CCTK_REAL :: tmpval
  CCTK_INT :: nanind
  CCTK_INT :: intp_npoints
  CCTK_INT :: number_of_killed_nans

  CCTK_REAL :: tmpx,tmpy,tmpz

  character(len=80) :: infoline

! ________________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO('Interpolating integrands onto sphere')

  if (do_nothing == 1) &
    return

  !print*,'current detector ',current_detector,detector_radius(current_detector)

  current_detector_radius = detector_radius(current_detector)

  ! Sanity Check
  if (current_detector_radius == zero) then
    do_nothing =1
    print*,"current_detector=",current_detector
    call CCTK_WARN(1,"This should never happen: The detector radius is 0!")
    return
  end if


  if (verbose>1) then
    write(infoline,'(a24,i4,a8,g16.8)') 'Analysing Detector No.: ',current_detector, &
                                        ' Radius ',current_detector_radius
    call CCTK_INFO(infoline)
  end if

  dthetainv = one / dtheta
  dphiinv = one / dphi

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

  call CCTK_CoordSystemHandle ( coord_system_handle, "cart3d" )
  if ( coord_system_handle .lt. 0) then
    call CCTK_WARN( 0, "Cannot get handle for cart3d coordinate system. Forgot to activate an implementation providing coordinates?" )
  endif

  ! make sure interp coordinates are initialized
  interp_x=zero
  interp_y=zero
  interp_z=zero

  interp_coords(1) = CCTK_PointerTo(interp_x)
  interp_coords(2) = CCTK_PointerTo(interp_y)
  interp_coords(3) = CCTK_PointerTo(interp_z)

  do i=1, nrdecompvars
    out_array(i)  = CCTK_PointerTo(interp_arrays(1,1,i))
  end do
  ! FIXME : output arrays should be dynamic
  num_out_arrays=nrdecompvars

  ! set all out_types(:) to type real
  out_types=CCTK_VARIABLE_REAL

  ! get coordinates
  rad=current_detector_radius
  do ip = 1, cnphi
    do it = 1, cntheta
      interp_x(it,ip) = origin_x + rad * sintheta(it,ip) * cosphi(it,ip)
      interp_y(it,ip) = origin_y + rad * sintheta(it,ip) * sinphi(it,ip)
      interp_z(it,ip) = origin_z + rad * costheta(it,ip)
    end do
  end do

  if (maxntheta.gt.cntheta) then
     if (verbose>1) then
        call CCTK_INFO("remove theta points out of range")
     end if
     
     tmpx=origin_x + rad
     tmpy=origin_y
     tmpz=origin_z
     do ip = 1, cnphi
        do it = cntheta+1, maxntheta
           interp_x(it,ip) = tmpx
           interp_y(it,ip) = tmpy
           interp_z(it,ip) = tmpz
        end do
     end do
  end if

  if (maxnphi.gt.cnphi) then
     if (verbose>1) then
        call CCTK_INFO("remove phi points out of range")
     end if

     tmpx=origin_x + rad
     tmpy=origin_y
     tmpz=origin_z
     do ip = cnphi+1, maxnphi
        do it = 1, cntheta
           interp_x(it,ip) = tmpx
           interp_y(it,ip) = tmpy
           interp_z(it,ip) = tmpz
        end do
     end do
  end if


  !print*,'inter_x',interp_x
  !print*,'inter_y',interp_y
  !print*,'inter_z',interp_z

  do i=1,nrdecompvars
    in_array(i)=varindices(i)
  end do
  num_in_arrays=nrdecompvars

  ! table entries for interpolator
  do i=1,nrdecompvars
    op_indices(i)=i-1
    op_codes(i)=0
  end do
  !print*,'op stuff',op_indices,op_codes

  call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_indices, "operand_indices" )
  if(status.lt. 0) then
    call CCTK_WARN(0, "Cannot set operand indices array in parameter table")
  end if
  call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_codes, "operation_codes" )
  if (status.lt. 0) then
    call CCTK_WARN(0, "Cannot set operation codes array in parameter table")
  end if

  ! XXX Check that this works with ntheta<maxntheta and intp_npoints=ntheta*nphi! that would be cool, but it would need the memory layout to be handled nicely...
  intp_npoints=maxntheta*maxnphi

  call CCTK_InterpGridArrays(ierr, cctkGH, 3, interp_handle, &
                             table_handle, coord_system_handle, &
                             intp_npoints, CCTK_VARIABLE_REAL, &
                             interp_coords, num_in_arrays, in_array, &
                             num_out_arrays, out_types, out_array )
  if (ierr.ne.0) then
    call CCTK_WARN(1,"Interpolation failed")
  end if

  ! clean up potential NaN's in grid functions
  !  print*,"get rid of nan's"
  number_of_killed_nans=0
  do vind=1, nrdecompvars
    do j=1, cnphi
      do i=1, cntheta
        tmpval=interp_arrays(i,j,vind)
        nanind=0
        call ylm_isnan(nanind, tmpval)
        if (nanind>0) then
           interp_arrays(i,j,vind) = 0.d0
           number_of_killed_nans=number_of_killed_nans+1
        end if
      end do
    end do
  end do
  if (number_of_killed_nans>0) then
    write(infoline,'(a13,i10,a5)') 'Have removed ',number_of_killed_nans,' NaNs'
    call CCTK_INFO(infoline)
  end if

  ! clean up out of range stuff
  if (maxntheta.gt.cntheta) then
    if (verbose>1) then
       call CCTK_INFO("wipe unused theta points with zero")
    end if
    do vind=1, nrdecompvars
      do j=1, cnphi
        do i=1+cntheta, maxntheta
           interp_arrays(i,j,vind) = 0.d0
        end do
      end do
    end do
  end if
  if (maxnphi.gt.cnphi) then
    if (verbose>1) then
       call CCTK_INFO("wipe unused phi points with zero")
    end if
    do vind=1, nrdecompvars
      do j=cnphi+1, maxnphi
        do i=1, cntheta
           interp_arrays(i,j,vind) = 0.d0
        end do
      end do
    end do
  end if

  !  print*,"get rid of memory"
  call Util_TableDestroy(status, table_handle)
  if(status.lt. 0) then
    call CCTK_WARN(0, "Cannot destroy interpolator table")
  end if

end subroutine Ylm_Project_Sphere

subroutine Ylm_Project_Dummy(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

!  CCTK_INT :: im, jm, it, ip
  CCTK_INT :: interp_handle, table_handle, coord_system_handle

  character(len=200) :: interp
  CCTK_INT :: interp_len
  character(len=7) :: interp_order

!  CCTK_REAL :: dthetainv, dphiinv
!  CCTK_REAL :: dxdth, dxdph, dydth, dydph, dzdth, dzdph
!  CCTK_REAL :: gtt, gtp, gpp
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER out_array(nrdecompvars)
  CCTK_INT in_array(nrdecompvars)
  CCTK_INT out_types(nrdecompvars)
  CCTK_INT :: i !, vind, myproc

! temp shortcuts
!  CCTK_REAL :: tgxx, dgxx, &
!               tgxy, dgxy, &
!               tgxz, dgxz, &
!               tgyy, dgyy, &
!               tgyz, dgyz, &
!               tgzz, dgzz

  CCTK_INT op_indices(nrdecompvars)
  CCTK_INT op_codes(nrdecompvars)

CCTK_INT :: ierr,num_in_arrays,num_out_arrays,status 
!,j
!  CCTK_REAL :: rad,r2,ct,st,ct2,st2,cp,cp2,sp,sp2

!  CCTK_REAL :: lapse,cor_angle
!  CCTK_REAL :: tmpval
!  CCTK_INT :: nanind
  CCTK_INT :: intp_npoints
!  CCTK_INT :: number_of_killed_nans

!  CCTK_REAL :: tmpx,tmpy,tmpz

  character(len=80) :: infoline

! ________________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO('Dummy interpolating')

  if (do_nothing == 1) &
    return

  !print*,'current detector ',current_detector,detector_radius(current_detector)

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

  call CCTK_CoordSystemHandle ( coord_system_handle, "cart3d" )
  if ( coord_system_handle .lt. 0) then
    call CCTK_WARN( 0, "Cannot get handle for cart3d coordinate system. Forgot to activate an implementation providing coordinates?" )
  endif

  interp_coords(1) = CCTK_NullPointer()
  interp_coords(2) = CCTK_NullPointer()
  interp_coords(3) = CCTK_NullPointer()

  do i=1, nrdecompvars
    out_array(i)  = CCTK_NullPointer()
  end do
  ! FIXME : output arrays should be dynamic
  num_out_arrays=nrdecompvars

  ! set all out_types(:) to type real
  out_types=CCTK_VARIABLE_REAL

  do i=1,nrdecompvars
    in_array(i)=varindices(i)
  end do
  num_in_arrays=nrdecompvars

  ! table entries for interpolator
  do i=1,nrdecompvars
    op_indices(i)=i-1
    op_codes(i)=0
  end do
  !print*,'op stuff',op_indices,op_codes

  call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_indices, "operand_indices" )
  if(status.lt. 0) then
    call CCTK_WARN(0, "Cannot set operand indices array in parameter table")
  end if
  call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_codes, "operation_codes" )
  if (status.lt. 0) then
    call CCTK_WARN(0, "Cannot set operation codes array in parameter table")
  end if

  intp_npoints=0

  call CCTK_InterpGridArrays(ierr, cctkGH, 3, interp_handle, &
                             table_handle, coord_system_handle, &
                             intp_npoints, CCTK_VARIABLE_REAL, &
                             interp_coords, num_in_arrays, in_array, &
                             num_out_arrays, out_types, out_array )
  if (ierr.ne.0) then
    call CCTK_WARN(1,"Dummy interpolation failed")
  end if

  !  print*,"get rid of memory"
  call Util_TableDestroy(status, table_handle)
  if(status.lt. 0) then
    call CCTK_WARN(0, "Cannot destroy interpolator table")
  end if

end subroutine Ylm_Project_Dummy
