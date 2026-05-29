! $Header: /numrelcvs/HerrmannCVS/WaveExtract/src/ProjectSphere.F90,v 1.9 2005/07/28 20:34:27 herrmann Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine WavExtr_ProjectSphere(CCTK_ARGUMENTS)

  use WavExtrConstants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: im, jm, it, ip, j
  CCTK_INT :: interp_handle, table_handle, coord_system_handle

  character(len=200) :: interp
  CCTK_INT :: interp_len
  character(len=7) :: interp_order

  CCTK_INT, dimension(2) :: lsh,lbnd

  CCTK_REAL :: dtheta, dphi, dthetainv, dphiinv
  CCTK_REAL :: dxdth, dxdph, dydth, dydph, dzdth, dzdph
!  CCTK_REAL :: gtt, gtp, gpp
  CCTK_POINTER, dimension(3) :: interp_coords
  CCTK_POINTER, dimension(28) :: out_array
  CCTK_INT, dimension(7) :: in_array
  CCTK_INT, dimension(28) :: out_types

! temp shortcuts
  CCTK_REAL :: tgxx, dgxx, &
               tgxy, dgxy, &
               tgxz, dgxz, &
               tgyy, dgyy, &
               tgyz, dgyz, &
               tgzz, dgzz

  CCTK_REAL :: tgxx_t1, tgxx_t2, &
               tgxy_t1, tgxy_t2, &
               tgxz_t1, tgxz_t2, &
               tgyy_t1, tgyy_t2, &
               tgyz_t1, tgyz_t2, &
               tgzz_t1, tgzz_t2

  CCTK_INT :: di

  CCTK_INT,dimension(28) :: op_indices
  CCTK_INT,dimension(28) :: op_codes

! FIXME : make this the real lapse !
  CCTK_REAL :: lapse

  CCTK_REAL :: dr_met, r_t1, r_t2
  CCTK_INT,dimension(7) :: op_indices_tmp, op_codes_tmp
  CCTK_INT,dimension(6) :: op_indices_ncf_tmp, op_codes_ncf_tmp

  CCTK_INT :: ierr,num_in_arrays,num_out_arrays,status
  CCTK_REAL :: rad,r2,ct,st,ct2,st2,cp,cp2,sp,sp2
  CCTK_REAL :: cor_angle,siO,coO,cosO,sinO
  character(len=80) :: infoline
  CCTK_INT :: metric_string_len, metric_group_index, metric_firstvar_index
  character(len=200) :: local_metric, warnline


  integer :: i,k

  CCTK_REAL :: g11,g12,g13,g21,g22,g23,g31,g32,g33
  CCTK_REAL :: g11_n,g12_n,g13_n,g21_n,g22_n,g23_n,g31_n,g32_n,g33_n
  CCTK_REAL :: d1_g11,d1_g12,d1_g13,d1_g21,d1_g22,d1_g23,d1_g31,d1_g32,d1_g33
  CCTK_REAL :: d2_g11,d2_g12,d2_g13,d2_g21,d2_g22,d2_g23,d2_g31,d2_g32,d2_g33
  CCTK_REAL :: d3_g11,d3_g12,d3_g13,d3_g21,d3_g22,d3_g23,d3_g31,d3_g32,d3_g33
  CCTK_REAL :: dx_g11,dx_g12,dx_g13,dx_g21,dx_g22,dx_g23,dx_g31,dx_g32,dx_g33
  CCTK_REAL :: dy_g11,dy_g12,dy_g13,dy_g21,dy_g22,dy_g23,dy_g31,dy_g32,dy_g33
  CCTK_REAL :: dz_g11,dz_g12,dz_g13,dz_g21,dz_g22,dz_g23,dz_g31,dz_g32,dz_g33

  CCTK_REAL :: jac(3,3), det_jac

  CCTK_REAL :: dxcdxs, dxcdys, dycdxs, dycdys

  CCTK_REAL :: xv,yv

  CCTK_REAL :: dt
! ________________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO('Interpolating metric and derivatives onto sphere')

  if (do_nothing == 1) &
    return

  if (cctk_iteration .ne. 0) then
    if (mod(cctk_iteration,out_every_det(current_detector)).ne.0) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    end if
  end if

  current_detector_radius = detector_radius(current_detector)

  di=current_detector

  ! Sanity Check
  if (current_detector_radius .lt. 1.d-10) then
    do_nothing =1
    if (verbose>0) then
       print*,"current_detector_radius=",current_detector_radius
       print*,"current_detector=",current_detector
       call CCTK_WARN(1,"This should never happen: The detector radius is 0!")
    end if
    return
  end if


  if (verbose>1) then
    write(infoline,'(a24,i4,a8,g16.8)') 'Analysing Detector No.: ',current_detector, &
                                        ' Radius ',current_detector_radius
    call CCTK_INFO(infoline)
  end if

  ! local shape of the 2D grid arrays
  call CCTK_GrouplshGN(ierr, cctkGH, 2, lsh, "WaveExtract::surface_arrays")
  if ( ierr .lt. 0 ) then
    call CCTK_WARN(0, "cannot get local size for surface arrays")
  end if

  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon.ne.0) then
    dphi = two*pi
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


  ! For observers we use it's own coord system
  if (observers.ne.0) then
    ! FIXME : add observers coordinate system handle !
    call CCTK_WARN(1,"IMPLEMENT COORDINATE SYSTEM HANDLE FOR OBSERVERS")
    call CCTK_CoordSystemHandle ( coord_system_handle, "cart3d" )
    if ( coord_system_handle .lt. 0) then
      call CCTK_WARN( 0, "Cannot get handle for observers coordinate system. Forgot to activate an implementation providing coordinates?" )
    end if
  else 
    call CCTK_CoordSystemHandle ( coord_system_handle, "cart3d" )
    if ( coord_system_handle .lt. 0) then
      call CCTK_WARN( 0, "Cannot get handle for cart3d coordinate system. Forgot to activate an implementation providing coordinates?" )
    end if
  end if

  interp_coords(1) = CCTK_PointerTo(interp_x)
  interp_coords(2) = CCTK_PointerTo(interp_y)
  interp_coords(3) = CCTK_PointerTo(interp_z)

  out_array(1)  = CCTK_PointerTo(gxxi)
  out_array(2)  = CCTK_PointerTo(dx_gxxi)
  out_array(3)  = CCTK_PointerTo(dy_gxxi)
  out_array(4)  = CCTK_PointerTo(dz_gxxi)
  out_array(5)  = CCTK_PointerTo(gxyi)
  out_array(6)  = CCTK_PointerTo(dx_gxyi)
  out_array(7)  = CCTK_PointerTo(dy_gxyi)
  out_array(8)  = CCTK_PointerTo(dz_gxyi)
  out_array(9)  = CCTK_PointerTo(gxzi)
  out_array(10) = CCTK_PointerTo(dx_gxzi)
  out_array(11) = CCTK_PointerTo(dy_gxzi)
  out_array(12) = CCTK_PointerTo(dz_gxzi)
  out_array(13) = CCTK_PointerTo(gyyi)
  out_array(14) = CCTK_PointerTo(dx_gyyi)
  out_array(15) = CCTK_PointerTo(dy_gyyi)
  out_array(16) = CCTK_PointerTo(dz_gyyi)
  out_array(17) = CCTK_PointerTo(gyzi)
  out_array(18) = CCTK_PointerTo(dx_gyzi)
  out_array(19) = CCTK_PointerTo(dy_gyzi)
  out_array(20) = CCTK_PointerTo(dz_gyzi)
  out_array(21) = CCTK_PointerTo(gzzi)
  out_array(22) = CCTK_PointerTo(dx_gzzi)
  out_array(23) = CCTK_PointerTo(dy_gzzi)
  out_array(24) = CCTK_PointerTo(dz_gzzi)
  num_out_arrays=24
  if ( (CCTK_EQUALS(metric_type,"static conformal")) .and. &
       (use_conformal_factor(current_detector) .ne. 0) ) then
    out_array(25) = CCTK_PointerTo(psii)
    out_array(26) = CCTK_PointerTo(dx_psii)
    out_array(27) = CCTK_PointerTo(dy_psii)
    out_array(28) = CCTK_PointerTo(dz_psii)
    num_out_arrays=28
  end if

  ! set all out_types(1:28) to type real
  out_types=CCTK_VARIABLE_REAL

  ! find the cartesian coordinates for the interpolation points
  if (use_spherical_surface.ne.0) then
    if(verbose>2) then
      call CCTK_INFO("go for spherical surface")
    end if

    ! Find out the lower bounds of the distributed integration grid arrays.
    call CCTK_GrouplbndGN(ierr, cctkGH,2,lbnd,"waveextract::surface_integrands")
    if ( ierr .lt. 0 ) then
      call CCTK_WARN(0, "cannot get lower bounds for surface integrands")
    end if

    if (ntheta(di).ne.sf_thorn_ntheta(surface_index(current_detector)+1)) then
      call CCTK_WARN(1,"different theta surface for spherical surface")
    end if
    if (nphi(di).ne.sf_thorn_nphi(surface_index(current_detector)+1)) then
      call CCTK_WARN(1,"different phi surface for spherical surface")
    end if

    do ip = 1, lsh(2)
      do it = 1, lsh(1)
        ! FIXME : check the index computation for parallel case
        rad=sf_radius(it+lbnd(1),ip+lbnd(2),surface_index(current_detector)+1)
        !print*,'FIXME rad is ',rad,it,ip,lbnd(1),lbnd(2),current_detector,surface_index(current_detector)+1
        if (it+lbnd(1)>ntheta(current_detector) .or. &
            ip+lbnd(2)>nphi(current_detector) ) then
          interp_x(it,ip)=zero
          interp_y(it,ip)=zero
          interp_z(it,ip)=zero
        else
          interp_x(it,ip) = origin_x + rad * sintheta(it,ip) * cosphi(it,ip)
          interp_y(it,ip) = origin_y + rad * sintheta(it,ip) * sinphi(it,ip)
          interp_z(it,ip) = origin_z + rad * costheta(it,ip)
        end if
      end do
    end do
  else
    rad=current_detector_radius
    do ip = 1, lsh(2)
      do it = 1, lsh(1)
        interp_x(it,ip) = origin_x + rad * sintheta(it,ip) * cosphi(it,ip)
        interp_y(it,ip) = origin_y + rad * sintheta(it,ip) * sinphi(it,ip)
        interp_z(it,ip) = origin_z + rad * costheta(it,ip)
      end do
    end do
  end if

  !print*,'interp_x',interp_x
  !print*,'interp_y',interp_y
  !print*,'interp_z',interp_z


  call CCTK_FortranString(metric_string_len, &
       detector_metric(current_detector), &
       local_metric)
  call CCTK_GroupIndex(metric_group_index, local_metric(1:metric_string_len))
  call CCTK_ActiveTimeLevelsGI(ierr, cctkGH, metric_group_index)
  if (ierr < 1) then
    write(warnline,*) "The metric ", &
         local_metric(1:metric_string_len), " for detector number", &
         current_detector, " does not have storage on"
    call CCTK_WARN(0, warnline)
  end if
  call CCTK_FirstVarIndexI(metric_firstvar_index, metric_group_index)

  if (verbose > 1) then
    write(infoline,*) "Extracting with respect to metric ",&
         local_metric(1:metric_string_len)
    call CCTK_INFO(infoline)
  end if

  in_array(1) = metric_firstvar_index + 0
  in_array(2) = metric_firstvar_index + 1
  in_array(3) = metric_firstvar_index + 2
  in_array(4) = metric_firstvar_index + 3
  in_array(5) = metric_firstvar_index + 4
  in_array(6) = metric_firstvar_index + 5

  num_in_arrays=6

  ! table entries for interpolator
  op_indices = (/ 0, 0, 0, 0, &
                  1, 1, 1, 1, &
                  2, 2, 2, 2, &
                  3, 3, 3, 3, &
                  4, 4, 4, 4, &
                  5, 5, 5, 5, &
                  6, 6, 6, 6 /)

  ! derivative operation codes
#define DERIV(x) x
  op_codes = (/ 0, DERIV(1), DERIV(2), DERIV(3), &
                0, DERIV(1), DERIV(2), DERIV(3), &
                0, DERIV(1), DERIV(2), DERIV(3), &
                0, DERIV(1), DERIV(2), DERIV(3), &
                0, DERIV(1), DERIV(2), DERIV(3), &
                0, DERIV(1), DERIV(2), DERIV(3), &
                0, DERIV(1), DERIV(2), DERIV(3) /)


  ! Do the actual interpolation depending on conformal factor
  if(CCTK_EQUALS(metric_type,"static conformal")) then
    if (observers.ne.0) then
      call CCTK_VarIndex(in_array(7), "observers::obs_psi")
    else
      call CCTK_VarIndex(in_array(7), "staticconformal::psi")
    end if
    num_in_arrays=7

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

    call CCTK_InterpGridArrays(ierr, cctkGH, 3, interp_handle, &
                                table_handle, coord_system_handle, &
                                lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                                interp_coords, num_in_arrays, in_array, &
                                num_out_arrays, out_types, out_array )
    if (ierr.ne.0) then
      call CCTK_WARN(1,"Interpolation of metric failed")
    end if
  else
    call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_indices(1:num_out_arrays), "operand_indices" )
    if(status.lt. 0) then
      call CCTK_WARN(0, "Cannot set operand indices array in parameter table")
    end if
    call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_codes(1:num_out_arrays), "operation_codes" )
    if (status.lt. 0) then
      call CCTK_WARN(0, "Cannot set operation codes array in parameter table")
    end if

    call CCTK_InterpGridArrays(ierr, cctkGH, 3, interp_handle, &
                                table_handle, coord_system_handle, &
                                lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                                interp_coords, &
                                num_in_arrays, &
                                    in_array(1:num_in_arrays), &
                                num_out_arrays, &
                                    out_types(1:num_out_arrays), &
                                    out_array(1:num_out_arrays))
    if (ierr.ne.0) then
      call CCTK_WARN(1,"Interpolation of metric failed")
    end if
  end if

  ! 2nd and 3rd interpolation if we compute the dr_g?? terms numerically
  if (CCTK_EQUALS(metric_derivatives,"dr_gtt")) then
    if(CCTK_EQUALS(metric_type,"static conformal")) then
      op_indices_tmp = (/ 0, 1, 2, 3, 4, 5, 6 /)
      op_codes_tmp   = (/ 0, 0, 0, 0, 0, 0, 0 /)
      num_in_arrays = 7
      num_out_arrays = 7
      call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_indices_tmp(1:num_out_arrays), "operand_indices" )
      if(status.lt. 0) then
        call CCTK_WARN(0, "Cannot set operand indices array in parameter table")
      end if
      call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_codes_tmp(1:num_out_arrays), "operation_codes" )
      if (status.lt. 0) then
        call CCTK_WARN(0, "Cannot set operation codes array in parameter table")
      end if
   else
      op_indices_ncf_tmp = (/ 0, 1, 2, 3, 4, 5 /)
      op_codes_ncf_tmp   = (/ 0, 0, 0, 0, 0, 0 /)
      num_in_arrays = 6
      num_out_arrays = 6
      call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                               op_indices_ncf_tmp(1:num_out_arrays), "operand_indices" )
      if(status.lt. 0) then
        call CCTK_WARN(0, "Cannot set operand indices array in parameter table")
      end if
      call Util_TableSetIntArray ( status, table_handle, num_out_arrays, &
                                 op_codes_ncf_tmp(1:num_out_arrays), "operation_codes" )
      if (status.lt. 0) then
        call CCTK_WARN(0, "Cannot set operation codes array in parameter table")
      end if
  end if

    ! input arrays are reused from above!

    ! output arrays for _t1
    out_array(1)  = CCTK_PointerTo(gxxi_t1)
    out_array(2)  = CCTK_PointerTo(gxyi_t1)
    out_array(3)  = CCTK_PointerTo(gxzi_t1)
    out_array(4)  = CCTK_PointerTo(gyyi_t1)
    out_array(5)  = CCTK_PointerTo(gyzi_t1)
    out_array(6)  = CCTK_PointerTo(gzzi_t1)
    if(CCTK_EQUALS(metric_type,"static conformal")) then
      out_array(7)  = CCTK_PointerTo(psii_t1)
    end if

    ! take the same resolution in r as in theta/phi
    dr_met = two*pi*rad/ntheta(current_detector)
    print*,'dr_met:',dr_met

    ! coords at rad-dr
    r_t1=current_detector_radius -dr_met

    do ip = 1, lsh(2)
      do it = 1, lsh(1)
        interp_x(it,ip) = origin_x + r_t1 * sintheta(it,ip) * cosphi(it,ip)
        interp_y(it,ip) = origin_y + r_t1 * sintheta(it,ip) * sinphi(it,ip)
        interp_z(it,ip) = origin_z + r_t1 * costheta(it,ip)
      end do
    end do

    call CCTK_InterpGridArrays(ierr, cctkGH, 3, interp_handle, &
                               table_handle, coord_system_handle, &
                               lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                               interp_coords, &
                               num_in_arrays, &
                                 in_array(1:num_in_arrays), &
                               num_out_arrays, &
                                 out_types(1:num_out_arrays), &
                                 out_array(1:num_out_arrays))
    if (ierr.ne.0) then
      call CCTK_WARN(1,"interpolator fucked up")
      do_nothing=1
    end if

    ! output arrays for _t2
    out_array(1)  = CCTK_PointerTo(gxxi_t2)
    out_array(2)  = CCTK_PointerTo(gxyi_t2)
    out_array(3)  = CCTK_PointerTo(gxzi_t2)
    out_array(4)  = CCTK_PointerTo(gyyi_t2)
    out_array(5)  = CCTK_PointerTo(gyzi_t2)
    out_array(6)  = CCTK_PointerTo(gzzi_t2)
    if(CCTK_EQUALS(metric_type,"static conformal")) then
      out_array(7)  = CCTK_PointerTo(psii_t2)
    end if

    ! coords at rad+dr
    r_t2=current_detector_radius +dr_met
    do ip = 1, lsh(2)
      do it = 1, lsh(1)
        interp_x(it,ip) = origin_x + r_t2 * sintheta(it,ip) * cosphi(it,ip)
        interp_y(it,ip) = origin_y + r_t2 * sintheta(it,ip) * sinphi(it,ip)
        interp_z(it,ip) = origin_z + r_t2 * costheta(it,ip)
      end do
    end do

    call CCTK_InterpGridArrays(ierr, cctkGH, 3, interp_handle, &
                               table_handle, coord_system_handle, &
                               lsh(1) * lsh(2), CCTK_VARIABLE_REAL, &
                               interp_coords, &
                               num_in_arrays, &
                                 in_array(1:num_in_arrays), &
                               num_out_arrays, &
                                 out_types(1:num_out_arrays), &
                                 out_array(1:num_out_arrays))
   if (ierr.ne.0) then
      call CCTK_WARN(1,"interpolator fucked up")
      do_nothing=1
    end if   
  end if


  ! compute radial derivatives
  ! note that rad=current_detector_radius which we know
  ! FIXME : radial derivative for sphericalsurface case!!
  if (CCTK_EQUALS(metric_derivatives,"di_gtt")) then
    do ip=1, lsh(2)
      do it=1, lsh(1)
        dr_gxxi(it,ip) = one/rad *( (interp_x(it,ip)-origin_x)*dx_gxxi(it,ip) &
                                   +(interp_y(it,ip)-origin_y)*dy_gxxi(it,ip) &
                                   +(interp_z(it,ip)-origin_z)*dz_gxxi(it,ip))
        dr_gxyi(it,ip) = one/rad *( (interp_x(it,ip)-origin_x)*dx_gxyi(it,ip) &
                                   +(interp_y(it,ip)-origin_y)*dy_gxyi(it,ip) &
                                   +(interp_z(it,ip)-origin_z)*dz_gxyi(it,ip))
        dr_gxzi(it,ip) = one/rad *( (interp_x(it,ip)-origin_x)*dx_gxzi(it,ip) &
                                   +(interp_y(it,ip)-origin_y)*dy_gxzi(it,ip) &
                                   +(interp_z(it,ip)-origin_z)*dz_gxzi(it,ip))
        dr_gyyi(it,ip) = one/rad *( (interp_x(it,ip)-origin_x)*dx_gyyi(it,ip) &
                                   +(interp_y(it,ip)-origin_y)*dy_gyyi(it,ip) &
                                   +(interp_z(it,ip)-origin_z)*dz_gyyi(it,ip))
        dr_gyzi(it,ip) = one/rad *( (interp_x(it,ip)-origin_x)*dx_gyzi(it,ip) &
                                   +(interp_y(it,ip)-origin_y)*dy_gyzi(it,ip) &
                                   +(interp_z(it,ip)-origin_z)*dz_gyzi(it,ip))
        dr_gzzi(it,ip) = one/rad *( (interp_x(it,ip)-origin_x)*dx_gzzi(it,ip) &
                                   +(interp_y(it,ip)-origin_y)*dy_gzzi(it,ip) &
                                   +(interp_z(it,ip)-origin_z)*dz_gzzi(it,ip))
        if(CCTK_EQUALS(metric_type,"static conformal")) then
          dr_psii(it,ip) = one/rad *( (interp_x(it,ip)-origin_x)*dx_psii(it,ip) &
                                     +(interp_y(it,ip)-origin_y)*dy_psii(it,ip) &
                                     +(interp_z(it,ip)-origin_z)*dz_psii(it,ip))
        end if
      end do
    end do
  end if

  ! Convert to spherical coord. system
  ! Note that these equations not take the conformal factor into
  ! account. That has to be done afterwards !!
  ! FIXME : WARN, THESE EQUATIONS PROBABLY CANT WORK FOR PHYSICAL METRIC
  r2=rad*rad
  do ip = 1, lsh(2)
    do it = 1, lsh(1)
      ct = costheta(it,ip) ; ct2 = ct*ct
      st = sintheta(it,ip) ; st2 = st*st
      cp = cosphi(it,ip) ; cp2 = cp*cp
      sp = sinphi(it,ip) ; sp2 = sp*sp

      tgxx = gxxi(it,ip) ; dgxx = dr_gxxi(it,ip)
      tgxy = gxyi(it,ip) ; dgxy = dr_gxyi(it,ip)
      tgxz = gxzi(it,ip) ; dgxz = dr_gxzi(it,ip)
      tgyy = gyyi(it,ip) ; dgyy = dr_gyyi(it,ip)
      tgyz = gyzi(it,ip) ; dgyz = dr_gyzi(it,ip)
      tgzz = gzzi(it,ip) ; dgzz = dr_gzzi(it,ip)

      grr(it,ip) = st2*cp2*tgxx +st2*sp2*tgyy +ct2*tgzz &
               +two*( st2*cp*sp*tgxy +st*cp*ct*tgxz +st*ct*sp*tgyz)

      grt(it,ip) = rad*(st*cp2*ct*tgxx +two*st*ct*sp*cp*tgxy &
               +cp*(ct2-st2)*tgxz +st*sp2*ct*tgyy &
               +sp*(ct2-st2)*tgyz-ct*st*tgzz)

      grp(it,ip) = rad*st*(-st*sp*cp*tgxx -st*(sp2-cp2)*tgxy &
               -sp*ct*tgxz +st*sp*cp*tgyy +ct*cp*tgyz)

      gtt(it,ip) = r2*(ct2*cp2*tgxx +two*ct2*sp*cp*tgxy &
               -two*st*ct*cp*tgxz +ct2*sp2*tgyy &
               -two*st*sp*ct*tgyz +st2*tgzz)

      gtp(it,ip) = r2*st*(-cp*sp*ct*tgxx -ct*(sp2-cp2)*tgxy &
               +st*sp*tgxz +cp*sp*ct*tgyy -st*cp*tgyz)

      gpp(it,ip) = r2*st2*(sp2*tgxx -two*cp*sp*tgxy +cp2*tgyy)

      if (compute_EPrad_direct_from_metric>0) then
         dt = cctk_time - last_time

         ! compute time derivatives
         if (cctk_iteration>0) then
           dt_grr(it,ip) = (grr(it,ip)-grr_p(it,ip,current_detector))/dt
           dt_grt(it,ip) = (grt(it,ip)-grt_p(it,ip,current_detector))/dt
           dt_grp(it,ip) = (grp(it,ip)-grp_p(it,ip,current_detector))/dt
           dt_gtt(it,ip) = (gtt(it,ip)-gtt_p(it,ip,current_detector))/dt
           dt_gtp(it,ip) = (gtp(it,ip)-gtp_p(it,ip,current_detector))/dt
           dt_gpp(it,ip) = (gpp(it,ip)-gpp_p(it,ip,current_detector))/dt
         else
           dt_grr(it,ip) = zero
           dt_grt(it,ip) = zero
           dt_grp(it,ip) = zero
           dt_gtt(it,ip) = zero
           dt_gtp(it,ip) = zero
           dt_gpp(it,ip) = zero
         end if

         ! swap time levels
         grr_p(it,ip,current_detector)=grr(it,ip)
         grt_p(it,ip,current_detector)=grt(it,ip)
         grp_p(it,ip,current_detector)=grp(it,ip)
         gtt_p(it,ip,current_detector)=gtt(it,ip)
         gtp_p(it,ip,current_detector)=gtp(it,ip)
         gpp_p(it,ip,current_detector)=gpp(it,ip)
      end if

      if (CCTK_EQUALS(metric_derivatives,"di_gtt")) then
        dr_gtt(it,ip) = two/rad*gtt(it,ip) +r2*(ct2*cp2*dgxx &
                       +two*ct2*sp*cp*dgxy -two*st*ct*cp*dgxz &
                       +ct2*sp2*dgyy -two*st*sp*ct*dgyz +st2*dgzz)

        dr_gtp(it,ip) = two/rad*gtp(it,ip) +r2*st*(-cp*sp*ct*dgxx &
                       -ct*(sp2-cp2)*dgxy +st*sp*dgxz &
                       +cp*sp*ct*dgyy -st*cp*dgyz)

        dr_gpp(it,ip) = two/rad*gpp(it,ip) +r2*st2*(sp2*dgxx &
                       -two*cp*sp*dgxy+cp2*dgyy)
      else if (CCTK_EQUALS(metric_derivatives,"dr_gtt")) then
        if(CCTK_EQUALS(metric_type,"static conformal")) then
          tgxx_t1 = gxxi_t1(it,ip)*psii_t1(it,ip)**4 ; tgxx_t2 = gxxi_t2(it,ip)*psii_t2(it,ip)**4
          tgxy_t1 = gxyi_t1(it,ip)*psii_t1(it,ip)**4 ; tgxy_t2 = gxyi_t2(it,ip)*psii_t2(it,ip)**4
          tgxz_t1 = gxzi_t1(it,ip)*psii_t1(it,ip)**4 ; tgxz_t2 = gxzi_t2(it,ip)*psii_t2(it,ip)**4
          tgyy_t1 = gyyi_t1(it,ip)*psii_t1(it,ip)**4 ; tgyy_t2 = gyyi_t2(it,ip)*psii_t2(it,ip)**4
          tgyz_t1 = gyzi_t1(it,ip)*psii_t1(it,ip)**4 ; tgyz_t2 = gyzi_t2(it,ip)*psii_t2(it,ip)**4
          tgzz_t1 = gzzi_t1(it,ip)*psii_t1(it,ip)**4 ; tgzz_t2 = gzzi_t2(it,ip)*psii_t2(it,ip)**4
        else
          tgxx_t1 = gxxi_t1(it,ip) ; tgxx_t2 = gxxi_t2(it,ip)
          tgxy_t1 = gxyi_t1(it,ip) ; tgxy_t2 = gxyi_t2(it,ip)
          tgxz_t1 = gxzi_t1(it,ip) ; tgxz_t2 = gxzi_t2(it,ip)
          tgyy_t1 = gyyi_t1(it,ip) ; tgyy_t2 = gyyi_t2(it,ip)
          tgyz_t1 = gyzi_t1(it,ip) ; tgyz_t2 = gyzi_t2(it,ip)
          tgzz_t1 = gzzi_t1(it,ip) ; tgzz_t2 = gzzi_t2(it,ip)
        end if
        gtt_t1(it,ip)=r_t1**2*(ct2*cp2*tgxx_t1 +two*ct2*sp*cp*tgxy_t1 &
                               -two*st*ct*cp*tgxz_t1 +ct2*sp2*tgyy_t1 &
                               -two*st*sp*ct*tgyz_t1 +st2*tgzz_t1)
        gtt_t2(it,ip)=r_t2**2*(ct2*cp2*tgxx_t2 +two*ct2*sp*cp*tgxy_t2 &
                               -two*st*ct*cp*tgxz_t2 +ct2*sp2*tgyy_t2 &
                               -two*st*sp*ct*tgyz_t2 +st2*tgzz_t2)
        gtp_t1(it,ip) = r_t1**2*st*(-cp*sp*ct*tgxx_t1 -ct*(sp2-cp2)*tgxy_t1 &
                                 +st*sp*tgxz_t1 +cp*sp*ct*tgyy_t1 -st*cp*tgyz_t1)
        gtp_t2(it,ip) = r_t2**2*st*(-cp*sp*ct*tgxx_t2 -ct*(sp2-cp2)*tgxy_t2 &
                                 +st*sp*tgxz_t2 +cp*sp*ct*tgyy_t2 -st*cp*tgyz_t2)
        gpp_t1(it,ip) = r_t1**2*st2*(sp2*tgxx_t1 -two*cp*sp*tgxy_t1 +cp2*tgyy_t1)
        gpp_t2(it,ip) = r_t2**2*st2*(sp2*tgxx_t2 -two*cp*sp*tgxy_t2 +cp2*tgyy_t2)

        ! finite diff approx of dr_g??
        dr_gtt(it,ip) = 1.d0/(2.d0*dr_met)*(gtt_t1(it,ip)-gtt_t2(it,ip))
        dr_gtp(it,ip) = 1.d0/(2.d0*dr_met)*(gtp_t1(it,ip)-gtp_t2(it,ip))
        dr_gpp(it,ip) = 1.d0/(2.d0*dr_met)*(gpp_t1(it,ip)-gpp_t2(it,ip))
      end if
    end do
  end do


  ! WARNING : psi transforms like tensor density. that means like
  !           det(jacobian)=r^2 sin(theta)
  ! FIXME : trafo of psi
  ! Check that other formulas don't contain this implicitely

  ! Convert to physical quantities
  if(CCTK_EQUALS(metric_type,"static conformal")) then

    ! numerical derivatives already contain the conformal factor!
    if (CCTK_EQUALS(metric_derivatives,"di_gtt")) then
      psi_ext_deriv = four*psii**3*dr_psii

      ! Note derivatives first, because they depend on g?? !!
      ! the equation assumes those are without conf. factor.
      dr_gtt = psii**4 *dr_gtt + psi_ext_deriv *gtt
      dr_gtp = psii**4 *dr_gtp + psi_ext_deriv *gtp
      dr_gpp = psii**4 *dr_gpp + psi_ext_deriv *gpp
    end if

    grr = psii**4 * grr
    grt = psii**4 * grt
    grp = psii**4 * grp
    gtt = psii**4 * gtt
    gtp = psii**4 * gtp
    gpp = psii**4 * gpp

  end if


! set additional entries to 0
  do it=1,lsh(1)
    do ip=1,lsh(2)
       if (it+lbnd(1)>ntheta(current_detector) .or. &
            ip+lbnd(2)>nphi(current_detector) ) then
!          print*,'FIXME : reset values at ',it,ip,lbnd,ntheta(current_detector),nphi(current_detector)
          grr(it,ip)=zero
          grt(it,ip)=zero
          grp(it,ip)=zero
          gtt(it,ip)=zero
          gtp(it,ip)=zero
          gpp(it,ip)=zero
          dr_gtt(it,ip)=zero
          dr_gtp(it,ip)=zero
          dr_gpp(it,ip)=zero
       end if
    end do
 end do
          
end subroutine WavExtr_ProjectSphere

