! $Header: /numrelcvs/AEIDevelopment/WaveExtract/src/ProjectSphere.F90,v 1.12 2008/02/19 04:35:46 schnetter Exp $

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


  integer, save :: interpolation_timer = -1

  CCTK_INT :: it, ip, j
  CCTK_INT :: interp_handle, table_handle, coord_system_handle

  character(len=200) :: interp
  CCTK_INT :: interp_len
  character(len=20) :: interp_order

  CCTK_INT, dimension(2) :: lsh,lbnd

  CCTK_REAL :: dtheta, dphi, dthetainv, dphiinv
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

  CCTK_INT :: di

  CCTK_INT,dimension(28) :: op_indices
  CCTK_INT,dimension(28) :: op_codes

! FIXME : make this the real lapse !
  CCTK_REAL :: lapse

  CCTK_INT :: ierr,num_in_arrays,num_out_arrays,status, istat
  CCTK_REAL :: rad,r2,ct,st,ct2,st2,cp,cp2,sp,sp2
  CCTK_REAL :: cor_angle,siO,coO,cosO,sinO
  character(len=200) :: infoline
  CCTK_INT :: metric_string_len, metric_group_index, metric_firstvar_index
  character(len=200) :: local_metric, warnline


  integer :: i,k
  integer :: interpolator_error_level
  integer :: sf_i, sf_j, sf_k

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
    call CCTK_WARN(1,"This should never happen: The detector radius is 0!")
    return
  end if

  if (calc_when_necessary .eq. 1) then
    if (cctk_time .lt. current_detector_radius-50) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    endif
    call CCTK_IsFunctionAliased(istat, "MergerHandler_WeHaveMerger")
    if (istat .eq. 1) then
      if (MergerHandler_WeHaveMerger() .eq. 1) then
        if (cctk_time .gt. MergerHandler_MergerTime()+current_detector_radius+ringdown_margin) then
          if (verbose>2) call CCTK_INFO("No time for this detector")
          return
        endif
      endif
    endiF
  end if

  if (verbose>1) then
    write(infoline,'(a,i4,a,g16.8)') 'Analysing Detector No.: ',current_detector, &
                                        ' Radius ',current_detector_radius
    call CCTK_INFO(infoline)
  end if

  if (interpolation_timer < 0) then
     call CCTK_TimerCreate(interpolation_timer, "Interpolation")
     if (interpolation_timer < 0) then
        call CCTK_WARN(1,"could not create Interpolation timer")
     end if
  end if


  ! local shape of the 2D grid arrays
  call CCTK_GrouplshGN(ierr, cctkGH, 2, lsh, "WaveExtract::surface_arrays")
  if ( ierr .lt. 0 ) then
    call CCTK_WARN(0, "cannot get local size for surface arrays")
  end if

  if (size(ctheta,1)<2) call CCTK_WARN (0, "internal error")
  if (size(cphi,2)<2) call CCTK_WARN (0, "internal error")
  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon .ne. 0) then
    dphi = two*pi
  end if

  dthetainv = one / dtheta
  dphiinv = one / dphi


! rotate phi coordinate if forced
  if (phicorotate .ne. 0) then
    call CCTK_INFO("rotate phi by adding an angle")
    lapse=one

    ! find index for psi.
    ! FIXME assumes spherical symmetry for psi at detector radius
    find_r: do i=1+cctk_nghostzones(3),cctk_lsh(1)-cctk_nghostzones(3)
      if (abs(x(i,1,1))<rad) then
        print*,'i criterium',i,x(i,1,1),rad
        exit find_r
      end if
    end do find_r

    cor_angle=cctk_time*rotation_omega/lapse/psi(i,1,1)**3 ! FIXME : general eq
    print*,'time',cctk_time,'rot_ome',rotation_omega,'lapse',lapse
    print*,'i',i,'x(i)',x(i,1,1),'psi',psi(i,1,1),'psi^3',psi(i,1,1)**3
    print*,'cor_angle is...',cor_angle

!   undo rotation on phi itself.
    do j = 1, lsh(2)
      cphi(:,j) = cphi(:,j) - cor_angle
    end do

    sinphi = sin(cphi)
    cosphi = cos(cphi)
 end if
 

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

  if (make_interpolator_warnings_fatal .ne. 0) then
    interpolator_error_level = 0
  else
    interpolator_error_level = 1
  end if

  ! For observers we use it's own coord system
  if (observers .ne. 0) then
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

  call CCTK_TimerStart(ierr,"Interpolation")

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

  ! Find out the lower bounds of the distributed integration grid arrays.
  call CCTK_GrouplbndGN(ierr, cctkGH,2,lbnd,"waveextract::surface_integrands")
  if ( ierr .lt. 0 ) then
    call CCTK_WARN(0, "cannot get lower bounds for surface integrands")
  end if

  ! find the cartesian coordinates for the interpolation points
  if (use_spherical_surface .ne. 0) then
    if(verbose>2) then
      call CCTK_INFO("go for spherical surface")
    end if


    !if (int_ntheta(di).ne.sf_ntheta(surface_index(current_detector)+1)) then
    !  call CCTK_WARN(1,"different theta surface for spherical surface")
    !end if

    !if (int_nphi(di).ne.sf_nphi(surface_index(current_detector)+1)) then
    !  call CCTK_WARN(1,"different phi surface for spherical surface")
    !end if
    
    
    do ip = 1, lsh(2)
      do it = 1, lsh(1)
        ! FIXME : check the index computation for parallel case
        
        sf_k = surface_index(current_detector)+1
        
        sf_i = it+lbnd(1) + nghoststheta(sf_k)
        sf_j = ip+lbnd(2) + nghostsphi(sf_k)
        
        ! get proper sf_radius-index by taking into account symmetries
        if (symmetric_z(sf_k)/=0) then
           if (sf_i > sf_ntheta(sf_k)-nghoststheta(sf_k)) then
              sf_i = sf_ntheta(sf_k)-nghoststheta(sf_k) - (sf_i - (sf_ntheta(sf_k)-nghoststheta(sf_k)))
           end if
        end if
        if (symmetric_x(sf_k)/=0) then
           if (symmetric_y(sf_k)/=0) then
              if (sf_j > sf_nphi(sf_k)-nghostsphi(sf_k) .and. sf_j < 2*(sf_nphi(sf_k)-nghostsphi(sf_k))) then
                 sf_j = sf_nphi(sf_k)-nghostsphi(sf_k) - (sf_j - (sf_nphi(sf_k)-nghostsphi(sf_k)))
              else if (sf_j > 2*(sf_nphi(sf_k)-nghostsphi(sf_k))) then
                 sf_j = sf_j - 2*(sf_nphi(sf_k)-nghostsphi(sf_k))
              end if
           else
              if (sf_j > sf_nphi(sf_k)-nghostsphi(sf_k)) then
                 sf_j = sf_nphi(sf_k)-nghostsphi(sf_k) - (sf_j - (sf_nphi(sf_k)-nghostsphi(sf_k)))
              end if
           end if
        end if
        if (symmetric_y(sf_k)/=0) then
           if (sf_j > sf_nphi(sf_k)-nghostsphi(sf_k)) then
              sf_j = sf_nphi(sf_k)-nghostsphi(sf_k) - (sf_j - (sf_nphi(sf_k)-nghostsphi(sf_k)))
           end if
        end if
        
        rad=sf_radius(sf_i,sf_j,sf_k)
        if ((it+lbnd(1) <= int_ntheta(current_detector) .and. ip+lbnd(2) <= int_nphi(current_detector)) .and. rad .eq. 0) then 
           !print*,'FIXME rad is ',rad,it,ip,lbnd(1),lbnd(2),sf_i,sf_j,int_ntheta(current_detector),int_nphi(current_detector),current_detector,surface_index(current_detector)+1
           call CCTK_WARN(1, "Detector radius as given by SphericalSurface is 0!")
        endif
        
        if (it+lbnd(1)>int_ntheta(current_detector) .or. &
            ip+lbnd(2)>int_nphi(current_detector) ) then
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
    rad=current_detector_radius  ! set to mean radius because we need that later on  FIXME: May want to set this to areal radius!
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

  if (corotate3d .ne. 0) then
    ! FIXME : take lapse into account
    call CCTK_INFO("corotation undone in 3D")
    lapse=one
    cor_angle=cctk_time*rotation_omega/lapse
    print*,'cor_angle',cor_angle

    coO = cos(cor_angle)
    siO = sin(cor_angle)

    ! FIXME : BUG ADD OTHER TERMS - CHECK THIS AGAINST COMP ALG !!
    ! FIXME : CHECK LAPSE COROTATION
    call CCTK_WARN (0, "internal error -- cannot do this in global mode")
    do k = 1,  cctk_lsh(3)
      do j = 1,  cctk_lsh(2)
        do i = 1,  cctk_lsh(1)
          lapse=alp(i,j,k)
          cor_angle=cctk_time*rotation_omega/lapse
          gxx_tmp(i,j,k)=coO*coO*gxx(i,j,k)+coO*siO*gxy(i,j,k) &
                     +siO*coO*gxy(i,j,k)+siO*siO*gyy(i,j,k)
          gyy_tmp(i,j,k)=siO*siO*gxx(i,j,k)-siO*coO*gxy(i,j,k) &
                     -siO*coO*gxy(i,j,k)+coO*coO*gyy(i,j,k)
!         gzz is left unchanged
          gxy_tmp(i,j,k)=-coO*siO*gxx(i,j,k)+coO*coO*gxy(i,j,k) &
                   -siO*siO*gxy(i,j,k)+siO*coO*gyy(i,j,k)
          gxz_tmp(i,j,k)=coO*gxz(i,j,k)+siO*gyz(i,j,k)
          gyz_tmp(i,j,k)=-siO*gxz(i,j,k)+coO*gyz(i,j,k)
        end do
      end do
    end do

    call CCTK_VarIndex(in_array(1), "WaveExtract::gxx_tmp")
    call CCTK_VarIndex(in_array(2), "WaveExtract::gxy_tmp")
    call CCTK_VarIndex(in_array(3), "WaveExtract::gxz_tmp")
    call CCTK_VarIndex(in_array(4), "WaveExtract::gyy_tmp")
    call CCTK_VarIndex(in_array(5), "WaveExtract::gyz_tmp")
    call CCTK_VarIndex(in_array(6), "WaveExtract::gzz_tmp")
    num_in_arrays=6

    ! debug io
    ! call CCTK_OutputVarByMethod(ierr,cctkGH,"WaveExtract::gxx_tmp","IOASCII_1D")

  end if





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
    if (observers .ne. 0) then
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
      call CCTK_WARN(interpolator_error_level,"Interpolation of metric failed")
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
      call CCTK_WARN(interpolator_error_level,"Interpolation of metric failed")
    end if
  end if

  call CCTK_TimerStop(ierr,"Interpolation")
  call Util_TableDestroy ( status, table_handle )


  ! take into account the corotation and rotate again.
  if (corotate .ne. 0) then
    if (verbose >4) then
      if (size(gxxi,1)<5 .or. size(gxxi,2)<5) call CCTK_WARN (0, "internal error")
      print*,'rotation trafo gxxi before corotate:',gxxi(5,5)
      print*,'rotation trafo dx_gxxi before corotate:',dx_gxxi(5,5)
      print*,'rotation trafo dy_gxxi before corotate:',dy_gxxi(5,5)
      print*,'rotation trafo dz_gxxi before corotate:',dz_gxxi(5,5)
    end if

    ! FIXME : take lapse into account
    call CCTK_INFO("undoing corotation")
    lapse=one
    cor_angle=cctk_time*rotation_omega/lapse
    print*,'cor_angle',cor_angle

    coO = cos(cor_angle)
    siO = sin(cor_angle)

    cosO=coO
    sinO=siO


    ! FIXME : BUG ADD OTHER TERMS - CHECK THIS AGAINST COMP ALG !!
    do ip=1, lsh(2)
      do it=1, lsh(1)
        xv=interp_x(it,ip)
        yv=interp_y(it,ip)

        g11=gxxi(it,ip)
        g12=gxyi(it,ip)
        g13=gxzi(it,ip)
        g21=gxyi(it,ip)
        g22=gyyi(it,ip)
        g23=gyzi(it,ip)
        g31=gxzi(it,ip)
        g32=gyzi(it,ip)
        g33=gzzi(it,ip)

        g11_n=cosO*cosO*g11+sinO*cosO*g21+cosO*sinO*g12+sinO*sinO*g22
        g12_n=-cosO*sinO*g11-sinO*sinO*g21+cosO*cosO*g12+sinO*cosO*g22
        g13_n=cosO*g13+sinO*g23
        g21_n=-cosO*sinO*g11+cosO*cosO*g21-sinO*sinO*g12+sinO*cosO*g22
        g22_n=sinO*sinO*g11-sinO*cosO*g21-cosO*sinO*g12+cosO*cosO*g22
        g23_n=-sinO*g13+cosO*g23
        g31_n=cosO*g31+sinO*g32
        g32_n=-sinO*g31+cosO*g32
        g33_n=g33

        gxxi(it,ip)=g11_n
        gxyi(it,ip)=g12_n
        gxzi(it,ip)=g13_n
        gyyi(it,ip)=g22_n
        gyzi(it,ip)=g23_n
        gzzi(it,ip)=g33_n


!       FIXME : DO NOT FORGET CONFORMAL FACTOR
        dxcdxs=cosO
        dxcdys=sinO
        dycdxs=-sinO
        dycdys=cosO
!        dxcdxs=one
!        dxcdys=one
!        dycdxs=one
!        dycdys=one

!       initialize jacobi matrix
        jac=0.d0
        jac(3,3)=1.d0
        jac(1,1)=dxcdxs
        jac(1,2)=dxcdys
        jac(2,1)=dycdxs
        jac(2,2)=dycdys
        det_jac = - jac(1,3) * jac(2,2) * jac(3,1) &
                 &+ jac(1,2) * jac(2,3) * jac(3,1) &
                 &+ jac(1,3) * jac(2,1) * jac(3,2) &
                 &- jac(1,1) * jac(2,3) * jac(3,2) &
                 &- jac(1,2) * jac(2,1) * jac(3,3) &
                 &+ jac(1,1) * jac(2,2) * jac(3,3)

        psii(it,ip)= sqrt(det_jac)*psii(it,ip) 

!       partial derivatives
        dx_g11=dx_gxxi(it,ip)
        dx_g12=dx_gxyi(it,ip)
        dx_g13=dx_gxzi(it,ip)
        dx_g22=dx_gyyi(it,ip)
        dx_g23=dx_gyzi(it,ip)
        dx_g33=dx_gzzi(it,ip)
        dx_g21=dx_g12
        dx_g31=dx_g31
        dx_g32=dx_g32

        dy_g11=dy_gxxi(it,ip)
        dy_g12=dy_gxyi(it,ip)
        dy_g13=dy_gxzi(it,ip)
        dy_g22=dy_gyyi(it,ip)
        dy_g23=dy_gyzi(it,ip)
        dy_g33=dy_gzzi(it,ip)
        dy_g21=dy_g12
        dy_g31=dy_g31
        dy_g32=dy_g32

        dz_g11=dz_gxxi(it,ip)
        dz_g12=dz_gxyi(it,ip)
        dz_g13=dz_gxzi(it,ip)
        dz_g22=dz_gyyi(it,ip)
        dz_g23=dz_gyzi(it,ip)
        dz_g33=dz_gzzi(it,ip)
        dz_g21=dz_g12
        dz_g31=dz_g31
        dz_g32=dz_g32

!       FIXME : out that coords have be trafo'd...
        d1_g11=cosO*cosO*dx_g11+sinO*cosO*dx_g21+cosO*sinO*dx_g12+sinO*sinO*dx_g22
        d2_g11=cosO*cosO*dy_g11+sinO*cosO*dy_g21+cosO*sinO*dy_g12+sinO*sinO*dy_g22
        d3_g11=cosO*cosO*dz_g11+sinO*cosO*dz_g21+cosO*sinO*dz_g12+sinO*sinO*dz_g22

        d1_g12=-cosO*sinO*dx_g11-sinO*sinO*dx_g21+cosO*cosO*dx_g12+sinO*cosO*dx_g22
        d2_g12=-cosO*sinO*dy_g11-sinO*sinO*dy_g21+cosO*cosO*dy_g12+sinO*cosO*dy_g22
        d3_g12=-cosO*sinO*dz_g11-sinO*sinO*dz_g21+cosO*cosO*dz_g12+sinO*cosO*dz_g22

        d1_g13=cosO*dx_g13+sinO*dx_g23
        d2_g13=cosO*dy_g13+sinO*dy_g23
        d3_g13=cosO*dz_g13+sinO*dz_g23

        d1_g21=-cosO*sinO*dx_g11+cosO*cosO*dx_g21-sinO*sinO*dx_g12+sinO*cosO*dx_g22
        d2_g21=-cosO*sinO*dy_g11+cosO*cosO*dy_g21-sinO*sinO*dy_g12+sinO*cosO*dy_g22
        d3_g21=-cosO*sinO*dz_g11+cosO*cosO*dz_g21-sinO*sinO*dz_g12+sinO*cosO*dz_g22

        d1_g22=sinO*sinO*dx_g11-sinO*cosO*dx_g21-cosO*sinO*dx_g12+cosO*cosO*dx_g22
        d2_g22=sinO*sinO*dy_g11-sinO*cosO*dy_g21-cosO*sinO*dy_g12+cosO*cosO*dy_g22
        d3_g22=sinO*sinO*dz_g11-sinO*cosO*dz_g21-cosO*sinO*dz_g12+cosO*cosO*dz_g22

        d1_g23=-sinO*dx_g13+cosO*dx_g23
        d2_g23=-sinO*dy_g13+cosO*dy_g23
        d3_g23=-sinO*dz_g13+cosO*dz_g23

        d1_g31=cosO*dx_g31+sinO*dx_g32
        d2_g31=cosO*dy_g31+sinO*dy_g32
        d3_g31=cosO*dz_g31+sinO*dz_g32

        d1_g32=-sinO*dx_g31+cosO*dx_g32
        d2_g32=-sinO*dy_g31+cosO*dy_g32
        d3_g32=-sinO*dz_g31+cosO*dz_g32

        d1_g33=dx_g33
        d2_g33=dy_g33
        d3_g33=dz_g33

        dx_gxxi(it,ip)=dxcdxs*d1_g11+dxcdys*d2_g11
        dx_gxyi(it,ip)=dxcdxs*d1_g12+dxcdys*d2_g12
        dx_gxzi(it,ip)=dxcdxs*d1_g13+dxcdys*d2_g13
        dx_gyyi(it,ip)=dxcdxs*d1_g22+dxcdys*d2_g22
        dx_gyzi(it,ip)=dxcdxs*d1_g23+dxcdys*d2_g23
        dx_gzzi(it,ip)=dxcdxs*d1_g33+dxcdys*d2_g33

        dy_gxxi(it,ip)=dycdxs*d1_g11+dycdys*d2_g11
        dy_gxyi(it,ip)=dycdxs*d1_g12+dycdys*d2_g12
        dy_gxzi(it,ip)=dycdxs*d1_g13+dycdys*d2_g13
        dy_gyyi(it,ip)=dycdxs*d1_g22+dycdys*d2_g22
        dy_gyzi(it,ip)=dycdxs*d1_g23+dycdys*d2_g23
        dy_gzzi(it,ip)=dycdxs*d1_g23+dycdys*d2_g33

        dz_gxxi(it,ip)=d3_g11
        dz_gxyi(it,ip)=d3_g12
        dz_gxzi(it,ip)=d3_g13
        dz_gyyi(it,ip)=d3_g22
        dz_gyzi(it,ip)=d3_g23
        dz_gzzi(it,ip)=d3_g33

        dx_psii(it,ip)=sqrt(det_jac)*dx_psii(it,ip)
        dy_psii(it,ip)=sqrt(det_jac)*dy_psii(it,ip)
        dz_psii(it,ip)=sqrt(det_jac)*dz_psii(it,ip)

      end do
    end do

  end if



  ! compute radial derivatives
  ! note that rad=current_detector_radius which we know
  ! FIXME : radial derivative for sphericalsurface case!!
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

  ! Convert to spherical coord. system
  ! Note that these equations do not take the conformal factor into
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


      dr_gtt(it,ip) = two/rad*gtt(it,ip) +r2*(ct2*cp2*dgxx &
                     +two*ct2*sp*cp*dgxy -two*st*ct*cp*dgxz &
                     +ct2*sp2*dgyy -two*st*sp*ct*dgyz +st2*dgzz)

      dr_gtp(it,ip) = two/rad*gtp(it,ip) +r2*st*(-cp*sp*ct*dgxx &
                     -ct*(sp2-cp2)*dgxy +st*sp*dgxz &
                     +cp*sp*ct*dgyy -st*cp*dgyz)

      dr_gpp(it,ip) = two/rad*gpp(it,ip) +r2*st2*(sp2*dgxx &
                     -two*cp*sp*dgxy+cp2*dgyy)

    end do
  end do


  ! WARNING : psi transforms like tensor density. that means like
  !           det(jacobian)=r^2 sin(theta)
  ! FIXME : trafo of psi
  ! Check that other formulas don't contain this implicitely

  ! Convert to physical quantities
  if(CCTK_EQUALS(metric_type,"static conformal")) then
    psi_ext_deriv = four*psii**3*dr_psii

    ! Note derivatives first, because they depend on g?? !!
    ! the equation assumes those are without conf. factor.
    dr_gtt = psii**4 *dr_gtt + psi_ext_deriv *gtt
    dr_gtp = psii**4 *dr_gtp + psi_ext_deriv *gtp
    dr_gpp = psii**4 *dr_gpp + psi_ext_deriv *gpp

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
       if (it+lbnd(1)>int_ntheta(current_detector) .or. &
            ip+lbnd(2)>int_nphi(current_detector) ) then
!          print*,'FIXME : reset values at ',it,ip,lbnd,int_ntheta(current_detector),int_nphi(current_detector)
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

