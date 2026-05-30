!/*@@
!  @file      Schwarzschild_Mass_Radius.F
!  @date      unknown
!  @author    unknown
!  @desc
!             Computes Schwarzschild Mass quantity and radius
!  @enddesc
!@@*/



#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


!/*@@
!  @routine    WavExtrL_Schw_Mass_Rad
!  @date       unknown
!  @author     unknown
!  @desc
!              Computes Schwarzschild Mass quantity and radius
!  @enddesc
!  @calls      spher_harm_combs
!  @@*/
  subroutine  WavExtrL_SchwarzMassRad(CCTK_ARGUMENTS)

  use WavExtrLConstants


  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS


  CCTK_REAL :: dtp,dtheta,dphi

  integer :: ierr, sumhandle, istat

  integer :: num_out_vals, num_in_fields, minus_one
  CCTK_REAL,dimension(6) :: out_vals, local_reduced_vals

  character(len=200) :: infoline

! _________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO("Calculate Schwarzschild Mass and Radius")

  if (do_nothing == 1) &
    return

  if (cctk_iteration .ne. 0) then
    if (mod(cctk_iteration,my_out_every_det(current_detector)).ne.0) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    end if
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

  ! Get sum reduction operator handle
  call CCTK_ReductionArrayHandle ( sum_handle, 'sum' )
  if ( sum_handle .lt. 0 ) then
    call CCTK_WARN(0,'Could not obtain a handle for sum reduction')
  end if


  if (size(ctheta,1)<2) call CCTK_WARN (0, "internal error")
  if (size(cphi,2)<2) call CCTK_WARN (0, "internal error")
  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon .ne. 0) then
    dphi=two*pi
  end if

  dtp= dtheta*dphi
  !print*,'dtheta,dphi,dtp',dtheta,dphi,dtp

  call CCTK_TimerStart(ierr,"Schwarzschild")

  ! spherical parts of the metric
  ! note we compute gpp/sin^2(theta)
  int_tmp1 = sym_factor*weights*sintheta    *dtp* grr
  int_tmp2 = sym_factor*weights*sintheta    *dtp* gtt
  int_tmp3 = sym_factor*weights*sintheta    *dtp* dr_gtt
  int_tmp4 = sym_factor*weights*one/sintheta*dtp* gpp


  ! Different ways to compute the radius of extraction.
  ! "aerial radius" is coordinate invariant and probably the "best"
  ! CactusEinstein/Extract uses "average Schwarzschild metric"

  ! we divide by sintheta and flag points which are in the range ntheta-maxntheta with
  ! zero
  where (sintheta .lt. 1.d-14)
     sintheta = one
  end where
  if (CCTK_EQUALS(rsch2_computation,"aerial radius")) then
    int_tmp5 = sym_factor*weights*dtp* sqrt(gtt*gpp -gtp**2)
  else if(CCTK_EQUALS(rsch2_computation,"average Schwarzschild metric")) then
    int_tmp5 = sym_factor*weights*sintheta*dtp* half*(gtt+gpp/sintheta**2)
  else if(CCTK_EQUALS(rsch2_computation,"Schwarzschild gtt")) then
    int_tmp5 = sym_factor*weights*dtp* gtt
  else if(CCTK_EQUALS(rsch2_computation,"Schwarzschild gpp")) then
    int_tmp5 = sym_factor*weights*dtp* gpp/sintheta**2
  end if

  ! FIXME : NO sintheta ??? check formula

  ! Compute the derivative of the schwarzschild radius
  ! with respect to the isotropic radius eta.
  ! It would be enough to just use dr_gtt. By adding in
  ! dr_gpp, the results become a bit better, but it still 
  ! assumes Schwarzschild coordinates.
  if (CCTK_EQUALS(drsch_dri_computation,"average dr_gtt dr_gpp")) then
    int_tmp6 = sym_factor*weights*sintheta*dtp* (dr_gtt+dr_gpp/sintheta**2)*half
  else if (CCTK_EQUALS(drsch_dri_computation,"dr_gtt")) then
    int_tmp6 = sym_factor*weights*sintheta*dtp* dr_gtt
  else if (CCTK_EQUALS(drsch_dri_computation,"dr_gpp")) then
    int_tmp6 = sym_factor*weights*sintheta*dtp* dr_gpp/sintheta**2
  end if

  num_out_vals =1
  num_in_fields=6
  sumhandle = sum_handle ! i.e., convert from CCTK_INT to integer
  minus_one = -1

  local_reduced_vals(1) = sum(int_tmp1,weights.gt.1.0e-15)
  local_reduced_vals(2) = sum(int_tmp2,weights.gt.1.0e-15)
  local_reduced_vals(3) = sum(int_tmp3,weights.gt.1.0e-15)
  local_reduced_vals(4) = sum(int_tmp4,weights.gt.1.0e-15)
  local_reduced_vals(5) = sum(int_tmp5,weights.gt.1.0e-15)
  local_reduced_vals(6) = sum(int_tmp6,weights.gt.1.0e-15)


  call CCTK_ReduceLocArrayToArray1D(ierr, cctkGH, minus_one,&
                        sumhandle, local_reduced_vals,&
                        out_vals, num_in_fields, CCTK_VARIABLE_REAL)

  if (ierr.ne.0) then
    call CCTK_WARN(1,"the reduction calculation of the Schwarzschild mass/radius/related failed")
  end if
  call CCTK_TimerStop(ierr,"Schwarzschild")

  sph_grr    =out_vals(1)
  sph_gtt    =out_vals(2)
  sph_dr_gtt =out_vals(3)
  sph_gpp    =out_vals(4)
  rsch2      =out_vals(5)
  drsch_dri  =out_vals(6)

  ! Normalizations
  ! FIXME DTAU_DT ADD THIS FUNCTIONALITY
  ! dtau_dt    = sqrt(one/(four*Pi)*dtau_dt)
  sph_grr    = sph_grr/(four*Pi)
  sph_gtt    = sph_gtt/(four*Pi)
  sph_dr_gtt = sph_dr_gtt/(four*Pi)
  sph_gpp    = sph_gpp/(four*Pi)
  rsch2      = rsch2/(four*Pi)

! try rsch2/=rsch2
! check for nan in computation of rsch2
  ierr=0
  call NaNChecker_CheckVarsForNaN(ierr, cctkGH, 1, "waveextractl::rsch2", &
                                                   "both","just warn")
  if (ierr /= 0) then
    call CCTK_WARN(1,"NaN in rsch2 - stopping all further computations")
    ! print*,'FIXME ',rsch2,ISNAN(rsch2)
    do_nothing=1
    return
  end if

  if (rsch2.lt.1.d-10) then
    call CCTK_WARN(1,"rsch2 < 10^-10 - stopping all further computations")
    ! print*,'FIXME rsch2',rsch2
    do_nothing=1
    return
  end if

  ! Schwarzschild radius
  rsch = sqrt( rsch2 )
  Schwarzschild_radius = rsch
  Schw_Radii(current_detector) = rsch

  ! dr_schwarzschild/dr_isotropic
  drsch_dri = one/(eight*Pi*rsch)*drsch_dri
  dri_drsch = one/drsch_dri

  ! Calculate the Schwarzschild mass parameter and S factor
  S_factor = drsch_dri**2/sph_grr
  Schwarzschild_Mass = rsch*(one-S_factor)/two
  Schw_Masses(current_detector) = Schwarzschild_Mass

  if (verbose > 3) then
    write(infoline,'(A25,G20.8)') '  rsch2                = ', rsch2
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  dtau_dt              = ',dtau_dt
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  sph_grr              = ',sph_grr
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  sph_gtt              = ',sph_gtt
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  sph_dr_gtt           = ',sph_dr_gtt
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  sph_gpp              = ',sph_gpp
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  drsch_dri            = ',drsch_dri
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  S_factor             = ',S_factor
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  Schwarzschild_radius = ', &
                                     Schwarzschild_radius
    call CCTK_INFO(infoline)
    write(infoline,'(A25,G20.8)') '  Schwarzschild_mass   = ', &
                                     Schwarzschild_Mass
    call CCTK_INFO(infoline)
  end if

end subroutine WavExtrL_SchwarzMassRad
