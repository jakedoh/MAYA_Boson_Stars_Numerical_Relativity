
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


  subroutine  WavExtr_MomMass(CCTK_ARGUMENTS)

  use WavExtrConstants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS


  CCTK_REAL :: dtp,dtheta,dphi
  CCTK_REAL :: aux,aux2

  CCTK_INT :: ierr

  CCTK_INT :: int_var1,int_var2,int_var3,int_var4,int_var5,int_var6

  CCTK_INT :: num_out_vals, num_in_fields
  CCTK_REAL,dimension(6) :: out_vals

  character(len=80) :: infoline

! _________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO("Calculate Mass and Momentum Stuff")

  if (do_nothing == 1) &
    return

  if (cctk_iteration .ne. 0) then
    if (mod(cctk_iteration,out_every_det(current_detector)).ne.0) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    end if
  end if

  call CCTK_VarIndex ( int_var1, "WaveExtract::int_tmp1" )
  call CCTK_VarIndex ( int_var2, "WaveExtract::int_tmp2" )
  call CCTK_VarIndex ( int_var3, "WaveExtract::int_tmp3" )
  call CCTK_VarIndex ( int_var4, "WaveExtract::int_tmp4" )
  call CCTK_VarIndex ( int_var5, "WaveExtract::int_tmp5" )
  call CCTK_VarIndex ( int_var6, "WaveExtract::int_tmp6" )
  if ( min(int_var1,int_var2,int_var3,int_var4,int_var5,int_var6) .lt. 0 ) then
    call CCTK_WARN(0, 'Could not get index to one of the int_tmpX grid arrays')
  end if


  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon.ne.0) then
    dphi=two*pi
  end if

  dtp= dtheta*dphi
! print*,'dtheta,dphi,dtp',dtheta,dphi,dtp


  ! spherical parts of the metric
  ! note we compute gpp/sin^2(theta)
  int_tmp1 = sym_factor*weights*sintheta    *dtp* grr
  int_tmp2 = sym_factor*weights*sintheta    *dtp* gtt
  int_tmp3 = sym_factor*weights*sintheta    *dtp* dr_gtt
  int_tmp4 = sym_factor*weights*one/sintheta*dtp* gpp


  ! Different ways to compute the radius of extraction.
  ! "areal radius" is coordinate invariant and probably the "best"
  ! CactusEinstein/Extract uses "average Schwarzschild metric"

  ! we divide by sintheta and flag points which are in the range ntheta-maxntheta with
  ! zero
  where (abs(sintheta) .lt. 1.d-14)
     sintheta = one
  end where
  if (CCTK_EQUALS(rsch2_computation,"areal radius")) then
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

  ! actual integrals
  call CCTK_Reduce(ierr, cctkGH, -1, sum_handle, &
                   num_out_vals, CCTK_VARIABLE_REAL, &
                   out_vals, num_in_fields, &
                   int_var1,int_var2,int_var3, &
                   int_var4,int_var5,int_var6)
  if (ierr.ne.0) then
    call CCTK_WARN(1,"the reduction calculation of the Schwarzschild mass/radius/related failed")
  end if

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
  call NaNChecker_CheckVarsForNaN(ierr, cctkGH, 1, "waveextract::rsch2", &
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

end subroutine WavExtr_SchwarzMassRad
