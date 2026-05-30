!/*@@
!  @file      Schwarzschild_Mass_Radius.F
!  @date      unknown
!  @author    unknown
!  @desc
!             Computes Schwarzschild Mass quantity and radius
!  @enddesc
!  @version   $Id: SchwarzschildMassRadius.F90,v 1.6 2005/08/09 14:45:06 herrmann Exp $
!@@*/



#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


!/*@@
!  @routine    WavExtr_Schw_Mass_Rad
!  @date       unknown
!  @author     unknown
!  @desc
!              Computes Schwarzschild Mass quantity and radius
!  @enddesc
!  @calls      spher_harm_combs
!  @@*/
  subroutine  WavExtr_SchwarzMassRad(CCTK_ARGUMENTS)

  use WavExtrConstants


  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS


  CCTK_REAL :: dtp,dtheta,dphi
  CCTK_REAL :: aux,aux2

  CCTK_INT :: ierr

  CCTK_INT :: int_var1,int_var2,int_var3,int_var4,int_var5,int_var6,int_var7,int_var8,int_var9,int_var10

  CCTK_INT :: num_out_vals, num_in_fields
  CCTK_REAL,dimension(10) :: out_vals
  CCTK_REAL,dimension(6) :: out_vals_6

  CCTK_INT :: cdt, create_file

  CCTK_REAL :: dt

  character(len=80) :: infoline

  CCTK_INT :: myproc

! _________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO("Calculate Schwarzschild Mass and Radius")

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
  call CCTK_VarIndex ( int_var7, "WaveExtract::int_tmp7" )
  call CCTK_VarIndex ( int_var8, "WaveExtract::int_tmp8" )
  call CCTK_VarIndex ( int_var9, "WaveExtract::int_tmp9" )
  call CCTK_VarIndex ( int_var10, "WaveExtract::int_tmp10" )

  if ( min(int_var1,int_var2,int_var3,int_var4,int_var5,int_var6,int_var7,int_var8,int_var9,int_var10) .lt. 0 ) then
    call CCTK_WARN(0, 'Could not get index to one of the int_tmpX grid arrays')
  end if

  dtheta = ctheta(2,1) - ctheta(1,1)
  dphi = cphi(1,2) - cphi(1,1)

  if (cartoon.ne.0) then
    dphi=two*pi
  end if

  dtp= dtheta*dphi

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

  if (compute_EPrad_direct_from_metric>0) then
    ! int_tmp7 = d^2 E /dtdOmega
    int_tmp7 = weights/(thirtytwo*pi)*((dt_gtp)**2+ one/four*(dt_gtt-dt_gpp)**2)*sintheta*dtp
    int_tmp8 = cosphi*sintheta* int_tmp7
    int_tmp9 = sinphi*sintheta* int_tmp7
    int_tmp10=        costheta* int_tmp7
    
    num_out_vals =1
    num_in_fields=10

    call CCTK_Reduce(ierr, cctkGH, -1, sum_handle, &
         num_out_vals, CCTK_VARIABLE_REAL, &
         out_vals, num_in_fields, &
         int_var1,int_var2,int_var3, &
         int_var4,int_var5,int_var6, &
         int_var7,int_var8,int_var9, &
         int_var10)
  else 
    num_out_vals =1
    num_in_fields=6

    call CCTK_Reduce(ierr, cctkGH, -1, sum_handle, &
         num_out_vals, CCTK_VARIABLE_REAL, &
         out_vals_6, num_in_fields, &
         int_var1,int_var2,int_var3, &
         int_var4,int_var5,int_var6)
    out_vals=zero
    out_vals(1)=out_vals_6(1)
    out_vals(2)=out_vals_6(2)
    out_vals(3)=out_vals_6(3)
    out_vals(4)=out_vals_6(4)
    out_vals(5)=out_vals_6(5)
    out_vals(6)=out_vals_6(6)
  end if

  if (ierr.ne.0) then
    call CCTK_WARN(1,"the reduction calculation of the Schwarzschild mass/radius/related failed")
  end if

  sph_grr    =out_vals(1)
  sph_gtt    =out_vals(2)
  sph_dr_gtt =out_vals(3)
  sph_gpp    =out_vals(4)
  rsch2      =out_vals(5)
  drsch_dri  =out_vals(6)

  if (compute_EPrad_direct_from_metric>0) then
     cdt=current_detector
     dEdt(cdt)    =out_vals(7)
     dPxdt(cdt)   =out_vals(8)
     dPydt(cdt)   =out_vals(9)
     dPzdt(cdt)   =out_vals(10)
     
     dt = cctk_time - last_time

     if (cctk_iteration>0) then
       Erad(cdt)       =Erad(cdt)+dt*dEdt(cdt)
       Pxrad(cdt)      =Pxrad(cdt)+dt*dPxdt(cdt)
       Pyrad(cdt)      =Pyrad(cdt)+dt*dPydt(cdt)
       Pzrad(cdt)      =Pzrad(cdt)+dt*dPzdt(cdt)
     else if (cctk_iteration.eq.0) then
       dEdt(cdt)    =zero
       dPxdt(cdt)   =zero
       dPydt(cdt)   =zero
       dPzdt(cdt)   =zero

       Erad(cdt)       =zero
       Pxrad(cdt)      =zero
       Pyrad(cdt)      =zero
       Pzrad(cdt)      =zero
     end if

     ! XXX useless output
     if (print_EPrad_to_stdout>0) then
        print*,"*************************************"
        print*,"cdt=",cdt,"rad=",current_detector_radius
        print*,"dt=",dt
        print*,"dEdt=",dEdt(cdt)
        print*,"dPxdt=",dPxdt(cdt)
        print*,"dPydt=",dPydt(cdt)
        print*,"dPzdt=",dPzdt(cdt)
        print*,"E=",Erad(cdt)
        print*,"Px=",Pxrad(cdt)
        print*,"Py=",Pyrad(cdt)
        print*,"Pz=",Pzrad(cdt)
        print*,"*************************************"
     end if
        
     ! last detector on first iteration
     if (cctk_iteration.eq.0) then
        create_file=1
     else
        create_file=0
     end if

     ! only write IO from first detector
     myproc = CCTK_MyProc(cctkGH)
     if (myproc.eq.0) then
        call write_radiated_energy_momentum(ierr,cdt,create_file,&
             current_detector_radius,cctk_time,&
             Erad(cdt),dEdt(cdt),&
             Pxrad(cdt),Pyrad(cdt),Pzrad(cdt),dPxdt(cdt),dPydt(cdt),dPzdt(cdt))
     end if
        
  end if

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

