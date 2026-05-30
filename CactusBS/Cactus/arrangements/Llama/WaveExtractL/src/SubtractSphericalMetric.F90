!/*@@
!  @file      Subtract_spherical_metric.F90
!  @date      unknown
!  @author    unknown
!  @desc
!             Subtract spherical background from metric
!  @enddesc
!  @version   $Id: SubtractSphericalMetric.F90 55 2008-10-03 19:53:10Z reisswig $
!  @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

!/*@@
!  @routine    WavExtrL_Subtr_spher_metric
!  @date       unknown
!  @author     unknown
!  @desc
!             Subtract spherical background from metric
!  @enddesc
!@@*/
subroutine WavExtrL_SubtrSpherMetric(CCTK_ARGUMENTS)

  implicit none

  CCTK_INT :: istat

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
! _________________________________________________________________\

  if (verbose>2) &
    call CCTK_INFO("Subtract Spherical Background")

  if (do_nothing==1) &
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

  grr = grr - sph_grr
  gtt = gtt - sph_gtt
  gpp = gpp - sph_gtt*sintheta**2
  dr_gtt = dr_gtt - sph_dr_gtt
  dr_gpp = dr_gpp - sph_dr_gtt*sintheta**2

end subroutine WavExtrL_SubtrSpherMetric
