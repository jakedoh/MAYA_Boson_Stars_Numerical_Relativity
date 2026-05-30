!/*@@
!  @file      Subtract_spherical_metric.F90
!  @date      unknown
!  @author    unknown
!  @desc
!             Subtract spherical background from metric
!  @enddesc
!  @version   $Id: SubtractSphericalMetric.F90,v 1.1.1.1 2004/06/10 15:08:38 herrmann Exp $
!  @@*/


#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

!/*@@
!  @routine    WavExtr_Subtr_spher_metric
!  @date       unknown
!  @author     unknown
!  @desc
!             Subtract spherical background from metric
!  @enddesc
!@@*/
subroutine WavExtr_SubtrSpherMetric(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
! _________________________________________________________________\

  if (verbose>2) &
    call CCTK_INFO("Subtract Spherical Background")

  if (do_nothing==1) &
    return

  if (cctk_iteration .ne. 0) then
    if (mod(cctk_iteration,out_every_det(current_detector)).ne.0) then
      if (verbose>2) call CCTK_INFO("No time for this detector")
      return
    end if
  end if

  grr = grr - sph_grr
  gtt = gtt - sph_gtt
  gpp = gpp - sph_gtt*sintheta**2
  dr_gtt = dr_gtt - sph_dr_gtt
  dr_gpp = dr_gpp - sph_dr_gtt*sintheta**2

end subroutine WavExtr_SubtrSpherMetric
