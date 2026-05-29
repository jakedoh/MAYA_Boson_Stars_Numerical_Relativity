#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine  Ylm_Setup_SphericalSurface(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT i,rind

  if (verbose>2) call CCTK_INFO("spherical surface setup")


  ! initial setup is done by spherical surface thorn.
  do i=1,number_of_detectors
    if(surface_index(i).eq.-1) then
       cycle
    end if
    detector_radius(i)=sf_mean_radius(surface_index(i)+1)
    !print*,'i',i,'detector_radius', detector_radius(i)
    if (detector_radius(i).eq.zero) then
      !print*,'warning that detector radius is useless'
      surface_index(i)=-two
      detector_radius(i)=-two
    end if
  end do
  number_of_detectors=i-1

  ! reshuffle detectors
  rind=1
  do i=1,number_of_detectors
    if (surface_index(i).eq.-two .or. surface_index(i).eq.-one) then
      cycle
    end if
    detector_radius(rind)=sf_mean_radius(surface_index(i)+1)
    surface_index(rind)=surface_index(i)
    rind=rind+1
  end do
  number_of_detectors=rind-1
  current_detector=rind-1
  if (number_of_detectors>0) then
    do_nothing=0
  else
    do_nothing=1
  end if

end subroutine  Ylm_Setup_SphericalSurface
