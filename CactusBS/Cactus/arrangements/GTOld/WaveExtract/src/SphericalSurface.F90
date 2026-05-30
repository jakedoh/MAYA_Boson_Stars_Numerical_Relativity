#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine  WaEx_Setup_SphericalSurface(CCTK_ARGUMENTS)

  use WavExtrConstants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT i,rind

  if (verbose>2) call CCTK_INFO("spherical surface setup")

  ! initial setup is done by spherical surface thorn.
  do i=1,number_of_detectors
    if(surface_index(i).eq.-1) then
      call CCTK_WARN(0,"this surface was not selected as a spherical surface")
      continue
    end if
    if(sf_valid(i)<=0) then
      print*,'invalid: i=',i,'surface_index',surface_index(i)
      call CCTK_WARN(1,"surface is invalid from sf_valid")
      continue
    end if
    detector_radius(i)=sf_mean_radius(surface_index(i)+1)
    ntheta(i)=sf_thorn_ntheta(surface_index(i)+1)
    nphi(i)=sf_thorn_nphi(surface_index(i)+1)
  end do
  number_of_detectors=i-1

  ! reshuffle detectors
  rind=1
  do i=1,number_of_detectors
    if (surface_index(i).eq.-two .or. surface_index(i).eq.-one) then
      continue
    end if
    detector_radius(rind)=sf_mean_radius(surface_index(i)+1)
    ntheta(rind)=sf_thorn_ntheta(surface_index(i)+1)
    nphi(rind)=sf_thorn_nphi(surface_index(i)+1)
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
  !print*,'max det no.', number_of_detectors
  !print*,'new detect arrays:',detector_radius
  !print*,'surf ind',surface_index

end subroutine  WaEx_Setup_SphericalSurface

