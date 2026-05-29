#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine WavExtrL_Setup_SphericalSurface(CCTK_ARGUMENTS)

  use WavExtrLConstants

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT i,rind

  CCTK_INT, save :: firstcall = 1

  if (verbose>2) call CCTK_INFO("spherical surface setup")

  !if (firstcall .eq. 0) then 
!     return
 ! end if
  
  firstcall = 0

  ! initial setup is done by spherical surface thorn.
  do i=1,maximum_detector_number
    if(surface_index(i).eq.-1) then
      call CCTK_WARN(0,"this surface was not selected as a spherical surface")
      continue
    end if
    if(sf_valid(surface_index(i)+1)<=0) then
      call CCTK_WARN(1,"surface is invalid from sf_valid")
      continue
    end if
    detector_radius(i)=sf_mean_radius(surface_index(i)+1)
    
    int_ntheta(i)=sf_ntheta(surface_index(i)+1) - 2*nghoststheta(surface_index(i)+1)
    int_nphi(i)=sf_nphi(surface_index(i)+1) - 2*nghostsphi(surface_index(i)+1)
    if (symmetric_z(surface_index(i)+1)/=0) then
      int_ntheta(i) = int_ntheta(i)*2
    endif
    if (symmetric_x(surface_index(i)+1)/=0) then
      int_nphi(i) = int_nphi(i)*2
    endif
    if (symmetric_y(surface_index(i)+1)/=0) then
      int_nphi(i) = int_nphi(i)*2
    endif
    
    !print*,'WaveExtractL: int_ntheta.', int_ntheta(i)
    !print*,'WaveExtractL: int_nphi.', int_nphi(i)
    !print*,'WaveExtractL: r.', detector_radius(i)
    
  end do
  maximum_detector_number=i-1

  ! reshuffle detectors
!   rind=1
!   do i=1,maximum_detector_number
!     if (surface_index(i).eq.-two .or. surface_index(i).eq.-one) then
!       continue
!     end if
!     detector_radius(rind)=sf_mean_radius(surface_index(i)+1)
!     
!     ntheta(i)=sf_ntheta(surface_index(i)+1) - 2*nghoststheta(surface_index(i)+1)
!     nphi(i)=sf_nphi(surface_index(i)+1) - 2*nghostsphi(surface_index(i)+1)
!     if (symmetric_z(surface_index(i)+1)/=0) then
!       ntheta(i) = ntheta(i)*2
!     endif
!     if (symmetric_x(surface_index(i)+1)/=0) then
!       nphi(i) = nphi(i)*2
!     endif
!     if (symmetric_y(surface_index(i)+1)/=0) then
!       nphi(i) = nphi(i)*2
!     endif
!     
! !    print*,'WaveExtractL: ntheta.', ntheta(i)
! !    print*,'WaveExtractL: nphi.', nphi(i)
!     
!     surface_index(rind)=surface_index(i)
!     rind=rind+1
!   end do

!   maximum_detector_number=rind-1

   current_detector=0  !maximum_detector_number  !rind-1

   if (maximum_detector_number>0) then
     do_nothing=0
   else
     do_nothing=1
   end if
  !print*,'max det no.', maximum_detector_number
  !print*,'new detect arrays:',detector_radius
  !print*,'surf ind',surface_index

end subroutine  WavExtrL_Setup_SphericalSurface
