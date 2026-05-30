#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


subroutine WavExtr_loop_over_detectors(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: detl,i

  if (verbose>1) call CCTK_INFO("in loop over detectors")

  if (number_of_detectors.eq.0 .and. verbose>1) &
       call CCTK_INFO("no detectors requested")

  first_detector_this_time=1

  if (first_detector_ever>0) then
     wavextr_files_exist=0
  else
     wavextr_files_exist=1
  end if

  do_nothing=0
  detl=number_of_detectors
  do i = 1,detl
     current_detector=i

     if (verbose>1) print*,'analyse detector=',i

     if (verbose>1) call CCTK_INFO("call setup sphere")
     call WavExtr_SetupSphere(CCTK_PASS_FTOF)

     if (verbose>1) call CCTK_INFO("call project sphere")
     call WavExtr_ProjectSphere(CCTK_PASS_FTOF)

     if (verbose>1) call CCTK_INFO("call Schwarzschild Mass Rad")
     call WavExtr_SchwarzMassRad(CCTK_PASS_FTOF)

     if (subtract_spherical_background>0) then
        if (verbose>1) call CCTK_INFO("call subtract spherical background")
        call WavExtr_SubtrSpherMetric(CCTK_PASS_FTOF)
     end if

     if (verbose>1) call CCTK_INFO("call MoncriefQ")
     call WavExtr_MoncriefQ(CCTK_PASS_FTOF)

     if (verbose>1) print*,'done with detector=',i

     ! global flags for IO
     if (write_single_file_per_mode>0) wavextr_files_exist=1
     first_detector_ever=0
     first_detector_this_time=0
  end do

end subroutine
