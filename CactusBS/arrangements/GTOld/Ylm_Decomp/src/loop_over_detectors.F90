#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


subroutine Ylm_loop_over_detectors(CCTK_ARGUMENTS)
  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: detl,i,myproc,nprocs
  CCTK_INT, dimension(number_of_detectors) :: live_detectors
  CCTK_INT :: number_of_live_detectors

  if (verbose>1) call CCTK_INFO("in loop over detectors")

  if (number_of_detectors.eq.0) then
    if (verbose>1) call CCTK_INFO("no detectors requested")
    return
  end if

  if (do_nothing .eq. 1) then
    if (verbose>1) call CCTK_INFO("not doing anything here")
    return
  end if

  first_detector_this_time=1

  if (first_detector_ever>0) then
     ylm_files_exist=0
  else
     ylm_files_exist=1
  end if

  ! only ever consider detectors for which output is requested
  ! collect live detectors in a new array 
  number_of_live_detectors = 0
  do i = 1,number_of_detectors
     if (out_every_det(i)<0) then
        print*,'out_every_det<0: i=',i,'out_every_det=',out_every_det(i),'detl=',i
        call CCTK_WARN(1,"detector out_every_det<0. shouldn't happen!")
     end if
     if (mod(cctk_iteration,out_every_det(i)).eq.0) then
       number_of_live_detectors = number_of_live_detectors+1
       live_detectors(number_of_live_detectors) = i
     else
       if (verbose>1) call CCTK_INFO("skip detector")
     end if
  end do

  ! no detector requires output, do nothing
  if (number_of_live_detectors .eq. 0) then
    return
  end if

  ! how many processors to use for the calculation?
  if (distribute_spheres.ne.0) then 
    nprocs=CCTK_NProcs(cctkGH)
  else
    nprocs=1
  end if
  myproc=CCTK_MyProc(cctkGH) ! counts from zero

  ! round number_of_live_detectors upwards to next multiple of nprocs
  if (mod(number_of_live_detectors,nprocs) .eq. 0) then ! already an integer multiple
    detl=number_of_live_detectors
  else                                                ! round up
    detl=number_of_live_detectors - mod(number_of_live_detectors,nprocs) + nprocs
  end if
  do i = 1,detl,nprocs ! we process nprocs spheres in parallel
     ! even if we do not do any calculation ourselves, we still have to take
     ! part in the global interpolation if any of the other processors do
     if (distribute_spheres.ne.0 .and. i+myproc.gt.number_of_live_detectors .or. &
         distribute_spheres.eq.0 .and. myproc .ne. 0) then
        if (verbose>1) call CCTK_INFO("dummy detector")
        call Ylm_Project_Dummy(CCTK_PASS_FTOF)
        cycle
     end if
    
     current_detector=live_detectors(i+myproc)

     if (verbose>1) print*,'analyse detector=',current_detector

     ! first set the current ntheta and nphi
     cntheta=maxntheta
     cnphi=maxnphi
     if (ntheta_det(current_detector).gt.maxntheta) then
        call CCTK_WARN(1,"ntheta_det is larger than maxntheta, using maxntheta instead")
        cntheta=maxntheta
     end if

     if (nphi_det(current_detector).gt.maxnphi) then
        call CCTK_WARN(1,"nphi_det is larger than maxnphi, using maxnphi instead")
        cnphi=maxnphi
     end if

     if (ntheta_det(current_detector).gt.0) then
        cntheta=ntheta_det(current_detector)
        if (verbose>1) then
           print*,'setting cntheta to',cntheta
        end if
     end if

     if (nphi_det(current_detector).gt.0) then
        cnphi=nphi_det(current_detector)
        if (verbose>1) then
           print*,'setting cnphi to',cnphi
        end if
     end if

     ! check that we can use the precomputed quantities
     if (cntheta.eq.maxntheta .and. cnphi.eq.maxnphi) then
        if (verbose>1) call CCTK_INFO("use precomputed values")
        use_precomputed_trigs=1
     else
        use_precomputed_trigs=0
     end if

     if (only_spwght_m2.eq.1 .and. cntheta.eq.maxntheta .and. cnphi.eq.maxnphi) then
        use_precomputed_spwght=1
     else 
        use_precomputed_spwght=0
     end if

     if (verbose>1) call CCTK_INFO("call Ylm_Setup_Sphere")
     call Ylm_Setup_Sphere(CCTK_PASS_FTOF)

     if (verbose>1) call CCTK_INFO("call project sphere")
     call Ylm_Project_Sphere(CCTK_PASS_FTOF)

     if (verbose>1) call CCTK_INFO("call integrate sphere")
     call Ylm_Integrate_Sphere(CCTK_PASS_FTOF)

     if (verbose>1) print*,'done with detector=',current_detector

     ! global flags for IO
     if (write_single_file_per_mode>0) ylm_files_exist=1     
     first_detector_ever=0
     first_detector_this_time=0

  end do

end subroutine
