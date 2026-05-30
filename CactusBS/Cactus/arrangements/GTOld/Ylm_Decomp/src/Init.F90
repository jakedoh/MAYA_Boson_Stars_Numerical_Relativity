! Init.F90
!
! Initialization of the integral weights
!
! $Id: Init.F90,v 1.6 2004/10/07 13:00:19 herrmann Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine Ylm_Init_tempGF(CCTK_ARGUMENTS)

  use Ylm_Constants
  implicit none
  DECLARE_CCTK_ARGUMENTS

  ! Initialize temporary 3d storage to get max radius
  temp3d = zero

end subroutine Ylm_Init_tempGF

subroutine Ylm_Init(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, im, jm
  CCTK_INT :: ierr

  !________________________________________________________________________

  if (verbose>2) call CCTK_INFO("Ylm_Init")

  ! files initially don't exist
  ylm_files_exist=0
  first_detector_ever=1
  first_detector_this_time=1

  ! setup indicator arrays
  cmplgf = -one
  imagIndex = -one

  ! store the maximum number of detectors.
  ! some arrays use this information and the original parameter gets
  ! overwritten/readjusted
  max_det_no_param=number_of_detectors

  ! Get sum reduction operator handle
  call CCTK_ReductionArrayHandle ( sum_handle, 'sum' )
  if ( sum_handle .lt. 0 ) then
    call CCTK_WARN(0,'Could not obtain a handle for sum reduction')
  end if

  if (out_offset.ne.0) then
     call CCTK_WARN(1,"out_offset is deprecated - assuming 0")
  end if

  if (lmax .gt. 0) then
     call CCTK_INFO("using specific range given by lmax")
     l_mode=lmax
  end if

  ! allocate local 2D and 3D arrays
  call Ylm_Allocate_Arrays(CCTK_PASS_FTOF)

  ! finally precompute values for maxntheta and maxnphi - saves time for the
  ! default case of maxntheta,maxnphi
  cntheta=maxntheta
  cnphi=maxnphi
  use_precomputed_trigs=0
  call Ylm_Setup_Sphere(CCTK_PASS_FTOF)
  ctheta_max=ctheta
  cphi_max=cphi
  sintheta_max=sintheta
  costheta_max=costheta
  sinphi_max=sinphi
  cosphi_max=cosphi

  ! also precompute an array which holds the sw=-2 l,m,th,ph information
  ! again, this saves time for the default case
  call Ylm_precompute_swm2(CCTK_PASS_FTOF)

  ! assume we only have spin weight -2 and check in Setup_Vars if that is true
  use_precomputed_spwght=1
  only_spwght_m2=1

end subroutine Ylm_Init


subroutine Ylm_Allocate_Arrays(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT mmodes
  
  CCTK_REAL arrsiz
  character(len=2048) :: infoline

  ! note that this only gets run on start and after recovery!
  ! the memory should be setup only once per run
  call CCTK_INFO("allocating storage for local arrays")

  ! spin weight -2 precomputed
  mmodes=2*l_mode+1
  if (.not.allocated(swm2_ylm_pre)) then
     arrsiz=two*dble(maxntheta*maxnphi*mmodes*l_mode)*eight/(1024d0**2)
     write(infoline,'(A30,F6.2,A3)') 'size of precomputation array: ',arrsiz,' MB'
     call CCTK_INFO(infoline)
     allocate(swm2_ylm_pre(2,maxntheta,maxnphi,mmodes,l_mode))
  else
     call CCTK_WARN(1,"VERY BAD: swm2_ylm_pre is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  ! interpolation coords: maxntheta * maxnphi
  if (.not.allocated(interp_x)) then
     arrsiz=dble(maxntheta*maxnphi)*eight/(1024d0**2)
     write(infoline,'(A33,F6.2,A3)') 'size of maxntheta,maxnphi array: ',arrsiz,' MB'
     call CCTK_INFO(infoline)
     allocate(interp_x(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: interp_x is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(interp_y)) then
     allocate(interp_y(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: interp_y is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(interp_z)) then
     allocate(interp_z(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: interp_z is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if


  ! theta, phi and trig functions: maxntheta * maxnphi
  if (.not.allocated(ctheta)) then
     allocate(ctheta(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: ctheta is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(cphi)) then
     allocate(cphi(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: cphi is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(sintheta)) then
     allocate(sintheta(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: sintheta is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(costheta)) then
     allocate(costheta(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: costheta is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(sinphi)) then
     allocate(sinphi(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: sinphi is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(cosphi)) then
     allocate(cosphi(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: cosphi is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(ctheta_max)) then
     allocate(ctheta_max(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: ctheta_max is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(cphi_max)) then
     allocate(cphi_max(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: cphi_max is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(sintheta_max)) then
     allocate(sintheta_max(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: sintheta_max is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(costheta_max)) then
     allocate(costheta_max(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: costheta_max is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(sinphi_max)) then
     allocate(sinphi_max(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: sinphi_max is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

  if (.not.allocated(cosphi_max)) then
     allocate(cosphi_max(maxntheta,maxnphi))
  else
     call CCTK_WARN(1,"VERY BAD: cosphi_max is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if


  ! 3D arrays: maxntheta, maxnphi, number_of_variables
  if (.not.allocated(interp_arrays)) then
     arrsiz=dble(maxntheta*maxnphi*max_nr_vars)*eight/(1024d0**2)
     write(infoline,'(A45,F6.2,A3)') 'size of maxntheta,maxnphi,max_nr_vars array: ',arrsiz,' MB'
     call CCTK_INFO(infoline)
     allocate(interp_arrays(maxntheta,maxnphi,max_nr_vars))
  else
     call CCTK_WARN(1,"VERY BAD: interp_arrays is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if
  
  if (.not.allocated(integrand)) then
     allocate(integrand(maxntheta,maxnphi,max_nr_vars))
  else
     call CCTK_WARN(1,"VERY BAD: integrand is already allocated - THIS SHOULD NOT HAPPEN - CROSS YOUR FINGERS!")
  end if

end subroutine Ylm_Allocate_Arrays
