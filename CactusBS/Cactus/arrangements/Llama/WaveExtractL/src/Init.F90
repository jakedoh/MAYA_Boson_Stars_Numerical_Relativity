! Init.F90
!
! Initialization of the integral weights

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"


! This routine sets up the weights for the Simpson rule integration
! over the surface.
! taken from Peters EHFinder
subroutine WavExtrL_Init(CCTK_ARGUMENTS)

  use WavExtrLConstants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j
  CCTK_INT :: ierr

  CCTK_INT, dimension(2) :: lsh, lbnd

  !________________________________________________________________________

  ! Find out the lower bounds of the distributed integration grid arrays.
  call CCTK_GrouplbndGN(ierr, cctkGH,2,lbnd,"WaveExtractL::surface_integrands")
  if ( ierr .lt. 0 ) then
    call CCTK_WARN(0, "cannot get lower bounds for surface integrands")
  end if

  ! Find out the local size of the distributed integration grid arrays
  call CCTK_GrouplshGN(ierr, cctkGH, 2, lsh, "WaveExtractL::surface_integrands")
  if ( ierr .lt. 0 ) then
    call CCTK_WARN (0, "cannot get local size for surface integrands")
  end if

  ! set output iterations
  if (out_every>0) then
    my_out_every_det=out_every
  else
    my_out_every_det=out_every_det
  end if

  if (verbose > 4) then
    print*,'out_every_det',my_out_every_det
  end if


  ! initialize
  do_nothing = 0
  gxx_tmp = zero
  gxy_tmp = zero
  gxz_tmp = zero
  gyy_tmp = zero
  gyz_tmp = zero
  gzz_tmp = zero

  ! set ntheta and nphi for detectors which are not explicitely given in par-file
  do i=1,maximum_detector_number
    if (ntheta(i).eq.0 .and. use_spherical_surface .eq. 0) then
      int_ntheta(i)=maxntheta
    end if
    if (nphi(i).eq.0 .and. use_spherical_surface .eq. 0) then
      int_nphi(i)=maxnphi
    end if
  end do


  ! store the maximum number of detectors.
  ! some arrays use this information and the original parameter gets
  ! overwritten/readjusted
  max_det_no_param=maximum_detector_number

  if (CCTK_EQUALS(integration_method,"extended midpoint rule")) then
    ! we have to stagger the points, hence we can use the extended midpoint
    ! rule which just assigns a weight of 1 to each point. see Numerical
    ! Recipes p. 135 for details.
    ! Note that this method is only accurate up to O(1/N^2). usually that
    ! is good enough.
    weights = one
    phiweights = one
    thetaweights = one
  else if (CCTK_EQUALS(integration_method,"open extended")) then
    call CCTK_WARN(1,"VERY BAD CHOICE 'open extended'. this code is broken")
    ! we stagger the origin, so we need an open end formula.
    ! but we have the points at half values.
    ! FIXME: check the weights - these are not correct.
    ! FIXME: it won't work for maxntheta!=ntheta !!
    ! FIXME: BUGBUG: This code is broken, but extended midpoint is accurate enough.
    ! Initialise the weight grid array for the 2D Simpsons rule integration.
    ! To do this I need to figure out the global location of the given point.
    ! There are 3 cases in the one dimensional case. If the point is on the
    !  boundary the weight is 1/3. If it is at an even position the weight
    ! is 4/3 and if it is at an odd position the weight is 2/3. 

    weights = one
    do j = 1, lsh(2)

      ! This is first done in the phi direction. Meaning that all points with
      ! the same theta coordinate are set to the same weight.
      if ( ( lbnd(2)+j .eq. 1 ) .or. ( lbnd(2)+j .eq. maxnphi ) ) then
        weights(:,j) = onethird
        phiweights(:,j) = onethird
      else if ( mod(lbnd(2)+j,2) .eq. 0 ) then
        weights(:,j) = fourthirds
        phiweights(:,j) = fourthirds
      else
        weights(:,j) = twothirds
        phiweights(:,j) = twothirds
      end if

      ! Then it is done in the theta direction with the one-dimensional
      ! weights beeing multiplied.
      do i = 1, lsh(1)
        if ( ( lbnd(1)+i .eq. 1 ) .or. ( lbnd(1)+i .eq. maxntheta ) ) then
          weights(i,j) = onethird * weights(i,j)
          thetaweights(i,j) = onethird
        else if ( mod(lbnd(1)+i,2) .eq. 0 ) then
          weights(i,j) = fourthirds * weights(i,j)
          thetaweights(i,j) = fourthirds
        else
          weights(i,j) = twothirds * weights(i,j)
          thetaweights(i,j) = twothirds
        end if
      end do
    end do

    ! FIXME: DESCRIPTION is not accurate for open end formula, taken directly from Peters EHFinder.
    ! The end result is a 2D array with the weights in the following pattern.
    !    ie 2D Simpson _ WARN : 
    !  1/9   4/9  2/9   4/9  2/9   4/9  1/9
    !  4/9  16/9  8/9  16/9  8/9  16/9  4/9
    !  2/9   8/9  4/9   8/9  4/9   8/9  2/9
    !  4/9  16/9  8/9  16/9  8/9  16/9  4/9
    !  2/9   8/9  4/9   8/9  4/9   8/9  2/9
    !  4/9  16/9  8/9  16/9  8/9  16/9  4/9
    !  1/9   4/9  2/9   4/9  2/9   4/9  1/9

  end if

  ! setup the sym_factor to account for symmetries in the integrals. 
  ! the integrals are multiplied by this.
  ! default is full mode
  sym_factor=one
  if (cartoon .ne. 0) then
    if (CCTK_EQUALS(domain,"bitant")) then
      sym_factor=four
    else if (CCTK_EQUALS(domain,"full")) then
      sym_factor=two
    end if
  else if (CCTK_EQUALS(domain,"bitant")) then
    sym_factor=two
  else if (CCTK_EQUALS(domain,"quadrant")) then
    sym_factor=four
  else if (CCTK_EQUALS(domain,"octant")) then
    sym_factor=eight
  else if (CCTK_EQUALS(domain,"quadrant")) then
    sym_factor=four
  end if

  ! Get sum reduction operator handle
  call CCTK_ReductionArrayHandle ( sum_handle, 'sum' )
  if ( sum_handle .lt. 0 ) then
    call CCTK_WARN(0,'Could not obtain a handle for sum reduction')
  end if


end subroutine WavExtrL_Init
