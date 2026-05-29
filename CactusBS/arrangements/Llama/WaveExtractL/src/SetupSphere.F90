#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine WavExtrL_SetupSphere(CCTK_ARGUMENTS)

  use WavExtrLConstants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, dimension(2) :: lsh, lbnd

  CCTK_INT :: i, j, l

  CCTK_INT :: di

  CCTK_INT :: status

  CCTK_REAL :: ltheta, utheta, lphi, uphi
  CCTK_REAL :: dtheta, dphi

  CCTK_REAL :: theta, sum1

  character(len=256) :: infoline

  ! _________________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO("Setup Sphere")

  if (do_nothing==1) &
    return

  ! get local shape of 2D grid arrays
  call CCTK_GrouplbndGN(status, cctkGH, 2, lbnd,"WaveExtractL::surface_arrays")
  if ( status .lt. 0 ) then
    call CCTK_WARN(0, "cannot get lower bounds for surface arrays")
  end if

  call CCTK_GrouplshGN(status, cctkGH, 2, lsh, "WaveExtractL::surface_arrays")
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, "cannot get local size for surface arrays" )
  end if

  ! shorthand for current detector index
  di=current_detector

  
  ! set ntheta and nphi for detectors which are not explicitely given in par-file
  !int_ntheta = 0
  !int_nphi = 0
  do i=1,maximum_detector_number
    if (ntheta(i).eq.0 .and. use_spherical_surface .eq. 0) then
      int_ntheta(i)=maxntheta
    else if (use_spherical_surface .eq. 0) then
      int_ntheta(i)=ntheta(i)
    end if
    if (nphi(i).eq.0 .and. use_spherical_surface .eq. 0) then
      int_nphi(i)=maxnphi
    else if (use_spherical_surface .eq. 0) then
      int_nphi(i)=nphi(i)
    end if
    !print*,'WaveExtractL: ntheta.',int_ntheta(i)
    !print*,'WaveExtractL: nphi.',int_nphi(i)
    !print*,'WaveExtractL: detector_radius.',detector_radius(i)
  end do

  

  ! setup weights

  if (CCTK_EQUALS(integration_method,"extended midpoint rule")) then
    ! we have to stagger the points, hence we can use the extended midpoint
    ! rule which just assigns a weight of 1 to each point. see Numerical
    ! Recipes p. 135 for details.
    ! Note that this method is only accurate up to O(1/N^2). usually that
    ! is good enough.
    weights = one
    phiweights = one
    thetaweights = one
  else if (CCTK_EQUALS(integration_method,"Gauss")) then
    weights = one
    phiweights = one
    thetaweights = one
    
    do i=1, lsh(1)
     theta = pi*((lbnd(1)+i)-0.5d0)/int_ntheta(di)
     sum1 = 0d0
     do l=0, (int_ntheta(di)-1)/2
        sum1 = sum1 + sin((2*l+1)*theta)/(2*l+1)
     end do
     weights(i,:) = 8d0/(2d0*pi)*sum1 !8d0/int_nphi(di)*pi/int_ntheta(di)*sum1/dtheta/dphi
   end do
    
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
      if ( ( lbnd(2)+j .eq. 1 ) .or. ( lbnd(2)+j .eq. int_nphi(di) ) ) then
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
        if ( ( lbnd(1)+i .eq. 1 ) .or. ( lbnd(1)+i .eq. int_ntheta(di) ) ) then
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

  ! set the weights to zero outside the grid requested on the current sphere
  do j = 1, lsh(2)
     do i = 1, lsh(1)
       if ((lbnd(2)+j.gt.int_nphi(di)) .or. (lbnd(1)+i.gt.int_ntheta(di))) then
!         write (*,*) 'setting weights to zero at ', i, j
          weights(i,j) = 0
          phiweights(i,j) = 0
          thetaweights(i,j) = 0
       end if
    end do
 end do



  ! Theta and phi setup.
  ! Full mode is the default
  ltheta = zero; utheta = pi
  lphi = zero; uphi = two * pi
  if(CCTK_EQUALS(domain,"bitant")) then
    if (CCTK_EQUALS(bitant_plane,"xy")) then
      ltheta = zero; utheta = half * pi
    else
      ltheta = zero; utheta = pi
    end if
  else if(CCTK_EQUALS(domain, 'quadrant')) then
    if(CCTK_EQUALS(quadrant_direction, 'x')) then
      ltheta = zero; utheta = half * pi
      lphi = zero; uphi = pi
    else if(CCTK_EQUALS(quadrant_direction, 'y')) then
      ltheta = zero; utheta = half * pi
      lphi = zero; uphi = pi
    else if(CCTK_EQUALS(quadrant_direction, 'z')) then
      ltheta = zero; utheta = pi
      lphi = zero; uphi = half * pi
    else
      call CCTK_WARN(1,"unknown quadrant_direction")
    end if
  else if(CCTK_EQUALS(domain, 'octant')) then
    ltheta = zero; utheta = half * pi
    lphi = zero; uphi = half*pi
  end if

  ! Find dtheta and dphi and initialise the theta and phi grid arrays.
  ! Here i + lbnd(1) - 1 is the global index for theta and
  !      j + lbnd(2) - 1 is the global index for phi.
  dtheta = ( utheta - ltheta ) / int_ntheta(di)
  dphi = ( uphi - lphi ) / int_nphi(di)
  if (cartoon .ne. 0) then
    utheta = pi
    ltheta = zero
! FIXME : OTHER MODES??
    if (CCTK_EQUALS(domain,"bitant")) then
      dtheta = pi/(two*int_ntheta(di))
    else
      dtheta = pi/int_ntheta(di)
    end if
    dphi = two*pi
    uphi = zero
    lphi = zero
  end if

  if (verbose>1) then
    write(infoline,'(A5,G20.8,A5,G20.8)') 'dphi=',dphi,'nphi=',int_nphi(di)
    call CCTK_INFO(infoline)
    write(infoline,'(A7,G20.8,A7,G20.8)') 'dtheta=',dtheta,'ntheta=',int_ntheta(di)
    call CCTK_INFO(infoline)
  end if


  ! We stagger the origin, because we divide by sin(theta) and a division by sin(0) is
  ! not healthy.
  ! this is also useful to be able to use the "extended midpoint rule" for�
  ! integration.
  do i = 1, lsh(1)
    ctheta(i,:) = ltheta + dtheta * ( dble(i) + lbnd(1) - half )
  end do
  do j = 1, lsh(2)
    cphi(:,j) = lphi + dphi * ( dble(j) +lbnd(2) - half )
  end do

  ! FIXME : unneccessary or at least combine with normal case
  if (cartoon .ne. 0) then
    cphi = zero
  end if 

  if (verbose>7) then
    print*,'ctheta',ctheta
    print*,'cphi',cphi
  endif

  ! Calculate sines and cosines and store them in grid arrays since they
  ! are expensive.
  sintheta = sin(ctheta)
  costheta = cos(ctheta)
  sinphi = sin(cphi)
  cosphi = cos(cphi)

  ! zero the parts of the array which are not used. the arrays are of size 
  ! (maxntheta,maxnphi), but we only use (int_ntheta,int_nphi)
  do j=1,lsh(2)
    do i=1,lsh(1)
       if (i+lbnd(1)>int_ntheta(current_detector) .or. &
            j+lbnd(2)>int_nphi(current_detector) ) then
           ctheta(i,j)=zero
           cphi(i,j)=zero
           sintheta(i,j)=zero
           costheta(i,j)=zero
           sinphi(i,j)=zero
           cosphi(i,j)=zero
        end if
     end do
  end do

end subroutine WavExtrL_SetupSphere
