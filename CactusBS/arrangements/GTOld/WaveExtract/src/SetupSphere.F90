! $Id: SetupSphere.F90,v 1.5 2005/01/14 17:36:19 herrmann Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine WavExtr_SetupSphere(CCTK_ARGUMENTS)

  use WavExtrConstants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT, dimension(2) :: lsh, lbnd

  CCTK_INT :: i, j

  CCTK_INT :: di

  CCTK_INT :: status

  CCTK_REAL :: ltheta, utheta, lphi, uphi
  CCTK_REAL :: dtheta, dphi

  character(len=256) :: infoline

  ! _________________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO("Setup Sphere")

  if (do_nothing==1) then
     if (verbose>1) call CCTK_INFO("do_nothing set for this detector")
     return
  end if

  ! get local shape of 2D grid arrays
  call CCTK_GrouplbndGN(status, cctkGH, 2, lbnd,"WaveExtract::surface_arrays")
  if ( status .lt. 0 ) then
    call CCTK_WARN(0, "cannot get lower bounds for surface arrays")
  end if

  call CCTK_GrouplshGN(status, cctkGH, 2, lsh, "WaveExtract::surface_arrays")
  if ( status .lt. 0 ) then
    call CCTK_WARN ( 0, "cannot get local size for surface arrays" )
  end if

  ! shorthand for current detector index
  di=current_detector

  ! Theta and phi setup. ltheta: lower theta, utheta: upper theta
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
  dtheta = ( utheta - ltheta ) / ntheta(di)
  dphi = ( uphi - lphi ) / nphi(di)
  if (cartoon.ne.0) then
    utheta = pi
    ltheta = zero
! FIXME : OTHER MODES??
    if (CCTK_EQUALS(domain,"bitant")) then
      dtheta = pi/(two*ntheta(di))
    else
      dtheta = pi/ntheta(di)
    end if
    dphi = two*pi
    uphi = zero
    lphi = zero
  end if

  if (verbose>1) then
    write(infoline,'(A5,G20.8,A5,G20.8)') 'dphi=',dphi,'nphi=',nphi(di)
    call CCTK_INFO(infoline)
    write(infoline,'(A7,G20.8,A7,G20.8)') 'dtheta=',dtheta,'ntheta=',ntheta(di)
    call CCTK_INFO(infoline)
  end if


  ! We stagger the origin, because we divide by sin(theta) and a division by sin(0) is
  ! not healthy.
  ! this is also useful to be able to use the "extended midpoint rule" for·
  ! integration.
  do i = 1, lsh(1)
    ctheta(i,:) = ltheta + dtheta * ( dble(i) + lbnd(1) - half )
  end do
  do j = 1, lsh(2)
    cphi(:,j) = lphi + dphi * ( dble(j) +lbnd(2) - half )
  end do

  ! FIXME : unneccessary or at least combine with normal case
  if (cartoon.ne.0) then
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
  ! (maxntheta,maxnphi), but we only use (ntheta,nphi)
  do j=1,lsh(2)
    do i=1,lsh(1)
       if (i+lbnd(1)>ntheta(current_detector) .or. &
            j+lbnd(2)>nphi(current_detector) ) then
           ctheta(i,j)=zero
           cphi(i,j)=zero
           sintheta(i,j)=zero
           costheta(i,j)=zero
           sinphi(i,j)=zero
           cosphi(i,j)=zero
        end if
     end do
  end do

end subroutine WavExtr_SetupSphere
