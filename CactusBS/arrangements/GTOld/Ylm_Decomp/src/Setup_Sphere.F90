! $Id: Setup_Sphere.F90,v 1.4 2004/02/17 16:27:40 herrmann Exp $

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine Ylm_Setup_Sphere(CCTK_ARGUMENTS)

  use Ylm_Constants

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j

  CCTK_INT :: status

  CCTK_REAL :: ltheta, utheta, lphi, uphi
  ! _________________________________________________________________________

  if (verbose>2) &
    call CCTK_INFO("Setup Sphere")

  if (do_nothing==1) &
    return

  ! recompute the values
  if (use_precomputed_trigs==0) then

     ! Theta and phi setup.
     ! Full mode is the default - symmetries are now done by special thorns
     ltheta = zero; utheta = pi
     lphi = zero; uphi = two * pi
     
     ! Find dtheta and dphi and initialise the theta and phi grid arrays.
     dtheta = ( utheta - ltheta ) / (cntheta -1)
     dphi = ( uphi - lphi ) / (cnphi -1)
     
     ! We stagger the z-axis, because we divide by sin(theta) and a
     ! division by sin(0) is not healthy.
     ctheta(1,:)=1d-2
     ctheta(cntheta,:)=1d-2
     do i = 2, cntheta-1
        ctheta(i,:) = ltheta + dtheta * ( dble(i) - one)
     end do
     do j = 1, cnphi
        cphi(:,j) = lphi + dphi * ( dble(j) - one)
     end do
     
     ! Calculate sines and cosines and store them in grid arrays since they
     ! are expensive.
     sintheta = sin(ctheta)
     costheta = cos(ctheta)
     sinphi = sin(cphi)
     cosphi = cos(cphi)

  else
     ! precomputed values for maxntheta, maxnphi

     ctheta=ctheta_max
     cphi=cphi_max
     sintheta=sintheta_max
     costheta=costheta_max
     sinphi=sinphi_max
     cosphi=cosphi_max

  end if


end subroutine Ylm_Setup_Sphere
