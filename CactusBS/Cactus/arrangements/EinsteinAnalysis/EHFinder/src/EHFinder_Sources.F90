! Calculation of the sources for the level set function.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_Sources(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l
  CCTK_REAL :: idx, idy, idz
  CCTK_REAL :: psito4
  CCTK_REAL :: idetg, alp2, tmp1, tmp2, tmp3
  CCTK_REAL :: cfactor, ssign
  CCTK_REAL, dimension(3) :: cdx, dfup
  CCTK_REAL :: al, ar, bl, br, cl, cr
  CCTK_REAL :: alminus, alplus, blminus, blplus, clminus, clplus
  CCTK_REAL :: arminus, arplus, brminus, brplus, crminus, crplus

#include "include/physical_part.h"

! calculate 1/(2*delta) in each direction
  idx = half / CCTK_DELTA_SPACE(1)
  idy = half / CCTK_DELTA_SPACE(2)
  idz = half / CCTK_DELTA_SPACE(3)

! Set the sign depending on the surface direction.
  if ( CCTK_EQUALS ( surface_direction, 'outward' ) ) ssign = one
  if ( CCTK_EQUALS ( surface_direction, 'inward' ) ) ssign = -one

  do l = 1, eh_number_level_sets

    do k = kzl, kzr
      do j = jyl, jyr
        do i = ixl, ixr
!         Calculate the inverse of the 3-metric
# include "include/metric.h"
        end do
      end do
    end do

!   Calculate the derivatives of the level set function using the intrinsic
!   scheme. Note, this should never be used and may disappear in later versions.
    if ( CCTK_EQUALS ( upwind_type, 'intrinsic' ) ) then
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
#include "include/upwind_second2.h"
          end do
        end do
      end do
    end if
    
!   Calculate the derivatives of the level set function using shift upwinding.
    if ( CCTK_EQUALS ( upwind_type, 'shift' ) ) then
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr
#include "include/upwind_shift_second2.h"
          end do
        end do
      end do
    end if

!   If the three metric is the static conformal metric we convert the inverse
!   three metric to the physical inverse three metric by multiplying with
!   psi^(-4).
    if ( CCTK_EQUALS ( metric_type, 'static conformal' ) ) then
      do k = kzl, kzr
        do j = jyl,jyr
          do i = ixl, ixr
            if ( eh_mask(i,j,k,l) .ge. 0 ) then
              psito4 = psi(i,j,k)**(-4)
              g3xx(i,j,k) = g3xx(i,j,k) * psito4
              g3xy(i,j,k) = g3xy(i,j,k) * psito4
              g3xz(i,j,k) = g3xz(i,j,k) * psito4
              g3yy(i,j,k) = g3yy(i,j,k) * psito4
              g3yz(i,j,k) = g3yz(i,j,k) * psito4
              g3zz(i,j,k) = g3zz(i,j,k) * psito4
            end if
          end do
        end do
      end do
    end if

!   Calculate the derivatives of the level set using characteristic upwinding.
    if ( CCTK_EQUALS ( upwind_type, 'characteristic' ) ) then
      do k = kzl, kzr
        do j = jyl, jyr
          do i = ixl, ixr

!           We use centered derivatives to figure out which direction
!           to upwind in.
#include "include/upwind_characteristic_second2.h"
          end do
        end do
      end do
    end if

    do k = kzl, kzr
      do j = jyl, jyr
        do i = ixl, ixr

!         If the current point is active ...
          if ( eh_mask(i,j,k,l) .ge. 0 ) then

!           Square the lapse.
            alp2 = alp(i,j,k)**2
  
!           Calculate beta^i df_i.
            tmp1 = betax(i,j,k) * dfx(i,j,k,l) + &
                   betay(i,j,k) * dfy(i,j,k,l) + &
                   betaz(i,j,k) * dfz(i,j,k,l)

!           Calculate gamma^ij df_i df_j.
            tmp2 = g3xx(i,j,k) * dfx(i,j,k,l)**2 + &
                   g3yy(i,j,k) * dfy(i,j,k,l)**2 + &
                   g3zz(i,j,k) * dfz(i,j,k,l)**2 + &
                   two * ( g3xy(i,j,k) * dfx(i,j,k,l) * dfy(i,j,k,l) + &
                           g3xz(i,j,k) * dfx(i,j,k,l) * dfz(i,j,k,l) + &
                           g3yz(i,j,k) * dfy(i,j,k,l) * dfz(i,j,k,l) )

!           If the metric is positive definite ...
            if ( tmp2 .ge. zero ) then

!             Calculate the right hand side.
              sf(i,j,k,l) = tmp1 - ssign * sqrt ( alp2 * tmp2 )

!             If the lapse is negative we change the sign of the right hand
!             side function. This is done to be able to handle for example
!             Schwarzschild in isotropic coordinates with the isotropic
!             lapse.          
              sf(i,j,k,l) = sf(i,j,k,l) * sign ( one, alp(i,j,k) )
            else

!             Otherwise print a level 0 warning.
              call CCTK_WARN ( 0, '3-metric not positive definite: Stopping' )
            end if
          else
            sf(i,j,k,l) = zero
          end if
        end do
      end do
    end do 
  end do

  return
end subroutine EHFinder_Sources
