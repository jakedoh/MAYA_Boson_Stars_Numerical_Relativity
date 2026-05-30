! Initialisation of the level set function and various other things.
! $Header$

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EHFinder_Init_F(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT :: i, j, k, l, status
  CCTK_REAL, dimension(3) :: xp, xpt
  CCTK_REAL, dimension(3,3) :: txyz
  CCTK_REAL :: cosa, sina, cosb, sinb, cosc, sinc
  CCTK_REAL :: last_time
  CCTK_REAL :: theta, dtheta, thetamin, thetamax, r_el
  CCTK_REAL :: phi, dphi, phimin, phimax
  CCTK_REAL :: costh, sinth, cosph, sinph
  CCTK_INT, dimension(1) :: lsh, lbnd
  CCTK_INT, dimension(2) :: lsh2, lbnd2

! Get the size of the local grid.
  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)

! Initialize the last_time varible. Note the parameters should be
! chosen to be consistent with the run producing the numerical data.
! F.ex. if dt was 0.1 but the data was only stored every 4 iterations, the
! parameters should be chosen so that dt now is 0.4.
  last_time = abs(cctk_delta_time) * ( last_iteration_number + &
                cheat_iterations ) / saved_iteration_every

  if ( cheat == 1 ) then
  last_time = abs(cctk_delta_time) * ( last_iteration_number + &
                cheat_iterations ) / saved_iteration_every
    cctk_iteration = -cheat_iterations
  else
    last_time = abs(cctk_delta_time) * last_iteration_number &
                                     / saved_iteration_every
  end if
  cctk_time = last_time

! Allocate the logical array containing the flag determining if the
! corresponding level set should be re-initialized.

  if ( allocated(re_init_this_level_set) ) then
    deallocate ( re_init_this_level_set )
  end if
  allocate ( re_init_this_level_set(eh_number_level_sets) )
  if ( allocated(re_initialize_undone) ) then
    deallocate ( re_initialize_undone )
  end if
  allocate ( re_initialize_undone(eh_number_level_sets) )

  if ( evolve_generators .gt. 0 ) then

    if ( CCTK_EQUALS( domain, 'full' ) ) then
      thetamin = zero; thetamax = pi
      phimin = zero; phimax = 2*pi
    else if ( CCTK_EQUALS( domain, 'bitant') ) then
      if ( CCTK_EQUALS( bitant_plane, 'xy' ) ) then
        thetamin = zero; thetamax = half * pi
        phimin = zero; phimax = 2*pi
      else  if ( CCTK_EQUALS( bitant_plane, 'xz' ) ) then
        thetamin = zero; thetamax = pi
        phimin = zero; phimax = pi
      else
        thetamin = zero; thetamax = pi
        phimin = -half * pi; phimax = half * pi
      end if
    else if ( CCTK_EQUALS( domain, 'quadrant' ) ) then
      if ( CCTK_EQUALS( quadrant_direction, 'x' ) ) then
        thetamin = zero; thetamax = half * pi
        phimin = zero; phimax = pi
      else if ( CCTK_EQUALS( quadrant_direction, 'y' ) ) then
        thetamin = zero; thetamax = half * pi
        phimin = -half * pi; phimax = half * pi
      else
        thetamin = zero; thetamax = pi
        phimin = zero; phimax = half * pi
      end if
    else if ( CCTK_EQUALS( domain, 'octant' ) ) then
      thetamin = zero; thetamax = half * pi
      phimin = zero; phimax = half * pi
    end if
                   
    if ( CCTK_EQUALS( generator_distribution, 'line' ) ) then

      call CCTK_GrouplbndGN ( status, cctkGH, 1, lbnd, 'ehfinder::xg' )
      if ( status .lt. 0 ) then
        call CCTK_WARN ( 0, 'cannot get lower bounds for generator arrays' )
      end if
      call CCTK_GrouplshGN ( status, cctkGH, 1, lsh, 'ehfinder::xg' )
      if ( status .lt. 0 ) then
        call CCTK_WARN ( 0, 'cannot get local size for generator arrays' )
      end if

      if ( number_of_generators .eq. 1 ) then
        theta = half * ( thetamax - thetamin ) + thetamin
      else
        dtheta = ( thetamax - thetamin ) / ( number_of_generators - 1 )
      end if

    else if ( CCTK_EQUALS( generator_distribution, '2D array' ) ) then

      call CCTK_GrouplbndGN ( status, cctkGH, 2, lbnd2, 'ehfinder::xg2' )
      if ( status .lt. 0 ) then
        call CCTK_WARN ( 0, 'cannot get lower bounds for generator arrays' )
      end if
      call CCTK_GrouplshGN ( status, cctkGH, 2, lsh2, 'ehfinder::xg2' )
      if ( status .lt. 0 ) then
        call CCTK_WARN ( 0, 'cannot get local size for generator arrays' )
      end if

      if ( number_of_generators_theta .eq. 1 ) then
        theta = half * ( thetamax - thetamin ) + thetamin
      else
        dtheta = ( thetamax - thetamin ) / ( number_of_generators_theta )
      end if

      if ( number_of_generators_phi .eq. 1 ) then
        phi = half * ( phimax - phimin ) + phimin
      else
        dphi = ( phimax - phimin ) / ( number_of_generators_phi )
      end if

    end if
  end if

  do l = 1, eh_number_level_sets

!   If a sphere is requested...
    if ( CCTK_EQUALS( initial_f(l), 'sphere' ) ) then

!     Set up a sphere of radius initial_rad and translated 
!     (translate_x,translate_y,translate_z) away from the origin.
      f(:,:,:,l) = sqrt( ( x - translate_x(l) )**2 + &
                         ( y - translate_y(l) )**2 + &
                         ( z - translate_z(l) )**2 ) - initial_rad(l)
      if ( evolve_generators .gt. 0 ) then
        if ( CCTK_EQUALS( generator_distribution, 'line' ) ) then
          do i = 1, lsh(1)
            theta = thetamin + dtheta * ( i + lbnd(1) - 1 ) + half * dtheta
            xg(i,l) = initial_rad(l) * sin(theta) + translate_x(l)
            yg(i,l) = translate_y(l)
            zg(i,l) = initial_rad(l) * cos(theta) + translate_z(l)
          end do
        else if ( CCTK_EQUALS( generator_distribution, '2D array' ) ) then
          do j = 1, lsh2(2)
            phi = phimin + dphi * ( j + lbnd2(2) - 1 ) + half * dphi
            cosph = cos(phi); sinph = sin(phi)
            do i = 1, lsh2(1)
              theta = thetamin + dtheta * ( i + lbnd2(1) - 1 )
              costh = cos(theta); sinth = sin(theta)
              xg2(i,j,l) = initial_rad(l) * sinth * cosph + translate_x(l)
              yg2(i,j,l) = initial_rad(l) * sinth * sinph + translate_x(l)
              zg2(i,j,l) = initial_rad(l) * costh + translate_z(l)
            end do
          end do
        end if
      end if
    end if

!   If an ellipsoid is requested...
    if ( CCTK_EQUALS( initial_f(l), 'ellipsoid' ) ) then

!     Calculate sines and cosines of the rotation parameters.
      cosa = cos(rotation_alpha(l))
      sina = sin(rotation_alpha(l))
      cosb = cos(rotation_beta(l))
      sinb = sin(rotation_beta(l))
      cosc = cos(rotation_gamma(l))
      sinc = sin(rotation_gamma(l))

!     Set up the rotation matrix. The order is alpha around the z-axis,
!     beta around the y-axis and finally gamma around the x-axis.
      txyz(1,1) = cosa * cosb
      txyz(1,2) = sina * cosb
      txyz(1,3) = -sinb
      txyz(2,1) = cosa * sinb * sinc - sina * cosc
      txyz(2,2) = sina * sinb * sinc + cosa * cosc
      txyz(2,3) = cosb * sinc
      txyz(3,1) = cosa * sinb * cosc + sina * sinc
      txyz(3,2) = sina * sinb * cosc - cosa * sinc
      txyz(3,3) = cosb * cosc

!     Apply the rotations and translation for all points on the grid.
!     Even though at first glance it looks like the translation is done
!     first, the opposite is actually true.
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            xp(1) = x(i,j,k) - translate_x(l)
            xp(2) = y(i,j,k) - translate_y(l)
            xp(3) = z(i,j,k) - translate_z(l)
            xpt = matmul ( txyz, xp )
            f(i,j,k,l) = sqrt( xpt(1)**2 / initial_a(l)**2 + &
                             xpt(2)**2 / initial_b(l)**2 + &
                             xpt(3)**2 / initial_c(l)**2) - 1.0
          end do
        end do
      end do

      if ( evolve_generators .gt. 0 ) then
        if ( CCTK_EQUALS( generator_distribution, 'line' ) ) then
          do i = 1, lsh(1)
            theta = thetamin + dtheta * ( i + lbnd(1) - 1 )
            costh = cos(theta); sinth = sin(theta)
            r_el = sqrt ( one / ( sinth**2 / initial_a(l)**2 + &
                           costh**2 / initial_c(l)**2 ) )
            xp(1) = r_el * sinth + translate_x(l)
            xp(2) = translate_y(l)
            xp(3) = r_el * costh + translate_z(l)
            xpt = matmul ( txyz, xp )
            xg(i,l) = xpt(1)
            yg(i,l) = xpt(2)
            zg(i,l) = xpt(3)
          end do
        else if ( CCTK_EQUALS( generator_distribution, '2D array' ) ) then
          do j = 1, lsh2(2)
            phi = phimin + dphi * ( j + lbnd2(2) - 1 ) + half * dphi
            cosph = cos(phi); sinph = sin(phi)
            do i = 1, lsh2(1)
              theta = thetamin + dtheta * ( i + lbnd2(1) - 1 ) + half * dtheta
              costh = cos(theta); sinth = sin(theta)
              r_el = sqrt ( one / &
                            ( sinth**2*cosph**2 / initial_a(l)**2 + &
                              sinth**2*sinph**2 / initial_b(l)**2 + &
                              costh**2 / initial_c(l)**2 ) )
              xp(1) = r_el * sinth * cosph + translate_x(l)
              xp(2) = r_el * sinth * sinph + translate_y(l)
              xp(3) = r_el * costh + translate_z(l)
              xpt = matmul ( txyz, xp )
              xg2(i,j,l) = xpt(1)
              yg2(i,j,l) = xpt(2)
              zg2(i,j,l) = xpt(3) 
            end do
          end do
!          print*,xg2
!          print*
!          print*,yg2
!          print*
!          print*,zg2
!          stop
        end if
      end if
    end if

!   if an ovaloid of Cassini is requested...
    if ( CCTK_EQUALS( initial_f(l), 'cassini' ) ) then
      f(:,:,:,l) = (x**2+y**2+z**2)**2 + cas_a(l)**4 - &
            2*cas_a(l)**2*(x**2 - (y**2+z**2)) - cas_b(l)**4
    end if
  end do

! Initialise the internal mask.
  eh_mask = 0

  return
end subroutine EHFinder_Init_F

subroutine EHFinder_Init(CCTK_ARGUMENTS)

  use EHFinder_mod

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS

! Find the maximal grid spacing.
  delta = max ( CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2), CCTK_DELTA_SPACE(3) )

! Set up the value used in interiour inactive cells.
  ex_value = - ( one + shell_width ) * delta 

! Get handles for various reduction operations.
  call CCTK_ReductionArrayHandle ( max_handle, 'maximum' )
  if ( max_handle .lt. 0 ) then
    call CCTK_WARN(0,'Could not obtain a handle for maximum reduction')
  end if
  call CCTK_ReductionArrayHandle ( min_handle, 'minimum' )
  if ( min_handle .lt. 0 ) then
    call CCTK_WARN(0,'Could not obtain a handle for minimum reduction')
  end if
  call CCTK_ReductionArrayHandle ( sum_handle, 'sum' )
  if ( sum_handle .lt. 0 ) then
    call CCTK_WARN(0,'Could not obtain a handle for sum reduction')
  end if

end subroutine EHFinder_Init
