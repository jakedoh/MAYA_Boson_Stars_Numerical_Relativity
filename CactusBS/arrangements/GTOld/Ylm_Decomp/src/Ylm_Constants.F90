! Various useful constants, taken from EHFinder.
! $Id: Ylm_Constants.F90,v 1.2 2003/11/24 14:48:55 herrmann Exp $

#include "cctk.h"

module Ylm_Constants
  implicit none

  ! keep all these variables around
  save
  
  CCTK_REAL, parameter :: zero = 0.0d0
  CCTK_REAL, parameter :: one = 1.0d0
  CCTK_REAL, parameter :: two = 2.0d0
  CCTK_REAL, parameter :: three = 3.0d0
  CCTK_REAL, parameter :: four = 4.0d0
  CCTK_REAL, parameter :: five = 5.0d0
  CCTK_REAL, parameter :: six = 6.0d0
  CCTK_REAL, parameter :: eight = 8.0d0
  CCTK_REAL, parameter :: ten = 10.0d0
  CCTK_REAL, parameter :: half = 0.5d0
  CCTK_REAL, parameter :: onethird = one / three
  CCTK_REAL, parameter :: twothirds = two / three
  CCTK_REAL, parameter :: fourthirds = four / three
  CCTK_REAL, parameter :: quarter = 0.25d0
  CCTK_REAL, parameter :: tenth = 0.1d0
  CCTK_REAL, parameter :: huge = 1d23
  CCTK_COMPLEX, parameter :: iunit=(0.d0,1.d0)
  CCTK_REAL, parameter :: pi = 3.1415926535897932385d0

  ! arrays:
  !   spin weight -2 pre computed values
  !   index : real/imag(1,2), theta, phi, m, l
  CCTK_REAL, dimension(:,:,:,:,:), allocatable :: swm2_ylm_pre

  !   2d surface arrays: theta,phi
  CCTK_REAL, dimension(:,:), allocatable ::   interp_x,interp_y,interp_z
  CCTK_REAL, dimension(:,:), allocatable ::   ctheta, cphi, sintheta, costheta, sinphi, cosphi
  CCTK_REAL, dimension(:,:), allocatable ::   ctheta_max, cphi_max, sintheta_max, costheta_max, sinphi_max, cosphi_max
  !   3d arrays: detector, theta, phi
  CCTK_REAL, dimension(:,:,:), allocatable :: interp_arrays, integrand                                 

end module Ylm_Constants
