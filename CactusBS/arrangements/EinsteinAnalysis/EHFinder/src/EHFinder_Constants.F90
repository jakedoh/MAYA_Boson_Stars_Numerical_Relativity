! Various useful constants.
! $Header$

#include "cctk.h"

module EHFinder_Constants
  CCTK_INT, parameter :: izero = 0
  integer, parameter :: ik = kind(izero)

  CCTK_REAL, parameter :: zero = 0.0d0
  CCTK_REAL, parameter :: one = 1.0d0
  CCTK_REAL, parameter :: two = 2.0d0
  CCTK_REAL, parameter :: three = 3.0d0
  CCTK_REAL, parameter :: four = 4.0d0
  CCTK_REAL, parameter :: eight = 8.0d0
  CCTK_REAL, parameter :: half = 0.5d0
  CCTK_REAL, parameter :: onethird = one / three
  CCTK_REAL, parameter :: twothirds = two / three
  CCTK_REAL, parameter :: fourthirds = four / three
  CCTK_REAL, parameter :: quarter = 0.25d0
  CCTK_REAL, parameter :: pi = 3.1415926535897932385d0
  CCTK_REAL, parameter :: small = 1.0d-10
  integer, parameter :: wp = kind(zero)
end module EHFinder_Constants
