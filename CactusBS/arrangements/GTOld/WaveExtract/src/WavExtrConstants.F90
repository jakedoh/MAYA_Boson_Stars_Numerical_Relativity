! Various useful constants, taken from EHFinder.
! $Id: WavExtrConstants.F90,v 1.1.1.1 2004/06/10 15:08:38 herrmann Exp $

#include "cctk.h"

module WavExtrConstants
  CCTK_REAL, parameter :: zero = 0.0d0
  CCTK_REAL, parameter :: one = 1.0d0
  CCTK_REAL, parameter :: two = 2.0d0
  CCTK_REAL, parameter :: three = 3.0d0
  CCTK_REAL, parameter :: four = 4.0d0
  CCTK_REAL, parameter :: eight = 8.0d0
  CCTK_REAL, parameter :: ten = 10.0d0
  CCTK_REAL, parameter :: sixteen = 16.0d0
  CCTK_REAL, parameter :: thirtytwo = 32.0d0
  CCTK_REAL, parameter :: half = 0.5d0
  CCTK_REAL, parameter :: onethird = one / three
  CCTK_REAL, parameter :: twothirds = two / three
  CCTK_REAL, parameter :: fourthirds = four / three
  CCTK_REAL, parameter :: quarter = 0.25d0
  CCTK_REAL, parameter :: tenth = 0.1d0
  CCTK_REAL, parameter :: huge = 1d23
  CCTK_REAL, parameter :: pi = 3.1415926535897932385d0
end module WavExtrConstants
