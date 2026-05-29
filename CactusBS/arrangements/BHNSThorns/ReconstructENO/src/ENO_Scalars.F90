 /*@@
   @file      ENO_Scalars.F90
   @date      Sat Apr  6 17:35:42 2002
   @author    Ian Hawke
   @desc 
   Module containing the coefficient array for ENO reconstruction.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 module ENO_Scalars

   implicit none

   CCTK_REAL, allocatable :: eno_coeffs(:, :)
   
   logical, save :: coeffs_allocated = .false.

 end module ENO_Scalars
