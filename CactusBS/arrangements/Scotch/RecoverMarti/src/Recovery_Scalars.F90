 /*@@
   @file      Recovery_Scalars.F90
   @date      Mon Feb 25 11:11:23 2002
   @author    
   @desc 
   Module containing various scalars to avoid having to use
   CCTK_EQUALS on keywords at every cell (slow).
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

 module Recovery_Scalars

   implicit none

   CCTK_REAL, dimension(3) :: radii, center_x, center_y, center_z 
   CCTK_REAL :: puncture_utau_min
   
 end module Recovery_Scalars
