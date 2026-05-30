 /*@@
   @file      PPM_Scalars.F90
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

 module PPM_Scalars

   implicit none

   LOGICAL :: PPM3, PPM4
   LOGICAL :: ANALYTICAL
   LOGICAL :: FAST

 end module PPM_Scalars
