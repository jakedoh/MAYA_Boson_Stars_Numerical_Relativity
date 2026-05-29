 /*@@
   @file      Marquina_Scalars.F90
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

 module Marquina_Scalars

   implicit none

   LOGICAL :: ANALYTICAL
   LOGICAL :: FAST

 end module Marquina_Scalars
