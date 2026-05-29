 /*@@
   @file      Whisky_Startup.F90
   @date      Sun Feb 10 00:02:52 2002
   @author    Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

 /*@@
   @routine    Whisky_Startup
   @date       Sun Feb 10 00:03:09 2002
   @author     Ian Hawke
   @desc 
   Startup banner.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

integer function Whisky_Startup(CCTK_ARGUMENTS)

  USE Whisky_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_ARGUMENTS

  CCTK_INT :: ierr
  
  call CCTK_RegisterBanner(ierr, "Maya Scotch via Whisky: relativistic hydrodynamics, no ice.")

  whisky_mhd_handle = 0
  if ( CCTK_EQUALS(evolution_method,"Scotch_QuasiMHD") ) then
     whisky_mhd_handle = 1 
  else if ( CCTK_EQUALS(evolution_method, "Scotch_MHD_Bvec") ) then
     whisky_mhd_handle = 2
  else if ( CCTK_EQUALS(evolution_method, "Scotch_MHD_Avec") ) then
     if ( CCTK_EQUALS(Avec_gauge, "algebraic" ) ) then
        whisky_mhd_handle = 3
     else
        whisky_mhd_handle = 4
     end if
  end if
  
  Whisky_Startup = 0

end function Whisky_Startup
