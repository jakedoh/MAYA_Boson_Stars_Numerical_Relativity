 /*@@
   @file      ReconstructPPM_PreLoop.F90
   @date      Mon Feb 25 11:43:36 2002
   @author    
   @desc 
   Sets up various scalars used for efficiency reasons.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

 /*@@
   @routine    PPM_Scalar_Setup
   @date       Mon Feb 25 11:25:27 2002
   @author     
   @desc 
   Sets up the logical scalars from the parameters. These are
   solely used for efficiency.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

subroutine PPM_Scalar_Setup(CCTK_ARGUMENTS)

  USE PPM_Scalars

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS



  PPM3=.false.
  PPM4=.false.

  if (CCTK_EQUALS(ppm_flatten,"stencil_3")) then
    PPM3=.true.
  else if (CCTK_EQUALS(ppm_flatten,"stencil_4")) then
    PPM4=.true.
  else
    call CCTK_WARN(0, "PPM Flattening Procedure not recognized!")
  end if

  FAST = .false.

  if (CCTK_EQUALS(numerical_viscosity,"fast")) then
    FAST = .true.
  else if (CCTK_EQUALS(numerical_viscosity,"normal")) then
    FAST = .false.
  else
    call CCTK_WARN(0, "Numerical Viscosity Mode not recognized!")
  end if

  ANALYTICAL = .false.

  if (CCTK_EQUALS(left_eigenvectors,"analytical")) then
    ANALYTICAL = .true.
  else if (CCTK_EQUALS(left_eigenvectors,"numerical")) then
    ANALYTICAL = .false.
  else 
    call CCTK_WARN(0, "Left Eigenvector Mode not recognized!")
  end if

end subroutine PPM_Scalar_Setup
