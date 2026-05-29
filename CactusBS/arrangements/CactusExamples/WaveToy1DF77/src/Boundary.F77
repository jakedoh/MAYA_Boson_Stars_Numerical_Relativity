 /*@@
   @file      Boundary.F77
   @date      July 15th 2000
   @author    Gabrielle Allen
   @desc 
   Apply boundary conditions for 1D Wave Equation
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToy1DF77_Boundary(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer ierr
      integer sw(1)
      CCTK_REAL zero,one
      
      sw(1) = 1
      zero = 0.0d0
      one = 1.0d0

      if (CCTK_EQUALS(bound,"zero")) then
        call BndScalarVN(ierr,cctkGH, sw, zero, "wavetoy1df77::phi")
      else if (CCTK_EQUALS(bound,"radiative")) then
	call BndRadiativeVN(ierr,cctkGH,sw,one,one,"wavetoy1df77::phi",
     &	                    "wavetoy1df77::phi")
      else if (CCTK_EQUALS(bound,"none")) then
	call CCTK_WARN(9,"Hopefully applying periodic boundary condition")
      else 
	call CCTK_INFO("No boundary condition applied")
      end if

      return
      end
