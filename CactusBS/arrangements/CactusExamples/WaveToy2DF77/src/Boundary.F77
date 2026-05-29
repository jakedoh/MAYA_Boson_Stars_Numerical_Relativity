#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

      subroutine WaveToy2DF77_Boundary(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS

      integer ierr
      integer sw(2)
      CCTK_REAL zero
      
      sw(1) = 1
      sw(2) = 1

      zero = 0.0d0 

      if (CCTK_EQUALS(bound,"zero")) then
        call BndScalarVN(ierr,cctkGH, sw, zero,"wavetoy2df77::phi")
      end if

      return
      end
