
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include "IncludeHeader2.h"

      subroutine PrintInclude2(CCTK_ARGUMENTS)

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS

      print *
      print *
      print *,"Testing Cactus Build Include Files from Interface2..."
      print *
      print *,"Source should ONLY be included from all active thorns"
      print *
      print *,"Headers will be included from all compiled thorns"

      if (CCTK_IsThornActive("Include1") .eq. 1) then
        print *
        print *,"Thorn Include1 is active"
      else
        print *
        print *,"Thorn Include1 is not active"
      end if

      if (CCTK_IsThornActive("Include2") .eq. 1) then
        print *
        print *,"Thorn Include2 is active"
      else
        print *
        print *,"Thorn Include2 is not active"
      end if

#include "IncludeSource2.h"

#ifdef ih1in1
      print *
      print *,"INCLUDED header from Interface1"
#endif

#ifdef ih1in2
      print *
      print *,"INCLUDED header from Interface1"
#endif

#ifdef ih2in1 
      print *
      print *,"INCLUDED header from Interface2"
#endif

#ifdef ih2in2
      print *
      print *,"INCLUDED header from Interface2"
#endif

      return
      end

  
