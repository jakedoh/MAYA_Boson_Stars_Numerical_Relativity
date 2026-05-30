 /*@@
   @file      Whisky_EMTensor_temps.inc
   @date      Sat Jan 26 02:03:56 2002
   @author    Ian Hawke
   @desc 
     The declaration of the temporaries used in Whisky_EMTensor.inc.
   @enddesc 
 @@*/
#include "SpaceMask.h"
#ifdef FCODE

      CCTK_REAL velxlow, velylow, velzlow
      CCTK_REAL betaxlow, betaylow, betazlow, beta2
      CCTK_REAL utlow, uxlow, uylow, uzlow
      CCTK_REAL rhoenthalpy
      CCTK_INT :: type_bits, atmosphere

#endif

#ifdef CCODE

      CCTK_REAL velxlow, velylow, velzlow;
      CCTK_REAL betaxlow, betaylow, betazlow, beta2;
      CCTK_REAL utlow, uxlow, uylow, uzlow;
      CCTK_REAL rhoenthalpy;
      CCTK_INT  type_bits, atmosphere;

#endif
