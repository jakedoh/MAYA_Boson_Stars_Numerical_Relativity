 /*@@
   @file      EOS_Polytrope_Scalars.F
   @date      Mon Aug 26 12:03:33 2002
   @author    Ian Hawke
   @desc 
   Constants for EOS routines
   @enddesc 
 @@*/

#include "cctk.h"

      module EOS_Polytrope_Scalars

      implicit none

      CCTK_REAL :: p_geom_factor, rho_geom_factor, 
     .     rho_geom_factor_inv, eos_k_cgs
      
      CCTK_REAL ::  m_solar_geom, m_solar_cgs, c_cgs, G_cgs
   
      end module EOS_Polytrope_Scalars
