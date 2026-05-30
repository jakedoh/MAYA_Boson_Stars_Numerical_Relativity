 /*@@
   @file      EOS_Polytrope.F
   @date      Wed Mar 20 14:56:35 2002
   @author    Ian Hawke
   @desc 
      Routines to calculate a polytropic EOS. This is a faked
      2d version that ALWAYS assumes you want the pressure from
      the density. As such the specific internal energy is not set. 
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    EOS_Polytrope
   @date       Mon Aug 26 15:12:39 2002
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

      CCTK_REAL function EOS_Polytrope_Pressure(rho, eps)

      USE EOS_Polytrope_Scalars

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rho
      CCTK_DECLARE(CCTK_REAL, eps, )

      EOS_Polytrope_Pressure = p_geom_factor * eos_k_cgs * 
     .     (rho * rho_geom_factor_inv) ** eos_gamma
      
      end function EOS_Polytrope_Pressure

 /*@@
   @routine    EOS_Polytrope
   @date       Mon Aug 26 15:12:39 2002
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

      CCTK_REAL function EOS_Polytrope_SpecificIE(rho, press)

      USE EOS_Polytrope_Scalars

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rho, press, EOS_Polytrope_Pressure

      EOS_Polytrope_SpecificIE = EOS_Polytrope_Pressure(rho, press) 
     .     / ( rho * (eos_gamma - 1.d0) )
      
      end function EOS_Polytrope_SpecificIE

 /*@@
   @routine    EOS_Polytrope
   @date       Mon Aug 26 15:12:39 2002
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

      CCTK_REAL function EOS_Polytrope_RestMassDens(eps, press)

      USE EOS_Polytrope_Scalars

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL eps, press

      EOS_Polytrope_RestMassDens = press / ((eos_gamma - 1.d0) * eps)
      
      end function EOS_Polytrope_RestMassDens

 /*@@
   @routine    EOS_Polytrope
   @date       Mon Aug 26 15:12:39 2002
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/
      
      CCTK_REAL function EOS_Polytrope_DPressByDRho(rho, eps)

      USE EOS_Polytrope_Scalars

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rho
      CCTK_DECLARE(CCTK_REAL, eps, )

      EOS_Polytrope_DPressByDRho = p_geom_factor * eos_k_cgs * 
     .     eos_gamma * rho_geom_factor_inv *  
     .     (rho*rho_geom_factor_inv) ** (eos_gamma - 1.d0)

      end function EOS_Polytrope_DPressByDRho

 /*@@
   @routine    EOS_Polytrope
   @date       Mon Aug 26 15:12:39 2002
   @author     Ian Hawke
   @desc 
   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

      CCTK_REAL function EOS_Polytrope_DPressByDEps(rho, eps)

      USE EOS_Polytrope_Scalars

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_DECLARE(CCTK_REAL, rho, )
      CCTK_DECLARE(CCTK_REAL, eps, )

      EOS_Polytrope_DPressByDEps = 0.d0

      end function EOS_Polytrope_DPressByDEps

