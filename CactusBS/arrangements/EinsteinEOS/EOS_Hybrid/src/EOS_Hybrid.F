 /*@@
   @file      EOS_Hybrid.F
   @date      Wed Mar 20 14:56:35 2002
   @author    Ian Hawke
   @desc 
      Routines to calculate the EOS used by 
      Dimmelmeier Novak Font Ibanez Mueller PRD71 064023 (2005)
      in supernova core collapse simulations.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"

      CCTK_REAL function EOS_Hybrid_Pressure(rho, eps)

      USE EOS_Polytrope_Scalars   
      USE EOS_Hybrid_Scalars

      implicit none

      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rho, eps, local_eos_gamma, p_poly, p_th, local_eos_k_cgs, zero

      zero = 0.d0

      if (rho > rho_nuc) then
        local_eos_gamma = eos_gamma_supernuclear
        local_eos_k_cgs = eos_k_supernuclear_cgs
      else
        local_eos_gamma = eos_gamma
        local_eos_k_cgs = eos_k_cgs
      end if

      p_poly = p_geom_factor * local_eos_k_cgs * 
     .     (rho * rho_geom_factor_inv)**local_eos_gamma

      p_th = - p_geom_factor * local_eos_k_cgs * (eos_gamma_th - 1.d0) / 
     .     (local_eos_gamma - 1.d0) * (rho * rho_geom_factor_inv)**local_eos_gamma + 
     .     (eos_gamma_th - 1.d0) * rho * eps - 
     .     (eos_gamma_th - 1.d0) * (local_eos_gamma - eos_gamma) / 
     .     (eos_gamma - 1.d0) / (eos_gamma_supernuclear - 1.d0) * 
     .     p_geom_factor * eos_k_cgs * rho_geom_factor_inv**eos_gamma * 
     .     rho_nuc**(eos_gamma - 1.d0) * rho

      p_th = max(zero, p_th)

      EOS_Hybrid_Pressure = p_poly + p_th

      end function EOS_Hybrid_Pressure

c     The specific internal energy isn''t correct yet

      CCTK_REAL function EOS_Hybrid_SpecificIE(rho, press)

      USE EOS_Polytrope_Scalars   
      USE EOS_Hybrid_Scalars

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_DECLARE(CCTK_REAL, rho, )
      CCTK_DECLARE(CCTK_REAL, press, )

c$$$      if (rho > rho_nuc) then
c$$$        EOS_Hybrid_SpecificIE = eos_k_cgs / (eos_gamma_2 - 1.d0) * 
c$$$     .       (rho / rho_nuc) ** eos_gamma_2 * rho_nuc ** eos_gamma_1 / 
c$$$     .       rho + 
c$$$     .       eos_k_cgs * (eos_gamma_2 - eos_gamma_1) / 
c$$$     .       (eos_gamma_2 - 1.d0) / (eos_gamma_1 - 1.d0) * 
c$$$     .       rho_nuc ** (eos_gamma_1 - 1.d0)
c$$$      else
c$$$        EOS_Hybrid_SpecificIE = eos_k_cgs * (rho ** (eos_gamma_1 - 1.d0)) /
c$$$     .       (eos_gamma_1 - 1.d0)
c$$$      end if

      EOS_Hybrid_SpecificIE = 2.34567890d0

      end function EOS_Hybrid_SpecificIE

c     The rest mass density isn''t correct yet

      CCTK_REAL function EOS_Hybrid_RestMassDens(eps, press)

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_DECLARE(CCTK_REAL, eps, )
      CCTK_DECLARE(CCTK_REAL, press, )

      EOS_Hybrid_RestMassDens = 1.23456789d0
      
      end function EOS_Hybrid_RestMassDens
      
      CCTK_REAL function EOS_Hybrid_DPressByDRho(rho, eps)

      USE EOS_Polytrope_Scalars   
      USE EOS_Hybrid_Scalars   

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rho, eps, local_eos_gamma, local_eos_k_cgs, d_p_poly, d_p_th_1, 
     .     d_p_th_2, zero

      zero = 0.d0

      if (rho > rho_nuc) then
        local_eos_gamma = eos_gamma_supernuclear
        local_eos_k_cgs = eos_k_supernuclear_cgs
      else
        local_eos_gamma = eos_gamma
        local_eos_k_cgs = eos_k_cgs
      end if

      d_p_poly = local_eos_gamma * p_geom_factor * local_eos_k_cgs * 
     .     rho**(local_eos_gamma - 1.d0) * rho_geom_factor_inv**local_eos_gamma

      d_p_th_1 = - local_eos_gamma * p_geom_factor * local_eos_k_cgs * 
     .     (eos_gamma_th - 1.d0) / (local_eos_gamma - 1.d0) * 
     .     rho**(local_eos_gamma - 1.d0) * rho_geom_factor_inv**local_eos_gamma

      d_p_th_2 = (eos_gamma_th - 1.d0) * eps
     .     - (eos_gamma_th - 1.d0) * (local_eos_gamma - eos_gamma) / 
     .     (eos_gamma - 1.d0) / (eos_gamma_supernuclear - 1.d0) * 
     .     p_geom_factor * eos_k_cgs * rho_geom_factor_inv**eos_gamma * 
     .     rho_nuc**(eos_gamma - 1.d0)

!      d_p_th_1 = max(d_p_th_1, zero)
!      d_p_th_2 = max(d_p_th_2, zero)
      
      EOS_Hybrid_DPressByDRho = d_p_poly + d_p_th_1 + d_p_th_2

      end function EOS_Hybrid_DPressByDRho

      CCTK_REAL function EOS_Hybrid_DPressByDEps(rho, eps)

      implicit none
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rho
      CCTK_DECLARE(CCTK_REAL, eps, )

      EOS_Hybrid_DPressByDEps = (eos_gamma_th - 1.d0) * rho

      end function EOS_Hybrid_DPressByDEps
