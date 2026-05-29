 /*@@
   @file      EOS_Hybrid_Analysis.F
   @date      Fri Apr 26 16:14:47 2002
   @author    Harry Dimmelmeier, Ian Hawke, Christian Ott
   @desc 
      Calculates the polytropic and thermal contributions to the 
      total pressure.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    Check_Poly_Contrib
   @date       Fri Apr 26 16:25:35 2002
   @author     Harry Dimmelmeier, Ian Hawke, Christian Ott
   @desc 
      The routine that calculates the contributions.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

      subroutine Check_Poly_Contrib(CCTK_ARGUMENTS)

      USE EOS_Polytrope_Scalars
      USE EOS_Hybrid_Scalars

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      CCTK_INT :: i,j,k


      CCTK_REAL local_eos_gamma, local_eos_k_cgs

      do k = 1, cctk_lsh(3)
        do j = 1, cctk_lsh(2)
          do i = 1, cctk_lsh(1)

            if (rho(i,j,k) > rho_nuc) then
              local_eos_gamma = eos_gamma_supernuclear
              local_eos_k_cgs = eos_k_supernuclear_cgs
            else
              local_eos_gamma = eos_gamma
              local_eos_k_cgs = eos_k_cgs
            end if

            pressure_poly(i,j,k) = p_geom_factor * local_eos_k_cgs * 
     .           (rho(i,j,k) * rho_geom_factor_inv)**local_eos_gamma
            
            pressure_th(i,j,k) = - p_geom_factor * local_eos_k_cgs * 
     .           (eos_gamma_th - 1.d0) / (local_eos_gamma - 1.d0) * 
     .           (rho(i,j,k) * rho_geom_factor_inv)**local_eos_gamma + 
     .           (eos_gamma_th - 1.d0) * rho(i,j,k) * eps(i,j,k) - 
     .           (eos_gamma_th - 1.d0) * (local_eos_gamma - eos_gamma) / 
     .           (eos_gamma - 1.d0) / (eos_gamma_supernuclear - 1.d0) * 
     .           p_geom_factor * eos_k_cgs * 
     .           rho_geom_factor_inv**eos_gamma * 
     .           rho_nuc**(eos_gamma - 1.d0) * rho(i,j,k)

c$$$            pressure_th(i,j,k) = dmax1(0.d0, pressure_th(i,j,k))

          end do
        end do
      end do

      end subroutine Check_Poly_Contrib
