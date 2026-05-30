 /*@@
   @file      EOS_Hybrid_Startup.F
   @date      Wed Mar 20 14:50:22 2002
   @author    Ian Hawke
   @desc 
   Startup for EOS_Hybrid. 
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    EOS_Hybrid_Startup
   @date       Wed Mar 20 14:52:41 2002
   @author     Ian Hawke
   @desc 
   Startup for EOS_Hybrid
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

      integer function EOS_Hybrid_Startup()

      USE EOS_Polytrope_Scalars
      USE EOS_Hybrid_Scalars

      implicit none

      DECLARE_CCTK_PARAMETERS

      integer handle, ierr
      
      external EOS_Hybrid_Pressure
      external EOS_Hybrid_SpecificIE
      external EOS_Hybrid_RestMassDens
      external EOS_Hybrid_DPressByDRho
      external EOS_Hybrid_DPressByDEps

      call EOS_RegisterMethod(handle,"EOS_Hybrid")

      if (handle .ge. 0) then

        call EOS_RegisterPressure(ierr, handle, EOS_Hybrid_Pressure)
        call EOS_RegisterSpecificIntEnergy(ierr, handle, 
     .                                     EOS_Hybrid_SpecificIE)
        call EOS_RegisterRestMassDens(ierr, handle, 
     .                                EOS_Hybrid_RestMassDens)
        call EOS_RegisterDPressByDRho(ierr, handle, 
     .                                EOS_Hybrid_DPressByDRho)
        call EOS_RegisterDPressByDEps(ierr, handle, 
     .                                EOS_Hybrid_DPressByDEps)

      else

        call CCTK_WARN(0, "Unable to register the EOS method!")

      end if

cccc      all these are now set by the polytropic EOS

cccc      These constants match the values in RNSID: consts.h

cccc      m_solar_cgs = 1.987d33
cccc      c_cgs = 2.9979d10
cccc      G_cgs = 6.6732d-8
cccc
cccc      m_solar_geom = G_cgs / c_cgs**2 * m_solar_cgs
cccc      p_geom_factor = G_cgs / c_cgs**4 * m_solar_geom**2
cccc      rho_geom_factor = p_geom_factor * c_cgs**2
cccc      rho_geom_factor_inv = 1.d0 / rho_geom_factor 
cccc
cccc      eos_k_cgs = eos_k * rho_geom_factor**gamma_ini / p_geom_factor 

      eos_k_supernuclear_cgs = eos_k_cgs * (rho_nuc * rho_geom_factor_inv)**
     .     (eos_gamma - eos_gamma_supernuclear)

      EOS_Hybrid_Startup = 0
        
      end function EOS_Hybrid_Startup
