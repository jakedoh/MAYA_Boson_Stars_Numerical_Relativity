 /*@@
   @file      EOS_Polytrope_Startup.F
   @date      Wed Mar 20 14:50:22 2002
   @author    Ian Hawke
   @desc 
   Startup for EOS_Polytrope. A 2d routine so that all the standard 
   functionality is retained. Based on EOS_Ideal_Fluid.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"

 /*@@
   @routine    EOS_Polytrope_Startup
   @date       Wed Mar 20 14:52:41 2002
   @author     Ian Hawke
   @desc 
   Startup for EOS_Polytrope
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
@@*/

      integer function EOS_Polytrope_Startup()

      USE EOS_Polytrope_Scalars

      implicit none

      DECLARE_CCTK_PARAMETERS

      integer handle, ierr
      
      external EOS_Polytrope_Pressure
      external EOS_Polytrope_SpecificIE
      external EOS_Polytrope_RestMassDens
      external EOS_Polytrope_DPressByDRho
      external EOS_Polytrope_DPressByDEps

      call EOS_RegisterMethod(handle,"2D_Polytrope")

      if (handle .ge. 0) then

        call EOS_RegisterPressure(ierr, handle, EOS_Polytrope_Pressure)
        call EOS_RegisterSpecificIntEnergy(ierr, handle, 
     .                                     EOS_Polytrope_SpecificIE)
        call EOS_RegisterRestMassDens(ierr, handle, 
     .                                EOS_Polytrope_RestMassDens)
        call EOS_RegisterDPressByDRho(ierr, handle, 
     .                                EOS_Polytrope_DPressByDRho)
        call EOS_RegisterDPressByDEps(ierr, handle, 
     .                                EOS_Polytrope_DPressByDEps)

      else

        call CCTK_WARN(0, "Unable to register the EOS method!")

      end if

      EOS_Polytrope_Startup = 0

c$$$      These constants match the values in RNSID: consts.h

      m_solar_cgs = 1.987d33
      c_cgs = 2.9979d10
      G_cgs = 6.6732d-8
      m_solar_geom = G_cgs / c_cgs**2 * m_solar_cgs

      if (use_cgs .ne. 0) then

        p_geom_factor = G_cgs / c_cgs**4 * m_solar_geom**2
        rho_geom_factor = p_geom_factor * c_cgs**2
        rho_geom_factor_inv = 1.d0 / rho_geom_factor 
        
        eos_k_cgs = eos_k * rho_geom_factor**gamma_ini / p_geom_factor
        
      else
        
        p_geom_factor = 1.d0
        rho_geom_factor = 1.d0
        rho_geom_factor_inv = 1.d0
        eos_k_cgs = eos_k

      end if
        
      end function EOS_Polytrope_Startup
