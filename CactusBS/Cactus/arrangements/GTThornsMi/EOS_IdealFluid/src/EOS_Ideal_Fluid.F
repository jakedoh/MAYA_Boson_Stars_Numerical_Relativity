c/*@@
c  @file      EOS_Ideal_Fluid.F
c  @date      December 1999
c  @author    Mark Miller
c  @desc 
c     Routines to calculate Ideal Fluid EOS through EOS_Base
c  @enddesc 
c@@*/

#include "cctk.h"
#include "cctk_Parameters.h"

c/*@@
c  @routine    EOS_Ideal_Fluid_Pressure
c  @date       December 1999
c  @author     Mark Miller
c  @desc 
c    Calculate the pressure  P=(gamma-1)*eps*rho
c  @enddesc 
c  @calls     
c  @calledby   
c  @history 
c
c  @endhistory 
c
c@@*/
      CCTK_REAL FUNCTION EOS_Ideal_Fluid_Pressure(rest_mass_density, 
     .                                            specific_internal_energy)

      IMPLICIT NONE
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rest_mass_density
      CCTK_REAL specific_internal_energy

      EOS_Ideal_Fluid_Pressure = (eos_ideal_fluid_gamma - 1.0d0) *
     .   rest_mass_density * specific_internal_energy

      END FUNCTION EOS_Ideal_Fluid_Pressure

c/*@@
c  @routine    EOS_Ideal_Fluid_SpecificIE
c  @date       December 1999
c  @author     Mark Miller
c  @desc 
c    Calculate the specific internal energy  
c       eps = press / ( (gamma-1)*rho)
c  @enddesc 
c  @calls     
c  @calledby   
c  @history 
c
c  @endhistory 
c
c@@*/
      CCTK_REAL FUNCTION EOS_Ideal_Fluid_SpecificIE(rest_mass_density,pressure)

      IMPLICIT NONE
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rest_mass_density
      CCTK_REAL pressure

      EOS_Ideal_Fluid_SpecificIE = pressure / ((eos_ideal_fluid_gamma - 1.0d0)*
     .   rest_mass_density)

      END FUNCTION EOS_Ideal_Fluid_SpecificIE

c/*@@
c  @routine    EOS_Ideal_Fluid_RestMassDens
c  @date       December 1999
c  @author     Mark Miller
c  @desc 
c    Calculate the rest mass density
c       rho = pressure / ((gamma-1)*eps)
c  @enddesc 
c  @calls     
c  @calledby   
c  @history 
c
c  @endhistory 
c
c@@*/
      CCTK_REAL FUNCTION EOS_Ideal_Fluid_RestMassDens(specific_internal_energy,
     .    pressure)

      IMPLICIT NONE
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL specific_internal_energy
      CCTK_REAL pressure

      EOS_Ideal_Fluid_RestMassDens = pressure / 
     .   ((eos_ideal_fluid_gamma - 1.0d0)*specific_internal_energy)

      END FUNCTION EOS_Ideal_Fluid_RestMassDens

c/*@@
c  @routine    EOS_Ideal_Fluid_DPressByDRho
c  @date       December 1999
c  @author     Mark Miller
c  @desc 
c    Calculate d(pressure)/d(rho), keeping eps fixed:
c       dp/drho = (gamma-1)*eps
c  @enddesc 
c  @calls     
c  @calledby   
c  @history 
c
c  @endhistory 
c
c@@*/
      CCTK_REAL FUNCTION EOS_Ideal_Fluid_DPressByDRho(rest_mass_density, 
     .      specific_internal_energy)

      IMPLICIT NONE
      DECLARE_CCTK_PARAMETERS

c     use CCTK_DECLARE to silence compiler warnings
      CCTK_DECLARE(CCTK_REAL, rest_mass_density, )
      CCTK_REAL specific_internal_energy

      EOS_Ideal_Fluid_DPressByDRho = (eos_ideal_fluid_gamma - 1.0d0) *
     .   specific_internal_energy

      END FUNCTION EOS_Ideal_Fluid_DPressByDRho

c/*@@
c  @routine    EOS_Ideal_Fluid_DPressByDEps
c  @date       December 1999
c  @author     Mark Miller
c  @desc 
c    Calculate d(pressure)/d(eps), keeping rho fixed
c      dp/deps = (gamma-1)*rho
c  @enddesc 
c  @calls     
c  @calledby   
c  @history 
c
c  @endhistory 
c
c@@*/
      CCTK_REAL FUNCTION EOS_Ideal_Fluid_DPressByDEps(rest_mass_density, 
     .   specific_internal_energy)

      IMPLICIT NONE
      DECLARE_CCTK_PARAMETERS

      CCTK_REAL rest_mass_density
c     use CCTK_DECLARE to silence compiler warnings
      CCTK_DECLARE(CCTK_REAL, specific_internal_energy, )

      EOS_Ideal_Fluid_DPressByDEps = (eos_ideal_fluid_gamma - 1.0d0) *
     .   rest_mass_density

      END FUNCTION EOS_Ideal_Fluid_DPressByDEps
