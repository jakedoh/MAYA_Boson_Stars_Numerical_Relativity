/*
@file      EOS_RealisticPP.F90
@date      Monday October 16 2017
@author    Christopher W. Evans
@contact   cevans216@gatech.edu
*/

#include "cctk.h"
#include "cctk_Parameters.h"
subroutine EOS_RealisticPP_LocalEOSParameters(rho, local_gamma, local_k, local_eps)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho
  real*8, intent(out) :: local_gamma, local_k, local_eps
  integer :: j

  local_gamma = piecewise_gamma(1)
  local_k = piecewise_k(1)
  local_eps = piecewise_eps(1)

  do j = 2, n_pieces
    if (rho .gt. piecewise_rho(j)) then
      local_gamma = piecewise_gamma(j)
      local_k = piecewise_k(j)
      local_eps = piecewise_eps(j)
    end if
  end do
end subroutine EOS_RealisticPP_LocalEOSParameters

CCTK_REAL function EOS_RealisticPP_Pressure(rho, eps)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, eps
  CCTK_REAL :: polytropic_pressure, thermal_pressure
  real*8 :: local_gamma, local_k, local_eps

  call EOS_RealisticPP_LocalEOSParameters(rho, local_gamma, local_k, local_eps)

  polytropic_pressure = local_k*rho**local_gamma
  thermal_pressure = rho*(gamma_thermal - 1.)*(eps - local_eps) - ((gamma_thermal - 1.)/(local_gamma - 1.))*polytropic_pressure

  thermal_pressure = max(0., thermal_pressure)

  EOS_RealisticPP_Pressure = polytropic_pressure + thermal_pressure
end function EOS_RealisticPP_Pressure

CCTK_REAL function EOS_RealisticPP_Cold_Pressure(rho, eps)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, eps
  real*8 :: local_gamma, local_k, local_eps

  call EOS_RealisticPP_LocalEOSParameters(rho, local_gamma, local_k, local_eps)

  EOS_RealisticPP_Cold_Pressure = local_k*rho**local_gamma
end function EOS_RealisticPP_Cold_Pressure

CCTK_REAL function EOS_RealisticPP_SpecificIE(rho, press)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, press
  real*8 :: local_gamma, local_k, local_eps

  call EOS_RealisticPP_LocalEOSParameters(rho, local_gamma, local_k, local_eps)

  EOS_RealisticPP_SpecificIE = local_eps + press/rho/(gamma_thermal - 1.) + local_k*(rho**(local_gamma - 1.))*(gamma_thermal - local_gamma)/(local_gamma - 1.)/(gamma_thermal - 1.)
end function EOS_RealisticPP_SpecificIE

CCTK_REAL function EOS_RealisticPP_Cold_SpecificIE(rho, press)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, press
  real*8 :: local_gamma, local_k, local_eps

  call EOS_RealisticPP_LocalEOSParameters(rho, local_gamma, local_k, local_eps)

  EOS_RealisticPP_Cold_SpecificIE = local_eps + local_k*rho**(local_gamma - 1.) /(local_gamma - 1.)
end function EOS_RealisticPP_Cold_SpecificIE


CCTK_REAL function EOS_RealisticPP_RestMassDens(eps, press)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: eps, press
  character(256) :: errorstring

  write(errorstring,*) "EOS_RealisticPP_RestMassDens() is not implemented!"
  call CCTK_ERROR(errorstring)

  EOS_RealisticPP_RestMassDens = 0.
end function EOS_RealisticPP_RestMassDens

CCTK_REAL function EOS_RealisticPP_Cold_RestMassDens(eps, press)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: eps, press
  character(256) :: errorstring

  write(errorstring,*) "EOS_RealisticPP_RestMassDens() is not implemented!"
  call CCTK_ERROR(errorstring)

  EOS_RealisticPP_Cold_RestMassDens = 0.
end function EOS_RealisticPP_Cold_RestMassDens

CCTK_REAL function EOS_RealisticPP_DPressByDRho(rho, eps)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, eps
  real*8 :: local_gamma, local_k, local_eps, polytropic_pressure, thermal_pressure

  call EOS_RealisticPP_LocalEOSParameters(rho, local_gamma, local_k, local_eps)

  polytropic_pressure = local_k*rho**local_gamma
  thermal_pressure = rho*(gamma_thermal - 1.)*(eps - local_eps) - ((gamma_thermal - 1.)/(local_gamma - 1.))*polytropic_pressure

  thermal_pressure = max(0., thermal_pressure)

  EOS_RealisticPP_DPressByDRho = ((local_gamma - gamma_thermal + 1.)*polytropic_pressure + thermal_pressure)/rho
end function EOS_RealisticPP_DPressByDRho

CCTK_REAL function EOS_RealisticPP_Cold_DPressByDRho(rho, eps)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, eps
  real*8 :: local_gamma, local_k, local_eps

  call EOS_RealisticPP_LocalEOSParameters(rho, local_gamma, local_k, local_eps)

  EOS_RealisticPP_Cold_DPressByDRho = local_gamma*local_k*rho**(local_gamma - 1.)
end function EOS_RealisticPP_Cold_DPressByDRho

CCTK_REAL function EOS_RealisticPP_DPressByDEps(rho, eps)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, eps

  EOS_RealisticPP_DPressByDEps = (gamma_thermal - 1.)*rho
end function EOS_RealisticPP_DPressByDEps

CCTK_REAL function EOS_RealisticPP_Cold_DPressByDEps(rho, eps)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL, intent(in) :: rho, eps

  EOS_RealisticPP_Cold_DPressByDEps = 0.
end function EOS_RealisticPP_Cold_DPressByDEps

CCTK_REAL function EOS_RealisticPP_HybridColdEps(rho)
  USE EOS_RealisticPP_Scalars

  implicit none
  DECLARE_CCTK_PARAMETERS


  CCTK_REAL, intent(in) :: rho

  CCTK_REAL, external :: EOS_RealisticPP_Cold_SpecificIE
  CCTK_REAL ZERO
  EOS_RealisticPP_HybridColdEps =  EOS_RealisticPP_Cold_SpecificIE(rho, 0.)
  
end function EOS_RealisticPP_HybridColdEps
