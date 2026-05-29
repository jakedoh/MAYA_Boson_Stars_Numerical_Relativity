/*
@file      EOS_RealisticPP_Scalars.F90
@date      Monday October 16 2017
@author    Christopher W. Evans
@contact   cevans216@gatech.edu
*/

#include "cctk.h"

module EOS_RealisticPP_Scalars
  implicit none

  CCTK_INT :: n_pieces = 7
  CCTK_REAL, DIMENSION(1:7) :: piecewise_gamma, piecewise_k, piecewise_eps, piecewise_rho
  CCTK_REAL :: press_unit = 1.79358D-39
  CCTK_REAL :: rho_unit = 1.62005D-18

end module EOS_RealisticPP_Scalars
