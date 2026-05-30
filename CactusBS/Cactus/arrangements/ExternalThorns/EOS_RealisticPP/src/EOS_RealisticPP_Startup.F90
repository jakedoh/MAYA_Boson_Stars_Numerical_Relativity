/*
@file      EOS_RealisticPP_Startup.F90
@date      Monday October 16 2017
@author    Christopher W. Evans
@contact   cevans216@gatech.edu
*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine EOS_RealisticPP_DetermineEOS(pone, rho_one, gamma_one, gamma_two, gamma_three)
  USE EOS_RealisticPP_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8 :: log_pone, log_rho_one
  real*8, intent(out) :: pone, rho_one, gamma_one, gamma_two, gamma_three
  integer :: j
  character(len=256) :: buffer
  character(len=32) :: eos_name_string

  if (CCTK_EQUALS(realistic_eos_name, 'PAL6')) then
    log_pone = 34.380
    log_rho_one = 1.
    gamma_one = 2.227
    gamma_two = 2.189
    gamma_three = 2.159
    eos_name_string = 'PAL6'
  else if (CCTK_EQUALS(realistic_eos_name, 'SLy')) then
    log_pone = 34.384
    log_rho_one = 14.165
    gamma_one = 3.005
    gamma_two = 2.988
    gamma_three = 2.851
    eos_name_string = 'SLy'
  else if (CCTK_EQUALS(realistic_eos_name, 'AP1')) then
    log_pone = 33.943
    log_rho_one = 1.
    gamma_one = 2.442
    gamma_two = 3.256
    gamma_three = 2.908
    eos_name_string = 'AP1'
  else if (CCTK_EQUALS(realistic_eos_name, 'AP2')) then
    log_pone = 34.126
    log_rho_one = 1.
    gamma_one = 2.643 
    gamma_two = 3.014
    gamma_three = 2.945
    eos_name_string = 'AP2'
  else if (CCTK_EQUALS(realistic_eos_name, 'AP3')) then
    log_pone = 34.392
    log_rho_one = 1.
    gamma_one = 3.166
    gamma_two = 3.573
    gamma_three = 3.281
    eos_name_string = 'AP3'
  else if (CCTK_EQUALS(realistic_eos_name, 'AP4')) then
    log_pone = 34.269
    log_rho_one = 1.
    gamma_one = 2.830
    gamma_two = 3.445
    gamma_three = 3.348
    eos_name_string = 'AP4'
  else if (CCTK_EQUALS(realistic_eos_name, 'FPS')) then
    log_pone = 34.283
    log_rho_one = 1.
    gamma_one = 2.985
    gamma_two = 2.863
    gamma_three = 2.600
    eos_name_string = 'FPS'
  else if (CCTK_EQUALS(realistic_eos_name, 'WFF1')) then
    log_pone = 34.031
    log_rho_one = 1.
    gamma_one = 2.519
    gamma_two = 3.791
    gamma_three = 3.660
    eos_name_string = 'WFF1'
  else if (CCTK_EQUALS(realistic_eos_name, 'WFF2')) then
    log_pone = 34.233
    log_rho_one = 1.
    gamma_one = 2.888
    gamma_two = 3.475
    gamma_three = 3.517
    eos_name_string = 'WFF2'
  else if (CCTK_EQUALS(realistic_eos_name, 'WFF3')) then
    log_pone = 34.283
    log_rho_one = 1.
    gamma_one = 3.329
    gamma_two = 2.952
    gamma_three = 2.589
    eos_name_string = 'WFF3'
  else if (CCTK_EQUALS(realistic_eos_name, 'BBB2')) then
    log_pone = 34.331
    log_rho_one = 1.
    gamma_one = 3.418
    gamma_two = 2.835
    gamma_three = 2.832
    eos_name_string = 'BBB2'
  else if (CCTK_EQUALS(realistic_eos_name, 'BPAL12')) then
    log_pone = 34.358
    log_rho_one = 1.
    gamma_one = 2.209
    gamma_two = 2.201
    gamma_three = 2.176
    eos_name_string = 'BPAL12'
  else if (CCTK_EQUALS(realistic_eos_name, 'ENG')) then
    log_pone = 34.437
    log_rho_one = 1.
    gamma_one = 3.514
    gamma_two = 3.130
    gamma_three = 3.168
    eos_name_string = 'ENG'
  else if (CCTK_EQUALS(realistic_eos_name, 'MPA1')) then
    log_pone = 34.495
    log_rho_one = 1.
    gamma_one = 3.446
    gamma_two = 3.572
    gamma_three = 2.887
    eos_name_string = 'MPA1'
  else if (CCTK_EQUALS(realistic_eos_name, 'MS1')) then
    log_pone = 34.858
    log_rho_one = 1.
    gamma_one = 3.224
    gamma_two = 3.033
    gamma_three = 1.325
    eos_name_string = 'MS1'
  else if (CCTK_EQUALS(realistic_eos_name, 'MS2')) then
    log_pone = 34.605
    log_rho_one = 1.
    gamma_one = 2.477
    gamma_two = 2.184
    gamma_three = 1.855
    eos_name_string = 'MS2'
  else if (CCTK_EQUALS(realistic_eos_name, 'MS1b')) then
    log_pone = 34.855
    log_rho_one = 14.05556938
    gamma_one = 3.456
    gamma_two = 3.011
    gamma_three = 1.425
    eos_name_string = 'MS1b'
  else if (CCTK_EQUALS(realistic_eos_name, 'PS')) then
    log_pone = 34.671
    log_rho_one = 1.
    gamma_one = 2.216
    gamma_two = 1.640
    gamma_three = 2.365
    eos_name_string = 'PS'
  else if (CCTK_EQUALS(realistic_eos_name, 'GS1')) then
    log_pone = 34.504
    log_rho_one = 1.
    gamma_one = 2.350
    gamma_two = 1.267
    gamma_three = 2.421
    eos_name_string = 'GS1'
  else if (CCTK_EQUALS(realistic_eos_name, 'GS2')) then
    log_pone = 34.642
    log_rho_one = 1.
    gamma_one = 2.519
    gamma_two = 1.571
    gamma_three = 2.314
    eos_name_string = 'GS2'
  else if (CCTK_EQUALS(realistic_eos_name, 'BGN1H1')) then
    log_pone = 34.623
    log_rho_one = 1.
    gamma_one = 3.258
    gamma_two = 1.472
    gamma_three = 2.464
    eos_name_string = 'BGN1H1'
  else if (CCTK_EQUALS(realistic_eos_name, 'GNH3')) then
    log_pone = 34.648
    log_rho_one = 1.
    gamma_one = 2.664
    gamma_two = 2.194
    gamma_three = 2.304
    eos_name_string = 'GNH3'
  else if (CCTK_EQUALS(realistic_eos_name, 'H1')) then
    log_pone = 34.564
    log_rho_one = 1.
    gamma_one = 2.595
    gamma_two = 1.845
    gamma_three = 1.897
    eos_name_string = 'H1'
  else if (CCTK_EQUALS(realistic_eos_name, 'H2')) then
    log_pone = 34.617
    log_rho_one = 1.
    gamma_one = 2.775
    gamma_two = 1.855
    gamma_three = 1.858
    eos_name_string = 'H2'
  else if (CCTK_EQUALS(realistic_eos_name, 'H3')) then
    log_pone = 34.646
    log_rho_one = 1.
    gamma_one = 2.787
    gamma_two = 1.951
    gamma_three = 1.901
    eos_name_string = 'H3'
  else if (CCTK_EQUALS(realistic_eos_name, 'H4')) then
    log_pone = 34.669
    log_rho_one = 1.
    gamma_one = 2.909
    gamma_two = 2.246
    gamma_three = 2.144
    eos_name_string = 'H4'
  else if (CCTK_EQUALS(realistic_eos_name, 'H5')) then
    log_pone = 34.609
    log_rho_one = 1.
    gamma_one = 2.793
    gamma_two = 1.974
    gamma_three = 1.915
    eos_name_string = 'H5'
  else if (CCTK_EQUALS(realistic_eos_name, 'H6')) then
    log_pone = 34.593
    log_rho_one = 1.
    gamma_one = 2.637
    gamma_two = 2.121
    gamma_three = 2.064
    eos_name_string = 'H6'
  else if (CCTK_EQUALS(realistic_eos_name, 'H7')) then
    log_pone = 34.559
    log_rho_one = 1.
    gamma_one = 2.621
    gamma_two = 2.048
    gamma_three = 2.006
    eos_name_string = 'H7'
  else if (CCTK_EQUALS(realistic_eos_name, 'PCL2')) then
    log_pone = 34.507
    log_rho_one = 1.
    gamma_one = 2.554
    gamma_two = 1.880
    gamma_three = 1.977
    eos_name_string = 'PCL2'
  else if (CCTK_EQUALS(realistic_eos_name, 'ALF1')) then
    log_pone = 34.055
    log_rho_one = 1.
    gamma_one = 2.013
    gamma_two = 3.389
    gamma_three = 2.033
    eos_name_string = 'ALF1'
  else if (CCTK_EQUALS(realistic_eos_name, 'ALF2')) then
    log_pone = 34.616
    log_rho_one = 1.
    gamma_one = 4.070
    gamma_two = 2.411
    gamma_three = 1.890
    eos_name_string = 'ALF2'
  else if (CCTK_EQUALS(realistic_eos_name, 'ALF3')) then
    log_pone = 34.283
    log_rho_one = 1.
    gamma_one = 2.883
    gamma_two = 2.653
    gamma_three = 1.952
    eos_name_string = 'ALF3'
  else if (CCTK_EQUALS(realistic_eos_name, 'ALF4')) then
    log_pone = 34.314
    log_rho_one = 1.
    gamma_one = 3.009
    gamma_two = 3.438
    gamma_three = 1.803
    eos_name_string = 'ALF4'
  end if

  write(buffer,*) 'EOS_RealisticPP is using the ',TRIM(eos_name_string),' equation of state.'
  call CCTK_INFO(buffer)

  pone = (10.**log_pone)*(press_unit*code_unit*code_unit)
  rho_one = (10.**log_rho_one)*(rho_unit*code_unit*code_unit)

end subroutine EOS_RealisticPP_DetermineEOS

subroutine EOS_RealisticPP_SetLowDensityEOS()
  USE EOS_RealisticPP_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: j

  piecewise_gamma(1) = 1.58425
  piecewise_gamma(2) = 1.28733
  piecewise_gamma(3) = 0.62223
  piecewise_gamma(4) = 1.35692

  piecewise_rho(1) = (1e-11)*code_unit*code_unit
  piecewise_rho(2) = (3.951156e-11)*code_unit*code_unit
  piecewise_rho(3) = (6.125960e-7)*code_unit*code_unit
  piecewise_rho(4) = (4.254672e-6)*code_unit*code_unit

  piecewise_k(1) = (k0)*code_unit**(2.*(1. - piecewise_gamma(1)))
  piecewise_eps(1) = 0.

  ! Also set rho separation values for EOS as they are also the same
  piecewise_rho(6) = (8.1114721e-4)*code_unit*code_unit
  piecewise_rho(7) = (1.619100e-3)*code_unit*code_unit

  ! Determine piecewise_k and piecewise_eps for these values
  do j=2,4
    piecewise_k(j) = piecewise_k(j-1)*piecewise_rho(j)**(piecewise_gamma(j-1) - piecewise_gamma(j))
    piecewise_eps(j) = piecewise_eps(j-1) + piecewise_k(j-1)*(piecewise_gamma(j) - piecewise_gamma(j-1))*(piecewise_rho(j)**(piecewise_gamma(j-1) - 1.))/(piecewise_gamma(j-1) - 1.)/(piecewise_gamma(j) - 1.)
  end do
end subroutine EOS_RealisticPP_SetLowDensityEOS

subroutine EOS_RealisticPP_SetFourPieceEOS(rho_one, gamma_one, gamma_two, gamma_three)
  USE EOS_RealisticPP_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8, intent(in) :: rho_one, gamma_one, gamma_two, gamma_three
  integer :: j
  character(len=256) :: buffer

  n_pieces = 4

  piecewise_rho(1) = (1e-11)*code_unit*code_unit
  piecewise_rho(2) = rho_one
  piecewise_rho(3) = (8.1114721e-4)*code_unit*code_unit
  piecewise_rho(4) = (1.619100e-3)*code_unit*code_unit

  piecewise_gamma(1) = gamma_zero
  piecewise_gamma(2) = gamma_one
  piecewise_gamma(3) = gamma_two
  piecewise_gamma(4) = gamma_three

  piecewise_k(1) = (k0)*code_unit**(2.*(1. - piecewise_gamma(1)))
  piecewise_eps(1) = 0.

  do j=2,n_pieces
    piecewise_k(j) = piecewise_k(j-1)*piecewise_rho(j)**(piecewise_gamma(j-1) - piecewise_gamma(j))
    piecewise_eps(j) = piecewise_eps(j-1) + piecewise_k(j-1)*(piecewise_gamma(j) - piecewise_gamma(j-1))*(piecewise_rho(j)**(piecewise_gamma(j-1) - 1.))/(piecewise_gamma(j-1) - 1.)/(piecewise_gamma(j) - 1.)
  end do

  do j=1,n_pieces
    write(buffer, *) 'piecewise_rho(',j,') = ',piecewise_rho(j)
    call CCTK_INFO(buffer)
  end do
  do j=1,n_pieces
    write(buffer, *) 'piecewise_gamma(',j,') = ',piecewise_gamma(j)
    call CCTK_INFO(buffer)
  end do
  do j=1,n_pieces
    write(buffer, *) 'piecewise_eps(',j,') = ',piecewise_eps(j)
    call CCTK_INFO(buffer)
  end do
  do j=1,n_pieces
    write(buffer, *) 'piecewise_k(',j,') = ',piecewise_k(j)
    call CCTK_INFO(buffer)
  end do
end subroutine EOS_RealisticPP_SetFourPieceEOS

subroutine EOS_RealisticPP_SetSevenPieceEOS(pone, gamma_one, gamma_two, gamma_three)
  USE EOS_RealisticPP_Scalars

  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  real*8, intent(in) :: pone, gamma_one, gamma_two, gamma_three
  integer :: j

  piecewise_gamma(5) = gamma_one
  piecewise_gamma(6) = gamma_two
  piecewise_gamma(7) = gamma_three

  ! Determine piecewise_rho(5) from EOS parameters
  piecewise_rho(5) = (pone/(piecewise_k(4)*piecewise_rho(6)**gamma_one))**(1./(piecewise_gamma(4) - gamma_one))

  ! Set the remainder of the piecewise_k and piecewise_eps values
  do j=5,n_pieces
    piecewise_k(j) = piecewise_k(j-1)*piecewise_rho(j)**(piecewise_gamma(j-1) - piecewise_gamma(j))
    piecewise_eps(j) = piecewise_eps(j-1) + piecewise_k(j-1)*(piecewise_gamma(j) - piecewise_gamma(j-1))*(piecewise_rho(j)**(piecewise_gamma(j-1) - 1.))/(piecewise_gamma(j-1) - 1.)/(piecewise_gamma(j) - 1.)
  end do
end subroutine EOS_RealisticPP_SetSevenPieceEOS

integer function EOS_RealisticPP_Startup()
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: handle, ierr
  real*8 :: pone, rho_one, gamma_one, gamma_two, gamma_three

  external EOS_RealisticPP_Pressure
  external EOS_RealisticPP_SpecificIE
  external EOS_RealisticPP_RestMassDens
  external EOS_RealisticPP_DPressByDRho
  external EOS_RealisticPP_DPressByDEps

  external EOS_RealisticPP_Cold_Pressure
  external EOS_RealisticPP_Cold_SpecificIE
  external EOS_RealisticPP_Cold_RestMassDens
  external EOS_RealisticPP_Cold_DPressByDRho
  external EOS_RealisticPP_Cold_DPressByDEps

  call EOS_RegisterMethod(handle, 'RealisticPP')

  if (handle .ge. 0) then
    if (gamma_thermal .eq. 1.0) then
      call EOS_RegisterPressure(ierr, handle, EOS_RealisticPP_Cold_Pressure)
      call EOS_RegisterSpecificIntEnergy(ierr, handle, EOS_RealisticPP_Cold_SpecificIE)
      call EOS_RegisterRestMassDens(ierr, handle, EOS_RealisticPP_Cold_RestMassDens)
      call EOS_RegisterDPressByDRho(ierr, handle, EOS_RealisticPP_Cold_DPressByDRho)
      call EOS_RegisterDPressByDEps(ierr, handle, EOS_RealisticPP_Cold_DPressByDEps)
    else
      call EOS_RegisterPressure(ierr, handle, EOS_RealisticPP_Pressure)
      call EOS_RegisterSpecificIntEnergy(ierr, handle, EOS_RealisticPP_SpecificIE)
      call EOS_RegisterRestMassDens(ierr, handle, EOS_RealisticPP_RestMassDens)
      call EOS_RegisterDPressByDRho(ierr, handle, EOS_RealisticPP_DPressByDRho)
      call EOS_RegisterDPressByDEps(ierr, handle, EOS_RealisticPP_DPressByDEps)
    end if
  else
    call CCTK_WARN(0, 'Unable to register the EOS method!')
  end if

  ! If not using a four piece approximant, set values for low-density EOS that closely matches SLy. 
  ! Same for all equations of state.
  if (four_piece .eq. 0) then
    call CCTK_INFO("Initializing low-density equation of state (using SLy).")
    call EOS_RealisticPP_SetLowDensityEOS()
  end if

  ! Set EOS specific values
  call EOS_RealisticPP_DetermineEOS(pone, rho_one, gamma_one, gamma_two, gamma_three)

  if (four_piece .ne. 0) then
    call EOS_RealisticPP_SetFourPieceEOS(rho_one, gamma_one, gamma_two, gamma_three)
    call CCTK_INFO("Initialized four piece equation of state.")
  else
    call EOS_RealisticPP_SetSevenPieceEOS(pone, gamma_one, gamma_two, gamma_three)
    call CCTK_INFO("Initialized full seven piece equation of state.")
  end if

  EOS_RealisticPP_Startup = 0
end function EOS_RealisticPP_Startup


