! $Header$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine CalcK_copy_to_prev (CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  gxx_prev = gxx
  gxy_prev = gxy
  gxz_prev = gxz
  gyy_prev = gyy
  gyz_prev = gyz
  gzz_prev = gzz
end subroutine CalcK_copy_to_prev

subroutine CalcK_copy_to_prev2 (CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  gxx_prev2 = gxx
  gxy_prev2 = gxy
  gxz_prev2 = gxz
  gyy_prev2 = gyy
  gyz_prev2 = gyz
  gzz_prev2 = gzz
end subroutine CalcK_copy_to_prev2

subroutine CalcK_copy_to_next (CCTK_ARGUMENTS)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  gxx_next = gxx
  gxy_next = gxy
  gxz_next = gxz
  gyy_next = gyy
  gyz_next = gyz
  gzz_next = gzz
end subroutine CalcK_copy_to_next
