 /*@@
   @file      IVP_TestData.F90
   @date      Sun Jul  7 15:32:33 2002
   @author    Ian Hawke
   @desc 
   Very simple test data for the IVP solver
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine IVP_TestData(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if (conformal_state .ne. 0) then
    psi = 1.d0
    psix = 0.d0
    psiy = 0.d0
    psiz = 0.d0
    psixx = 0.d0
    psixy = 0.d0
    psixz = 0.d0
    psiyy = 0.d0
    psiyz = 0.d0
    psizz = 0.d0
  endif

  gxx = 1.d0
  gxy = 0.d0
  gxz = 0.d0
  gyy = 1.d0
  gyz = 0.d0
  gzz = 1.d0

  kxx = 0.d0
  kxy = 0.d0
  kxz = 0.d0
  kyy = 0.d0
  kyz = 0.d0
  kzz = 0.d0

  rho = 0.d0
  eps = 0.d0
  
end subroutine IVP_TestData

