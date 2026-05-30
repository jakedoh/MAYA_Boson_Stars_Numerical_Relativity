#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
      
subroutine TestAllTypes_Fortran(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  
  vint = 1
  vint1 = 1
  vint2 = 1
  vint4 = 1
  vint8 = 1
  vint16 = 1
  
  vreal = 1.0d0
  vreal4 = 1.0d0
  vreal8 = 1.0d0
  vreal16 = 1.0d0
  
  vcomplex = (1.0d0, 1.0d0)
  vcomplex8 = (1.0d0, 1.0d0)
  vcomplex16 = (1.0d0, 1.0d0)
  vcomplex32 = (1.0d0, 1.0d0)
end subroutine TestAllTypes_Fortran
