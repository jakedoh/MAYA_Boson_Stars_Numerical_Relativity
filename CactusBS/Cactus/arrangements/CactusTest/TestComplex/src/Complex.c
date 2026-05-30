#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Id$";

CCTK_FILEVERSION(CactusTest_TestComplex_Complex_c);

void TestComplexPower(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_COMPLEX num, tmp;
  
  num = CCTK_Cmplx(real_value, imaginary_value);
  
  tmp = CCTK_CmplxPow(num, power_value);
  
  *real_part      = CCTK_CmplxReal(tmp);
  *imaginary_part = CCTK_CmplxImag(tmp); 
}
