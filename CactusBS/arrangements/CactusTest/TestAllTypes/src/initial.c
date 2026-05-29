#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <complex.h>

void TestAllTypes_C(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  *vint = 1;
  *vint1 = 1;
  *vint2 = 1;
  *vint4 = 1;
  *vint8 = 1;
  *vint16 = 1;
  
  *vreal = 1.0;
  *vreal4 = 1.0;
  *vreal8 = 1.0;
  *vreal16 = 1.0;
  
  *vcomplex = 1.0 + 1.0*I;
  *vcomplex8 = 1.0 + 1.0*I;
  *vcomplex16 = 1.0 + 1.0*I;
  *vcomplex32 = 1.0 + 1.0*I;
}
