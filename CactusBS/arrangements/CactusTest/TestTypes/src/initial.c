 /*@@
   @file      initial.c
   @date      Thu Apr  1 09:20:09 2004
   @author    Erik Schnetter
   @desc 
   
   @enddesc 
   @version $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestTypes_TestTypes_c);

void TestTypes_C (CCTK_ARGUMENTS);

void TestTypes_C (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  *vint = 1;
  *vint1 = 1;
  *vint2 = 1;
  *vint4 = 1;
  *vint8 = 1;
  
  *vreal = 1.0;
  *vreal4 = 1.0;
  *vreal8 = 1.0;
/*   *vreal16 = 1.0; */
  
  *vcomplex = CCTK_Cmplx (1.0, 1.0);
  *vcomplex8 = CCTK_Cmplx8 (1.0, 1.0);
  *vcomplex16 = CCTK_Cmplx16 (1.0, 1.0);
/*   *vcomplex32 = CCTK_Cmplx32 (1.0, 1.0); */
}
