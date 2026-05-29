#include "cctk.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_IDBrillData_Startup_c)

void BrillData_InitSymBound(CCTK_ARGUMENTS);

void BrillData_InitSymBound(CCTK_ARGUMENTS)
{ 
  DECLARE_CCTK_ARGUMENTS
  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  
  SetCartSymVN(cctkGH,sym,"IDBrillData::brillpsi");
  SetCartSymVN(cctkGH,sym,"IDBrillData::brillMlinear");
  SetCartSymVN(cctkGH,sym,"IDBrillData::brillNsource");

  return;
} 
