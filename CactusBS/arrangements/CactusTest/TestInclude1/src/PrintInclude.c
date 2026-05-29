
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusTest_TestInclude1_PrintInclude_c)

#include "IncludeHeader1.h"

void PrintInclude1(CCTK_ARGUMENTS);

void PrintInclude1(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS

  printf("\n\nTesting Cactus Build Include Files from Interface1...\n");

  printf("\nSource should ONLY be included from all active thorns\n");
  printf("\nHeaders will be included from all compiled thorns\n");

  if (CCTK_IsThornActive("Include1"))
  {
    printf("\nThorn Include1 is active\n");
  }
  else
  {
    printf("\nThorn Include1 is not active\n");
  }

  if (CCTK_IsThornActive("Include2"))
  {
    printf("\nThorn Include2 is active\n");
  }
  else
  { 
    printf("\nThorn Include2 is not active\n");
  }

#include "IncludeSource1.h"

#ifdef ih1in1
  printf("\nINCLUDED header from Interface1\n");
#endif

#ifdef ih1in2
  printf("\nINCLUDED header from Interface1\n");
#endif

#ifdef ih2in1
  printf("\nINCLUDED header from Interface2\n");
#endif

#ifdef ih2in2
  printf("\nINCLUDED header from Interface2\n");
#endif

}

  
