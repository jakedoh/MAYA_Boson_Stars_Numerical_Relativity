#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void TestPar_PrintStrings (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TestPar_PrintStrings;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_VInfo (CCTK_THORNSTRING, "String 1: \"%s\"", string1);
  CCTK_VInfo (CCTK_THORNSTRING, "String 2: \"%s\"", string2);
  CCTK_VInfo (CCTK_THORNSTRING, "String 3: \"%s\"", string3);
  CCTK_VInfo (CCTK_THORNSTRING, "String 4: \"%s\"", string4);
  CCTK_VInfo (CCTK_THORNSTRING, "String 5: \"%s\"", string5);
  CCTK_VInfo (CCTK_THORNSTRING, "String 6: \"%s\"", string6);
  CCTK_VInfo (CCTK_THORNSTRING, "String 7: \"%s\"", string7);
  CCTK_VInfo (CCTK_THORNSTRING, "String 8: \"%s\"", string8);
  CCTK_VInfo (CCTK_THORNSTRING, "String 9: \"%s\"", string9);
}
