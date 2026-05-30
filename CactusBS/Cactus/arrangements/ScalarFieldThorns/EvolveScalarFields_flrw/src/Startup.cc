/*  File produced by Kranc */

#include "cctk.h"

extern "C" int EvolveScalarFields_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "EvolveScalarFields";
  CCTK_RegisterBanner(banner);
  return 0;
}
