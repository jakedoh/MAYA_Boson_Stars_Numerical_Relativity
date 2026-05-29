/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_WaveToy_CL_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_WaveToy_CL";
  CCTK_RegisterBanner(banner);
  return 0;
}
