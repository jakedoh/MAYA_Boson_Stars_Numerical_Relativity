/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_WaveToy_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_WaveToy";
  CCTK_RegisterBanner(banner);
  return 0;
}
