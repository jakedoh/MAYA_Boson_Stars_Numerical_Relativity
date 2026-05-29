/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_ADMConstraints_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_ADMConstraints";
  CCTK_RegisterBanner(banner);
  return 0;
}
