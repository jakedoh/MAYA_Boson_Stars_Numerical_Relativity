/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ML_CCZ4_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ML_CCZ4";
  CCTK_RegisterBanner(banner);
  return 0;
}
