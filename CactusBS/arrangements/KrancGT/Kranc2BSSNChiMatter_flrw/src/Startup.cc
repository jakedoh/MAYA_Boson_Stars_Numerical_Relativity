/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Kranc2BSSNChiMatter_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "Kranc2BSSNChiMatter";
  CCTK_RegisterBanner(banner);
  return 0;
}
