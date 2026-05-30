/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Kranc2BSSNMatter_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "Kranc2BSSNMatter";
  CCTK_RegisterBanner(banner);
  return 0;
}
