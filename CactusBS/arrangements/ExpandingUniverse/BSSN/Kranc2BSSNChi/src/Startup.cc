/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Kranc2BSSNChi_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "Kranc2BSSNChi";
  CCTK_RegisterBanner(banner);
  return 0;
}
