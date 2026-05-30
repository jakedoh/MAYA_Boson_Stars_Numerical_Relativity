/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Kranc2BSSN_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "Kranc2BSSN";
  CCTK_RegisterBanner(banner);
  return 0;
}
